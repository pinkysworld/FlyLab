#!/usr/bin/env python3
"""Recompute VNC-wide target-gene detection from Allen et al. GEO matrices.

This is a cell-level transcript-detection summary, not a cell-type annotation.
It deliberately does not feed these pooled VNC rates into FlyLab's
MaleCNS-superclass gains. The GEO archive contains four 30,000-barcode DGE
matrices; the published nGene, nUMI and mitochondrial cutoffs are applied
before any detection fraction is calculated.

    python scripts/reanalyze_vnc_expression.py /path/to/GSE141807_RAW.tar
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import tarfile
from pathlib import Path
from typing import Any, BinaryIO

import numpy as np
import yaml

ARCHIVE_SHA256 = "bda31982a53a432bb4712ebdaf11ac7d8b9e7be580f72c2acaf7a58912e64daf"
GEO_ACCESSION = "GSE141807"
EXPECTED_MEMBERS = {
    "GSM4213596_female_rep1_dge.txt.gz": ("GSM4213596", "female", 2590),
    "GSM4213597_female_rep2_dge.txt.gz": ("GSM4213597", "female", 9060),
    "GSM4213598_male_rep1_dge.txt.gz": ("GSM4213598", "male", 6522),
    "GSM4213599_male_rep2_dge.txt.gz": ("GSM4213599", "male", 8596),
}

# These are gene-name rules aligned to the expression layer's coarse target
# keys. The Ace group is an acetylcholinesterase target, not a receptor.
TARGET_PREFIXES: dict[str, tuple[re.Pattern[str], ...]] = {
    "insect_nAChR": (re.compile(r"^nAChR", re.I),),
    "insect_RDL": (re.compile(r"^Rdl$", re.I),),
    "insect_GluCl": (re.compile(r"^GluCl", re.I),),
    "insect_OctR": (re.compile(r"^(OctR|Oct)", re.I),),
    "insect_AChE": (re.compile(r"^Ace", re.I),),
}

PAPER_QC = {
    "nGene_gt": 200,
    "nUMI_gt": 1200,
    "nUMI_lte": 10000,
    "prop_mito_lte": 0.15,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def analyze_matrix(stream: BinaryIO, member_name: str) -> dict[str, Any]:
    """Analyze one tab-separated genes-by-barcodes matrix from GEO's tar."""
    try:
        sample_id, sex, published_cells = EXPECTED_MEMBERS[member_name]
    except KeyError as exc:
        raise ValueError(f"unexpected GSE141807 matrix member: {member_name}") from exc

    header = stream.readline().rstrip(b"\r\n").split(b"\t")
    if not header or header[0] != b"GENE":
        raise ValueError(f"{member_name}: expected GENE as the first matrix header")
    n_cells = len(header) - 1
    if n_cells <= 0:
        raise ValueError(f"{member_name}: matrix contains no barcodes")

    n_umi = np.zeros(n_cells, dtype=np.int64)
    n_gene = np.zeros(n_cells, dtype=np.int32)
    n_mito = np.zeros(n_cells, dtype=np.int64)
    family_detected = {key: np.zeros(n_cells, dtype=bool) for key in TARGET_PREFIXES}
    gene_detected: dict[str, np.ndarray] = {}
    mitochondrial_genes: set[str] = set()
    matrix_genes: set[str] = set()
    n_rows = 0

    for line_number, line in enumerate(stream, start=2):
        line = line.rstrip(b"\r\n")
        gene_bytes, separator, counts_bytes = line.partition(b"\t")
        if not separator:
            raise ValueError(f"{member_name}:{line_number}: row has no count columns")
        gene = gene_bytes.decode("utf-8", "replace")
        if gene in matrix_genes:
            raise ValueError(f"{member_name}:{line_number}: duplicate gene symbol {gene!r}")
        matrix_genes.add(gene)
        counts = np.fromstring(counts_bytes, sep="\t", dtype=np.int32)
        if counts.size != n_cells:
            raise ValueError(
                f"{member_name}:{line_number}: expected {n_cells} counts, found {counts.size}"
            )
        if np.any(counts < 0):
            raise ValueError(f"{member_name}:{line_number}: negative UMI count")

        n_rows += 1
        n_umi += counts
        detected = counts > 0
        n_gene += detected
        if gene.startswith("mt:"):
            n_mito += counts
            mitochondrial_genes.add(gene)

        for key, patterns in TARGET_PREFIXES.items():
            if any(pattern.match(gene) for pattern in patterns):
                family_detected[key] |= detected
                if gene not in gene_detected:
                    gene_detected[gene] = detected.copy()
                else:  # keep duplicate gene rows conservative if source changes
                    gene_detected[gene] |= detected

    prop_mito = n_mito / np.maximum(n_umi, 1)
    keep = (
        (n_gene > PAPER_QC["nGene_gt"])
        & (n_umi > PAPER_QC["nUMI_gt"])
        & (n_umi <= PAPER_QC["nUMI_lte"])
        & (prop_mito <= PAPER_QC["prop_mito_lte"])
    )
    kept = int(keep.sum())
    target_summary: dict[str, Any] = {}
    for key in TARGET_PREFIXES:
        genes = {
            gene: {
                "n_detected": int(mask[keep].sum()),
                "fraction_detected": float(mask[keep].mean()) if kept else None,
            }
            for gene, mask in sorted(gene_detected.items())
            if any(pattern.match(gene) for pattern in TARGET_PREFIXES[key])
        }
        target_summary[key] = {
            "n_genes_matched": len(genes),
            "genes": genes,
            "n_cells_with_any_matched_gene_detected": int(family_detected[key][keep].sum()),
            "fraction_cells_with_any_matched_gene_detected": (
                float(family_detected[key][keep].mean()) if kept else None
            ),
        }

    return {
        "sample_id": sample_id,
        "sex": sex,
        "cells_in_matrix": n_cells,
        "genes_in_matrix": n_rows,
        "cells_after_published_qc": kept,
        "published_cells_after_qc": published_cells,
        "cell_count_difference": kept - published_cells,
        "mitochondrial_genes": sorted(mitochondrial_genes),
        "targets": target_summary,
    }


def analyze_archive(archive: Path) -> dict[str, Any]:
    """Validate the source archive and emit its four per-replicate summaries."""
    checksum = _sha256(archive)
    if checksum != ARCHIVE_SHA256:
        raise ValueError(
            f"unexpected GEO archive SHA-256 {checksum}; expected {ARCHIVE_SHA256}"
        )

    results: dict[str, Any] = {}
    with tarfile.open(archive, mode="r:") as tar:
        members = {Path(member.name).name: member for member in tar.getmembers() if member.isfile()}
        if set(members) != set(EXPECTED_MEMBERS):
            raise ValueError(
                "unexpected archive contents: " + ", ".join(sorted(members))
            )
        for name in EXPECTED_MEMBERS:
            member = members[name]
            raw = tar.extractfile(member)
            if raw is None:
                raise ValueError(f"cannot read {name} from GEO archive")
            with gzip.GzipFile(fileobj=raw, mode="rb") as matrix:
                results[EXPECTED_MEMBERS[name][0]] = analyze_matrix(matrix, name)

    all_cells = sum(row["cells_after_published_qc"] for row in results.values())
    all_published = sum(row["published_cells_after_qc"] for row in results.values())
    return {
        "schema_version": 1,
        "source": {
            "accession": GEO_ACCESSION,
            "archive": "GSE141807_RAW.tar",
            "archive_sha256": checksum,
            "geo_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141807",
            "paper_url": "https://elifesciences.org/articles/54074",
        },
        "method": {
            "scope": "adult VNC-wide transcript detection; cells are not assigned to cell classes",
            "qc_thresholds_from_paper": PAPER_QC,
            "detection_rule": "one or more UMI for the gene in a QC-passing barcode",
            "dropout_note": "a zero count is not evidence that a cell does not express the gene",
            "not_used_for": "MaleCNS superclass gains, receptor-complex inference, or MN9 expression",
        },
        "cells_after_qc_total": all_cells,
        "published_cells_after_qc_total": all_published,
        "all_replicate_cell_counts_match_paper": all(
            row["cell_count_difference"] == 0 for row in results.values()
        ),
        "cell_count_reproduction_note": (
            "Three replicate counts match the paper. GSM4213598 retains 6,524 cells versus "
            "6,522 reported (difference +2) under the cited QC thresholds; this small "
            "discrepancy is unresolved and was not adjusted away."
        ),
        "replicates": results,
        "interpretation": (
            "Descriptive VNC-wide RNA detection only. Detection fractions are not receptor "
            "complex prevalence or cell-type-specific expression fractions."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="GSE141807_RAW.tar downloaded from NCBI GEO")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("data/literature/vnc_receptor_detection_GSE141807.yaml"),
        help="output summary YAML (default: %(default)s)",
    )
    args = parser.parse_args(argv)
    result = analyze_archive(args.archive)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(yaml.safe_dump(result, sort_keys=False, allow_unicode=True))
    print(f"wrote {args.out}")
    for row in result["replicates"].values():
        print(
            f"  {row['sample_id']} ({row['sex']}): "
            f"{row['cells_after_published_qc']:,} retained; "
            f"published {row['published_cells_after_qc']:,}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
