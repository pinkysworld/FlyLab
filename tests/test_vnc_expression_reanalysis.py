from __future__ import annotations

from io import BytesIO

from scripts.reanalyze_vnc_expression import analyze_matrix


def test_vnc_matrix_qc_and_gene_detection_use_published_cutoffs():
    barcodes = ["pass", "too_few_genes", "at_umi_floor", "high_mito"]
    rows = ["GENE\t" + "\t".join(barcodes)]
    for index in range(200):
        rows.append(f"gene{index:03d}\t1\t1\t1\t1")
    # The nicotinic row brings barcode 1 above 1,200 UMIs. Barcode 2 has only
    # 200 genes; barcode 3 sits exactly at the excluded 1,200-UMI floor.
    rows.append("nAChRalpha1\t1100\t0\t1000\t1001")
    rows.append("mt:ND1\t0\t0\t0\t300")
    matrix = BytesIO(("\n".join(rows) + "\n").encode())

    result = analyze_matrix(matrix, "GSM4213596_female_rep1_dge.txt.gz")

    assert result["cells_in_matrix"] == 4
    assert result["cells_after_published_qc"] == 1
    assert result["published_cells_after_qc"] == 2590
    assert result["cell_count_difference"] == -2589
    assert result["targets"]["insect_nAChR"]["genes"]["nAChRalpha1"] == {
        "n_detected": 1,
        "fraction_detected": 1.0,
    }
    assert result["targets"]["insect_RDL"]["n_genes_matched"] == 0
