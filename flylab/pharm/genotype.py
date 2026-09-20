"""Resistance alleles as a genotype toggle on the compound library.

Data source: ``data/literature/resistance_alleles.yaml`` (nine published
target-site alleles, compiled with an ``evidence_tier`` and an ``assay_type``
on every number).  Nothing here invents a fold-shift: an allele/compound pair
with no sourced number is left **unchanged** and reported in ``unshifted``.

The one rule this module exists to enforce
------------------------------------------
**Shifts are per compound, never uniform across a receptor.**  Two published
results make a receptor-wide toggle indefensible:

* *Rdl* A301S (*N. lugens*) shifts ethiprole (5.1x) and GABA itself (33.9x) but
  **not** fipronil (1.1x, pIC50 5.74 -> 5.70; Garrood et al. 2017).
* *para* M918T (super-kdr) shifts deltamethrin 100x but leaves DDT essentially
  untouched (fold 1.0; Vais et al. 2000, Usherwood et al. 2005).

So :func:`apply_genotype` multiplies the EC50 of *only* the compounds the
literature actually measured for that allele, and records every other compound
at that receptor in ``unshifted``.

Whole-animal resistance ratios
------------------------------
Rows whose ``assay_type`` is a whole-animal bioassay (``use_as_ec50_shift:
false`` in the YAML) carry penetration, metabolism and excretion inside the
same number.  They are still applied when explicitly requested, but the shifted
row is marked ``shift_kind: "whole_animal_RR"`` and the result carries a
warning that such a ratio overstates the target-site component.

Efficacy
--------
The YAML carries no ``efficacy_scale`` field today.  If one is ever added it is
honoured (it multiplies the row's ``efficacy``); until then, alleles whose notes
record an efficacy effect - super-kdr M918T reduces the number of deltamethrin
binding sites per channel from two to one (Vais et al. 2000), which is a change
in efficacy, not only in affinity - emit a warning saying the pure EC50 shift
under-represents the allele.
"""

from __future__ import annotations

import contextlib
import copy
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable

import yaml

from flylab.pharm.occupancy import compare_compound, load_library

__all__ = [
    "ALLELES_PATH",
    "WILD_TYPE",
    "GENOTYPE_WARNING",
    "load_alleles",
    "list_genotypes",
    "genotype_info",
    "apply_genotype",
    "genotype_occupancy",
    "genotype_assay",
    "genotype_panel",
]

#: Search order for the literature YAML (repo checkout, then CWD).
_DATA_DIRS = [
    Path(__file__).resolve().parents[2] / "data" / "literature",
    Path("data/literature"),
]
ALLELES_PATH = next(
    (d / "resistance_alleles.yaml" for d in _DATA_DIRS if (d / "resistance_alleles.yaml").exists()),
    _DATA_DIRS[0] / "resistance_alleles.yaml",
)

#: The identity genotype.
WILD_TYPE = "wt"

GENOTYPE_WARNING = (
    "Genotype toggle: EC50 shifts are applied per compound from "
    "data/literature/resistance_alleles.yaml. Compounds with no sourced shift "
    "for this allele are left unchanged (see 'unshifted'); FlyLab never "
    "extrapolates a shift across a receptor."
)

WHOLE_ANIMAL_WARNING = (
    "Whole-animal resistance ratio used as a receptor shift: the published "
    "number contains penetration, metabolism and excretion as well as "
    "target-site insensitivity, so it OVERSTATES the target-site change. "
    "Treat the shifted EC50 as an upper bound on the receptor effect."
)

#: assay_type strings that mean "this is not a receptor measurement".
_WHOLE_ANIMAL_ASSAYS = {"whole_animal_bioassay"}

_FOLD_VARIANT_RE = re.compile(r"^fold_(.+)$")


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
@lru_cache(maxsize=4)
def _read_alleles(path_str: str) -> dict[str, Any]:
    p = Path(path_str)
    if not p.exists():  # pragma: no cover - only when the dataset is absent
        return {"alleles": [], "not_sourced": [], "_missing": str(p)}
    data = yaml.safe_load(p.read_text()) or {}
    if not isinstance(data, dict):  # pragma: no cover - defensive
        return {"alleles": [], "not_sourced": []}
    data.setdefault("alleles", [])
    return data


def load_alleles(path: Path | None = None) -> dict[str, Any]:
    """Deep copy of the resistance-allele dataset (defensive: never raises)."""
    return copy.deepcopy(_read_alleles(str(path or ALLELES_PATH)))


def _as_float(value: Any) -> float | None:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    if f <= 0 or f != f:
        return None
    return f


def _norm_compound(name: Any) -> str | None:
    """Library-style key for a compound string, or None if it is not a key."""
    if not isinstance(name, str):
        return None
    key = name.strip().lower()
    if not key or "(" in key or " " in key:
        # e.g. "combination allele VAYA (I161V+G265A+...)" is a genotype label,
        # and "fipronil sulfone" is a metabolite, not a library key.
        return None
    return key.replace("-", "_")


def _shift_rows(allele: dict[str, Any]) -> list[dict[str, Any]]:
    """Flatten one allele's ``shifts`` into ``{variant, compound, fold, ...}``.

    Handles every shape the YAML actually uses: a plain ``fold``; an
    ``allele_variant`` sibling row; a ``per_allele_fold`` mapping (Ace point
    mutations); a ``per_compound_fold`` mapping (the VAYA combination allele);
    and ``fold_<VARIANT>`` keys (the spider-mite GluCl rows).
    """
    rows: list[dict[str, Any]] = []
    allele_assay = allele.get("assay_type")
    allele_use = allele.get("use_as_ec50_shift")
    for shift in allele.get("shifts") or []:
        if not isinstance(shift, dict):  # pragma: no cover - defensive
            continue
        base = {
            "metric": shift.get("metric"),
            "assay_type": shift.get("assay_type") or allele_assay,
            "use_as_ec50_shift": shift.get("use_as_ec50_shift", allele_use),
            "evidence_tier": shift.get("evidence_tier"),
            "source": shift.get("source") or allele.get("source"),
            "notes": shift.get("notes"),
            "quote": shift.get("quote"),
            "preparation": shift.get("preparation") or allele.get("preparation"),
            "efficacy_scale": _as_float(shift.get("efficacy_scale")),
        }
        compound = shift.get("compound")

        per_allele = shift.get("per_allele_fold")
        per_compound = shift.get("per_compound_fold")
        fold_variants = {
            m.group(1): shift[k]
            for k in shift
            if (m := _FOLD_VARIANT_RE.match(str(k))) and k not in ("fold_with_PBO",)
        }

        if isinstance(per_allele, dict):
            for variant, fold in per_allele.items():
                rows.append({**base, "variant": str(variant), "compound": compound, "fold": _as_float(fold)})
        elif isinstance(per_compound, dict):
            variant = _variant_label(compound)
            for comp, fold in per_compound.items():
                rows.append({**base, "variant": variant, "compound": comp, "fold": _as_float(fold)})
        elif fold_variants:
            for variant, fold in fold_variants.items():
                rows.append({**base, "variant": str(variant), "compound": compound, "fold": _as_float(fold)})
        else:
            rows.append(
                {
                    **base,
                    "variant": shift.get("allele_variant"),
                    "compound": compound,
                    "fold": _as_float(shift.get("fold")),
                }
            )
    return rows


def _variant_label(compound: Any) -> str:
    """``"combination allele VAYA (…)"`` -> ``"VAYA"``."""
    text = str(compound or "variant")
    m = re.search(r"\b([A-Z]{2,})\b", text)
    return m.group(1) if m else text.split("(")[0].strip().replace(" ", "_")


def _genotype_index(path: Path | None = None) -> dict[str, dict[str, Any]]:
    """``genotype_id -> {allele metadata, rows}`` including every variant id."""
    data = load_alleles(path)
    out: dict[str, dict[str, Any]] = {}
    for allele in data.get("alleles") or []:
        if not isinstance(allele, dict) or not allele.get("id"):
            continue
        allele_id = str(allele["id"])
        rows = _shift_rows(allele)
        by_variant: dict[str | None, list[dict[str, Any]]] = {}
        for row in rows:
            by_variant.setdefault(row.get("variant"), []).append(row)
        meta = {
            "gene": allele.get("gene"),
            "allele": allele.get("allele"),
            "alternative_numbering": allele.get("alternative_numbering"),
            "common_name": allele.get("common_name"),
            "species": allele.get("species"),
            "target_receptor": allele.get("target_receptor"),
            "genotype_note": allele.get("genotype_note"),
            "allele_notes": _allele_notes(allele),
        }
        base_rows = by_variant.pop(None, [])
        if base_rows:
            out[allele_id] = {"id": allele_id, "variant": None, **meta, "rows": base_rows}
        for variant, vrows in by_variant.items():
            vid = f"{allele_id}:{variant}"
            out[vid] = {
                "id": vid,
                "variant": variant,
                **meta,
                "allele": f"{meta.get('allele')} [{variant}]" if base_rows else str(variant),
                "rows": vrows,
            }
    return out


def _allele_notes(allele: dict[str, Any]) -> str:
    """All free text attached to an allele, for caveat detection."""
    chunks: list[str] = []
    for key in ("notes", "genotype_note", "quote"):
        if allele.get(key):
            chunks.append(str(allele[key]))
    for shift in allele.get("shifts") or []:
        if isinstance(shift, dict):
            for key in ("notes", "quote"):
                if shift.get(key):
                    chunks.append(str(shift[key]))
    return " ".join(chunks)


# ---------------------------------------------------------------------------
# public listing
# ---------------------------------------------------------------------------
def _library_keys(library: dict[str, Any] | None = None) -> set[str]:
    lib = library or load_library()
    return set(lib.get("compounds") or {})


def list_genotypes(library: dict[str, Any] | None = None, path: Path | None = None) -> list[dict[str, Any]]:
    """Every selectable genotype, wild type first.

    Each row: ``id``, ``gene``, ``allele``, ``species``, ``target_receptor``,
    ``compounds_affected`` (library keys only), ``n_rows`` (sourced, usable
    fold-shifts for library compounds) and ``caveats``.
    """
    keys = _library_keys(library)
    out: list[dict[str, Any]] = [
        {
            "id": WILD_TYPE,
            "gene": None,
            "allele": "wild type",
            "species": "Drosophila melanogaster (library default)",
            "target_receptor": None,
            "compounds_affected": [],
            "n_rows": 0,
            "caveats": ["Identity genotype: the library is returned unchanged."],
        }
    ]
    for gid, entry in _genotype_index(path).items():
        usable = [
            row
            for row, _alt in _preferred_rows(entry["rows"])
            if row["fold"] is not None and _norm_compound(row["compound"]) in keys
        ]
        affected = sorted({_norm_compound(r["compound"]) for r in usable if _norm_compound(r["compound"])})
        out.append(
            {
                "id": gid,
                "gene": entry.get("gene"),
                "allele": entry.get("allele"),
                "species": entry.get("species"),
                "target_receptor": entry.get("target_receptor"),
                "compounds_affected": affected,
                "n_rows": len(usable),
                "caveats": _caveats(entry, keys),
            }
        )
    return out


def _caveats(entry: dict[str, Any], keys: set[str]) -> list[str]:
    caveats: list[str] = []
    rows = entry["rows"]
    if any(_is_whole_animal(r) and r["fold"] is not None for r in rows):
        caveats.append(WHOLE_ANIMAL_WARNING)
    null_rows = [str(r["compound"]) for r in rows if r["fold"] is None]
    if null_rows:
        caveats.append(
            "No fold-shift number was sourced for: "
            + ", ".join(sorted(set(null_rows)))
            + " (direction may be known; the number is deliberately null)."
        )
    off_library = sorted({str(r["compound"]) for r in rows if _norm_compound(r["compound"]) not in keys})
    if off_library:
        caveats.append("Measured on compounds outside the FlyLab library: " + ", ".join(off_library) + ".")
    species = str(entry.get("species") or "")
    if species and "Drosophila" not in species:
        caveats.append(f"Cross-species transfer: the shift was measured in {species}.")
    efficacy = _efficacy_caveat(entry)
    if efficacy:
        caveats.append(efficacy)
    return caveats


def _efficacy_caveat(entry: dict[str, Any]) -> str | None:
    notes = str(entry.get("allele_notes") or "")
    if "efficacy" not in notes.lower() and "binding sites" not in notes.lower():
        return None
    if any(r.get("efficacy_scale") for r in entry["rows"]):
        return None
    return (
        "Efficacy limitation: the cited work reports that this allele changes "
        "efficacy (for super-kdr M918T, Vais et al. 2000 report the number of "
        "deltamethrin binding sites per channel falling from two to one), not "
        "only affinity. The YAML carries no efficacy_scale field, so FlyLab "
        "applies the EC50 shift alone and therefore UNDER-represents the "
        "allele; the maximal Nav effect is not capped."
    )


def _preferred_rows(
    rows: Iterable[dict[str, Any]],
) -> list[tuple[dict[str, Any], list[dict[str, Any]]]]:
    """One shift row per compound, plus the rows that were not used.

    A compound can appear more than once for the same allele (the aphid R81T
    entry carries both a receptor binding shift, 50x, and a whole-animal
    resistance ratio, 1679x).  Applying both would multiply them, which is
    nonsense, so the receptor-grade row wins; if only whole-animal rows exist
    the largest-evidence one is used and flagged.  Never more than one shift
    per compound.
    """
    grouped: dict[str, list[dict[str, Any]]] = {}
    order: list[str] = []
    for row in rows:
        key = str(row.get("compound"))
        if key not in grouped:
            grouped[key] = []
            order.append(key)
        grouped[key].append(row)
    out: list[tuple[dict[str, Any], list[dict[str, Any]]]] = []
    for key in order:
        candidates = grouped[key]
        usable = [r for r in candidates if r.get("fold") is not None]
        receptor_grade = [r for r in usable if not _is_whole_animal(r)]
        chosen = (receptor_grade or usable or candidates)[0]
        out.append((chosen, [r for r in candidates if r is not chosen]))
    return out


def _is_whole_animal(row: dict[str, Any]) -> bool:
    if row.get("use_as_ec50_shift") is False:
        return True
    if str(row.get("assay_type") or "") in _WHOLE_ANIMAL_ASSAYS:
        return True
    return "whole_animal" in str(row.get("metric") or "")


def genotype_info(genotype_id: str, path: Path | None = None) -> dict[str, Any]:
    """Metadata + flattened shift rows for one genotype id."""
    if genotype_id == WILD_TYPE:
        return {"id": WILD_TYPE, "allele": "wild type", "rows": [], "target_receptor": None}
    index = _genotype_index(path)
    if genotype_id not in index:
        raise KeyError(
            f"unknown genotype {genotype_id!r}. known: {WILD_TYPE}, " + ", ".join(sorted(index))
        )
    return copy.deepcopy(index[genotype_id])


# ---------------------------------------------------------------------------
# apply
# ---------------------------------------------------------------------------
def _spec_value(spec: dict[str, Any]) -> float | None:
    """Sourced value of a library row in M, or None for a placeholder (v3)."""
    from flylab.pharm.occupancy import spec_value_M

    return spec_value_M(spec)


def apply_genotype(
    library: dict[str, Any] | None,
    genotype_id: str,
    path: Path | None = None,
) -> dict[str, Any]:
    """A NEW library dict with per-compound EC50 shifts applied.

    The returned dict is a normal FlyLab library (pass it straight to
    ``library=`` anywhere) with one extra top-level block, ``genotype``::

        {"id", "allele", "gene", "species", "target_receptor",
         "shifted":   [{compound, receptor, fold, ec50_M_wt, ec50_M_genotype,
                        shift_kind, assay_type, source, evidence_tier}],
         "unshifted": [{compound, receptor, reason}],
         "not_in_library": [...], "warnings": [...]}

    ``library["unshifted"]`` is an alias of ``library["genotype"]["unshifted"]``.

    Only compounds with a sourced fold for this allele move.  Every other
    compound carrying a row at the allele's target receptor is listed in
    ``unshifted`` with the reason - FlyLab does not extrapolate.
    """
    lib = copy.deepcopy(library) if library is not None else load_library()
    compounds = lib.get("compounds") or {}
    if genotype_id == WILD_TYPE:
        lib["genotype"] = {
            "id": WILD_TYPE,
            "allele": "wild type",
            "gene": None,
            "species": None,
            "target_receptor": None,
            "shifted": [],
            "unshifted": [],
            "not_in_library": [],
            "warnings": ["Wild type: the library is returned unchanged (identity genotype)."],
        }
        lib["unshifted"] = []
        return lib

    entry = genotype_info(genotype_id, path)
    receptor = entry.get("target_receptor")
    warnings: list[str] = [GENOTYPE_WARNING]
    shifted: list[dict[str, Any]] = []
    unshifted: list[dict[str, Any]] = []
    not_in_library: list[dict[str, Any]] = []
    seen: set[str] = set()

    for row, alternatives in _preferred_rows(entry["rows"]):
        key = _norm_compound(row.get("compound"))
        if key is None or key not in compounds:
            not_in_library.append(
                {
                    "compound": row.get("compound"),
                    "fold": row.get("fold"),
                    "reason": "compound is not a FlyLab library key; shift recorded but not applied",
                }
            )
            continue
        if row["fold"] is None:
            unshifted.append(
                {
                    "compound": key,
                    "receptor": receptor,
                    "reason": "no fold-shift number sourced for this allele/compound pair "
                    f"({row.get('notes') or row.get('metric') or 'value is null in the literature dataset'})",
                }
            )
            continue
        spec = (compounds[key].get("receptors") or {}).get(receptor)
        if not spec:
            unshifted.append(
                {
                    "compound": key,
                    "receptor": receptor,
                    "reason": f"{key} has no {receptor} row in the library, so the allele cannot act on it",
                }
            )
            continue
        # Schema v3: a placeholder row carries no number, so there is nothing to
        # multiply. Missing evidence must not become a shifted number.
        if _spec_value(spec) is None:
            unshifted.append(
                {
                    "compound": key,
                    "receptor": receptor,
                    "reason": f"{key}:{receptor} is a placeholder row (no sourced value), so the "
                    "allele cannot shift it",
                }
            )
            continue

        fold = float(row["fold"])
        ec50_wt = float(_spec_value(spec))
        spec["value_M"] = ec50_wt * fold
        spec["ec50_M"] = spec["value_M"]
        whole_animal = _is_whole_animal(row)
        shift_kind = "whole_animal_RR" if whole_animal else "receptor_shift"
        shift_block = {
            "allele": entry.get("allele"),
            "genotype_id": genotype_id,
            "fold": fold,
            "source": row.get("source"),
            "assay_type": row.get("assay_type"),
            "shift_kind": shift_kind,
            "metric": row.get("metric"),
            "species": entry.get("species"),
            "evidence_tier_of_shift": row.get("evidence_tier"),
        }
        if whole_animal:
            shift_block["warning"] = WHOLE_ANIMAL_WARNING
        if alternatives:
            shift_block["alternatives_not_applied"] = [
                {
                    "fold": a.get("fold"),
                    "metric": a.get("metric"),
                    "assay_type": a.get("assay_type"),
                    "source": a.get("source"),
                }
                for a in alternatives
            ]
        eff_scale = row.get("efficacy_scale")
        if eff_scale is not None and spec.get("efficacy") is not None:
            spec["efficacy"] = float(spec["efficacy"]) * float(eff_scale)
            shift_block["efficacy_scale"] = float(eff_scale)
        spec["genotype_shift"] = shift_block
        spec["evidence_tier"] = "literature_order"
        spec["source"] = (
            f"{spec.get('source', '')} | genotype {genotype_id} ({entry.get('allele')}): "
            f"EC50 x{fold:g} [{shift_kind}] from {row.get('source')}"
        ).strip(" |")
        shifted.append(
            {
                "compound": key,
                "receptor": receptor,
                "fold": fold,
                "ec50_M_wt": ec50_wt,
                "ec50_M_genotype": spec["value_M"],
                "shift_kind": shift_kind,
                "assay_type": row.get("assay_type"),
                "evidence_tier": row.get("evidence_tier"),
                "source": row.get("source"),
            }
        )
        seen.add(key)
        if whole_animal:
            warnings.append(
                f"{key}: {WHOLE_ANIMAL_WARNING} (allele {entry.get('allele')}, "
                f"metric {row.get('metric')}, fold {fold:g})"
            )

    # Everything else that touches the same receptor is explicitly untouched.
    for key, entry_c in compounds.items():
        if key in seen:
            continue
        if receptor and receptor in (entry_c.get("receptors") or {}):
            spec = entry_c["receptors"][receptor]
            if str(spec.get("direction", "none")) == "none":
                continue
            if any(u["compound"] == key for u in unshifted):
                continue
            unshifted.append(
                {
                    "compound": key,
                    "receptor": receptor,
                    "reason": "no shift measured for this compound with this allele; "
                    "left unchanged (FlyLab never extrapolates a shift across a receptor)",
                }
            )

    if unshifted:
        warnings.append(
            f"{len(unshifted)} compound(s) at {receptor} were left unchanged for {genotype_id}: "
            + ", ".join(sorted({u['compound'] for u in unshifted}))
            + ". Rdl A301S shifts ethiprole and GABA but not fipronil; para M918T shifts "
            "pyrethroids but not DDT - a receptor-wide toggle would be wrong."
        )
    for caveat in _caveats(entry, set(compounds)):
        if caveat not in warnings:
            warnings.append(caveat)

    lib["genotype"] = {
        "id": genotype_id,
        "allele": entry.get("allele"),
        "gene": entry.get("gene"),
        "species": entry.get("species"),
        "target_receptor": receptor,
        "shifted": shifted,
        "unshifted": unshifted,
        "not_in_library": not_in_library,
        "warnings": warnings,
    }
    lib["unshifted"] = unshifted
    return lib


def _genotype_block(lib: dict[str, Any]) -> dict[str, Any]:
    return copy.deepcopy(lib.get("genotype") or {})


# ---------------------------------------------------------------------------
# occupancy / assays
# ---------------------------------------------------------------------------
def genotype_occupancy(
    compound: str,
    conc_M: float,
    genotype_id: str = WILD_TYPE,
    library: dict[str, Any] | None = None,
    path: Path | None = None,
) -> dict[str, Any]:
    """:func:`compare_compound` against the genotype-shifted library.

    Adds a ``genotype`` block (only the rows that concern ``compound``) and
    ``warnings``.  The wild-type result is identical to ``compare_compound``.
    """
    shifted = apply_genotype(library, genotype_id, path=path)
    result = compare_compound(compound, conc_M, library=shifted)
    key = compound.lower().strip()
    block = _genotype_block(shifted)
    block["shifted"] = [r for r in block.get("shifted", []) if r["compound"] == key]
    block["unshifted"] = [r for r in block.get("unshifted", []) if r["compound"] == key]
    block["applies_to_this_compound"] = bool(block["shifted"])
    warnings = list(block.get("warnings", []))
    if genotype_id != WILD_TYPE and not block["shifted"]:
        warnings.insert(
            0,
            f"{key} is NOT shifted by {genotype_id}: no fold-shift for this compound is "
            "sourced for this allele, so its EC50 is the wild-type value.",
        )
    result["genotype"] = block
    result["warnings"] = warnings
    if genotype_id != WILD_TYPE:
        wt = compare_compound(compound, conc_M, library=library)
        by_wt = {r["receptor"]: r for r in wt["receptors"]}
        for row in result["receptors"]:
            ref = by_wt.get(row["receptor"])
            if ref:
                row["engagement_wt"] = ref["engagement"]
                row["occupancy_wt"] = ref["occupancy"]
                if row["engagement"] is not None and ref["engagement"] is not None:
                    delta = row["engagement"] - ref["engagement"]
                else:
                    delta = None  # one side is not modelled
                row["delta_engagement_vs_wt"] = delta
                row["delta_occupancy_vs_wt"] = delta
    return result


#: assays that accept ``library=`` directly.
_LIBRARY_ASSAYS = {"subgraph", "spiking", "taste_map"}
#: assays that read the library through ``compare_compound`` in their module.
_PATCHED_ASSAYS = {"taste": "flylab.assays.taste", "cns": "flylab.assays.wholens"}
ASSAYS = tuple(sorted(_LIBRARY_ASSAYS | set(_PATCHED_ASSAYS)))


@contextlib.contextmanager
def _patched_library(module_name: str, library: dict[str, Any]):
    """Temporarily make one assay module see the shifted library.

    ``assays/taste.py`` and ``assays/wholens.py`` call ``compare_compound``
    without a ``library=`` hook, so - exactly as
    ``flylab.assays.ensemble._taste_library`` already does for its resampled
    libraries - the module-level symbol is swapped for the duration of the run
    and restored in ``finally``.  No file outside this module is edited.
    """
    import importlib

    module = importlib.import_module(module_name)
    original = module.compare_compound

    def _shifted(name, conc, library_kw=None, **kw):
        return compare_compound(name, conc, library=library)

    module.compare_compound = _shifted
    try:
        yield
    finally:
        module.compare_compound = original


def genotype_assay(
    assay: str,
    compound: str | None,
    conc_M: float,
    genotype_id: str = WILD_TYPE,
    library: dict[str, Any] | None = None,
    path: Path | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Run a named assay with the genotype-shifted library.

    ``assay`` is one of ``"subgraph"``, ``"spiking"``, ``"taste_map"``,
    ``"taste"``, ``"cns"``.  The first three take ``library=`` directly; the
    last two are run under :func:`_patched_library`.  The returned notebook
    gains a ``genotype`` block and the genotype warnings.
    """
    if assay not in ASSAYS:
        raise ValueError(f"unknown assay {assay!r}; expected one of {ASSAYS}")
    shifted = apply_genotype(library, genotype_id, path=path)
    block = _genotype_block(shifted)

    if assay == "subgraph":
        from flylab.assays.subgraph import run_subgraph_assay

        nb = run_subgraph_assay(compound, conc_M, library=shifted, **kw)
    elif assay == "spiking":
        from flylab.assays.spiking import run_spiking_assay

        nb = run_spiking_assay(compound, conc_M, library=shifted, **kw)
    elif assay == "taste_map":
        from flylab.assays.taste_map import run_taste_map_assay

        nb = run_taste_map_assay(compound, conc_M, library=shifted, **kw)
    elif assay == "taste":
        from flylab.assays.taste import run_taste_assay

        with _patched_library(_PATCHED_ASSAYS["taste"], shifted):
            nb = run_taste_assay(compound, conc_M, **kw)
    else:  # cns
        from flylab.assays.wholens import run_wholens_assay

        with _patched_library(_PATCHED_ASSAYS["cns"], shifted):
            nb = run_wholens_assay(compound, conc_M, **kw)

    key = (compound or "").lower().strip()
    block["shifted_for_this_compound"] = [r for r in block.get("shifted", []) if r["compound"] == key]
    block["applies_to_this_compound"] = bool(block["shifted_for_this_compound"])
    nb["genotype"] = block
    nb.setdefault("warnings", [])
    for w in block.get("warnings", []):
        if w not in nb["warnings"]:
            nb["warnings"].append(w)
    if genotype_id != WILD_TYPE and not block["applies_to_this_compound"]:
        nb["warnings"].insert(
            0,
            f"{key or 'vehicle'} is NOT shifted by {genotype_id}: this run is "
            "numerically identical to wild type.",
        )
    return nb


# ---------------------------------------------------------------------------
# panel
# ---------------------------------------------------------------------------
def _relevant_genotypes(compound: str, library: dict[str, Any] | None, path: Path | None) -> list[str]:
    key = compound.lower().strip()
    lib = library or load_library()
    entry = (lib.get("compounds") or {}).get(key, {})
    active = {
        name
        for name, spec in (entry.get("receptors") or {}).items()
        if str(spec.get("direction", "none")) != "none"
    }
    affecting: list[str] = []
    same_receptor: list[str] = []
    for row in list_genotypes(library, path):
        if row["id"] == WILD_TYPE:
            continue
        if key in row["compounds_affected"]:
            affecting.append(row["id"])
        elif row["target_receptor"] in active:
            # Kept deliberately: an allele at this compound's receptor with NO
            # sourced shift for it is the result worth showing (Rdl A302S has
            # no fipronil number), not a row to hide.
            same_receptor.append(row["id"])
    return affecting + same_receptor


def _readout_of(assay: str | None, nb: dict[str, Any]) -> dict[str, Any]:
    if assay is None or not nb:
        return {}
    r = nb.get("readouts", {}) or {}
    if assay == "taste":
        return {"mn9_sugar_hz": r.get("mn9_sugar_hz"), "g_ach": r.get("g_ach")}
    if assay == "cns":
        return {
            "cns_excitation_index": (r.get("treated") or {}).get("cns_excitation_index"),
            "delta_excitation": r.get("delta_excitation"),
        }
    return {
        "mn9_hz": r.get("mn9_hz"),
        "dnp01_hz": r.get("dnp01_hz"),
        "mean_hz": r.get("mean_hz"),
    }


def genotype_panel(
    compound: str,
    conc_M: float,
    assay: str | None = "subgraph",
    genotypes: Iterable[str] | None = None,
    library: dict[str, Any] | None = None,
    path: Path | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Wild type vs each relevant allele: fold, occupancy, circuit effect.

    ``assay=None`` gives the occupancy-only panel (fast).  Rows carry the
    fold-shift, the target-receptor occupancy, the circuit readouts and their
    delta against wild type, which is the table the paper needs.
    """
    key = compound.lower().strip()
    ids = list(genotypes) if genotypes is not None else _relevant_genotypes(key, library, path)
    ids = [WILD_TYPE] + [g for g in ids if g != WILD_TYPE]

    rows: list[dict[str, Any]] = []
    warnings: list[str] = [GENOTYPE_WARNING]
    wt_ref: dict[str, Any] | None = None
    for gid in ids:
        occ = genotype_occupancy(key, conc_M, gid, library=library, path=path)
        block = occ["genotype"]
        receptor = block.get("target_receptor")
        by = {r["receptor"]: r for r in occ["receptors"]}
        if receptor not in by:
            # wild type has no target receptor: report the most-occupied insect row
            insect = [
                r for r in occ["receptors"]
                if r["receptor"].startswith("insect_") and r["direction"] != "none"
                and r.get("engagement", r.get("occupancy")) is not None
            ]
            receptor = max(insect, key=lambda r: r["engagement"])["receptor"] if insect else None
        target = by.get(receptor)
        nb = genotype_assay(assay, key, conc_M, gid, library=library, path=path, **kw) if assay else None
        readouts = _readout_of(assay, nb or {})
        shifted_rows = block.get("shifted", [])
        row = {
            "genotype": gid,
            "allele": block.get("allele") or ("wild type" if gid == WILD_TYPE else None),
            "gene": block.get("gene"),
            "species": block.get("species"),
            "receptor": receptor,
            "fold_shift": shifted_rows[0]["fold"] if shifted_rows else (1.0 if gid == WILD_TYPE else None),
            "shift_kind": shifted_rows[0]["shift_kind"] if shifted_rows else None,
            "ec50_M": target["ec50_M"] if target else None,
            "occupancy": target["occupancy"] if target else None,
            "readouts": readouts,
            "shifted": bool(shifted_rows),
            "assay": assay,
        }
        if gid == WILD_TYPE:
            wt_ref = row
        elif wt_ref:
            if row["occupancy"] is not None and wt_ref["occupancy"] is not None:
                row["delta_occupancy_vs_wt"] = row["occupancy"] - wt_ref["occupancy"]
            row["circuit_effect_vs_wt"] = {
                k: (None if readouts.get(k) is None or wt_ref["readouts"].get(k) is None
                    else float(readouts[k]) - float(wt_ref["readouts"][k]))
                for k in readouts
            }
            if not row["shifted"]:
                row["note"] = (
                    f"no sourced shift for {key} with this allele - identical to wild type "
                    "(reported, not extrapolated)"
                )
        for w in occ.get("warnings", []):
            if w not in warnings:
                warnings.append(w)
        rows.append(row)

    return {
        "compound": key,
        "concentration_M": conc_M,
        "assay": assay,
        "genotypes": ids,
        "rows": rows,
        "warnings": warnings,
        "disclaimer": (
            "Genotype panel. Fold-shifts are published per compound; a blank "
            "fold means the literature has no number for that allele/compound "
            "pair and the wild-type EC50 was kept."
        ),
    }
