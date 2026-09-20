"""Hill occupancy over the sourced compound library (schema v2).

Backward compatible with v0.4: :func:`hill_occupancy`, :func:`load_library` and
:func:`compare_compound` keep their signatures and their existing output keys.
v0.5 adds ``evidence_tier``/``efficacy`` per row, a ``selectivity`` block, an
occupancy-vs-concentration curve and a library hash for notebook provenance.
"""

from __future__ import annotations

import copy
import hashlib
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

LIBRARY_PATH = Path(__file__).with_name("library.yaml")

#: Row used when a compound does not list a receptor at all.
PLACEHOLDER_ROW: dict[str, Any] = {
    "ec50_M": 0.01,
    "n": 1.0,
    "direction": "none",
    "efficacy": 0.0,
    "source": "class-order placeholder: receptor not listed for this compound",
    "evidence_tier": "class_placeholder",
}

ALLOWED_DIRECTIONS = (
    "agonist",
    "partial_agonist",
    "antagonist",
    "positive_modulator",
    "negative_modulator",
    "inhibitor",
    "none",
)
ALLOWED_TIERS = ("literature_order", "class_placeholder", "measured_fit")

#: Fallback pairs if a library predates the ``selectivity_pairs`` section.
DEFAULT_SELECTIVITY_PAIRS: dict[str, dict[str, str]] = {
    "nAChR": {"insect": "insect_nAChR", "vertebrate": "vertebrate_nAChR_a4b2"},
    "GABA_A": {"insect": "insect_RDL", "vertebrate": "vertebrate_GABA_A"},
    "GluCl": {"insect": "insect_GluCl", "vertebrate": "vertebrate_GlyR"},
    "AChE": {"insect": "insect_AChE", "vertebrate": "vertebrate_AChE"},
    "Nav": {"insect": "insect_Nav", "vertebrate": "vertebrate_Nav1_x"},
}

#: Half-decade grid used by :func:`occupancy_curve` (1e-11 .. 1e-3 M).
DEFAULT_CURVE_CONCS: tuple[float, ...] = tuple(
    round(10 ** (-11 + 0.5 * i), 15) for i in range(17)
)


def hill_occupancy(conc_M: float, ec50_M: float, n: float = 1.0) -> float:
    """Fractional occupancy from the Hill equation."""
    if conc_M < 0 or ec50_M <= 0:
        raise ValueError("concentration must be >= 0 and EC50 > 0")
    if conc_M == 0:
        return 0.0
    cn = conc_M ** n
    return cn / (ec50_M ** n + cn)


@lru_cache(maxsize=4)
def _read_library(path_str: str) -> dict[str, Any]:
    with Path(path_str).open() as fh:
        return yaml.safe_load(fh)


def load_library(path: Path | None = None) -> dict[str, Any]:
    """Load (and cache) the YAML library. Returns a private deep copy."""
    return copy.deepcopy(_read_library(str(path or LIBRARY_PATH)))


def library_sha256(path: Path | None = None) -> str:
    """SHA-256 of the library file on disk, for notebook provenance."""
    return hashlib.sha256((path or LIBRARY_PATH).read_bytes()).hexdigest()


def receptor_table(library: dict[str, Any] | None = None) -> dict[str, Any]:
    """The ``receptors`` section (organism / transmitter / family per target)."""
    lib = library or load_library()
    return lib.get("receptors", {})


def selectivity_pairs(library: dict[str, Any] | None = None) -> dict[str, dict[str, str]]:
    """Insect/vertebrate counterpart pairs used by the scorecard."""
    lib = library or load_library()
    return lib.get("selectivity_pairs") or copy.deepcopy(DEFAULT_SELECTIVITY_PAIRS)


def list_compounds(library: dict[str, Any] | None = None) -> list[str]:
    """Sorted library keys."""
    lib = library or load_library()
    return sorted(lib["compounds"])


def _compound(name: str, lib: dict[str, Any]) -> dict[str, Any]:
    key = name.lower().strip()
    if key not in lib["compounds"]:
        raise KeyError(f"unknown compound {name!r}. known: {', '.join(sorted(lib['compounds']))}")
    return lib["compounds"][key]


def receptor_spec(compound: dict[str, Any], receptor: str) -> dict[str, Any]:
    """Declared spec for ``receptor``, or the inert placeholder row."""
    spec = compound.get("receptors", {}).get(receptor)
    return dict(spec) if spec else dict(PLACEHOLDER_ROW)


def _row(receptor: str, spec: dict[str, Any], conc_M: float) -> dict[str, Any]:
    occ = hill_occupancy(conc_M, float(spec["ec50_M"]), float(spec.get("n", 1.0)))
    row = {
        "receptor": receptor,
        "occupancy": occ,
        "ec50_M": float(spec["ec50_M"]),
        "direction": spec.get("direction", "unknown"),
        "source": spec.get("source", ""),
        "evidence_tier": spec.get("evidence_tier", "class_placeholder"),
        "n": float(spec.get("n", 1.0)),
    }
    if spec.get("efficacy") is not None:
        row["efficacy"] = float(spec["efficacy"])
    return row


def _selectivity(compound: dict[str, Any], conc_M: float, pairs: dict[str, dict[str, str]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for pair_name, pair in pairs.items():
        ins_name, vert_name = pair["insect"], pair["vertebrate"]
        ins = receptor_spec(compound, ins_name)
        vert = receptor_spec(compound, vert_name)
        ins_row = _row(ins_name, ins, conc_M)
        vert_row = _row(vert_name, vert, conc_M)
        ratio = vert_row["ec50_M"] / ins_row["ec50_M"]
        tiers = {ins_row["evidence_tier"], vert_row["evidence_tier"]}
        placeholder = "none" in (ins_row["direction"], vert_row["direction"])
        out[pair_name] = {
            "insect_receptor": ins_name,
            "vertebrate_receptor": vert_name,
            "insect_ec50_M": ins_row["ec50_M"],
            "vertebrate_ec50_M": vert_row["ec50_M"],
            # >1 means the insect target is the more potent one.
            "ec50_ratio_vert_over_insect": ratio,
            "log10_ec50_ratio_vert_over_insect": _log10(ratio),
            "insect_occupancy": ins_row["occupancy"],
            "vertebrate_occupancy": vert_row["occupancy"],
            "occupancy_difference": ins_row["occupancy"] - vert_row["occupancy"],
            "evidence_tier": "class_placeholder" if "class_placeholder" in tiers else sorted(tiers)[0],
            "placeholder": placeholder,
        }
    return out


def _log10(x: float) -> float:
    from math import log10

    return log10(x) if x > 0 else float("-inf")


def compare_compound(name: str, conc_M: float, library: dict[str, Any] | None = None) -> dict[str, Any]:
    """Occupancy of every listed receptor, plus the selectivity scorecard.

    Output keys from v0.4 (``compound``, ``class``, ``concentration_M``,
    ``receptors``, ``disclaimer``) are unchanged; ``receptors`` rows gain
    ``evidence_tier``/``n``/``efficacy`` and the result gains ``key``,
    ``cas``, ``selectivity`` and ``library_version``.
    """
    lib = library or load_library()
    key = name.lower().strip()
    compound = _compound(key, lib)
    rows = [_row(receptor, spec, conc_M) for receptor, spec in compound["receptors"].items()]
    return {
        "compound": compound["name"],
        "key": key,
        "class": compound.get("class"),
        "cas": compound.get("cas"),
        "concentration_M": conc_M,
        "receptors": rows,
        "selectivity": _selectivity(compound, conc_M, selectivity_pairs(lib)),
        "library_version": lib.get("library_version"),
        "notes": compound.get("notes"),
        "disclaimer": "EC50 values are literature-order teaching defaults. Circuit injection currently uses reduced_taste_v0, not a FlyWire/MaleCNS extract.",
    }


def occupancy_curve(
    compound: str,
    concs: list[float] | None = None,
    library: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Per-receptor occupancy across concentrations (default 1e-11..1e-3 M).

    Returns ``[{"conc_M": float, "receptors": {name: occupancy}}, ...]``.
    """
    lib = library or load_library()
    entry = _compound(compound, lib)
    grid = list(concs) if concs else list(DEFAULT_CURVE_CONCS)
    curve = []
    for c in grid:
        occ = {
            receptor: hill_occupancy(c, float(spec["ec50_M"]), float(spec.get("n", 1.0)))
            for receptor, spec in entry["receptors"].items()
        }
        curve.append({"conc_M": float(c), "receptors": occ})
    return curve
