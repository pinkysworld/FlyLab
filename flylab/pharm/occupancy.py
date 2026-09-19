from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

LIBRARY_PATH = Path(__file__).with_name("library.yaml")


def hill_occupancy(conc_M: float, ec50_M: float, n: float = 1.0) -> float:
    """Fractional occupancy. conc and EC50 in mol/L."""
    if conc_M < 0 or ec50_M <= 0:
        raise ValueError("concentration must be >= 0 and EC50 > 0")
    if conc_M == 0:
        return 0.0
    cn = conc_M**n
    return cn / (ec50_M**n + cn)


def load_library(path: Path | None = None) -> dict[str, Any]:
    with (path or LIBRARY_PATH).open() as fh:
        return yaml.safe_load(fh)


def compare_compound(name: str, conc_M: float, library: dict[str, Any] | None = None) -> dict[str, Any]:
    lib = library or load_library()
    key = name.lower().strip()
    if key not in lib["compounds"]:
        known = ", ".join(sorted(lib["compounds"]))
        raise KeyError(f"unknown compound {name!r}. known: {known}")
    compound = lib["compounds"][key]
    rows = []
    for receptor, spec in compound["receptors"].items():
        occ = hill_occupancy(conc_M, float(spec["ec50_M"]), float(spec.get("n", 1.0)))
        rows.append(
            {
                "receptor": receptor,
                "occupancy": occ,
                "ec50_M": float(spec["ec50_M"]),
                "direction": spec.get("direction", "unknown"),
                "source": spec.get("source", ""),
            }
        )
    return {
        "compound": compound["name"],
        "class": compound.get("class"),
        "concentration_M": conc_M,
        "receptors": rows,
        "disclaimer": (
            "Occupancy only. Circuit injection is Phase 1. "
            "EC50 values are teaching defaults with citations, not fitted FlyLab data."
        ),
    }
