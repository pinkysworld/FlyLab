"""Monte-Carlo uncertainty on log10 EC50.

Library EC50s are order-of-magnitude values, so the honest error bar is on
their logarithm. ``sd_log10=0.3`` means roughly a factor of two either way
(68% interval), which is the spread typically seen between preparations for
the same compound and receptor.
"""

from __future__ import annotations

import copy
from typing import Any

import numpy as np

from flylab.pharm.occupancy import hill_occupancy, load_library

__all__ = ["sample_library", "occupancy_ci"]

#: Rows with these tiers are inert placeholders and are never perturbed.
_FIXED_TIERS = {"class_placeholder"}


def _is_estimate(spec: dict[str, Any]) -> bool:
    return (
        spec.get("evidence_tier") not in _FIXED_TIERS
        and str(spec.get("direction", "none")) != "none"
    )


def sample_library(rng: np.random.Generator, sd_log10: float = 0.3,
                   library: dict[str, Any] | None = None) -> dict[str, Any]:
    """Return a copy of the library with log-normally perturbed EC50s.

    Only real estimates move: ``class_placeholder`` rows and ``direction:
    none`` rows keep their inert values, because perturbing a placeholder
    would invent evidence. Other agents (circuit ensembles) can call this to
    build a library ensemble that matches the occupancy CI.
    """
    if sd_log10 < 0:
        raise ValueError("sd_log10 must be >= 0")
    lib = copy.deepcopy(library or load_library())
    for entry in lib["compounds"].values():
        for spec in entry["receptors"].values():
            if not _is_estimate(spec):
                continue
            shift = float(rng.normal(0.0, sd_log10))
            spec["ec50_M"] = float(spec["ec50_M"]) * (10.0 ** shift)
            spec["ec50_shift_log10"] = shift
    return lib


def occupancy_ci(
    compound: str,
    conc_M: float,
    sd_log10: float = 0.3,
    n: int = 500,
    seed: int = 0,
    library: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Monte-Carlo credible interval on occupancy at one concentration.

    Draws ``n`` log10-normal perturbations of each receptor's EC50 and reports
    the point estimate plus mean and 2.5 / 97.5 percentiles of the resulting
    occupancy. Placeholder rows are reported with a zero-width interval.
    """
    if n < 2:
        raise ValueError("n must be >= 2")
    lib = library or load_library()
    key = compound.lower().strip()
    if key not in lib["compounds"]:
        raise KeyError(f"unknown compound {compound!r}")
    entry = lib["compounds"][key]
    rng = np.random.default_rng(seed)

    receptors: dict[str, Any] = {}
    for receptor, spec in entry["receptors"].items():
        ec50 = float(spec["ec50_M"])
        hill = float(spec.get("n", 1.0))
        point = hill_occupancy(conc_M, ec50, hill)
        if _is_estimate(spec) and sd_log10 > 0:
            draws = ec50 * 10.0 ** rng.normal(0.0, sd_log10, size=n)
            samples = np.array([hill_occupancy(conc_M, float(e), hill) for e in draws])
        else:
            samples = np.full(n, point)
        lo, hi = np.percentile(samples, [2.5, 97.5])
        receptors[receptor] = {
            "point": float(point),
            "mean": float(samples.mean()),
            "p2_5": float(lo),
            "p97_5": float(hi),
            "evidence_tier": spec.get("evidence_tier", "class_placeholder"),
            "perturbed": bool(_is_estimate(spec) and sd_log10 > 0),
        }

    return {
        "compound": entry["name"],
        "key": key,
        "concentration_M": conc_M,
        "n": int(n),
        "sd_log10": float(sd_log10),
        "seed": int(seed),
        "receptors": receptors,
        "method": f"Monte-Carlo over log10 EC50, {n} draws, sd_log10={sd_log10}, seed={seed}",
        "warnings": [
            "Uncertainty is on the library EC50 only; the circuit patch rules are treated as exact.",
            "sd_log10 is a stated assumption, not a fitted standard error.",
        ],
    }
