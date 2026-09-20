"""Monte-Carlo uncertainty on the log10 of a library potency/affinity value.

Library values are order-of-magnitude numbers, so the honest error bar is on
their logarithm. ``sd_log10=0.3`` means roughly a factor of two either way
(68% interval), which is the spread typically seen between preparations for
the same compound and receptor.

Schema v3: a row with no sourced value (``param_type: unknown``) is **not
modelled**. Perturbing it would invent evidence, so it reports ``None``
everywhere rather than a point value -- the peer review's "Missing evidence
should not become a small quantitative response" applied to the CI as well.
"""

from __future__ import annotations

import copy
from typing import Any

import numpy as np

from flylab.pharm.occupancy import engagement, load_library, spec_value_M

__all__ = ["sample_library", "occupancy_ci"]

#: Rows with these tiers are inert placeholders and are never perturbed.
_FIXED_TIERS = {"class_placeholder"}


def _is_estimate(spec: dict[str, Any]) -> bool:
    return (
        spec.get("evidence_tier") not in _FIXED_TIERS
        and str(spec.get("direction", "none")) != "none"
        and spec_value_M(spec) is not None
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
    # Draw the per-row shifts from a child stream seeded by ONE value taken from
    # ``rng``. The caller's stream then advances by a fixed amount however many
    # rows the library has, so adding a row (schema v3 added the subunit-resolved
    # nicotinic rows) cannot perturb whatever the caller draws next.
    rng = np.random.default_rng(int(rng.integers(0, 2**62)))
    for entry in lib["compounds"].values():
        for spec in entry["receptors"].values():
            if not _is_estimate(spec):
                continue
            shift = float(rng.normal(0.0, sd_log10))
            value = float(spec_value_M(spec)) * (10.0 ** shift)
            spec["value_M"] = value
            spec["ec50_M"] = value  # deprecated alias, kept in step
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
    engagement. Placeholder rows are reported as ``None`` (not modelled), not
    as a point value with a zero-width interval.
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
        value = spec_value_M(spec)
        hill = float(spec.get("n", 1.0))
        param_type = spec.get("param_type")
        if value is None:
            # Not modelled: no point estimate, so no interval either.
            receptors[receptor] = {
                "point": None,
                "mean": None,
                "p2_5": None,
                "p97_5": None,
                "param_type": param_type,
                "engagement_model": "not_modelled",
                "evidence_tier": spec.get("evidence_tier", "class_placeholder"),
                "perturbed": False,
                "not_modelled_reason": "placeholder: no sourced value at this receptor",
            }
            continue
        point = engagement(conc_M, value, hill, param_type=param_type,
                           relation=spec.get("relation"))
        if _is_estimate(spec) and sd_log10 > 0:
            draws = value * 10.0 ** rng.normal(0.0, sd_log10, size=n)
            samples = np.array([
                engagement(conc_M, float(e), hill, param_type=param_type,
                           relation=spec.get("relation"))
                for e in draws
            ])
        else:
            samples = np.full(n, point)
        lo, hi = np.percentile(samples, [2.5, 97.5])
        receptors[receptor] = {
            "point": float(point),
            "mean": float(samples.mean()),
            "p2_5": float(lo),
            "p97_5": float(hi),
            "param_type": param_type,
            "engagement_model": "binding_occupancy" if str(param_type) in ("Kd", "Ki")
            else "functional_engagement",
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
        "method": f"Monte-Carlo over log10 of the sourced value, {n} draws, "
                  f"sd_log10={sd_log10}, seed={seed}",
        "warnings": [
            "Uncertainty is on the library value only; the circuit patch rules are treated as exact.",
            "Rows with no sourced value are reported as None (not modelled), never as a number.",
            "sd_log10 is a stated assumption, not a fitted standard error.",
        ],
    }
