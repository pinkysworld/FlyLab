"""Classical receptor-binding maths: Gaddum, Schild, Black and Leff.

Pure functions, no library access, numpy-friendly (scalars or arrays in, same
shape out). References:

* Gaddum, J.H. (1937) The quantitative effects of antagonistic drugs.
  J Physiol 89:7P. -- competitive shift of the occupancy curve.
* Schild, H.O. (1947) pA, a new scale for the measurement of drug antagonism.
  Br J Pharmacol Chemother 2:189. -- dose ratio r = 1 + [B]/Kb.
* Black, J.W. and Leff, P. (1983) Operational models of pharmacological
  agonism. Proc R Soc Lond B 220:141. -- response from occupancy and efficacy.
"""

from __future__ import annotations

from typing import Any, Iterable, Mapping

import numpy as np

__all__ = [
    "competitive_occupancy",
    "schild_shift",
    "operational_response",
    "effective_ach_gain",
]

ArrayLike = Any


def _positive(x: ArrayLike, what: str) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    if np.any(arr <= 0):
        raise ValueError(f"{what} must be > 0")
    return arr


def competitive_occupancy(
    conc_M: ArrayLike,
    kd_M: float,
    competitor_M: ArrayLike = 0.0,
    competitor_kd_M: float = 1.0,
    n: float = 1.0,
) -> np.ndarray:
    """Fractional occupancy of the agonist in the presence of a competitor.

    Gaddum (1937): the competitor multiplies the apparent Kd by the Schild
    dose ratio, ``Kd_app = Kd * (1 + [B]/Kb)``; occupancy is then the Hill
    expression at that shifted Kd.
    """
    a = np.asarray(conc_M, dtype=float)
    if np.any(a < 0):
        raise ValueError("concentration must be >= 0")
    kd = _positive(kd_M, "kd_M")
    kb = _positive(competitor_kd_M, "competitor_kd_M")
    b = np.asarray(competitor_M, dtype=float)
    if np.any(b < 0):
        raise ValueError("competitor concentration must be >= 0")
    kd_app = kd * (1.0 + b / kb)
    an = np.power(a, n)
    return an / (np.power(kd_app, n) + an)


def schild_shift(antagonist_M: ArrayLike, kb_M: float) -> np.ndarray:
    """Dose ratio ``r = 1 + [B]/Kb`` for a competitive antagonist (Schild 1947).

    ``r`` is the factor by which the agonist EC50 moves right; ``log10(r-1)``
    against ``log10[B]`` is the Schild plot, slope 1 for simple competition.
    """
    kb = _positive(kb_M, "kb_M")
    b = np.asarray(antagonist_M, dtype=float)
    if np.any(b < 0):
        raise ValueError("antagonist concentration must be >= 0")
    return 1.0 + b / kb


def operational_response(
    conc_M: ArrayLike,
    ka_M: float,
    tau: float,
    n: float = 1.0,
) -> np.ndarray:
    """Fractional response of the Black and Leff (1983) operational model.

    ``E/Em = (tau*[A])^n / (([A] + Ka)^n + (tau*[A])^n)``, where ``Ka`` is the
    agonist dissociation constant and ``tau`` the efficacy (receptor reserve)
    parameter: ``tau >> 1`` is a full agonist with spare receptors, ``tau < 1``
    a partial agonist. Returns 0 for ``tau <= 0``.
    """
    a = np.asarray(conc_M, dtype=float)
    if np.any(a < 0):
        raise ValueError("concentration must be >= 0")
    ka = _positive(ka_M, "ka_M")
    if tau <= 0:
        return np.zeros_like(a)
    num = np.power(tau * a, n)
    return num / (np.power(a + ka, n) + num)


def effective_ach_gain(rows: Iterable[Mapping[str, Any]], ach_tone: float = 1.0) -> float:
    """Cholinergic gain from occupancy rows, scaled by an external ACh tone.

    Convenience for callers (e.g. an exposure or ensemble run) that already
    know the ambient acetylcholine level, for instance from a separate
    acetylcholinesterase model. ``ach_tone`` multiplies the mechanism gain and
    the result keeps the 0.05 floor used everywhere else.
    """
    from flylab.pharm.mechanisms import gains_from_occupancy

    if ach_tone <= 0:
        raise ValueError("ach_tone must be > 0")
    g = gains_from_occupancy(list(rows))
    return max(0.05, g["g_ach"] * float(ach_tone))
