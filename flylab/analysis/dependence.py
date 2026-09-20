"""Connectome-Dependence Analysis (CDA): *which* level of network information
does a pharmacological prediction actually need?

The null-model panel in :mod:`flylab.analysis.nullmodels` asks one question per
degradation ("is the real graph distinguishable from this shuffle?").  This module turns that
into the general method the paper leads with.  For one
``compound x concentration x circuit x readout`` it runs the real graph against
progressively degraded versions of itself and reports the **connectome-
dependence profile**

    D = (z_sign, z_weight, z_degree, z_ER)

together with the empirical two-sided permutation probability for each mode -
which is the statistic to read first, since ``z`` on a degraded graph is a
standardised distance from a usually non-normal null, not a probability.

The information ladder
----------------------
Each null keeps a different amount of the real connectome, so a drug effect
that is distinguishable from a given null tells you that *that* level of
information is necessary:

======  ============================  ================================================
rank    graph model (null mode)       what the model still knows
======  ============================  ================================================
0       ``erdos_renyi``               N, E, the weight histogram, the transmitter
                                      census.  Nothing about who connects to whom.
1       ``rewire_degree_preserving``  + every node's in/out degree
2       ``weight_permute``            + the real edge list (weights reshuffled)
3       ``sign_permute``              + the real weights on the real edges
                                      (only which cell carries which transmitter
                                      is destroyed)
4       ``real_connectome``           the MaleCNS cut itself
======  ============================  ================================================

The ladder is ordered by how much of the real cut each model retains; ranks 2
and 3 are not strictly nested (one keeps transmitter identity and loses the
weight-topology pairing, the other the reverse), and that is stated wherever
the level is reported.  :func:`necessary_information_level` returns the
**weakest** model in this ladder that is **not distinguishable from** the real
graph, together with whether that model is *equivalent within tolerance* or
merely *indeterminate*.  It deliberately does not say "reproduces": failing to
reject a null is not evidence that the null is right.

Three-way verdict per mode (the only thing a permutation test can support)
--------------------------------------------------------------------------
For each mode the module reports one of three verdicts, never two:

``distinguishable``
    the empirical two-sided permutation probability is at or below ``alpha``.
    The real graph and this graph model give different drug effects.
``equivalent_within_tolerance``
    the test did not reject **and** the absolute gap between the real effect
    and the null ensemble's median is smaller than a prespecified margin
    ``delta``.  Only this verdict licenses the phrase "this graph model
    reproduces the effect".
``indeterminate``
    the test did not reject and the gap is not smaller than ``delta``: the
    data are consistent with a difference this study cannot resolve.  This is
    the honest home of imidacloprid's sign / weight / degree modes
    (p = 0.275 / 0.586 / 0.472 at n = 1000), which earlier versions of this
    module reported as "reproduces the effect".

``delta`` is prespecified, not fitted.  Its default is
``DEFAULT_DELTA_FRAC`` (5 %) of the **vehicle readout** of the real graph --
the untreated baseline the drug contrast is measured against, which is the one
scale in the cell that is independent of the drug, of the null and of the
sample size.  Five per cent of baseline is one tenth of the 50 % relative
change :mod:`flylab.analysis.selectivity` treats as "the circuit responded", so
the margin is an order of magnitude below the smallest change this pipeline
anywhere calls a circuit effect.  The absolute value actually used is recorded
in every result (``delta``, ``delta_frac``, ``delta_scale``); pass ``delta=``
to override it with an absolute value in readout units.  Where the vehicle
readout is undefined or zero no margin exists, equivalence cannot be claimed,
and every non-rejected mode is ``indeterminate``.

Multiplicity
------------
A single prespecified profile (:func:`dependence_profile` for imidacloprid or
fipronil at ``PAPER_N``) is a **confirmatory** test and needs no correction.
The 84-cell landscape is not: it runs 21 compounds x 4 concentrations x 4
modes, and turning the resulting count into a claim without correction inflates
it.  :func:`dependence_landscape` therefore applies Benjamini-Hochberg FDR
across the structural tests of the run (:data:`STRUCTURAL_MODES` x cells),
reports ``p_adjusted`` / ``q_value`` / ``fdr_alpha`` per cell, classifies on
the adjusted values and keeps the raw ones beside them, and reports **both**
counts so that the change in the headline number is visible.  A landscape is
labelled ``exploratory`` unless its permutation resolution ``1/(n+1)`` is fine
enough for the adjusted threshold ``fdr_alpha / m`` to be reachable by a single
test; when the resolution cannot support any rejection at all, the result says
so instead of silently reporting zero topology-dependent cells.

Classification rule (a label, not a hypothesis test)
----------------------------------------------------
For a cell with real effect ``e`` and per-mode verdicts, with ``alpha``
(default 0.05) and an ``effect_floor`` below which the drug simply does nothing
on this readout:

* ``undefined`` - the readout does not exist on the real graph (e.g. MN9 is
  silenced, so a ratio has no denominator).  No dependence statement is made.
* ``no-effect`` - ``|e| <= effect_floor``.
* ``topology-dependent`` - at least one *structural* null
  (``weight_permute`` or ``rewire_degree_preserving``) is distinguishable from
  the real graph: the prediction needs the wiring pattern, not just the
  graph's size or census.
* ``mixed`` - no structural null is distinguishable but ``sign_permute`` is:
  the prediction needs to know *which* cells carry which transmitter, but not
  the wiring pattern itself.
* ``composition-dominated`` - neither, whether or not ``erdos_renyi`` is
  distinguishable: no structure-preserving degradation could be told apart
  from the real graph on this readout.  Whether that is *equivalence* or
  merely *no resolution* is carried by the per-mode verdicts and summarised in
  ``equivalence``; the class name is not evidence of either.
  ``er_distinguishable = False`` in addition means even a random graph of the
  same size could not be told apart: the prediction is network-insensitive on
  this readout.

The *class* is a label built on top of the verdicts, and within one landscape
it is computed from the FDR-adjusted probabilities.  A cell can only ever be
as trustworthy as its permutation resolution ``1/(n+1)`` allows; cells whose p
had not stabilised are flagged.

Everything here is ``model_derived``: it describes FlyLab's simulation of one
hops-limited MaleCNS cut under a teaching EC50 library, never a measured drug
effect in a fly.
"""
from __future__ import annotations

import time
from typing import Any, Iterable, Sequence

import numpy as np

from flylab.analysis.nullmodels import (
    BASE_WARNINGS,
    DEFAULT_ALPHA,
    DEFAULT_CHECKPOINTS,
    DEFAULT_GRAPH,
    DEFAULT_READOUT,
    DEFAULT_TOL,
    MODES,
    _fast_effect,
    _fast_supported,
    _plan,
    _shuffle_state,
    _state_for,
    convergence_report,
    empirical_p,
    null_distribution,
    p_resolution,
)

__all__ = [
    "MODE_INFORMATION",
    "INFORMATION_LADDER",
    "STRUCTURAL_MODES",
    "CLASSES",
    "VERDICTS",
    "DEFAULT_CONCS",
    "DEFAULT_EFFECT_FLOOR",
    "DEFAULT_DELTA_FRAC",
    "CONFIRMATORY_COMPOUNDS",
    "FAST_N",
    "DEFAULT_N",
    "PAPER_N",
    "Z_KEYS",
    "benjamini_hochberg",
    "equivalence_margin",
    "mode_verdict",
    "dependence_profile",
    "dependence_landscape",
    "necessary_information_level",
    "classify",
    "estimate_landscape_runtime",
]

#: what each null model still knows about the real cut, weakest first
MODE_INFORMATION: dict[str, dict[str, Any]] = {
    "erdos_renyi": {
        "rank": 0,
        "level": "size_and_composition",
        "keeps": "N, E, the weight histogram and the transmitter census",
        "destroys": "every trace of who connects to whom",
    },
    "rewire_degree_preserving": {
        "rank": 1,
        "level": "degree_sequence",
        "keeps": "every node's in/out degree and its transmitter",
        "destroys": "the wiring pattern",
    },
    "weight_permute": {
        "rank": 2,
        "level": "topology_without_weight_pairing",
        "keeps": "the real edge list and the transmitters",
        "destroys": "which synapse count sits on which edge",
    },
    "sign_permute": {
        "rank": 3,
        "level": "wiring_without_transmitter_identity",
        "keeps": "the real edges, the real weights and the NT histogram",
        "destroys": "which cell carries which transmitter",
    },
}

#: modes ordered weakest -> richest, plus the real graph as the top rung
INFORMATION_LADDER: tuple[str, ...] = (
    "erdos_renyi",
    "rewire_degree_preserving",
    "weight_permute",
    "sign_permute",
    "real_connectome",
)

#: nulls that keep the graph's size, degree and census but destroy wiring
STRUCTURAL_MODES: tuple[str, ...] = ("weight_permute", "rewire_degree_preserving")

#: the labels :func:`classify` can return
CLASSES: tuple[str, ...] = (
    "topology-dependent",
    "mixed",
    "composition-dominated",
    "no-effect",
    "undefined",
)

#: the three verdicts a permutation test plus a prespecified margin can
#: support for one mode.  "not significant" is NOT one of them, and neither is
#: anything that means "the shuffle reproduces the effect" on its own.
VERDICTS: tuple[str, ...] = (
    "distinguishable",
    "equivalent_within_tolerance",
    "indeterminate",
)

#: D is reported in this order, as the advisor wrote it
Z_KEYS: dict[str, str] = {
    "sign_permute": "z_sign",
    "weight_permute": "z_weight",
    "rewire_degree_preserving": "z_degree",
    "erdos_renyi": "z_ER",
}

#: default concentration ladder for a landscape (four decades)
DEFAULT_CONCS: tuple[float, ...] = (1e-8, 1e-7, 1e-6, 1e-5)

#: |effect| at or below this counts as "the drug does nothing on this readout"
DEFAULT_EFFECT_FLOOR = 1e-6

#: Prespecified equivalence margin, as a fraction of the **vehicle readout**
#: of the real graph (see the module docstring for the justification): a gap
#: between the real drug effect and the null ensemble's median smaller than
#: 5 % of untreated baseline activity is the largest gap this module is
#: willing to call "the same effect".  It is one tenth of the 50 % relative
#: change :mod:`flylab.analysis.selectivity` treats as a circuit response.
DEFAULT_DELTA_FRAC = 0.05

#: The compounds whose profiles the paper prespecifies as confirmatory tests.
#: A profile is ``confirmatory`` only for one of these, at ``n >= PAPER_N``;
#: everything else (in particular the landscape) is exploratory and is
#: multiplicity-corrected instead.
CONFIRMATORY_COMPOUNDS: tuple[str, ...] = ("imidacloprid", "fipronil")

#: ``--fast``-style n ladder, so one landscape function serves three callers.
#: ``FAST_N`` sits just above the smallest n at which anything can be
#: significant at alpha = 0.05 (n = 19 gives exactly 0.05; FAST_N gives 0.048)
#: and is meant for the browser;
#: ``DEFAULT_N`` is the screening default; ``PAPER_N`` is what a headline claim
#: needs.  Runtimes for the 84-cell default landscape on the ``named`` cut are
#: about 11 s, 55 s and 9 min in process (see
#: :func:`estimate_landscape_runtime`), or roughly a quarter of that with
#: ``n_jobs = 4``.
FAST_N = 20
DEFAULT_N = 100
PAPER_N = 1000

CDA_WARNINGS = [
    "The dependence class is a documented label built on the per-mode "
    "verdicts; within one cell no correction is applied across the four null "
    "modes, and the alpha threshold is a reporting convention.",
    "'Not distinguishable' is NOT 'equivalent': only a mode whose gap from "
    "the null median is below the prespecified margin delta is reported as "
    "equivalent_within_tolerance. Everything else that fails to reject is "
    "indeterminate and is bounded by the permutation resolution 1/(n+1).",
]


def _mode_key(mode: str) -> str:
    return Z_KEYS.get(mode, f"z_{mode}")


# --------------------------------------------------------------------------
# multiplicity: Benjamini-Hochberg over the structural tests of one run
# --------------------------------------------------------------------------
def benjamini_hochberg(
    pvalues: Sequence[float | None],
    alpha: float = DEFAULT_ALPHA,
    resolution: float | None = None,
) -> dict[str, Any]:
    """Benjamini-Hochberg step-up adjusted probabilities (q-values).

    ``adjusted[i]`` is the smallest FDR level at which test ``i`` would be
    rejected, enforced monotone from the largest p downwards, and capped at 1.
    ``None`` entries pass through as ``None`` and are excluded from ``m``.

    The extra diagnostics exist because an *empirical permutation* p cannot go
    below its resolution ``1/(n+1)``, so a coarse run can make rejection
    arithmetically impossible however strong the effect is:

    ``min_rejections``
        the smallest number of tests that must sit at the resolution floor
        before BH can reject anything at all, ``ceil(m * resolution / alpha)``.
    ``can_reject``
        False when that number exceeds ``m``: at this resolution the procedure
        cannot reject any test, and a count of zero would be an artefact of the
        permutation budget rather than a result.
    ``resolution_supports_single_rejection``
        True when ``resolution <= alpha / m``, i.e. one isolated test could be
        rejected on its own.  This is the condition a *confirmatory* landscape
        has to meet.
    """
    import math

    a = float(alpha)
    idx = [i for i, p in enumerate(pvalues) if p is not None]
    m = len(idx)
    adjusted: list[float | None] = [None] * len(pvalues)
    if m:
        order = sorted(idx, key=lambda i: float(pvalues[i]))  # type: ignore[arg-type]
        running = 1.0
        for rank in range(m, 0, -1):
            i = order[rank - 1]
            q = min(running, float(pvalues[i]) * m / rank)  # type: ignore[arg-type]
            running = q
            adjusted[i] = float(min(1.0, q))
    rejected = [(q is not None and q <= a) for q in adjusted]
    res = None if resolution is None else float(resolution)
    min_rej = None
    if res is not None and m:
        min_rej = int(math.ceil(m * res / a)) if a > 0 else None
    return {
        "adjusted": adjusted,
        "rejected": rejected,
        "m": m,
        "alpha": a,
        "n_rejected": int(sum(rejected)),
        "resolution": res,
        "min_rejections": min_rej,
        "can_reject": (True if min_rej is None else bool(min_rej <= m)),
        "resolution_supports_single_rejection": (
            None if res is None or not m else bool(res <= a / m)
        ),
        "method": "benjamini-hochberg (step-up, independent or positively dependent tests)",
    }


# --------------------------------------------------------------------------
# the three-way verdict: distinguishable / equivalent / indeterminate
# --------------------------------------------------------------------------
def equivalence_margin(
    vehicle: float | None,
    delta_frac: float = DEFAULT_DELTA_FRAC,
    delta: float | None = None,
) -> dict[str, Any]:
    """The prespecified equivalence margin for one cell, and its provenance.

    The margin is ``delta_frac * |vehicle readout|`` of the *real* graph unless
    an absolute ``delta`` is passed.  The vehicle readout is the untreated
    baseline the drug contrast is measured against: it is in the units of the
    readout, it does not depend on the drug, on the null ensemble or on the
    number of shuffles, and it is the same scale for every mode of a cell, so
    an equivalence claim made against it is comparable across the ladder.

    Returns ``delta=None`` when no margin can be formed (no vehicle value, or a
    silenced baseline).  A cell with no margin can never report
    ``equivalent_within_tolerance``; its non-rejected modes are
    ``indeterminate``, which is the correct answer rather than a missing one.
    """
    if delta is not None:
        return {
            "delta": abs(float(delta)),
            "delta_frac": None,
            "delta_scale": "absolute (caller-supplied, in readout units)",
            "vehicle": (None if vehicle is None else float(vehicle)),
            "prespecified": True,
        }
    v = None if vehicle is None else abs(float(vehicle))
    if v is None or v <= 0.0:
        return {
            "delta": None,
            "delta_frac": float(delta_frac),
            "delta_scale": "|vehicle readout| (unavailable: baseline is undefined or zero)",
            "vehicle": (None if vehicle is None else float(vehicle)),
            "prespecified": True,
        }
    return {
        "delta": float(delta_frac) * v,
        "delta_frac": float(delta_frac),
        "delta_scale": f"{float(delta_frac):.3g} x |vehicle readout| = {float(delta_frac) * v:.6g}",
        "vehicle": float(vehicle),
        "prespecified": True,
    }


def mode_verdict(
    real_effect: float | None,
    null_median: float | None,
    p: float | None,
    alpha: float = DEFAULT_ALPHA,
    delta: float | None = None,
) -> dict[str, Any]:
    """One mode's three-way verdict - the only claim this design supports.

    ``distinguishable`` when the permutation test rejects; otherwise
    ``equivalent_within_tolerance`` when ``|real - median(null)| < delta``;
    otherwise ``indeterminate``.  A rejected test that *also* sits inside the
    margin is still reported as ``distinguishable`` (the test is the decision
    rule), but ``within_tolerance`` stays True so that "statistically
    distinguishable, practically small" is visible rather than hidden.
    """
    gap = (
        None
        if (real_effect is None or null_median is None)
        else abs(float(real_effect) - float(null_median))
    )
    within = None if (gap is None or delta is None) else bool(gap < float(delta))
    distinguishable = bool(p is not None and float(p) <= float(alpha))
    if distinguishable:
        verdict = "distinguishable"
        why = (
            f"the empirical two-sided permutation p = {float(p):.4g} is at or below "
            f"alpha = {float(alpha):.3g}: the real graph and this graph model give "
            "different drug effects."
        )
    elif within:
        verdict = "equivalent_within_tolerance"
        why = (
            f"the test did not reject (p = {'n/a' if p is None else format(float(p), '.4g')}) "
            f"and the gap between the real effect and the null median, {gap:.4g}, is "
            f"below the prespecified margin delta = {float(delta):.4g}: this graph "
            "model gives the same effect to within the margin."
        )
    else:
        verdict = "indeterminate"
        if gap is None or delta is None:
            why = (
                f"the test did not reject (p = {'n/a' if p is None else format(float(p), '.4g')}) "
                "and no equivalence margin could be formed, so nothing stronger than "
                "'not distinguishable' is supported."
            )
        else:
            why = (
                f"the test did not reject (p = {float(p):.4g}) but the gap between the "
                f"real effect and the null median, {gap:.4g}, is not below the "
                f"prespecified margin delta = {float(delta):.4g}: failure to reject is "
                "not evidence of equivalence."
            )
    return {
        "verdict": verdict,
        "verdict_reason": why,
        "distinguishable": distinguishable,
        "within_tolerance": within,
        "abs_gap_from_null_median": gap,
        "delta": (None if delta is None else float(delta)),
        "p_used": (None if p is None else float(p)),
        "alpha": float(alpha),
    }


def _decision_p(row: dict[str, Any]) -> float | None:
    """The probability a verdict is taken on: FDR-adjusted where one exists."""
    q = row.get("p_adjusted")
    return q if q is not None else row.get("p_two_sided")


def _attach_verdict(
    row: dict[str, Any], alpha: float, delta: float | None
) -> dict[str, Any]:
    """Write the three-way verdict of one mode row back into that row."""
    v = mode_verdict(
        row.get("real_effect"), row.get("null_median"), _decision_p(row), alpha, delta
    )
    row.update(v)
    row["p_basis"] = "fdr_adjusted" if row.get("p_adjusted") is not None else "raw"
    # kept for the CLI / bridge, but it now means exactly "distinguishable"
    row["beats_null"] = v["distinguishable"]
    return row


# --------------------------------------------------------------------------
# classification and the information level
# --------------------------------------------------------------------------
def classify(
    real_effect: float | None,
    per_mode: dict[str, dict[str, Any]],
    alpha: float = DEFAULT_ALPHA,
    effect_floor: float = DEFAULT_EFFECT_FLOOR,
    delta: float | None = None,
) -> dict[str, Any]:
    """Label one cell from its real effect and its per-mode verdicts.

    ``per_mode`` maps a mode name onto a row carrying ``p_two_sided`` (and,
    inside a landscape, the FDR-adjusted ``p_adjusted``, which is what the
    decision is then taken on).  Rows that already carry a ``verdict`` from
    :func:`mode_verdict` are used as they stand; bare ``{"p_two_sided": ...}``
    rows are given a verdict here, with ``delta`` (which may be ``None``) as
    the equivalence margin.

    The rule is the one in the module docstring, applied verbatim so the paper
    can cite this function.  Note what it deliberately does not do: a class of
    ``composition-dominated`` is a statement that no structure-preserving
    degradation could be *told apart* from the real graph, never a statement
    that those degradations reproduce the effect.  Which of them are
    equivalent within tolerance and which are merely indeterminate is in
    ``mode_verdicts`` and summarised in ``equivalence``.
    """
    verdicts: dict[str, str] = {}
    reasons: dict[str, str] = {}
    within: dict[str, bool | None] = {}
    for m, row in per_mode.items():
        if row.get("verdict") in VERDICTS:
            verdicts[m] = str(row["verdict"])
            reasons[m] = str(row.get("verdict_reason") or "")
            within[m] = row.get("within_tolerance")
        else:
            v = mode_verdict(
                row.get("real_effect", real_effect),
                row.get("null_median"),
                _decision_p(row),
                alpha,
                delta,
            )
            verdicts[m] = v["verdict"]
            reasons[m] = v["verdict_reason"]
            within[m] = v["within_tolerance"]

    distinguishable = {m: (v == "distinguishable") for m, v in verdicts.items()}
    structural = [m for m in STRUCTURAL_MODES if distinguishable.get(m)]
    sign_distinguishable = bool(distinguishable.get("sign_permute"))
    er_distinguishable = bool(distinguishable.get("erdos_renyi"))
    equivalent = sorted(m for m, v in verdicts.items() if v == "equivalent_within_tolerance")
    indeterminate = sorted(m for m, v in verdicts.items() if v == "indeterminate")

    if real_effect is None:
        label = "undefined"
        why = (
            "the readout is undefined on the real graph (e.g. the motor neuron "
            "is silenced), so no dependence statement is possible."
        )
    elif abs(float(real_effect)) <= float(effect_floor):
        label = "no-effect"
        why = (
            f"|real effect| {abs(float(real_effect)):.3g} is at or below the "
            f"effect floor {float(effect_floor):.3g}: the compound does not move "
            "this readout, so there is nothing for the wiring to explain."
        )
    elif structural:
        label = "topology-dependent"
        why = (
            "the effect is distinguishable from "
            + " and ".join(structural)
            + ", so it needs the wiring pattern and not only the graph's size, "
            "degree sequence and transmitter census."
        )
    elif sign_distinguishable:
        label = "mixed"
        why = (
            "no structure-preserving null is distinguishable from the real "
            "graph but sign_permute is: the prediction needs to know which "
            "cells carry which transmitter, not the wiring pattern itself."
        )
    else:
        label = "composition-dominated"
        detail = []
        if equivalent:
            detail.append(
                "equivalent within the prespecified margin for "
                + ", ".join(equivalent)
            )
        if indeterminate:
            detail.append(
                "merely indeterminate (not distinguishable, not shown equivalent) for "
                + ", ".join(indeterminate)
            )
        why = (
            "no structure-preserving degradation could be distinguished from "
            "the real graph on this readout, so nothing here requires its "
            "topology"
            + ("; " + "; ".join(detail) if detail else "")
            + (
                ". A same-size random graph IS distinguishable (erdos_renyi), "
                "so the prediction is not wholly network-insensitive."
                if er_distinguishable
                else ". Even a same-size random graph is not distinguishable, "
                "so the prediction is network-insensitive on this readout."
            )
        )
    return {
        "class": label,
        "reason": why,
        "alpha": float(alpha),
        "effect_floor": float(effect_floor),
        "delta": (None if delta is None else float(delta)),
        "mode_verdicts": verdicts,
        "mode_verdict_reasons": reasons,
        "within_tolerance": within,
        "distinguishable": distinguishable,
        "structural_distinguishable": structural,
        "equivalent_modes": equivalent,
        "indeterminate_modes": indeterminate,
        "n_indeterminate": len(indeterminate),
        "er_distinguishable": er_distinguishable,
        "equivalence": (
            "equivalence shown within tolerance for " + ", ".join(equivalent)
            if equivalent
            else "no mode was shown equivalent within tolerance"
        ),
        "network_insensitive": bool(
            label in ("composition-dominated", "no-effect") and not er_distinguishable
        ),
    }


def necessary_information_level(
    profile: dict[str, Any] | Iterable[dict[str, Any]],
    alpha: float | None = None,
    delta: float | None = None,
) -> dict[str, Any]:
    """The weakest graph model that is **not distinguishable** from the real graph.

    ``profile`` is a :func:`dependence_profile` result (or just its ``modes``
    rows).  Walking :data:`INFORMATION_LADDER` from the weakest model upwards,
    the answer is the first model the permutation test could not distinguish
    from the MaleCNS cut.  It is reported together with ``verdict``:

    * ``equivalent_within_tolerance`` - the gap between the real effect and
      that model's null median is below the prespecified margin ``delta``, so
      "this level of information suffices" is supported;
    * ``indeterminate`` - the test simply did not reject.  That is a statement
      about this study's resolution, not about the graph model, and the
      returned ``description`` says so.

    This function used to report "the weakest graph model that already
    reproduces the effect".  It no longer does, and no string it returns says
    "reproduces": failing to reject a null is not evidence of equivalence.
    When every null is distinguishable the answer is ``real_connectome``:
    nothing less than the MaleCNS cut is indistinguishable from it.
    """
    if isinstance(profile, dict):
        rows = list(profile.get("modes") or [])
        a = float(alpha if alpha is not None else profile.get("alpha", DEFAULT_ALPHA))
        real_effect = profile.get("real_effect")
        cls = (profile.get("classification") or {}).get("class")
        d = delta if delta is not None else profile.get("delta")
    else:
        rows = list(profile)
        a = float(alpha if alpha is not None else DEFAULT_ALPHA)
        real_effect = rows[0].get("real_effect") if rows else None
        cls = None
        d = delta
    by_mode = {r["mode"]: r for r in rows}
    warnings: list[str] = []
    floor = float(
        profile.get("effect_floor", DEFAULT_EFFECT_FLOOR)
        if isinstance(profile, dict)
        else DEFAULT_EFFECT_FLOOR
    )

    def verdict_of(row: dict[str, Any]) -> dict[str, Any]:
        if row.get("verdict") in VERDICTS:
            return row
        return {**row, **mode_verdict(
            row.get("real_effect", real_effect), row.get("null_median"),
            _decision_p(row), a, d,
        )}

    if real_effect is not None and abs(float(real_effect)) <= floor:
        return {
            "level": "none_required",
            "mode": None,
            "rank": -1,
            "keeps": "nothing: the compound does not move this readout",
            "destroys": "n/a",
            "p_two_sided": None,
            "verdict": None,
            "equivalent_within_tolerance": None,
            "delta": (None if d is None else float(d)),
            "description": (
                "there is no effect to reproduce at this concentration, so no "
                "level of network information is necessary."
            ),
            "class": cls,
            "warnings": warnings,
        }
    if real_effect is None:
        return {
            "level": None,
            "mode": None,
            "rank": None,
            "verdict": None,
            "equivalent_within_tolerance": None,
            "delta": (None if d is None else float(d)),
            "description": "no effect to compare: the readout is undefined on the real graph.",
            "class": cls,
            "warnings": warnings,
        }

    for mode in INFORMATION_LADDER[:-1]:
        row = by_mode.get(mode)
        if row is None:
            continue
        v = verdict_of(row)
        if v["verdict"] == "distinguishable":
            continue
        info = MODE_INFORMATION[mode]
        # ranks 1-3 are not strictly nested, so a RICHER model can be
        # distinguishable where this one is not. Say so instead of hiding it.
        richer_distinguishable = [
            m
            for m in INFORMATION_LADDER[INFORMATION_LADDER.index(mode) + 1 : -1]
            if m in by_mode and verdict_of(by_mode[m])["verdict"] == "distinguishable"
        ]
        if richer_distinguishable:
            warnings.append(
                f"the ladder is not monotone here: {mode} cannot be "
                "distinguished from the real graph while the richer model(s) "
                + ", ".join(richer_distinguishable)
                + " can, because they destroy different things (weights vs "
                "transmitter identity). Read the level as 'the cheapest graph "
                "model this test cannot tell apart from the real one', not as "
                "'everything above it is indistinguishable too'."
            )
        if row.get("stabilised") is False:
            warnings.append(
                f"{mode}: the permutation p had not stabilised at n="
                f"{row.get('n')}, so this level is provisional."
            )
        if v.get("p_used") is None:
            warnings.append(
                f"{mode}: no permutation p could be formed, so this level rests "
                "on no test at all."
            )
        equivalent = v["verdict"] == "equivalent_within_tolerance"
        if equivalent:
            description = (
                f"a graph model that keeps only {info['keeps']} is equivalent to "
                "the real cut on this readout within the prespecified margin "
                f"delta = {float(v['delta']):.4g} (gap "
                f"{float(v['abs_gap_from_null_median']):.4g}), so nothing more "
                "detailed is necessary for this prediction."
            )
        else:
            description = (
                f"a graph model that keeps only {info['keeps']} could not be "
                "distinguished from the real cut on this readout, but the gap "
                + (
                    f"({float(v['abs_gap_from_null_median']):.4g}) is not below the "
                    f"prespecified margin delta = {float(v['delta']):.4g}"
                    if v.get("abs_gap_from_null_median") is not None and v.get("delta") is not None
                    else "could not be compared against a prespecified margin"
                )
                + ", so this is the weakest level this study can rule out needing, "
                "not a demonstration that the level suffices."
            )
            warnings.append(
                f"{mode} is INDETERMINATE, not equivalent: the necessary "
                "information level is bounded by what n shuffles could resolve. "
                "Do not report it as 'this graph model reproduces the effect'."
            )
        return {
            "level": info["level"],
            "mode": mode,
            "rank": info["rank"],
            "keeps": info["keeps"],
            "destroys": info["destroys"],
            "p_two_sided": row.get("p_two_sided"),
            "p_used": v.get("p_used"),
            "p_basis": row.get("p_basis", "raw"),
            "verdict": v["verdict"],
            "equivalent_within_tolerance": bool(equivalent),
            "abs_gap_from_null_median": v.get("abs_gap_from_null_median"),
            "delta": v.get("delta"),
            "non_monotone": bool(richer_distinguishable),
            "richer_models_distinguishable": richer_distinguishable,
            "description": description,
            "ladder_note": (
                "ranks 2 and 3 are not strictly nested: weight_permute keeps "
                "transmitter identity and loses the weight-topology pairing, "
                "sign_permute the reverse."
            ),
            "class": cls,
            "warnings": warnings,
        }
    missing = [m for m in INFORMATION_LADDER[:-1] if m not in by_mode]
    if missing:
        warnings.append(
            "modes " + ", ".join(missing) + " were not run, so the ladder is incomplete."
        )
    return {
        "level": "real_connectome",
        "mode": "real_connectome",
        "rank": 4,
        "keeps": "the MaleCNS cut itself",
        "destroys": "nothing",
        "p_two_sided": None,
        "verdict": "distinguishable",
        "equivalent_within_tolerance": False,
        "delta": (None if d is None else float(d)),
        "description": (
            "every degradation of the graph was distinguishable from the real "
            "one, so this prediction needs the real MaleCNS wiring."
        ),
        "class": cls,
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# one cell: the dependence profile
# --------------------------------------------------------------------------
def dependence_profile(
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    assay: str = "subgraph",
    readout: str = "auto",
    graph: str | None = None,
    modes: Sequence[str] = MODES,
    n: int = 200,
    seed: int = 0,
    alpha: float = DEFAULT_ALPHA,
    effect_floor: float = DEFAULT_EFFECT_FLOOR,
    delta_frac: float = DEFAULT_DELTA_FRAC,
    delta: float | None = None,
    confirmatory: bool | None = None,
    tol: float = DEFAULT_TOL,
    checkpoints: Sequence[float] = DEFAULT_CHECKPOINTS,
    n_jobs: int = 1,
    **kw: Any,
) -> dict[str, Any]:
    """Connectome-dependence profile for one compound x conc x circuit x readout.

    Runs the real graph against every mode in ``modes`` and returns the profile
    ``D = (z_sign, z_weight, z_degree, z_ER)``, the empirical two-sided
    permutation probability for each mode with its resolution ``1/(n+1)`` and
    stabilisation flag, the **three-way verdict** per mode
    (:data:`VERDICTS`), the classification and the necessary information
    level.

    ``delta_frac`` (or an absolute ``delta``) sets the equivalence margin used
    to separate ``equivalent_within_tolerance`` from ``indeterminate``; see
    :func:`equivalence_margin` and the module docstring.  The absolute value
    used is always recorded in ``delta``.

    ``confirmatory`` marks a prespecified test.  Left at ``None`` it is
    inferred: True for one of :data:`CONFIRMATORY_COMPOUNDS` at
    ``n >= PAPER_N``, False otherwise.  A confirmatory profile is a single
    planned comparison and carries no multiplicity correction; an exploratory
    one says so in its warnings, and the corrected form of the exploratory
    question is :func:`dependence_landscape`.

    Runtime is ``len(modes) * n`` shuffles; on the 1126-node ``named`` cut that
    is about 4.6 s per 100 shuffles across all four modes (the degree-
    preserving rewire is ~90 % of it), so a 1000-shuffle profile takes ~45 s in
    process and ~15 s with ``n_jobs=4``.
    """
    t0 = time.perf_counter()
    if readout == "auto":
        readout = DEFAULT_READOUT[assay]
    graph_name = graph if graph is not None else DEFAULT_GRAPH[assay]
    modes = tuple(modes)
    rows: list[dict[str, Any]] = []
    warnings = list(BASE_WARNINGS) + list(CDA_WARNINGS)
    real_effect: float | None = None
    real_vehicle: float | None = None

    for mode in modes:
        res = null_distribution(
            assay,
            compound,
            conc_M,
            mode,
            n=n,
            seed=seed,
            readout=readout,
            graph=graph_name,
            n_jobs=n_jobs,
            alpha=alpha,
            tol=tol,
            checkpoints=checkpoints,
            **kw,
        )
        real_effect = res["real_effect"]
        real_vehicle = res.get("real_vehicle")
        rows.append(_row_from_null(res))
        for w in res["warnings"]:
            if w not in warnings:
                warnings.append(w)

    if confirmatory is None:
        confirmatory = bool(
            str(compound) in CONFIRMATORY_COMPOUNDS and int(n) >= PAPER_N
        )
    return _assemble_profile(
        compound=compound,
        conc_M=conc_M,
        assay=assay,
        graph_name=graph_name,
        readout=readout,
        n=n,
        seed=seed,
        alpha=alpha,
        effect_floor=effect_floor,
        delta_frac=delta_frac,
        delta=delta,
        real_vehicle=real_vehicle,
        real_effect=real_effect,
        rows=rows,
        warnings=warnings,
        runtime_s=time.perf_counter() - t0,
        confirmatory=bool(confirmatory),
    )


def _row_from_null(res: dict[str, Any]) -> dict[str, Any]:
    """One mode row of a profile, permutation p first."""
    mode = res["mode"]
    return {
        "mode": mode,
        "information_kept": MODE_INFORMATION[mode]["keeps"],
        "rank": MODE_INFORMATION[mode]["rank"],
        "real_effect": res["real_effect"],
        "p_two_sided": res["p_two_sided"],
        "p_resolution": res["p_resolution"],
        "beats_null": res["beats_null"],
        "stabilised": res["stabilised"],
        "n_stabilised": res["n_stabilised"],
        "resolution_limited": res.get("resolution_limited"),
        "verdict_stable": res.get("verdict_stable"),
        "z": res["z"],
        "null_mean": res["null_mean"],
        "null_median": res["null_median"],
        "null_sd": res["null_sd"],
        "n": res["n"],
        "n_ok": res["n_ok"],
        "runtime_s": res.get("runtime_s"),
    }


def _assemble_profile(
    *,
    compound: str,
    conc_M: float,
    assay: str,
    graph_name: str,
    readout: str,
    n: int,
    seed: int,
    alpha: float,
    effect_floor: float,
    real_effect: float | None,
    rows: list[dict[str, Any]],
    warnings: list[str],
    runtime_s: float,
) -> dict[str, Any]:
    by_mode = {r["mode"]: r for r in rows}
    D = {_mode_key(m): (by_mode[m]["z"] if m in by_mode else None) for m in MODES}
    P = {
        _mode_key(m).replace("z_", "p_"): (by_mode[m]["p_two_sided"] if m in by_mode else None)
        for m in MODES
    }
    cls = classify(real_effect, by_mode, alpha=alpha, effect_floor=effect_floor)
    res = p_resolution(min((r["n_ok"] for r in rows), default=0))
    if res > float(alpha):
        warnings.append(
            f"n = {int(n)} gives a permutation resolution of {res:.4f}, which is "
            f"coarser than alpha = {float(alpha)}: NO mode can be beaten at this "
            "n, so the class cannot be 'topology-dependent' whatever the wiring "
            f"does. Use at least n = {int(round(1.0 / float(alpha))) - 1} "
            "(and 500-1000 for a headline claim)."
        )
    profile = {
        "compound": compound,
        "conc_M": float(conc_M),
        "assay": assay,
        "graph": graph_name,
        "readout": readout,
        "n": int(n),
        "seed": int(seed),
        "alpha": float(alpha),
        "effect_floor": float(effect_floor),
        "real_effect": real_effect,
        # permutation probabilities first, then the z profile D
        "p": P,
        "p_resolution": res,
        "resolution_coarser_than_alpha": bool(res > float(alpha)),
        "D": D,
        "D_vector": [D[_mode_key(m)] for m in MODES],
        "D_keys": [_mode_key(m) for m in MODES],
        "modes": rows,
        "classification": cls,
        "class": cls["class"],
        "stabilised": all(bool(r["stabilised"]) for r in rows) if rows else False,
        "verdict_stable": (
            all(r.get("verdict_stable") is not False for r in rows) if rows else False
        ),
        "n_stabilised": (
            max((r["n_stabilised"] or 0) for r in rows) if rows else None
        )
        or None,
        "runtime_s": float(runtime_s),
        "label": "model_derived",
        "warnings": warnings,
    }
    profile["necessary_information_level"] = necessary_information_level(profile)
    lvl = profile["necessary_information_level"]
    for w in lvl.get("warnings", []):
        if w not in profile["warnings"]:
            profile["warnings"].append(w)
    return profile


# --------------------------------------------------------------------------
# the landscape
# --------------------------------------------------------------------------
def estimate_landscape_runtime(
    n_cells: int,
    n: int,
    modes: Sequence[str] = MODES,
    graph: str = "named",
) -> dict[str, Any]:
    """Rough serial runtime for a landscape, from measured per-shuffle costs.

    Measured on the committed cuts (one core): building a shuffle costs about
    0.002 s (sign / weight / ER) or 0.035 s (degree-preserving rewire) on
    ``named``, and 0.02 / 0.35 s on ``taste_motor``; evaluating one cell on an
    existing shuffle costs about 0.0015 s (``named``) or 0.02 s
    (``taste_motor``).  Because :func:`dependence_landscape` reuses one shuffle
    stream across every cell, the cost is
    ``n * sum_modes(shuffle_cost + n_cells * cell_cost)``.
    """
    shuffle_cost = {"named": {"rewire_degree_preserving": 0.035}, "taste_motor": {"rewire_degree_preserving": 0.35}}
    base_shuffle = 0.002 if graph == "named" else 0.02
    cell_cost = 0.0015 if graph == "named" else 0.02
    total = 0.0
    for mode in modes:
        s = shuffle_cost.get(graph, {}).get(mode, base_shuffle)
        total += int(n) * (s + int(n_cells) * cell_cost)
    return {
        "n_cells": int(n_cells),
        "n": int(n),
        "modes": list(modes),
        "graph": graph,
        "estimate_s": float(total),
        "note": (
            "order-of-magnitude only; divide by n_jobs for the parallel path. "
            "Scale n down for an interactive (browser) landscape and up for the "
            "paper's."
        ),
    }


def dependence_landscape(
    compounds: Sequence[str] | None = None,
    concs_M: Sequence[float] = DEFAULT_CONCS,
    assay: str = "subgraph",
    graph: str | None = None,
    readout: str = "auto",
    modes: Sequence[str] = MODES,
    n: int = 100,
    seed: int = 0,
    alpha: float = DEFAULT_ALPHA,
    effect_floor: float = DEFAULT_EFFECT_FLOOR,
    tol: float = DEFAULT_TOL,
    checkpoints: Sequence[float] = DEFAULT_CHECKPOINTS,
    n_jobs: int = 1,
    **kw: Any,
) -> dict[str, Any]:
    """The compound x concentration dependence landscape.

    One cell per ``(compound, conc_M)``: its real effect, the per-mode
    permutation p and z, the dependence class and the necessary information
    level.  ``compounds=None`` uses the whole library.

    The same ``n`` shuffled graphs are reused by every cell (a paired design,
    and the reason this is affordable): shuffle *i* depends only on
    ``(seed, i)``, so each cell's numbers are identical to running
    :func:`dependence_profile` for it on its own.

    Runtime is bounded and reported, and every result carries its own
    ``runtime_s``.  For the default 4-concentration ladder over the
    21-compound library (84 cells) on the ``named`` cut,
    :func:`estimate_landscape_runtime` predicts 11 s at ``n = FAST_N``, 55 s at
    ``n = DEFAULT_N``, 272 s at ``n = 500`` and 545 s at ``n = PAPER_N`` for a
    single idle core.  Measured: the ``n = 500`` landscape took 473 s wall
    clock with ``n_jobs = 2`` on a busy four-core machine.  Scale ``n`` down
    for the browser and up for the paper; ``n`` below 19 cannot resolve
    ``alpha = 0.05`` at all and the result says so.
    """
    t0 = time.perf_counter()
    from flylab.circuit.rate import compute_gains

    if readout == "auto":
        readout = DEFAULT_READOUT[assay]
    graph_name = graph if graph is not None else DEFAULT_GRAPH[assay]
    modes = tuple(modes)
    if compounds is None:
        from flylab.pharm.occupancy import list_compounds

        compounds = list(list_compounds())
    compounds = [str(c) for c in compounds]
    concs = [float(c) for c in concs_M]
    warnings = list(BASE_WARNINGS) + list(CDA_WARNINGS)

    ok_fast, why = _fast_supported(assay, readout, kw)
    cells: list[dict[str, Any]] = []
    if not ok_fast:
        warnings.append(
            f"rate engine unavailable ({why}); the landscape falls back to one "
            "independent dependence_profile per cell, which is much slower."
        )
        for compound in compounds:
            for conc in concs:
                cells.append(
                    dependence_profile(
                        compound,
                        conc,
                        assay=assay,
                        readout=readout,
                        graph=graph_name,
                        modes=modes,
                        n=n,
                        seed=seed,
                        alpha=alpha,
                        effect_floor=effect_floor,
                        tol=tol,
                        checkpoints=checkpoints,
                        n_jobs=n_jobs,
                        **kw,
                    )
                )
        return _landscape_result(
            cells, compounds, concs, assay, graph_name, readout, modes, n, seed,
            alpha, effect_floor, warnings, time.perf_counter() - t0, n_jobs,
        )

    base = _state_for(graph_name)
    specs: list[dict[str, Any]] = []
    for compound in compounds:
        for conc in concs:
            gains, _occ = compute_gains(
                compound, conc, library=kw.get("library"), rule_overrides=kw.get("rule_overrides")
            )
            specs.append({"compound": compound, "conc_M": conc, "gains": dict(gains)})
    plans = [_plan(base, assay, readout, s["gains"], kw) for s in specs]
    reals = [_fast_effect(base, p)["effect"] for p in plans]

    nulls: dict[str, list[list[float | None]]] = {}
    for mode in modes:
        nulls[mode] = _landscape_nulls(
            base=base,
            specs=specs,
            plans=plans,
            mode=mode,
            n=int(n),
            seed=int(seed),
            assay=assay,
            readout=readout,
            graph_name=graph_name,
            kw=kw,
            n_jobs=int(n_jobs),
            warnings=warnings,
        )

    for idx, spec in enumerate(specs):
        rows: list[dict[str, Any]] = []
        for mode in modes:
            draws = nulls[mode][idx]
            rows.append(
                _row_from_draws(
                    mode, reals[idx], draws, n=int(n), alpha=alpha, tol=tol, checkpoints=checkpoints
                )
            )
        cells.append(
            _assemble_profile(
                compound=spec["compound"],
                conc_M=spec["conc_M"],
                assay=assay,
                graph_name=graph_name,
                readout=readout,
                n=int(n),
                seed=int(seed),
                alpha=alpha,
                effect_floor=effect_floor,
                real_effect=reals[idx],
                rows=rows,
                warnings=[BASE_WARNINGS[1], BASE_WARNINGS[3], CDA_WARNINGS[0]],
                runtime_s=0.0,
            )
        )
    return _landscape_result(
        cells, compounds, concs, assay, graph_name, readout, modes, n, seed,
        alpha, effect_floor, warnings, time.perf_counter() - t0, n_jobs,
    )


def _row_from_draws(
    mode: str,
    real: float | None,
    draws: Sequence[float | None],
    n: int,
    alpha: float,
    tol: float,
    checkpoints: Sequence[float],
) -> dict[str, Any]:
    ok = np.asarray([v for v in draws if v is not None], dtype=float)
    n_ok = int(ok.size)
    mean = float(ok.mean()) if n_ok else None
    sd = float(ok.std(ddof=1)) if n_ok > 1 else 0.0
    pp = empirical_p(real, ok, centre=mean) if n_ok else {"p": None, "k": None}
    conv = convergence_report(real, list(draws), tol=tol, alpha=alpha, checkpoints=checkpoints)
    return {
        "mode": mode,
        "information_kept": MODE_INFORMATION[mode]["keeps"],
        "rank": MODE_INFORMATION[mode]["rank"],
        "real_effect": real,
        "p_two_sided": pp["p"],
        "p_resolution": p_resolution(n_ok),
        "beats_null": (pp["p"] is not None and pp["p"] <= float(alpha)),
        "stabilised": bool(conv["stabilised"]),
        "n_stabilised": conv["n_stabilised"],
        "resolution_limited": bool(conv.get("resolution_limited")),
        "verdict_stable": conv.get("verdict_stable"),
        "z": (
            float((real - mean) / sd)
            if (real is not None and mean is not None and sd > 0)
            else None
        ),
        "null_mean": mean,
        "null_median": (float(np.median(ok)) if n_ok else None),
        "null_sd": sd,
        "n": int(n),
        "n_ok": n_ok,
    }


def _landscape_block(payload: tuple[dict[str, Any], list[int]]) -> list[list[float | None]]:
    """Module-level worker: effects of every cell on a block of shuffles."""
    job, indices = payload
    base = _state_for(job["graph"])
    plans = [
        _plan(base, job["assay"], job["readout"], s["gains"], job["kw"]) for s in job["specs"]
    ]
    out: list[list[float | None]] = [[] for _ in plans]
    for i in indices:
        st = _shuffle_state(base, job["mode"], np.random.default_rng([int(job["seed"]), int(i)]))
        for c, plan in enumerate(plans):
            out[c].append(_fast_effect(st, plan)["effect"])
    return out


def _landscape_nulls(
    *,
    base: Any,
    specs: list[dict[str, Any]],
    plans: list[dict[str, Any]],
    mode: str,
    n: int,
    seed: int,
    assay: str,
    readout: str,
    graph_name: str,
    kw: dict[str, Any],
    n_jobs: int,
    warnings: list[str],
) -> list[list[float | None]]:
    """``n`` shuffles, every cell evaluated on each one."""
    if n_jobs <= 1 or n < 2:
        out: list[list[float | None]] = [[] for _ in plans]
        for i in range(n):
            st = _shuffle_state(base, mode, np.random.default_rng([seed, i]))
            for c, plan in enumerate(plans):
                out[c].append(_fast_effect(st, plan)["effect"])
        return out

    job = {
        "graph": graph_name,
        "assay": assay,
        "readout": readout,
        "mode": mode,
        "seed": seed,
        "kw": dict(kw),
        "specs": [{"gains": s["gains"]} for s in specs],
    }
    indices = list(range(n))
    jobs = min(n_jobs, n)
    size = (n + jobs - 1) // jobs
    blocks = [indices[i : i + size] for i in range(0, n, size)]
    try:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=jobs) as pool:
            chunks = list(pool.map(_landscape_block, [(job, b) for b in blocks]))
        out = [[] for _ in plans]
        for chunk in chunks:
            for c, vals in enumerate(chunk):
                out[c].extend(vals)
        note = (
            f"landscape shuffles were run on {jobs} processes; shuffle i is "
            "seeded by (seed, i) alone, so the draws are identical to serial."
        )
        if note not in warnings:
            warnings.append(note)
        return out
    except Exception as exc:  # pragma: no cover - platform dependent
        warnings.append(f"parallel execution failed ({exc!r}); fell back to serial.")
        return _landscape_nulls(
            base=base, specs=specs, plans=plans, mode=mode, n=n, seed=seed,
            assay=assay, readout=readout, graph_name=graph_name, kw=kw,
            n_jobs=1, warnings=warnings,
        )


def _landscape_result(
    cells: list[dict[str, Any]],
    compounds: Sequence[str],
    concs: Sequence[float],
    assay: str,
    graph_name: str,
    readout: str,
    modes: Sequence[str],
    n: int,
    seed: int,
    alpha: float,
    effect_floor: float,
    warnings: list[str],
    runtime_s: float,
    n_jobs: int,
) -> dict[str, Any]:
    counts: dict[str, int] = {c: 0 for c in CLASSES}
    levels: dict[str, int] = {}
    table: list[dict[str, Any]] = []
    for cell in cells:
        counts[cell["class"]] = counts.get(cell["class"], 0) + 1
        lvl = cell["necessary_information_level"]["level"] or "undefined"
        levels[lvl] = levels.get(lvl, 0) + 1
        row = {
            "compound": cell["compound"],
            "conc_M": cell["conc_M"],
            "real_effect": cell["real_effect"],
            "class": cell["class"],
            "necessary_information_level": lvl,
            "stabilised": cell["stabilised"],
        }
        for mode in modes:
            key = _mode_key(mode)
            row[key.replace("z_", "p_")] = cell["p"][key.replace("z_", "p_")]
            row[key] = cell["D"][key]
        table.append(row)
    res = p_resolution(int(n))
    if res > float(alpha):
        warnings.append(
            f"n = {int(n)} gives a permutation resolution of {res:.4f}, coarser "
            f"than alpha = {float(alpha)}: no cell in this landscape can be "
            "classified as topology-dependent at this n. It is a preview, not a "
            "result; the paper's landscape is run at n = 500-1000."
        )
    unstable = [r for r in table if not r["stabilised"]]
    if unstable:
        warnings.append(
            f"{len(unstable)} of {len(table)} cells had at least one mode whose "
            f"permutation p had not stabilised at n = {int(n)}."
        )
    topo = sorted({r["compound"] for r in table if r["class"] == "topology-dependent"})
    comp = sorted({r["compound"] for r in table if r["class"] == "composition-dominated"})
    return {
        "assay": assay,
        "graph": graph_name,
        "readout": readout,
        "compounds": list(compounds),
        "concs_M": list(concs),
        "modes": list(modes),
        "n": int(n),
        "n_jobs": int(n_jobs),
        "seed": int(seed),
        "alpha": float(alpha),
        "effect_floor": float(effect_floor),
        "p_resolution": res,
        "resolution_coarser_than_alpha": bool(res > float(alpha)),
        "shape": [len(compounds), len(concs)],
        "n_cells": len(cells),
        "cells": cells,
        "table": table,
        "summary": {
            "class_counts": counts,
            "level_counts": levels,
            "topology_dependent_compounds": topo,
            "composition_dominated_compounds": comp,
            "n_unstabilised_cells": len(unstable),
        },
        "runtime_s": float(runtime_s),
        "label": "model_derived",
        "warnings": warnings,
    }
