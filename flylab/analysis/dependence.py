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
    module reported as though the shuffle had given the same effect.

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

What a negative result here is worth
-----------------------------------
Two things bound it, and both are computed rather than assumed.

*The null has to destroy only what it claims to.*  ``sign_permute`` permutes
the node label array, which preserves the transmitter histogram by **count**
but not by **weight**: out-strength is heterogeneous, so the share of total
outgoing synaptic weight carried by cholinergic cells is 0.619 on the
``named`` cut against 0.545 +- 0.026 over 1000 permutations, putting the real
graph outside its own null.  Since a gain patch acts through transmitter
identity, that null moves the target set and the sign matrix together.
:func:`balance_report` measures this, and
``sign_permute_weight_matched`` - a constrained shuffle that holds each
transmitter's weighted share to within ``DEFAULT_BALANCE_TOL`` while still
relabelling most of the nodes - is the rank-3 rung of the ladder instead.

*The graph has to be rewirable at all.*  :func:`cut_census` reports the
statistics that decide this.  On the ``named`` cut 1154 of 1360 edges (85 %)
terminate on the four seed cells, only 184 of 1126 nodes have any input, and
only 25 % of edges leave a node that itself receives input: it is an in-star,
and a degree-preserving rewire of an in-star is close to the identity.  A
negative topology result there is weak evidence, and the census says so in the
result.  Repeat it on ``taste_motor`` (1841 nodes, 19 066 edges, 97 % of edges
leaving a node with input) before generalising.

Does the instrument work?
-------------------------
A ladder that never detects anything would produce exactly the negatives this
module is used to report, so the ladder is tested against ground truth.
:func:`synthetic_cut` plants a recurrent cholinergic loop of known strength in
a synthetic graph; :func:`ladder_recovery` checks that the ladder calls the
planted effect topology-dependent and the unplanted control not; and
:func:`ladder_power` maps detection rate against planted effect size and
permutation count, with ``loop_strength = 0`` giving the empirical
false-positive rate.

Everything here is ``model_derived``: it describes FlyLab's simulation of one
hops-limited MaleCNS cut under a teaching EC50 library, never a measured drug
effect in a fly.  The synthetic experiments describe the instrument, not a
connectome.
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
    "DEFAULT_MODES",
    "LOCAL_MODES",
    "LADDER_FALLBACKS",
    "BALANCED_TRANSMITTERS",
    "weighted_transmitter_shares",
    "balance_report",
    "DEFAULT_CONCS",
    "DEFAULT_EFFECT_FLOOR",
    "DEFAULT_EFFECT_FLOOR_FRAC",
    "effective_effect_floor",
    "DEFAULT_DELTA_FRAC",
    "CONFIRMATORY_COMPOUNDS",
    "FAST_N",
    "DEFAULT_N",
    "PAPER_N",
    "Z_KEYS",
    "benjamini_hochberg",
    "cut_census",
    "synthetic_cut",
    "ladder_recovery",
    "ladder_power",
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
        "keeps": "the real edges, the real weights and the unweighted NT histogram",
        "destroys": (
            "which cell carries which transmitter AND the pairing between "
            "transmitter and out-strength, so the weighted excitation/"
            "inhibition balance moves too"
        ),
        "joint": True,
        "caveat": (
            "JOINT null: permuting the node label array preserves the NT "
            "histogram by count, not by weight. Out-strength is heterogeneous, "
            "so the share of total outgoing synaptic weight carried by "
            "cholinergic cells moves: on the named cut it is 0.619 on the real "
            "graph against 0.545 +- 0.026 over 1000 permutations, which puts "
            "the real graph outside the whole null. Because a gain patch acts "
            "THROUGH transmitter identity, this null moves the target set and "
            "the sign matrix at once, and the extra variance makes it "
            "conservative. Use sign_permute_weight_matched for a claim about "
            "transmitter identity alone."
        ),
    },
    "sign_permute_weight_matched": {
        "rank": 3,
        "level": "wiring_with_weighted_transmitter_balance",
        "keeps": (
            "the real edges, the real weights, the NT histogram AND the share "
            "of total outgoing synaptic weight carried by each transmitter"
        ),
        "destroys": "which cell carries which transmitter, at matched out-strength",
    },
}

#: modes ordered weakest -> richest, plus the real graph as the top rung.
#: Rank 3 is the **weight-matched** transmitter null: the plain
#: ``sign_permute`` moves the weighted excitation/inhibition balance as well as
#: transmitter identity (see its ``caveat``), so it is a joint null and is not
#: a rung of this ladder.
INFORMATION_LADDER: tuple[str, ...] = (
    "erdos_renyi",
    "rewire_degree_preserving",
    "weight_permute",
    "sign_permute_weight_matched",
    "real_connectome",
)

#: If a rung was not run, the ladder may fall back on these modes for it, with
#: a warning naming the substitution.  ``sign_permute`` is a *weaker* stand-in
#: for the weight-matched null, never an equal one.
LADDER_FALLBACKS: dict[str, tuple[str, ...]] = {
    "sign_permute_weight_matched": ("sign_permute",),
}

#: Null modes this module implements itself, on top of
#: :data:`flylab.analysis.nullmodels.MODES`.
LOCAL_MODES: tuple[str, ...] = ("sign_permute_weight_matched",)

#: the modes a profile or landscape runs unless told otherwise: the four
#: nullmodels degradations plus the weight-matched transmitter null that
#: replaces ``sign_permute`` on the ladder.  ``sign_permute`` is kept in the
#: default set because the contrast between the two is itself informative.
DEFAULT_MODES: tuple[str, ...] = tuple(MODES) + LOCAL_MODES

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
    "sign_permute_weight_matched": "z_sign_wm",
    "weight_permute": "z_weight",
    "rewire_degree_preserving": "z_degree",
    "erdos_renyi": "z_ER",
}

#: default concentration ladder for a landscape (four decades)
DEFAULT_CONCS: tuple[float, ...] = (1e-8, 1e-7, 1e-6, 1e-5)

#: |effect| at or below this counts as "the drug does nothing on this readout".
#: An absolute floor alone is not enough: 1e-6 Hz against a ~6.6 Hz vehicle
#: baseline admits effects at numerical-noise level, and a cell moving the
#: network by 2.7e-05 Hz was being classified and counted like any other.  The
#: effective floor is therefore ``max(DEFAULT_EFFECT_FLOOR,
#: DEFAULT_EFFECT_FLOOR_FRAC * |vehicle readout|)``.
DEFAULT_EFFECT_FLOOR = 1e-6

#: the prespecified **relative** floor: an effect smaller than 1 % of the
#: untreated baseline is not a drug effect for this purpose.  Landscapes report
#: their class counts both with it and with the absolute floor alone, so the
#: cost of the convention is visible.
DEFAULT_EFFECT_FLOOR_FRAC = 0.01

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
# the weight-matched transmitter null (this module's own shuffle)
# --------------------------------------------------------------------------
#: transmitters whose share of total outgoing synaptic weight the weight-
#: matched null holds fixed.  These are the labels the gain patch acts on.
BALANCED_TRANSMITTERS: tuple[str, ...] = ("acetylcholine", "gaba", "glutamate")

#: each held share may drift by at most this much (absolute) from the real
#: graph's.  0.002 = two tenths of a percentage point of total outgoing
#: weight, against the 15.7 points by which a plain ``sign_permute`` moves the
#: cholinergic share on the ``named`` cut.  Measured on that cut, tightening
#: the band from 0.01 to 0.002 costs almost no mixing (56.7 % of nodes
#: relabelled against 50.1 %), and the band matters: the drug effect runs
#: through the weighted cholinergic share, so one percentage point of drift is
#: worth about 0.1 Hz on a readout whose equivalence margin is 0.33 Hz.
DEFAULT_BALANCE_TOL = 0.002

#: accepted label swaps per node.  3 sweeps relabel roughly 57 % of the nodes
#: on the ``named`` cut, which is enough mixing for a transmitter-identity
#: null while the weighted balance never leaves the tolerance band.
DEFAULT_BALANCE_SWEEPS = 3.0


def _out_strength(state: Any) -> "np.ndarray":
    return np.bincount(state.pre, weights=state.w, minlength=state.n).astype(float)


def weighted_transmitter_shares(
    state: Any, transmitters: Sequence[str] = BALANCED_TRANSMITTERS
) -> dict[str, float]:
    """Share of total *outgoing synaptic weight* carried by each transmitter.

    This is the quantity a transmitter null has to preserve if it is to be a
    null about transmitter *identity*: a gain patch scales every synapse made
    by a cell of that transmitter, so what the drug sees is weight, not the
    node count.
    """
    so = _out_strength(state)
    total = float(so.sum()) or 1.0
    labels = list(state.nt.tolist())
    return {
        t: float(sum(so[i] for i, lab in enumerate(labels) if lab == t) / total)
        for t in transmitters
    }


def _weight_matched_sign_permute(
    state: Any,
    rng: Any,
    tol: float = DEFAULT_BALANCE_TOL,
    sweeps: float = DEFAULT_BALANCE_SWEEPS,
    transmitters: Sequence[str] = BALANCED_TRANSMITTERS,
) -> Any:
    """Permute transmitter labels while holding the weighted E/I balance fixed.

    A constrained shuffle rather than a free permutation: starting from the
    real labelling (which satisfies the constraint exactly), random pairs of
    differently-labelled nodes are swapped and a swap is kept only if every
    tracked transmitter's share of total outgoing weight stays within ``tol``
    of the real graph's.  The node count of each label is preserved exactly by
    construction, the weighted share to within ``tol``, and the walk mixes
    because a swap between nodes of similar out-strength is always admissible.

    This is the null ``sign_permute`` was meant to be.  ``sign_permute``
    permutes the label array outright, which also destroys the pairing between
    transmitter and out-strength; the two differ by exactly the amount of
    weighted-balance variance that pairing carries.
    """
    r = np.random.default_rng(rng) if not isinstance(rng, np.random.Generator) else rng
    so = _out_strength(state)
    total = float(so.sum()) or 1.0
    labels = list(state.nt.tolist())
    tracked = {t: i for i, t in enumerate(transmitters)}
    real = [0.0] * len(transmitters)
    for i, lab in enumerate(labels):
        j = tracked.get(lab)
        if j is not None:
            real[j] += so[i]
    real = [v / total for v in real]
    cur = list(real)
    n = int(state.n)
    m = max(1, int(float(sweeps) * n))
    I = r.integers(0, n, m)
    J = r.integers(0, n, m)
    accepted = 0
    for k in range(m):
        i = int(I[k])
        j = int(J[k])
        a = labels[i]
        b = labels[j]
        if a == b:
            continue
        d = (so[j] - so[i]) / total
        ia = tracked.get(a)
        ib = tracked.get(b)
        if ia is not None and abs(cur[ia] + d - real[ia]) > tol:
            continue
        if ib is not None and abs(cur[ib] - d - real[ib]) > tol:
            continue
        labels[i], labels[j] = b, a
        if ia is not None:
            cur[ia] += d
        if ib is not None:
            cur[ib] -= d
        accepted += 1
    nt = np.empty(n, dtype=object)
    for i, lab in enumerate(labels):
        nt[i] = lab
    relabelled = int(sum(1 for i, lab in enumerate(labels) if lab != state.nt[i]))
    return state.replace(
        nt=nt,
        meta={
            "mode": "sign_permute_weight_matched",
            "swaps_accepted": accepted,
            "swaps_proposed": m,
            "nodes_relabelled": relabelled,
            "balance_tol": float(tol),
            "weighted_shares": {t: cur[k] for t, k in tracked.items()},
        },
    )


def _shuffle_state_any(base: Any, mode: str, rng: Any, **kw: Any) -> Any:
    """:func:`flylab.analysis.nullmodels._shuffle_state` plus this module's modes."""
    if mode == "sign_permute_weight_matched":
        return _weight_matched_sign_permute(
            base,
            rng,
            tol=float(kw.get("balance_tol", DEFAULT_BALANCE_TOL)),
            sweeps=float(kw.get("balance_sweeps", DEFAULT_BALANCE_SWEEPS)),
        )
    return _shuffle_state(base, mode, rng)


def balance_report(
    graph: str = "named",
    n: int = 200,
    seed: int = 0,
    transmitters: Sequence[str] = BALANCED_TRANSMITTERS,
    tol: float = DEFAULT_BALANCE_TOL,
) -> dict[str, Any]:
    """How far each transmitter null moves the weighted E/I balance.

    The evidence behind :data:`MODE_INFORMATION`'s ``caveat`` on
    ``sign_permute``, computed rather than asserted, and the check that the
    weight-matched null does what it claims.
    """
    base = _state_for(graph)
    real = weighted_transmitter_shares(base, transmitters)
    out: dict[str, Any] = {
        "graph": graph,
        "n": int(n),
        "transmitters": list(transmitters),
        "real": real,
        "tol": float(tol),
        "modes": {},
    }
    for mode in ("sign_permute", "sign_permute_weight_matched"):
        draws = {t: [] for t in transmitters}
        for i in range(int(n)):
            st = _shuffle_state_any(base, mode, np.random.default_rng([int(seed), i]))
            for t, v in weighted_transmitter_shares(st, transmitters).items():
                draws[t].append(v)
        rows = {}
        for t in transmitters:
            arr = np.asarray(draws[t], dtype=float)
            rows[t] = {
                "real": real[t],
                "null_mean": float(arr.mean()),
                "null_sd": float(arr.std(ddof=1)) if arr.size > 1 else 0.0,
                "null_min": float(arr.min()),
                "null_max": float(arr.max()),
                "max_abs_deviation": float(np.max(np.abs(arr - real[t]))),
                "percentile_of_real": float(100.0 * np.mean(arr < real[t])),
                "within_tol": bool(np.max(np.abs(arr - real[t])) <= float(tol)),
            }
        out["modes"][mode] = rows
    out["preserves_weighted_balance"] = {
        mode: all(r["within_tol"] for r in rows.values())
        for mode, rows in out["modes"].items()
    }
    out["label"] = "model_derived"
    return out


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


def effective_effect_floor(
    vehicle: float | None,
    effect_floor: float = DEFAULT_EFFECT_FLOOR,
    effect_floor_frac: float = DEFAULT_EFFECT_FLOOR_FRAC,
) -> dict[str, Any]:
    """The floor below which a drug effect is not an effect, and why.

    ``max(absolute, frac * |vehicle readout|)``.  The relative term is what
    keeps numerical noise out of the counts; the absolute term keeps a floor
    when the baseline itself is tiny or undefined.
    """
    rel = (
        None
        if (vehicle is None or float(effect_floor_frac) <= 0.0)
        else float(effect_floor_frac) * abs(float(vehicle))
    )
    eff = float(effect_floor) if rel is None else max(float(effect_floor), rel)
    return {
        "effect_floor": eff,
        "effect_floor_absolute": float(effect_floor),
        "effect_floor_frac": float(effect_floor_frac),
        "effect_floor_relative": rel,
        "binding": (
            "relative" if (rel is not None and rel >= float(effect_floor)) else "absolute"
        ),
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

    substitutions: list[str] = []
    for rung in INFORMATION_LADDER[:-1]:
        mode = rung
        row = by_mode.get(rung)
        if row is None:
            for alt in LADDER_FALLBACKS.get(rung, ()):
                if alt in by_mode:
                    mode, row = alt, by_mode[alt]
                    substitutions.append(f"{rung} -> {alt}")
                    warnings.append(
                        f"the ladder's {rung} rung was not run; {alt} was used "
                        "in its place. "
                        + str(MODE_INFORMATION.get(alt, {}).get("caveat") or "")
                    )
                    break
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
            for m in INFORMATION_LADDER[INFORMATION_LADDER.index(rung) + 1 : -1]
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
                "Report it as 'not distinguishable from the real graph', never "
                "as 'this graph model gives the same effect'."
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
            "ladder_substitutions": substitutions,
            "description": description,
            "ladder_note": (
                "ranks 2 and 3 are not strictly nested: weight_permute keeps "
                "transmitter identity and loses the weight-topology pairing, "
                "sign_permute the reverse."
            ),
            "class": cls,
            "warnings": warnings,
        }
    missing = [
        m
        for m in INFORMATION_LADDER[:-1]
        if m not in by_mode and not any(a in by_mode for a in LADDER_FALLBACKS.get(m, ()))
    ]
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
    modes: Sequence[str] = DEFAULT_MODES,
    n: int = 200,
    seed: int = 0,
    alpha: float = DEFAULT_ALPHA,
    effect_floor: float = DEFAULT_EFFECT_FLOOR,
    effect_floor_frac: float = DEFAULT_EFFECT_FLOOR_FRAC,
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

    local = [m for m in modes if m in LOCAL_MODES]
    if local:
        ok_fast, why = _fast_supported(assay, readout, kw)
        if not ok_fast:
            warnings.append(
                "the rate engine is unavailable here (" + why + "), so "
                + ", ".join(local)
                + " could not be run; the ladder is missing its rank-3 rung."
            )
            modes = tuple(m for m in modes if m not in LOCAL_MODES)

    for mode in modes:
        if mode in LOCAL_MODES:
            draws, real = _local_null_draws(
                assay=assay,
                compound=compound,
                conc_M=conc_M,
                readout=readout,
                graph_name=graph_name,
                mode=mode,
                n=int(n),
                seed=int(seed),
                kw=kw,
            )
            real_effect = real["effect"]
            real_vehicle = real["vehicle"]
            rows.append(
                _row_from_draws(
                    mode, real_effect, draws, n=int(n), alpha=alpha, tol=tol,
                    checkpoints=checkpoints,
                )
            )
            continue
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
    rows.sort(key=lambda r: [m for m in modes].index(r["mode"]))

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
        effect_floor_frac=effect_floor_frac,
        delta_frac=delta_frac,
        delta=delta,
        real_vehicle=real_vehicle,
        real_effect=real_effect,
        rows=rows,
        warnings=warnings,
        runtime_s=time.perf_counter() - t0,
        confirmatory=bool(confirmatory),
    )


def _local_null_draws(
    *,
    assay: str,
    compound: str,
    conc_M: float,
    readout: str,
    graph_name: str,
    mode: str,
    n: int,
    seed: int,
    kw: dict[str, Any],
) -> tuple[list[float | None], dict[str, Any]]:
    """``n`` draws for one of this module's own null modes, on the fast engine.

    ``flylab.analysis.nullmodels.null_distribution`` refuses a mode it does
    not know, so the modes this module adds are run here instead, through the
    same ``_state_for`` / ``_plan`` / ``_fast_effect`` path and with the same
    ``(seed, i)`` per-shuffle seeding, so a draw is identical to the one the
    landscape would take.
    """
    from flylab.circuit.rate import compute_gains

    base = _state_for(graph_name)
    gains = kw.get("gains")
    if gains is None:
        gains, _occ = compute_gains(
            compound, conc_M, library=kw.get("library"),
            rule_overrides=kw.get("rule_overrides"),
        )
    plan = _plan(base, assay, readout, dict(gains), kw)
    real = _fast_effect(base, plan)
    draws: list[float | None] = []
    for i in range(int(n)):
        st = _shuffle_state_any(base, mode, np.random.default_rng([int(seed), i]), **kw)
        draws.append(_fast_effect(st, plan)["effect"])
    return draws, real


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
    effect_floor_frac: float = DEFAULT_EFFECT_FLOOR_FRAC,
    delta_frac: float = DEFAULT_DELTA_FRAC,
    delta: float | None = None,
    real_vehicle: float | None = None,
    confirmatory: bool = False,
) -> dict[str, Any]:
    margin = equivalence_margin(real_vehicle, delta_frac=delta_frac, delta=delta)
    floors = effective_effect_floor(real_vehicle, effect_floor, effect_floor_frac)
    effect_floor = floors["effect_floor"]
    for row in rows:
        _attach_verdict(row, alpha, margin["delta"])
    by_mode = {r["mode"]: r for r in rows}
    # the four nullmodels degradations are always keyed, run or not, so a
    # consumer can rely on p_sign / p_weight / p_degree / p_ER; any further
    # mode that was run is keyed beside them
    keyed = list(MODES) + [m for m in by_mode if m not in MODES]
    D = {_mode_key(m): (by_mode[m]["z"] if m in by_mode else None) for m in keyed}
    P = {
        _mode_key(m).replace("z_", "p_"): (by_mode[m]["p_two_sided"] if m in by_mode else None)
        for m in keyed
    }
    Q = {
        _mode_key(m).replace("z_", "q_"): (by_mode[m].get("p_adjusted") if m in by_mode else None)
        for m in keyed
    }
    cls = classify(
        real_effect, by_mode, alpha=alpha, effect_floor=effect_floor, delta=margin["delta"]
    )
    res = p_resolution(min((r["n_ok"] for r in rows), default=0))
    if res > float(alpha):
        warnings.append(
            f"n = {int(n)} gives a permutation resolution of {res:.4f}, which is "
            f"coarser than alpha = {float(alpha)}: NO mode can be distinguished "
            "from the real graph at this n, so the class cannot be "
            "'topology-dependent' whatever the wiring does. Use at least n = "
            f"{int(round(1.0 / float(alpha))) - 1} (and 500-1000 for a headline "
            "claim)."
        )
    if margin["delta"] is None:
        warnings.append(
            "no equivalence margin could be formed (the vehicle readout is "
            "undefined or zero), so no mode can be reported as equivalent "
            "within tolerance; every mode that fails to reject is "
            "'indeterminate'."
        )
    if cls["indeterminate_modes"]:
        warnings.append(
            "indeterminate (neither distinguishable nor equivalent within "
            f"delta = {margin['delta'] if margin['delta'] is None else format(margin['delta'], '.4g')}): "
            + ", ".join(cls["indeterminate_modes"])
            + ". These modes support 'not distinguishable from this null "
            "ensemble' and nothing stronger."
        )
    if not confirmatory:
        warnings.append(
            "exploratory: this profile is not one of the prespecified "
            "confirmatory tests (" + ", ".join(CONFIRMATORY_COMPOUNDS)
            + f" at n >= {PAPER_N}) and carries no multiplicity correction of "
            "its own. Counting verdicts over many such profiles needs the "
            "FDR-controlled dependence_landscape instead."
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
        "effect_floor_absolute": floors["effect_floor_absolute"],
        "effect_floor_frac": floors["effect_floor_frac"],
        "effect_floor_relative": floors["effect_floor_relative"],
        "effect_floor_binding": floors["binding"],
        "real_effect": real_effect,
        "real_vehicle": real_vehicle,
        # the prespecified equivalence margin, as an absolute readout value
        "delta": margin["delta"],
        "delta_frac": margin["delta_frac"],
        "delta_scale": margin["delta_scale"],
        "confirmatory": bool(confirmatory),
        "multiplicity_correction": None,
        "design": "confirmatory (prespecified)" if confirmatory else "exploratory",
        # permutation probabilities first, then the z profile D
        "p": P,
        "q": Q,
        "p_resolution": res,
        "resolution_coarser_than_alpha": bool(res > float(alpha)),
        "D": D,
        "D_vector": [D[_mode_key(m)] for m in MODES],
        "D_keys": [_mode_key(m) for m in MODES],
        "modes_run": [r["mode"] for r in rows],
        "modes": rows,
        "verdicts": cls["mode_verdicts"],
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


def _refinalise(
    cell: dict[str, Any], alpha: float, effect_floor: float, delta: float | None
) -> dict[str, Any]:
    """Recompute one cell's verdicts, class and level on the p now in its rows.

    Used by :func:`dependence_landscape` after Benjamini-Hochberg has written
    ``p_adjusted`` into the structural rows: the decision is then taken on the
    adjusted probabilities, while the raw ones stay in ``p_two_sided``.
    """
    for row in cell["modes"]:
        _attach_verdict(row, alpha, delta)
    by_mode = {r["mode"]: r for r in cell["modes"]}
    cls = classify(
        cell["real_effect"], by_mode, alpha=alpha, effect_floor=effect_floor, delta=delta
    )
    cell["classification"] = cls
    cell["class"] = cls["class"]
    cell["verdicts"] = cls["mode_verdicts"]
    cell["q"] = {
        _mode_key(m).replace("z_", "q_"): (by_mode[m].get("p_adjusted") if m in by_mode else None)
        for m in list(MODES) + [m for m in by_mode if m not in MODES]
    }
    cell["necessary_information_level"] = necessary_information_level(cell)
    return cell


# --------------------------------------------------------------------------
# does the instrument work? ground-truth recovery on a synthetic cut
# --------------------------------------------------------------------------
#: transmitter composition of the synthetic cut, close to the ``named`` cut's
SYNTHETIC_COMPOSITION: dict[str, float] = {
    "acetylcholine": 0.62,
    "gaba": 0.30,
    "glutamate": 0.08,
}


def synthetic_cut(
    n_nodes: int = 250,
    n_edges: int = 1500,
    n_seeds: int = 4,
    loop_size: int = 10,
    loop_strength: float = 1.0,
    seed: int = 0,
    composition: dict[str, float] | None = None,
) -> dict[str, Any]:
    """A graph with a **known** topology-dependent effect planted in it.

    Background: ``n_edges`` uniformly random directed edges over ``n_nodes``
    with log-normal weights and a transmitter composition close to the real
    cut's.  Planted structure: a directed cycle of ``loop_size`` cholinergic
    nodes that passes through the driven seed nodes, carrying
    ``loop_strength`` times the mean background weight.

    The cycle is recurrent cholinergic gain, so a patch that collapses
    ``g_ach`` removes an amplification that only exists while the cycle is
    intact.  Degree-preserving rewiring redistributes those edges and destroys
    it, so the drug *contrast* -- not merely the rate -- depends on the wiring.
    ``loop_strength = 0`` plants nothing and is the negative control: the same
    graph, the same composition, no topology dependence to find.

    Nothing here touches the connectome or the compound library; the
    experiment is about the instrument, not about flies.
    """
    rng = np.random.default_rng(int(seed))
    comp = dict(composition or SYNTHETIC_COMPOSITION)
    names = list(comp)
    probs = np.asarray([comp[k] for k in names], dtype=float)
    probs = probs / probs.sum()
    labels = rng.choice(names, size=int(n_nodes), p=probs)
    nodes = [
        {"bodyId": 1_000_000 + i, "consensus_nt": str(labels[i])}
        for i in range(int(n_nodes))
    ]
    pre = rng.integers(0, int(n_nodes), int(n_edges))
    post = rng.integers(0, int(n_nodes), int(n_edges))
    w = np.exp(rng.normal(1.2, 0.8, int(n_edges)))
    keep = pre != post
    pre, post, w = pre[keep], post[keep], w[keep]
    mean_w = float(w.mean())
    edges = [
        {"pre": 1_000_000 + int(a), "post": 1_000_000 + int(b), "weight": float(v)}
        for a, b, v in zip(pre, post, w)
    ]
    seed_idx = list(range(int(n_seeds)))
    if float(loop_strength) > 0 and int(loop_size) >= 2:
        ring = seed_idx + [
            int(i) for i in rng.choice(
                np.arange(int(n_seeds), int(n_nodes)),
                size=max(0, int(loop_size) - len(seed_idx)),
                replace=False,
            )
        ]
        for i in ring:
            nodes[i]["consensus_nt"] = "acetylcholine"
        for a, b in zip(ring, ring[1:] + ring[:1]):
            edges.append(
                {
                    "pre": 1_000_000 + int(a),
                    "post": 1_000_000 + int(b),
                    "weight": float(loop_strength) * mean_w * 6.0,
                }
            )
    else:
        ring = []
    seeds = {
        "MN9": [1_000_000 + seed_idx[0]],
        "DNp01": [1_000_000 + i for i in seed_idx[1:]] or [1_000_000 + seed_idx[0]],
    }
    return {
        "nodes": nodes,
        "edges": edges,
        "seeds": seeds,
        "meta": {
            "synthetic": True,
            "loop_strength": float(loop_strength),
            "loop_nodes": [1_000_000 + int(i) for i in ring],
            "n_nodes": int(n_nodes),
            "n_edges": len(edges),
            "seed": int(seed),
        },
    }


def _ladder_on_state(
    state: Any,
    gains: dict[str, float],
    modes: Sequence[str] = DEFAULT_MODES,
    n: int = 200,
    seed: int = 0,
    alpha: float = DEFAULT_ALPHA,
    effect_floor: float = DEFAULT_EFFECT_FLOOR,
    effect_floor_frac: float = DEFAULT_EFFECT_FLOOR_FRAC,
    delta_frac: float = DEFAULT_DELTA_FRAC,
    readout: str = "mean_hz",
    label: str = "synthetic",
    **kw: Any,
) -> dict[str, Any]:
    """Run the whole ladder on an arbitrary graph state and one gain patch.

    The same shuffles, the same empirical p, the same three-way verdict and
    the same classification as :func:`dependence_profile`; only the source of
    the graph and of the gain patch differ, which is what lets the instrument
    be pointed at a graph whose answer is known in advance.
    """
    t0 = time.perf_counter()
    plan = _plan(state, "subgraph", readout, dict(gains), kw)
    real = _fast_effect(state, plan)
    rows: list[dict[str, Any]] = []
    for mode in modes:
        draws: list[float | None] = []
        for i in range(int(n)):
            st = _shuffle_state_any(state, mode, np.random.default_rng([int(seed), i]), **kw)
            draws.append(_fast_effect(st, plan)["effect"])
        rows.append(
            _row_from_draws(
                mode, real["effect"], draws, n=int(n), alpha=alpha,
                tol=DEFAULT_TOL, checkpoints=DEFAULT_CHECKPOINTS,
            )
        )
    return _assemble_profile(
        compound=label,
        conc_M=0.0,
        assay="subgraph",
        graph_name=label,
        readout=readout,
        n=int(n),
        seed=int(seed),
        alpha=alpha,
        effect_floor=effect_floor,
        effect_floor_frac=effect_floor_frac,
        delta_frac=delta_frac,
        real_vehicle=real["vehicle"],
        real_effect=real["effect"],
        rows=rows,
        warnings=[
            "synthetic ground truth: this graph is generated, not measured, "
            "and exists to test the instrument rather than to describe a fly.",
        ],
        runtime_s=time.perf_counter() - t0,
        confirmatory=True,
    )


def ladder_recovery(
    strengths: Sequence[float] = (0.0, 0.5, 1.0, 2.0),
    n: int = 200,
    seed: int = 0,
    gains: dict[str, float] | None = None,
    modes: Sequence[str] = DEFAULT_MODES,
    n_nodes: int = 250,
    n_edges: int = 1500,
    loop_size: int = 10,
    **kw: Any,
) -> dict[str, Any]:
    """Can the ladder find topology dependence that is known to be there?

    Without this, a negative result from the ladder is indistinguishable from
    an under-powered test.  For each planted ``loop_strength`` the full ladder
    is run on a synthetic cut and asked for its class: ``0.0`` is the negative
    control and must NOT come back topology-dependent; the positive strengths
    must.

    Returns one row per strength with the class, the necessary information
    level, the per-mode probabilities and the real effect, plus a ``recovered``
    verdict for the whole experiment.
    """
    from flylab.analysis.nullmodels import _GraphState

    patch = dict(gains or {"g_ach": 0.05})
    t0 = time.perf_counter()
    rows: list[dict[str, Any]] = []
    for strength in strengths:
        graph = synthetic_cut(
            n_nodes=n_nodes, n_edges=n_edges, loop_size=loop_size,
            loop_strength=float(strength), seed=int(seed),
        )
        state = _GraphState.from_graph(graph)
        prof = _ladder_on_state(
            state, patch, modes=modes, n=int(n), seed=int(seed),
            label=f"planted@{float(strength):g}", **kw
        )
        rows.append(
            {
                "loop_strength": float(strength),
                "planted": bool(float(strength) > 0),
                "class": prof["class"],
                "topology_dependent": prof["class"] == "topology-dependent",
                "necessary_information_level": prof["necessary_information_level"]["level"],
                "real_effect": prof["real_effect"],
                "real_vehicle": prof["real_vehicle"],
                "p": dict(prof["p"]),
                "verdicts": dict(prof["verdicts"]),
                "n": int(n),
            }
        )
    negatives = [r for r in rows if not r["planted"]]
    positives = [r for r in rows if r["planted"]]
    false_positive = any(r["topology_dependent"] for r in negatives)
    detected = [r for r in positives if r["topology_dependent"]]
    return {
        "design": (
            "planted recurrent cholinergic loop on a synthetic cut; the gain "
            "patch collapses g_ach, so the drug contrast depends on the loop "
            "surviving the shuffle"
        ),
        "gains": patch,
        "n": int(n),
        "seed": int(seed),
        "n_nodes": int(n_nodes),
        "n_edges": int(n_edges),
        "loop_size": int(loop_size),
        "modes": list(modes),
        "rows": rows,
        "n_positive": len(positives),
        "n_detected": len(detected),
        "false_positive_on_control": bool(false_positive),
        "recovered": bool(detected and not false_positive),
        "smallest_detected_strength": (
            min(r["loop_strength"] for r in detected) if detected else None
        ),
        "statement": (
            f"the ladder recovered {len(detected)} of {len(positives)} planted "
            f"topology-dependent effects at n = {int(n)} and returned "
            + ("a false positive" if false_positive else "no false positive")
            + " on the unplanted control."
        ),
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
    }


def ladder_power(
    strengths: Sequence[float] = (0.0, 0.25, 0.5, 1.0),
    ns: Sequence[int] = (50, 200),
    replicates: int = 5,
    seed: int = 0,
    gains: dict[str, float] | None = None,
    modes: Sequence[str] = STRUCTURAL_MODES,
    n_nodes: int = 250,
    n_edges: int = 1500,
    loop_size: int = 10,
    **kw: Any,
) -> dict[str, Any]:
    """Detection rate against planted effect size and permutation count.

    ``replicates`` independently generated graphs per cell; the reported rate
    is the share classified topology-dependent.  At ``loop_strength = 0`` the
    rate is the empirical false-positive rate and should sit near ``alpha``.

    Only the structural modes are run by default: they are the ones the class
    turns on, and the experiment is about power, not about the ladder rung.
    """
    from flylab.analysis.nullmodels import _GraphState

    patch = dict(gains or {"g_ach": 0.05})
    t0 = time.perf_counter()
    grid: list[dict[str, Any]] = []
    for strength in strengths:
        for n in ns:
            hits = 0
            done = 0
            for rep in range(int(replicates)):
                graph = synthetic_cut(
                    n_nodes=n_nodes, n_edges=n_edges, loop_size=loop_size,
                    loop_strength=float(strength), seed=int(seed) + rep,
                )
                prof = _ladder_on_state(
                    _GraphState.from_graph(graph), patch, modes=modes,
                    n=int(n), seed=int(seed) + rep,
                    label=f"planted@{float(strength):g}", **kw
                )
                done += 1
                hits += int(prof["class"] == "topology-dependent")
            grid.append(
                {
                    "loop_strength": float(strength),
                    "n": int(n),
                    "replicates": done,
                    "n_detected": hits,
                    "detection_rate": float(hits / done) if done else None,
                    "is_control": bool(float(strength) == 0.0),
                }
            )
    controls = [g for g in grid if g["is_control"]]
    return {
        "grid": grid,
        "strengths": [float(x) for x in strengths],
        "ns": [int(x) for x in ns],
        "replicates": int(replicates),
        "modes": list(modes),
        "false_positive_rate": (
            float(sum(g["n_detected"] for g in controls))
            / float(sum(g["replicates"] for g in controls))
            if controls
            else None
        ),
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
    }


# --------------------------------------------------------------------------
# what the graph under test actually looks like
# --------------------------------------------------------------------------
def cut_census(graph: str = "named", seed_types: Sequence[str] | None = None) -> dict[str, Any]:
    """The structural facts that decide how much a rewiring null can do.

    A degree-preserving rewire of an **in-star** is close to the identity: if
    almost every edge terminates on a handful of hub cells and almost no node
    has any input, there is very little wiring pattern left to destroy, so a
    negative topology result on such a cut is weak evidence that the
    prediction does not need wiring.  The paper has to be able to state these
    numbers, so they are computed here rather than asserted:

    ``edges_onto_seeds`` / ``share_onto_seeds``
        how much of the cut is "everything points at the seed cells".
    ``n_nodes_with_in_degree``
        how many nodes receive anything at all.
    ``share_edges_from_nodes_with_input``
        the share of edges whose source itself receives input - the only edges
        that can carry a path longer than one hop, i.e. the recurrence budget.
    ``top_in_degrees``
        the hubs, by name where the cut names them.
    """
    base = _state_for(graph)
    n = int(base.n)
    pre = np.asarray(base.pre)
    post = np.asarray(base.post)
    n_edges = int(pre.size)
    in_deg = np.bincount(post, minlength=n)
    out_deg = np.bincount(pre, minlength=n)
    has_input = in_deg > 0
    seeds = {k: list(v) for k, v in (base.seeds or {}).items()}
    if seed_types is None:
        seed_types = list(seeds)
    seed_idx = sorted(set(base.seed_indices(list(seed_types))))
    onto_seeds = int(np.isin(post, seed_idx).sum()) if seed_idx else 0
    order = np.argsort(-in_deg)[: max(6, len(seed_idx))]
    by_index: dict[int, str] = {}
    for name, bodies in seeds.items():
        for b in bodies:
            i = base.index.get(int(b))
            if i is not None:
                by_index[i] = name
    return {
        "graph": graph,
        "n_nodes": n,
        "n_edges": n_edges,
        "mean_degree": float(n_edges / n) if n else 0.0,
        "density": float(n_edges / (n * (n - 1))) if n > 1 else 0.0,
        "seed_types": list(seed_types),
        "n_seed_nodes": len(seed_idx),
        "edges_onto_seeds": onto_seeds,
        "share_onto_seeds": float(onto_seeds / n_edges) if n_edges else 0.0,
        "n_nodes_with_in_degree": int(has_input.sum()),
        "share_nodes_with_in_degree": float(has_input.mean()) if n else 0.0,
        "n_nodes_with_out_degree": int((out_deg > 0).sum()),
        "share_edges_from_nodes_with_input": (
            float(has_input[pre].mean()) if n_edges else 0.0
        ),
        "max_in_degree": int(in_deg.max()) if n else 0,
        "top_in_degrees": [
            {
                "index": int(i),
                "seed_type": by_index.get(int(i)),
                "in_degree": int(in_deg[i]),
                "out_degree": int(out_deg[i]),
            }
            for i in order.tolist()
        ],
        "in_star": bool(
            n_edges and (onto_seeds / n_edges) > 0.5 and has_input.mean() < 0.5
        ),
        "note": (
            "A cut where most edges terminate on a few seed cells and most "
            "nodes have no input is an in-star: a degree-preserving rewire of "
            "it is nearly the identity, so 'not distinguishable from the "
            "rewiring null' carries little information there. Repeat any "
            "negative topology result on a denser cut before generalising it."
        ),
        "label": "model_derived",
    }


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
    modes: Sequence[str] = DEFAULT_MODES,
    n: int = 100,
    seed: int = 0,
    alpha: float = DEFAULT_ALPHA,
    effect_floor: float = DEFAULT_EFFECT_FLOOR,
    effect_floor_frac: float = DEFAULT_EFFECT_FLOOR_FRAC,
    delta_frac: float = DEFAULT_DELTA_FRAC,
    delta: float | None = None,
    fdr_alpha: float | None = None,
    fdr_modes: Sequence[str] = STRUCTURAL_MODES,
    gains_by_compound: dict[str, dict[str, float]] | None = None,
    tol: float = DEFAULT_TOL,
    checkpoints: Sequence[float] = DEFAULT_CHECKPOINTS,
    n_jobs: int = 1,
    **kw: Any,
) -> dict[str, Any]:
    """The compound x concentration dependence landscape, FDR-controlled.

    One cell per ``(compound, conc_M)``: its real effect, the per-mode
    permutation p and z, the three-way verdict, the dependence class and the
    necessary information level.  ``compounds=None`` uses the whole library.

    **Multiplicity.** A landscape is a screen, not a planned comparison: the
    default grid is 21 compounds x 4 concentrations x 4 modes, and the count
    of topology-dependent cells is the number the paper quotes.  Benjamini-
    Hochberg FDR is therefore applied across the structural tests of the run
    (``fdr_modes`` x cells, ``fdr_alpha`` defaulting to ``alpha``), every cell
    carries ``p_adjusted`` / ``q_value`` / ``fdr_alpha``, and the class is
    decided on the adjusted probabilities while the raw ones stay beside them.
    The summary reports **both** counts (``n_topology_dependent`` and
    ``n_topology_dependent_raw``) so the effect of the correction on the
    headline number is visible rather than silent.

    Because an empirical permutation p cannot go below ``1/(n+1)``, a coarse
    run can make FDR rejection arithmetically impossible.  The result reports
    ``fdr.min_rejections`` (how many tests must sit at the resolution floor
    before anything can be rejected) and ``fdr.can_reject``; when nothing can
    be rejected the summary says so instead of reporting zero
    topology-dependent cells as if that were a finding.  The landscape is
    ``confirmatory`` only when its resolution is fine enough for one isolated
    test to clear ``fdr_alpha / m``; otherwise it is labelled ``exploratory``.

    ``gains_by_compound`` replaces the gain patch this landscape would
    otherwise compute from the library for the named compounds.  It exists so
    that a caller re-deriving the analysis under an *alternative* mechanism
    specification (:mod:`flylab.analysis.robustness`) can push that
    specification all the way into the permutation engine: the engine resolves
    gains through ``flylab.circuit.rate.compute_gains``, which
    ``robustness.mechanism_spec`` does not rebind, so passing the gains
    explicitly is the only way a specification can actually change a null-model
    result.  The vehicle arm always uses the default gains, which is correct:
    every member of that family reduces to no change at zero engagement.

    The same ``n`` shuffled graphs are reused by every cell (a paired design,
    and the reason this is affordable): shuffle *i* depends only on
    ``(seed, i)``, so each cell's numbers are identical to running
    :func:`dependence_profile` for it on its own.

    Runtime is bounded and reported, and every result carries its own
    ``runtime_s``.  :func:`estimate_landscape_runtime` is an upper bound from
    per-shuffle costs measured on an older, slower null path and now
    overestimates by roughly an order of magnitude.  Measured on the 84-cell
    default grid over the ``named`` cut with ``n_jobs = 4``: about 5.5 s at
    ``n = FAST_N``, and about 5 minutes at ``n = PAPER_N`` with the five
    default modes.  The denser ``taste_motor`` cut costs roughly 15x that per
    shuffle.  Scale ``n`` down for the browser and up for the paper; ``n``
    below 19 cannot resolve ``alpha = 0.05`` at all, and a *corrected*
    landscape needs enough resolution for ``fdr_alpha / m`` as well, which the
    result reports.
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
    if gains_by_compound and not ok_fast:
        raise ValueError(
            "gains_by_compound needs the in-module rate engine, which is not "
            f"available here ({why}); the notebook-assay fallback resolves its "
            "own gains and would silently ignore the override."
        )
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
                        effect_floor_frac=effect_floor_frac,
                        delta_frac=delta_frac,
                        delta=delta,
                        confirmatory=False,
                        tol=tol,
                        checkpoints=checkpoints,
                        n_jobs=n_jobs,
                        **kw,
                    )
                )
        return _landscape_result(
            cells, compounds, concs, assay, graph_name, readout, modes, n, seed,
            alpha, effect_floor, warnings, time.perf_counter() - t0, n_jobs,
            fdr_alpha=fdr_alpha, fdr_modes=fdr_modes, delta_frac=delta_frac,
            delta=delta, effect_floor_frac=effect_floor_frac,
        )

    base = _state_for(graph_name)
    specs: list[dict[str, Any]] = []
    overrides = dict(gains_by_compound or {})
    for compound in compounds:
        for conc in concs:
            if compound in overrides:
                gains = dict(overrides[compound])
            else:
                gains, _occ = compute_gains(
                    compound, conc, library=kw.get("library"), rule_overrides=kw.get("rule_overrides")
                )
            specs.append({"compound": compound, "conc_M": conc, "gains": dict(gains)})
    if overrides:
        warnings.append(
            "gains_by_compound overrides the library-derived gain patch for "
            + ", ".join(sorted(overrides))
            + "; these cells describe the caller's specification, not the "
            "shipped mechanism rules."
        )
    plans = [_plan(base, assay, readout, s["gains"], kw) for s in specs]
    real_cells = [_fast_effect(base, p) for p in plans]
    reals = [r["effect"] for r in real_cells]
    vehicles = [r["vehicle"] for r in real_cells]

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
                effect_floor_frac=effect_floor_frac,
                delta_frac=delta_frac,
                delta=delta,
                real_vehicle=vehicles[idx],
                real_effect=reals[idx],
                rows=rows,
                warnings=[BASE_WARNINGS[1], BASE_WARNINGS[3], CDA_WARNINGS[0], CDA_WARNINGS[1]],
                runtime_s=0.0,
                confirmatory=False,
            )
        )
    return _landscape_result(
        cells, compounds, concs, assay, graph_name, readout, modes, n, seed,
        alpha, effect_floor, warnings, time.perf_counter() - t0, n_jobs,
        fdr_alpha=fdr_alpha, fdr_modes=fdr_modes, delta_frac=delta_frac,
        delta=delta, effect_floor_frac=effect_floor_frac,
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
        st = _shuffle_state_any(
            base, job["mode"], np.random.default_rng([int(job["seed"]), int(i)]), **job["kw"]
        )
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
            st = _shuffle_state_any(base, mode, np.random.default_rng([seed, i]), **kw)
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


def _apply_fdr(
    cells: list[dict[str, Any]],
    alpha: float,
    effect_floor: float,
    fdr_alpha: float,
    fdr_modes: Sequence[str],
    resolution: float,
) -> dict[str, Any]:
    """Benjamini-Hochberg across the structural tests of a whole landscape.

    Writes ``p_adjusted`` / ``q_value`` / ``fdr_alpha`` into every row of the
    family, records each cell's own ``q_value`` (its smallest adjusted
    structural probability), keeps the uncorrected class as ``class_raw`` and
    re-decides ``class`` on the adjusted values.
    """
    family = [m for m in fdr_modes if m in MODES]
    index: list[tuple[int, str]] = []
    pvals: list[float | None] = []
    for i, cell in enumerate(cells):
        cell["class_raw"] = cell["class"]
        cell["classification_raw"] = cell["classification"]
        cell["necessary_information_level_raw"] = cell["necessary_information_level"]
        by_mode = {r["mode"]: r for r in cell["modes"]}
        for mode in family:
            row = by_mode.get(mode)
            if row is None:
                continue
            index.append((i, mode))
            pvals.append(row.get("p_two_sided"))
    bh = benjamini_hochberg(pvals, alpha=fdr_alpha, resolution=resolution)
    for (i, mode), q in zip(index, bh["adjusted"]):
        by_mode = {r["mode"]: r for r in cells[i]["modes"]}
        row = by_mode[mode]
        row["p_adjusted"] = q
        row["q_value"] = q
        row["fdr_alpha"] = float(fdr_alpha)
        row["fdr_family_size"] = bh["m"]
    for cell in cells:
        qs = [
            r.get("q_value")
            for r in cell["modes"]
            if r["mode"] in family and r.get("q_value") is not None
        ]
        cell["q_value"] = (min(qs) if qs else None)
        cell["fdr_alpha"] = float(fdr_alpha)
        cell["fdr_family_size"] = bh["m"]
        cell["fdr_modes"] = list(family)
        cell["multiplicity_correction"] = bh["method"]
        # each cell keeps its own equivalence margin (delta_frac x |vehicle|)
        # and its own effective effect floor (max(absolute, frac x |vehicle|))
        _refinalise(
            cell,
            alpha,
            float(cell.get("effect_floor", effect_floor)),
            cell.get("delta"),
        )
    bh["modes"] = list(family)
    bh.pop("adjusted", None)
    bh.pop("rejected", None)
    return bh


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
    *,
    fdr_alpha: float | None = None,
    fdr_modes: Sequence[str] = STRUCTURAL_MODES,
    delta_frac: float = DEFAULT_DELTA_FRAC,
    delta: float | None = None,
    effect_floor_frac: float = DEFAULT_EFFECT_FLOOR_FRAC,
) -> dict[str, Any]:
    res = p_resolution(int(n))
    q_alpha = float(alpha if fdr_alpha is None else fdr_alpha)
    deltas = [c.get("delta") for c in cells if c.get("delta") is not None]
    bh = _apply_fdr(
        cells,
        alpha=float(alpha),
        effect_floor=float(effect_floor),
        fdr_alpha=q_alpha,
        fdr_modes=fdr_modes,
        resolution=res,
    )

    counts: dict[str, int] = {c: 0 for c in CLASSES}
    counts_raw: dict[str, int] = {c: 0 for c in CLASSES}
    counts_absolute_floor: dict[str, int] = {c: 0 for c in CLASSES}
    levels: dict[str, int] = {}
    verdict_counts: dict[str, int] = {v: 0 for v in VERDICTS}
    non_monotone: list[dict[str, Any]] = []
    below_relative_floor: list[dict[str, Any]] = []
    table: list[dict[str, Any]] = []
    for cell in cells:
        # the same cell classified with the absolute floor only, so the cost of
        # the 1 %-of-baseline convention is reported rather than assumed
        by_mode_c = {r["mode"]: r for r in cell["modes"]}
        abs_cls = classify(
            cell["real_effect"], by_mode_c, alpha=float(alpha),
            effect_floor=cell.get("effect_floor_absolute", float(effect_floor)),
            delta=cell.get("delta"),
        )["class"]
        cell["class_absolute_floor"] = abs_cls
        counts_absolute_floor[abs_cls] = counts_absolute_floor.get(abs_cls, 0) + 1
        if (
            cell["real_effect"] is not None
            and cell.get("effect_floor_relative")
            and abs(float(cell["real_effect"])) <= float(cell["effect_floor_relative"])
            and abs_cls not in ("no-effect", "undefined")
        ):
            below_relative_floor.append(
                {
                    "compound": cell["compound"],
                    "conc_M": cell["conc_M"],
                    "real_effect": cell["real_effect"],
                    "class_absolute_floor": abs_cls,
                }
            )
        lvl_row = cell["necessary_information_level"]
        if lvl_row.get("non_monotone"):
            non_monotone.append(
                {
                    "compound": cell["compound"],
                    "conc_M": cell["conc_M"],
                    "class": cell["class"],
                    "necessary_information_level": lvl_row.get("level"),
                    "richer_models_distinguishable": lvl_row.get(
                        "richer_models_distinguishable"
                    ),
                }
            )
        counts[cell["class"]] = counts.get(cell["class"], 0) + 1
        counts_raw[cell["class_raw"]] = counts_raw.get(cell["class_raw"], 0) + 1
        lvl = cell["necessary_information_level"]["level"] or "undefined"
        levels[lvl] = levels.get(lvl, 0) + 1
        for v in cell["verdicts"].values():
            verdict_counts[v] = verdict_counts.get(v, 0) + 1
        row = {
            "compound": cell["compound"],
            "conc_M": cell["conc_M"],
            "real_effect": cell["real_effect"],
            "class": cell["class"],
            "class_raw": cell["class_raw"],
            "q_value": cell["q_value"],
            "fdr_alpha": cell["fdr_alpha"],
            "delta": cell["delta"],
            "necessary_information_level": lvl,
            "necessary_level_verdict": cell["necessary_information_level"].get("verdict"),
            "stabilised": cell["stabilised"],
        }
        by_mode = {r["mode"]: r for r in cell["modes"]}
        for mode in modes:
            key = _mode_key(mode)
            row[key.replace("z_", "p_")] = cell["p"][key.replace("z_", "p_")]
            row[key.replace("z_", "q_")] = (by_mode.get(mode) or {}).get("q_value")
            row[key] = cell["D"][key]
            row[f"verdict_{mode}"] = cell["verdicts"].get(mode)
        table.append(row)

    if res > float(alpha):
        warnings.append(
            f"n = {int(n)} gives a permutation resolution of {res:.4f}, coarser "
            f"than alpha = {float(alpha)}: no cell in this landscape can be "
            "classified as topology-dependent at this n. It is a preview, not a "
            "result; the paper's landscape is run at n = 500-1000."
        )
    if not bh["can_reject"]:
        warnings.append(
            f"FDR CANNOT REJECT ANYTHING AT THIS RESOLUTION: with {bh['m']} "
            f"structural tests, a permutation resolution of {res:.5f} and "
            f"fdr_alpha = {q_alpha}, at least {bh['min_rejections']} tests would "
            "have to sit at the resolution floor simultaneously before "
            "Benjamini-Hochberg could reject any of them. A count of zero "
            "topology-dependent cells here is a property of the permutation "
            "budget, NOT a finding. Raise n to at least "
            f"{int(round(bh['m'] / q_alpha)) - 1} for a single isolated test to "
            "be rejectable."
        )
    elif bh["min_rejections"] and bh["min_rejections"] > 1:
        warnings.append(
            f"the FDR-adjusted threshold is resolution-limited: with {bh['m']} "
            f"structural tests at resolution {res:.5f}, at least "
            f"{bh['min_rejections']} tests must sit at the floor together before "
            "any is rejected, so an isolated strong cell cannot survive the "
            "correction at this n."
        )
    confirmatory = bool(bh.get("resolution_supports_single_rejection"))
    if not confirmatory:
        warnings.append(
            "EXPLORATORY: this landscape's permutation resolution "
            f"({res:.5f}) does not reach the adjusted threshold fdr_alpha/m = "
            f"{q_alpha / max(bh['m'], 1):.6f}, so its counts are a screen to be "
            "confirmed, not a confirmatory result. The prespecified "
            "confirmatory tests are the "
            + "/".join(CONFIRMATORY_COMPOUNDS)
            + f" profiles at n >= {PAPER_N}."
        )
    unstable = [r for r in table if not r["stabilised"]]
    if unstable:
        warnings.append(
            f"{len(unstable)} of {len(table)} cells had at least one mode whose "
            f"permutation p had not stabilised at n = {int(n)}."
        )
    if non_monotone:
        warnings.append(
            f"{len(non_monotone)} of {len(cells)} cells have a NON-MONOTONE "
            "ladder: a weaker graph model could not be distinguished from the "
            "real cut while a richer one could, so the reported necessary "
            "information level sits below a rung that was rejected. These "
            "cells' levels must be quoted with that caveat: "
            + ", ".join(
                f"{r['compound']}@{r['conc_M']:g}" for r in non_monotone[:12]
            )
            + ("..." if len(non_monotone) > 12 else "")
        )
    if below_relative_floor:
        warnings.append(
            f"{len(below_relative_floor)} cells move the readout by less than "
            "the relative effect floor (1 % of the vehicle baseline) and are "
            "counted as no-effect here, but would be classified on the "
            "absolute floor alone: "
            + ", ".join(
                f"{r['compound']}@{r['conc_M']:g} ({r['real_effect']:.3g})"
                for r in below_relative_floor[:12]
            )
            + ("..." if len(below_relative_floor) > 12 else "")
        )
    n_topo = counts.get("topology-dependent", 0)
    n_topo_raw = counts_raw.get("topology-dependent", 0)
    if n_topo != n_topo_raw:
        warnings.append(
            f"multiplicity changes the headline count: {n_topo_raw} cells are "
            f"topology-dependent on the raw permutation p and {n_topo} survive "
            f"Benjamini-Hochberg at fdr_alpha = {q_alpha} over {bh['m']} "
            "structural tests. The corrected count is the one to quote."
        )
    topo = sorted({r["compound"] for r in table if r["class"] == "topology-dependent"})
    topo_raw = sorted({r["compound"] for r in table if r["class_raw"] == "topology-dependent"})
    comp = sorted({r["compound"] for r in table if r["class"] == "composition-dominated"})
    statement = (
        "exploratory screen, FDR-controlled: "
        f"{n_topo} of {len(cells)} cells are topology-dependent after "
        f"Benjamini-Hochberg at fdr_alpha = {q_alpha} across {bh['m']} structural "
        f"tests ({n_topo_raw} before correction)."
        if not confirmatory
        else (
            f"confirmatory-resolution landscape: {n_topo} of {len(cells)} cells "
            f"are topology-dependent after Benjamini-Hochberg at fdr_alpha = "
            f"{q_alpha} across {bh['m']} structural tests ({n_topo_raw} before "
            "correction); the permutation resolution supports the adjusted "
            "threshold for a single isolated test."
        )
    )
    if not bh["can_reject"]:
        statement += (
            " No cell CAN be rejected at this permutation resolution, so the "
            "corrected count is uninformative rather than zero."
        )
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
        "fdr_alpha": q_alpha,
        "fdr_modes": list(bh["modes"]),
        "effect_floor": float(effect_floor),
        "effect_floor_frac": float(effect_floor_frac),
        "delta_frac": (None if delta is not None else float(delta_frac)),
        "delta_absolute": (None if delta is None else abs(float(delta))),
        "delta_range": ([min(deltas), max(deltas)] if deltas else None),
        "p_resolution": res,
        "resolution_coarser_than_alpha": bool(res > float(alpha)),
        "confirmatory": confirmatory,
        "design": "confirmatory" if confirmatory else "exploratory",
        "design_statement": statement,
        "multiplicity_correction": bh["method"],
        "fdr": bh,
        "shape": [len(compounds), len(concs)],
        "n_cells": len(cells),
        "cells": cells,
        "table": table,
        "summary": {
            "design": "confirmatory" if confirmatory else "exploratory",
            "statement": statement,
            "class_counts": counts,
            "class_counts_raw": counts_raw,
            "class_counts_absolute_floor": counts_absolute_floor,
            "n_topology_dependent_absolute_floor": counts_absolute_floor.get(
                "topology-dependent", 0
            ),
            "n_non_monotone_ladders": len(non_monotone),
            "non_monotone_cells": non_monotone,
            "n_below_relative_effect_floor": len(below_relative_floor),
            "below_relative_effect_floor": below_relative_floor,
            "n_topology_dependent": n_topo,
            "n_topology_dependent_raw": n_topo_raw,
            "n_structural_tests": bh["m"],
            "n_structural_rejected_fdr": bh["n_rejected"],
            "fdr_alpha": q_alpha,
            "fdr_can_reject": bh["can_reject"],
            "fdr_min_rejections": bh["min_rejections"],
            "verdict_counts": verdict_counts,
            "level_counts": levels,
            "topology_dependent_compounds": topo,
            "topology_dependent_compounds_raw": topo_raw,
            "composition_dominated_compounds": comp,
            "n_unstabilised_cells": len(unstable),
        },
        "runtime_s": float(runtime_s),
        "label": "model_derived",
        "warnings": warnings,
    }
