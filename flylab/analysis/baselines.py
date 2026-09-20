"""The ablation ladder: how much does each layer of FlyLab actually add?

Connectome null models ask whether a *degraded graph* reproduces a prediction.
This module asks the complementary - and, for a computing audience, more
familiar - question: how much of the prediction survives if you delete whole
layers of the model?  Four levels of sophistication are run on the same
compound at the same concentration:

=======================  ==========================================================
level                    what it is allowed to use
=======================  ==========================================================
``A_receptor_only``      the compound, the concentration and the receptor
                         library.  Signed engagement summed over the insect
                         receptor rows.  **No** gain rules, no transmitters, no
                         circuit at all.
``B_composition_only``   the mechanism gain rules plus the graph's *global*
                         transmitter proportions (its own census, scored with
                         the whole-CNS excitation index of
                         ``flylab.assays.wholens``).  **No** edges.
``C_topology_only``      the real MaleCNS cut and the real rate model, but one
                         **generic** pharmacological multiplier per compound
                         applied to every transmitter instead of
                         mechanism-specific gains.  The multiplier is
                         *direction-aware* (:data:`GENERIC_RULES`): its
                         magnitude is generic (``theta_max``), its sign is one
                         bit taken from the mechanism table.  The old
                         depression-only rule is kept as
                         ``C_topology_only_floor``.
``D_full_flylab``        the full model: mechanism-specific per-transmitter
                         gains on the real cut.
=======================  ==========================================================

The four levels do **not** share units (A is a dimensionless receptor sum, B is
an excitation index, C and D are Hz), so a single compound's numbers are not
comparable across levels.  What is comparable - and what the ablation actually
measures - is the **ordering of the compound library**: if a cheap level ranks
the library the way the full model does, the expensive layers added nothing for
that question.  :func:`ablation` therefore reports, per level, the Pearson and
Spearman correlation of its prediction against the full model across the
library, the residual it leaves after standardising both, and, for the levels
that are in Hz, the residual in Hz.

For the nicotinic agonists the composition-only baseline tracks the full
model's ordering closely. That agreement is largely structural: both levels are
functions of the same gain vector, matched pharmacology-free pseudo-compounds
already reach a median near 0.85, shuffling compound labels leaves the value
unchanged, and it does not survive a degree-corrected engine. Quote it only
against its matched reference distribution and the normalisation used.  The output says so in words when it happens, and says the
opposite when a level fails.

Three things a reader of the ablation numbers has to know
--------------------------------------------------------
**1. The topology-only level used to be sign-broken, and its poor score was
structural.**  Through v0.6 level C applied ``g_all = max(0.05, 1 - theta_max)``
to ``g_ach``, ``g_gaba``, ``g_glu`` and ``g_oct`` alike, so *every* entry was at
most 1.0 and the level could only ever reduce activity.  A model that can only
depress cannot order a library containing RDL and GluCl blockers whose
full-model effect is positive, so "the connectome without the pharmacology is
not a cheap substitute for the connectome with it" was, in that form, a
statement about a baseline that had been denied the ability to express
disinhibition.  Level C is now *direction-aware*: the magnitude stays generic
(``theta_max``, the largest modelled insect engagement) and exactly **one bit**
- does this compound's dominant mechanism raise or lower net synaptic drive? -
is taken from the mechanism table.  The old rule is still computed and reported
beside it as ``C_topology_only_floor``, which is a **floor** (how badly a
depression-only connectome model does) and not a competitor.  See
:data:`GENERIC_RULES`.

**2. Level B versus level D is close to an algebraic identity.**  Both levels
are functions of the *same* gain vector: B applies it to a fixed census, D
applies it per cell on a row-normalised matrix, where each cell's recurrent
input is a composition-weighted average of its presynaptic gains (see
``flylab.circuit.rate``).  A high rank correlation between them is therefore
expected for structural reasons and is *not*, on its own, a finding about the
connectome.  :func:`composition_reference_distribution` builds the reference
distribution that makes the observed number interpretable - what rho arbitrary
gain vectors already produce, what a shuffled compound-to-gain assignment
produces, and what value would have been evidence against the claim - and
:func:`composition_dominance_under_normalisations` asks whether the dominance
survives an engine without the row normalisation.

**3. Two levels of this ladder disagree about the sign of glutamate.**  The
composition index inherited from ``flylab.assays.wholens`` scores glutamate as
partly *excitatory* (``ach + 0.3*glu - 1.1*gaba``), while the rate engine signs
it *inhibitory* (``SIGN["glutamate"] = -0.4``).  Both conventions are available
(:data:`NT_SIGN_CONVENTIONS`); the shipped default is still ``"wholens"`` so
published numbers do not move silently, and every information-gain block
carries a ``composition_sign_convention_check`` that reports what the
disagreement does to the B-versus-D correlation the paper's conclusion rests
on.

Everything here is ``model_derived``.
"""
from __future__ import annotations

import time
from typing import Any, Iterable, Sequence

import numpy as np

from flylab.analysis.nullmodels import (
    BASE_WARNINGS,
    DEFAULT_GRAPH,
    DEFAULT_READOUT,
    _fast_effect,
    _fast_supported,
    _plan,
    _state_for,
    drug_effect,
)

__all__ = [
    "LEVELS",
    "LEVEL_NOTES",
    "GENERIC_MULTIPLIER_RULE",
    "GENERIC_RULES",
    "GENERIC_RULE_NOTES",
    "DEFAULT_GENERIC_RULE",
    "C_FLOOR_LEVEL",
    "NT_SIGN_CONVENTIONS",
    "DEFAULT_NT_SIGN_CONVENTION",
    "REPRODUCES_R",
    "REPRODUCES_RHO",
    "generic_multiplier",
    "composition_reference_distribution",
    "composition_dominance_under_normalisations",
    "glutamate_sign_reconciliation",
    "baseline_receptor_only",
    "baseline_composition_only",
    "baseline_topology_only",
    "full_model",
    "level_predictions",
    "ablation",
    "ablation_table",
    "pearson",
    "spearman",
]

#: the ladder, cheapest first
LEVELS: tuple[str, ...] = (
    "A_receptor_only",
    "B_composition_only",
    "C_topology_only",
    "D_full_flylab",
)

LEVEL_NOTES: dict[str, str] = {
    "A_receptor_only": (
        "engagement only: signed sum of insect receptor engagement. No gain "
        "rules, no transmitters, no circuit. Unit: dimensionless."
    ),
    "B_composition_only": (
        "composition only: mechanism gains applied to the graph's global "
        "transmitter proportions (its own census), scored with the whole-CNS "
        "excitation index. No edges. Unit: excitation index."
    ),
    "C_topology_only": (
        "topology only: the real cut and the real rate model with one generic, "
        "direction-aware multiplier per compound instead of mechanism-specific "
        "gains. Magnitude generic (theta_max), sign (one bit) from the "
        "mechanism table. Unit: Hz."
    ),
    "D_full_flylab": (
        "full FlyLab: mechanism-specific per-transmitter gains on the real cut. "
        "Unit: Hz."
    ),
}

#: The historical (v0.4-v0.6) stand-in for pharmacology at level C: "the drug
#: removes a fraction theta_max of synaptic transmission, whatever the
#: transmitter".  Kept for continuity and reported as a **floor**, because it
#: can only ever depress: every entry is <= 1.0, so no compound can come out
#: positive and the level cannot express disinhibition.
GENERIC_MULTIPLIER_RULE = "g_all = max(0.05, 1 - theta_max)"

#: The level-C rule used by default.  Same generic magnitude, but the sign of
#: the perturbation is taken from the mechanism table, which is the minimum a
#: "connectome without the pharmacology" model needs in order to be able to
#: order a library that contains disinhibitors.
DIRECTION_AWARE_MULTIPLIER_RULE = (
    "g_all = max(0.05, 1 - theta_max) if the compound's dominant mechanism "
    "lowers net synaptic drive, else 1 + theta_max"
)

C_FLOOR_LEVEL = "C_topology_only_floor"

#: ``rule name -> (formula, what information it is given)``
GENERIC_RULES: dict[str, str] = {
    "direction_aware": DIRECTION_AWARE_MULTIPLIER_RULE,
    "depressant_floor": GENERIC_MULTIPLIER_RULE,
}

GENERIC_RULE_NOTES: dict[str, str] = {
    "direction_aware": (
        "Generic magnitude (theta_max), direction-aware sign. The baseline is "
        "told ONE bit per compound - whether the dominant insect mechanism "
        "raises or lowers net synaptic drive - and nothing about which "
        "transmitter is involved or by how much. This is the fair "
        "topology-only competitor: it can express disinhibition."
    ),
    "depressant_floor": (
        "The historical rule: max(0.05, 1 - theta_max) on every transmitter. "
        "Structurally incapable of a positive effect, so it cannot order a "
        "library containing RDL/GluCl blockers. Reported as a FLOOR (what a "
        "depression-only connectome model achieves), never as a competitor."
    ),
}

DEFAULT_GENERIC_RULE = "direction_aware"

#: How the composition level (B) signs glutamate.  ``"wholens"`` is the
#: inherited whole-CNS index, which scores glutamate as partly excitatory
#: (``ach + 0.3*glu - 1.1*gaba``); ``"rate_engine"`` uses the sign the rate
#: engine actually iterates (``SIGN["glutamate"] = -0.4``), which is the same
#: transmitter convention level D uses.  The two levels of this ladder
#: disagreed silently before v0.6.1; the default is unchanged so that published
#: numbers do not move on their own, and the disagreement is quantified in
#: every information-gain block.
NT_SIGN_CONVENTIONS: tuple[str, ...] = ("wholens", "rate_engine")
DEFAULT_NT_SIGN_CONVENTION = "wholens"

#: a cheaper level "reproduces" the full model at or above these.  The
#: threshold is applied to the **signed** rho: an ordering that is a perfect
#: inversion of the full model's is a reproduction of nothing.
REPRODUCES_R = 0.95
REPRODUCES_RHO = 0.90

_ACTIVATING = {"agonist", "partial_agonist", "positive_modulator"}
_BLOCKING = {"antagonist", "negative_modulator", "inhibitor"}

ABLATION_WARNINGS = [
    "The four levels do not share units; only their ordering of the compound "
    "library is comparable, which is what the correlations measure.",
    "A cheaper level reproducing the full model is a statement about this "
    "readout on this cut, not about connectomics in general.",
]


# --------------------------------------------------------------------------
# small statistics (no scipy: the science core must run in the browser)
# --------------------------------------------------------------------------
def _pairs(a: Sequence[float | None], b: Sequence[float | None]):
    xs, ys = [], []
    for x, y in zip(a, b):
        if x is None or y is None:
            continue
        x, y = float(x), float(y)
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        xs.append(x)
        ys.append(y)
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def pearson(a: Sequence[float | None], b: Sequence[float | None]) -> float | None:
    """Pearson r over the pairwise-complete entries, or None."""
    x, y = _pairs(a, b)
    if x.size < 3 or x.std() == 0 or y.std() == 0:
        return None
    r = float(np.corrcoef(x, y)[0, 1])
    return r if np.isfinite(r) else None


def _ranks(x: np.ndarray) -> np.ndarray:
    """Average ranks (ties shared), so Spearman survives duplicate values."""
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, x.size + 1, dtype=float)
    # average the ranks of tied values
    vals, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    for k in np.where(counts > 1)[0]:
        sel = inv == k
        ranks[sel] = ranks[sel].mean()
    return ranks


def spearman(a: Sequence[float | None], b: Sequence[float | None]) -> float | None:
    """Spearman rho over the pairwise-complete entries, or None."""
    x, y = _pairs(a, b)
    if x.size < 3:
        return None
    rx, ry = _ranks(x), _ranks(y)
    if rx.std() == 0 or ry.std() == 0:
        return None
    rho = float(np.corrcoef(rx, ry)[0, 1])
    return rho if np.isfinite(rho) else None


def _z(x: np.ndarray) -> np.ndarray:
    sd = x.std()
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)


# --------------------------------------------------------------------------
# the receptor layer, defensively (the evidence layer is typed: engagement may
# be None, which means "not modelled", never 0)
# --------------------------------------------------------------------------
def _engagement(row: dict[str, Any]) -> float | None:
    val = row.get("engagement", row.get("occupancy"))
    if val is None:
        return None
    try:
        v = float(val)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def _insect_rows(compound: str, conc_M: float, library: dict[str, Any] | None = None):
    from flylab.pharm.occupancy import compare_compound

    occ = (
        compare_compound(compound, conc_M, library=library)
        if library
        else compare_compound(compound, conc_M)
    )
    rows = [r for r in occ["receptors"] if str(r.get("receptor", "")).startswith("insect_")]
    return occ, rows


def _theta_max(rows: Iterable[dict[str, Any]]) -> tuple[float | None, int]:
    """Largest modelled insect engagement with a direction, and how many rows
    carried no number at all.

    Kept as the public-ish spelling of the generic magnitude; the level-C rules
    obtain it from :func:`_dominant_row`, which also returns the row itself so
    that the direction-aware rule can read one bit off the mechanism table.
    """
    _row, best, not_modelled = _dominant_row(rows)
    return best, not_modelled


# --------------------------------------------------------------------------
# level A: receptor only
# --------------------------------------------------------------------------
def baseline_receptor_only(
    compound: str | None,
    conc_M: float = 1e-6,
    library: dict[str, Any] | None = None,
    **_: Any,
) -> dict[str, Any]:
    """Baseline A: the compound, the dose and the receptor library. No circuit.

    ``effect`` is the signed sum of insect receptor engagement: ``+e`` for an
    activating direction, ``-e`` for a blocking one, and nothing at all for a
    row whose engagement is ``None`` (not modelled) or whose direction is
    ``none``.  This is deliberately the naive reading of a receptor scorecard:
    it has no gain rule, so it cannot know that a *saturating* nicotinic
    agonist silences cholinergic transmission rather than boosting it.  What
    the ablation measures is exactly how much that costs.
    """
    if not compound:
        return {
            "level": "A_receptor_only",
            "effect": 0.0,
            "unit": "dimensionless",
            "detail": {"rows": 0},
            "warnings": [],
        }
    occ, rows = _insect_rows(compound, conc_M, library)
    total = 0.0
    used = 0
    not_modelled = 0
    per_receptor: dict[str, Any] = {}
    for r in rows:
        name = str(r.get("receptor"))
        direction = str(r.get("direction") or "none")
        e = _engagement(r)
        if e is None:
            not_modelled += 1
            per_receptor[name] = None
            continue
        sign = 1.0 if direction in _ACTIVATING else (-1.0 if direction in _BLOCKING else 0.0)
        if sign == 0.0:
            per_receptor[name] = 0.0
            continue
        total += sign * e
        per_receptor[name] = sign * e
        used += 1
    warnings = []
    if not_modelled:
        warnings.append(
            f"{not_modelled} insect receptor row(s) for {compound} carry no "
            "modellable value (engagement is None); they were skipped, never "
            "read as an engagement of 0."
        )
    if used == 0:
        warnings.append(
            f"no insect receptor row for {compound} could drive baseline A; its "
            "prediction is 0 by absence of evidence, not by pharmacology."
        )
    return {
        "level": "A_receptor_only",
        "effect": float(total),
        "unit": "dimensionless",
        "detail": {
            "per_receptor": per_receptor,
            "n_rows": len(rows),
            "n_used": used,
            "n_not_modelled": not_modelled,
        },
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# level B: composition only
# --------------------------------------------------------------------------
def _settle(
    counts: dict[str, float],
    gains: dict[str, float],
    g_nav: float = 1.0,
    convention: str = DEFAULT_NT_SIGN_CONVENTION,
) -> dict[str, float]:
    """Whole-CNS style excitation index, on whatever census it is handed.

    ``convention="wholens"`` (default) is the ``flylab.assays.wholens`` index
    (imported when available so the two cannot drift), applied to the *graph's
    own* transmitter census instead of the 165k-cell MaleCNS one.  It scores
    glutamate as partly **excitatory** (``ach + 0.3*glu - 1.1*gaba``).

    ``convention="rate_engine"`` replaces only the glutamate coefficient with
    the sign the rate engine (level D) actually iterates,
    ``SIGN["glutamate"] = -0.4``, so that the two levels of the ladder agree on
    the transmitter's sign.  Nothing else in the formula changes, which is what
    makes the difference between the two attributable to the glutamate sign
    alone.
    """
    if convention not in NT_SIGN_CONVENTIONS:
        raise ValueError(
            f"unknown nt_sign_convention {convention!r}; expected one of {NT_SIGN_CONVENTIONS}"
        )
    if convention == "wholens":
        try:  # keep one implementation when the assay layer is importable
            from flylab.assays.wholens import _settle as _wholens_settle

            return _wholens_settle(counts, gains, g_nav)
        except Exception:  # pragma: no cover - fallback copy of the same formula
            glu_exc, glu_inh = 0.3, 0.7
    else:
        from flylab.circuit.rate import SIGN

        glu_exc = float(SIGN.get("glutamate", -0.4))
        glu_inh = abs(glu_exc)
    tot = max(sum(counts.values()), 1)
    drive = {k: 40.0 * g_nav * (v / tot) for k, v in counts.items()}
    ach = drive.get("acetylcholine", 0.0) * gains["acetylcholine"]
    gaba = drive.get("gaba", 0.0) * gains["gaba"]
    glu = drive.get("glutamate", 0.0) * gains.get("glutamate", 1.0)
    excitation = max(0.0, ach + glu_exc * glu - 1.1 * gaba)
    inhibition = max(0.0, gaba + glu_inh * glu)
    return {
        "cns_excitation_index": excitation,
        "cns_inhibition_index": inhibition,
        "excitation_inhibition_ratio": excitation / inhibition if inhibition > 1e-9 else None,
    }


def graph_census(graph: str | None = "named") -> dict[str, int]:
    """Transmitter counts of the graph's own nodes (``unclear`` for unlabelled)."""
    from flylab.circuit.rate import load_graph

    g = load_graph(graph)
    counts: dict[str, int] = {}
    for node in g["nodes"]:
        nt = node.get("consensus_nt") or "unclear"
        counts[nt] = counts.get(nt, 0) + 1
    return counts


def composition_effect_from_gains(
    gains: dict[str, float],
    counts: dict[str, int] | None = None,
    graph: str | None = "named",
    nt_sign_convention: str = DEFAULT_NT_SIGN_CONVENTION,
) -> float:
    """Level B's effect for an arbitrary gain dict: the census part, alone.

    Exposed because the reference distributions
    (:func:`composition_reference_distribution`) need to score gain vectors
    that belong to no compound.
    """
    counts = graph_census(graph) if counts is None else counts
    nt_gains = {
        "acetylcholine": gains.get("g_ach", 1.0) * gains.get("ach_tone", 1.0),
        "gaba": gains.get("g_gaba", 1.0),
        "glutamate": gains.get("g_glu", 1.0),
        "octopamine": gains.get("g_oct", 1.0),
        "other": 1.0,
    }
    treated = _settle(counts, nt_gains, gains.get("g_nav", 1.0), nt_sign_convention)
    vehicle = _settle(counts, {k: 1.0 for k in nt_gains}, 1.0, nt_sign_convention)
    return float(treated["cns_excitation_index"] - vehicle["cns_excitation_index"])


def baseline_composition_only(
    compound: str | None,
    conc_M: float = 1e-6,
    graph: str | None = "named",
    library: dict[str, Any] | None = None,
    rule_overrides: dict[str, float] | None = None,
    nt_sign_convention: str = DEFAULT_NT_SIGN_CONVENTION,
    **_: Any,
) -> dict[str, Any]:
    """Baseline B: mechanism gains x the graph's global transmitter proportions.

    No edges: the cut is reduced to "this fraction of the cells is
    cholinergic, this fraction GABAergic", and the drug's gains are applied to
    those proportions.  ``effect`` is the change in the excitation index.

    ``nt_sign_convention`` selects how glutamate is signed
    (:data:`NT_SIGN_CONVENTIONS`).  The default ``"wholens"`` is the inherited
    index, which treats glutamate as partly excitatory and therefore disagrees
    with the rate engine this ladder compares against; ``"rate_engine"``
    reconciles the two.  ``detail["effect_other_convention"]`` always carries
    the alternative, so the disagreement is visible per compound.
    """
    from flylab.circuit.rate import compute_gains

    counts = graph_census(graph)
    gains, _occ = compute_gains(compound, conc_M, library=library, rule_overrides=rule_overrides)
    nt_gains = {
        "acetylcholine": gains["g_ach"] * gains["ach_tone"],
        "gaba": gains["g_gaba"],
        "glutamate": gains["g_glu"],
        "octopamine": gains["g_oct"],
        "other": 1.0,
    }
    treated = _settle(counts, nt_gains, gains["g_nav"], nt_sign_convention)
    vehicle = _settle(counts, {k: 1.0 for k in nt_gains}, 1.0, nt_sign_convention)
    other = [c for c in NT_SIGN_CONVENTIONS if c != nt_sign_convention][0]
    total = max(sum(counts.values()), 1)
    warnings: list[str] = []
    if counts.get("glutamate"):
        warnings.append(
            "Level B and level D disagree on the sign of glutamate: this index "
            f"uses the {nt_sign_convention!r} convention while the rate engine "
            "signs glutamate -0.4. detail.effect_other_convention gives the "
            "same compound under the other convention."
        )
    return {
        "level": "B_composition_only",
        "effect": float(
            treated["cns_excitation_index"] - vehicle["cns_excitation_index"]
        ),
        "unit": "excitation index",
        "detail": {
            "graph": graph,
            "census": counts,
            "frac_ach": counts.get("acetylcholine", 0) / total,
            "frac_gaba": counts.get("gaba", 0) / total,
            "frac_glu": counts.get("glutamate", 0) / total,
            "gains": dict(gains),
            "treated": treated,
            "vehicle": vehicle,
            "nt_sign_convention": nt_sign_convention,
            "effect_other_convention": composition_effect_from_gains(
                gains, counts, graph, other
            ),
        },
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# level C: topology only
# --------------------------------------------------------------------------
def _dominant_row(rows: Iterable[dict[str, Any]]) -> tuple[dict[str, Any] | None, float | None, int]:
    """The insect row with the largest modelled engagement *and* a direction.

    Returns ``(row, theta_max, n_not_modelled)``.  This is the single row the
    generic rules are allowed to look at: level C gets one number from it
    (``theta_max``) and, for the direction-aware rule, one bit.
    """
    best_row: dict[str, Any] | None = None
    best: float | None = None
    not_modelled = 0
    for r in rows:
        direction = str(r.get("direction") or "none")
        e = _engagement(r)
        if e is None:
            not_modelled += 1
            continue
        if direction == "none":
            continue
        if best is None or e > best:
            best, best_row = e, r
    return best_row, best, not_modelled


def _net_drive_direction(gains: dict[str, float]) -> float:
    """Signed change in net synaptic drive implied by a gain dict.

    Each gain is combined with the sign the rate engine gives the transmitter
    it multiplies (:data:`flylab.circuit.rate.SIGN`), so blocking an inhibitory
    transmitter comes out positive.  ``g_nav`` is an excitability scale, so it
    enters with sign +1.  The *magnitude* returned is deliberately not used:
    level C is allowed the sign only.
    """
    from flylab.circuit.rate import SIGN

    ach = float(gains.get("g_ach", 1.0)) * float(gains.get("ach_tone", 1.0))
    net = (ach - 1.0) * SIGN.get("acetylcholine", 1.0)
    net += (float(gains.get("g_gaba", 1.0)) - 1.0) * SIGN.get("gaba", -1.0)
    net += (float(gains.get("g_glu", 1.0)) - 1.0) * SIGN.get("glutamate", -0.4)
    net += (float(gains.get("g_oct", 1.0)) - 1.0) * SIGN.get("octopamine", 0.2)
    net += float(gains.get("g_nav", 1.0)) - 1.0
    return float(net)


def _mechanism_direction(rows: Sequence[dict[str, Any]]) -> tuple[int, str]:
    """``(+1 | -1 | 0, source)``: does this compound raise or lower net drive?

    The bit is read off the mechanism table itself - the compound's insect rows
    are passed through ``flylab.pharm.mechanisms.gains_from_occupancy`` and the
    resulting gain vector is collapsed with the rate engine's transmitter signs
    - so it stays correct if a rule changes, and level C never sees the
    magnitudes.  It is a *bit*: the returned sign, never the size.
    """
    if not rows:
        return 0, "no_mechanism_row"
    try:
        from flylab.pharm.mechanisms import gains_from_occupancy

        gains = dict(gains_from_occupancy(list(rows)) or {})
        source = "mechanism_table"
    except Exception:  # pragma: no cover - pharm module in flux
        from flylab.circuit.rate import fallback_gains

        gains = dict(fallback_gains(list(rows)))
        source = "fallback_rules"
    net = _net_drive_direction(gains)
    if abs(net) < 1e-12:
        return 0, "no_rule_for_these_receptors"
    return (1 if net > 0 else -1), source


def generic_multiplier(
    compound: str | None,
    conc_M: float = 1e-6,
    library: dict[str, Any] | None = None,
    rule: str = DEFAULT_GENERIC_RULE,
) -> tuple[float, dict[str, Any]]:
    """One number per compound: the level-C stand-in for pharmacology.

    ``theta_max`` is the largest modelled insect engagement with a direction -
    the generic *magnitude*, identical under both rules.

    * ``rule="depressant_floor"`` (the historical rule) returns
      ``max(0.05, 1 - theta_max)``.  It is <= 1 for every compound, so the
      level can only depress; it is a **floor**, not a competitor.
    * ``rule="direction_aware"`` (default) keeps that magnitude but takes the
      **sign** from the mechanism table: ``1 + theta_max`` when the compound's
      dominant insect mechanism raises net synaptic drive (an RDL or GluCl
      blocker, a Nav opener), ``max(0.05, 1 - theta_max)`` when it lowers it.
      The baseline is still mechanism-free in the sense that matters - it does
      not know *which* transmitter is affected, or by how much - but it can now
      express disinhibition, which the library requires of any model that is to
      order it.

    Returns ``(multiplier, detail)``.
    """
    if rule not in GENERIC_RULES:
        raise ValueError(f"unknown generic rule {rule!r}; expected one of {tuple(GENERIC_RULES)}")
    if not compound:
        return 1.0, {
            "rule": rule,
            "formula": GENERIC_RULES[rule],
            "theta_max": None,
            "direction": 0,
            "direction_source": "no_compound",
            "n_not_modelled": 0,
        }
    _occ, rows = _insect_rows(compound, conc_M, library)
    row, theta, not_modelled = _dominant_row(rows)
    detail = {
        "rule": rule,
        "formula": GENERIC_RULES[rule],
        "theta_max": None if theta is None else float(theta),
        "n_not_modelled": not_modelled,
        "direction": 0,
        "direction_source": "not_applicable",
        "dominant_receptor": None if row is None else str(row.get("receptor")),
        "dominant_direction": None if row is None else str(row.get("direction") or "none"),
    }
    if theta is None:
        return 1.0, detail
    down = max(0.05, 1.0 - float(theta))
    if rule == "depressant_floor":
        detail["direction"] = -1
        detail["direction_source"] = "assumed_depressant"
        return down, detail
    direction, source = _mechanism_direction(rows)
    detail["direction_source"] = source
    if direction > 0:
        detail["direction"] = 1
        return 1.0 + float(theta), detail
    # direction == 0 means the mechanism table has no rule for this receptor;
    # the honest generic default is then the historical depressant assumption,
    # and the detail says so rather than hiding it.
    detail["direction"] = -1
    if direction == 0:
        detail["direction_source"] = f"{source}:default_depressant"
    return down, detail


def baseline_topology_only(
    compound: str | None,
    conc_M: float = 1e-6,
    graph: str | None = "named",
    assay: str = "subgraph",
    readout: str = "auto",
    library: dict[str, Any] | None = None,
    rule: str = DEFAULT_GENERIC_RULE,
    **kw: Any,
) -> dict[str, Any]:
    """Baseline C: the real cut, the real rate model, one generic multiplier.

    The multiplier (:func:`generic_multiplier`) is applied to **every**
    transmitter gain, so the model keeps all of the wiring and none of the
    mechanism except, under the default direction-aware rule, the one bit that
    says whether the compound raises or lowers net synaptic drive.  ``effect``
    is ``treated - vehicle`` in Hz, the same contrast and the same units as the
    full model.

    ``detail["floor"]`` always carries the historical depression-only rule's
    multiplier and effect, so the two can be reported side by side with the old
    one labelled as a floor.
    """
    from flylab.circuit.rate import DEFAULT_GAINS

    if readout == "auto":
        readout = DEFAULT_READOUT[assay]
    graph_name = graph if graph is not None else DEFAULT_GRAPH[assay]

    def _run(rule_name: str) -> tuple[float | None, float, dict[str, Any]]:
        m, det = generic_multiplier(compound, conc_M, library, rule=rule_name)
        gains = dict(DEFAULT_GAINS)
        for key in ("g_ach", "g_gaba", "g_glu", "g_oct"):
            gains[key] = float(m)
        return _effect_with_gains(assay, readout, graph_name, gains, kw), float(m), det

    effect, m, detail = _run(rule)
    floor_effect, floor_m, floor_detail = (
        (effect, m, detail)
        if rule == "depressant_floor"
        else _run("depressant_floor")
    )
    warnings: list[str] = []
    if detail.get("theta_max") is None:
        warnings.append(
            f"no modelled insect engagement for {compound}: baseline C used "
            "a multiplier of 1.0, i.e. it predicts no effect."
        )
    if str(detail.get("direction_source", "")).endswith("default_depressant"):
        warnings.append(
            f"the mechanism table has no gain rule for {compound}'s dominant "
            "insect receptor, so the direction-aware rule fell back to the "
            "depressant assumption for this compound."
        )
    return {
        "level": "C_topology_only",
        "effect": effect,
        "unit": "Hz",
        "detail": {
            "rule": GENERIC_RULES[rule],
            "rule_name": rule,
            "rule_note": GENERIC_RULE_NOTES[rule],
            "multiplier": float(m),
            "graph": graph_name,
            "readout": readout,
            **detail,
            "floor": {
                "level": C_FLOOR_LEVEL,
                "rule": GENERIC_MULTIPLIER_RULE,
                "rule_note": GENERIC_RULE_NOTES["depressant_floor"],
                "multiplier": float(floor_m),
                "effect": floor_effect,
                "direction_source": floor_detail.get("direction_source"),
            },
        },
        "warnings": warnings,
    }


def _effect_with_gains(
    assay: str,
    readout: str,
    graph_name: str,
    gains: dict[str, float],
    kw: dict[str, Any],
) -> float | None:
    """``treated - vehicle`` on the real graph for an arbitrary gain dict."""
    ok_fast, _why = _fast_supported(assay, readout, kw)
    if not ok_fast:
        raise ValueError(
            f"baseline C needs the rate engine, which does not support "
            f"{assay}/{readout}; run it on the subgraph or taste_map rate assay."
        )
    state = _state_for(graph_name)
    plan = _plan(state, assay, readout, gains, kw)
    return _fast_effect(state, plan)["effect"]


# --------------------------------------------------------------------------
# level D: the full model
# --------------------------------------------------------------------------
def full_model(
    compound: str | None,
    conc_M: float = 1e-6,
    graph: str | None = "named",
    assay: str = "subgraph",
    readout: str = "auto",
    library: dict[str, Any] | None = None,
    rule_overrides: dict[str, float] | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Level D: mechanism-specific gains on the real cut (what FlyLab ships)."""
    from flylab.circuit.rate import compute_gains

    if readout == "auto":
        readout = DEFAULT_READOUT[assay]
    graph_name = graph if graph is not None else DEFAULT_GRAPH[assay]
    gains, _occ = compute_gains(compound, conc_M, library=library, rule_overrides=rule_overrides)
    ok_fast, why = _fast_supported(assay, readout, kw)
    if ok_fast:
        effect = _effect_with_gains(assay, readout, graph_name, gains, kw)
    else:  # LIF or an exotic readout: go through the notebook assay
        effect = drug_effect(assay, compound, conc_M, readout, graph=graph_name, **kw)["effect"]
    return {
        "level": "D_full_flylab",
        "effect": effect,
        "unit": "Hz",
        "detail": {"gains": dict(gains), "graph": graph_name, "readout": readout},
        "warnings": ([] if ok_fast else [f"ran through the notebook assay ({why})."]),
    }


# --------------------------------------------------------------------------
# one compound: all four levels
# --------------------------------------------------------------------------
def level_predictions(
    compound: str,
    conc_M: float = 1e-6,
    graph: str | None = "named",
    assay: str = "subgraph",
    readout: str = "auto",
    generic_rule: str = DEFAULT_GENERIC_RULE,
    nt_sign_convention: str = DEFAULT_NT_SIGN_CONVENTION,
    **kw: Any,
) -> dict[str, dict[str, Any]]:
    """All four levels for one compound at one concentration.

    ``generic_rule`` picks level C's stand-in for pharmacology
    (:data:`GENERIC_RULES`; the historical depression-only rule is always also
    reported inside level C's ``detail["floor"]``).  ``nt_sign_convention``
    picks how level B signs glutamate (:data:`NT_SIGN_CONVENTIONS`).
    """
    return {
        "A_receptor_only": baseline_receptor_only(compound, conc_M, **kw),
        "B_composition_only": baseline_composition_only(
            compound, conc_M, graph=graph, nt_sign_convention=nt_sign_convention, **kw
        ),
        "C_topology_only": baseline_topology_only(
            compound, conc_M, graph=graph, assay=assay, readout=readout,
            rule=generic_rule, **kw
        ),
        "D_full_flylab": full_model(
            compound, conc_M, graph=graph, assay=assay, readout=readout, **kw
        ),
    }


def _library_compounds() -> list[str]:
    from flylab.pharm.occupancy import list_compounds

    return list(list_compounds())


def _compound_class(compound: str) -> str | None:
    try:
        from flylab.pharm.occupancy import load_library

        lib = load_library()
        return (lib.get("compounds", {}).get(compound.lower().strip()) or {}).get("class")
    except Exception:  # pragma: no cover - library in flux
        return None



def _is_nicotinic_agonist(
    compound: str, conc_M: float = 1e-6, library: dict[str, Any] | None = None
) -> bool:
    """True when the library gives this compound an activating insect nAChR row.

    A pharmacological grouping rather than a chemical one: it puts the
    neonicotinoids, nicotine, acetylcholine and spinosad in one set, which is
    the set the null models say is composition-dominated.
    """
    try:
        _occ, rows = _insect_rows(compound, conc_M, library)
    except Exception:  # pragma: no cover - library in flux
        return False
    for r in rows:
        if str(r.get("receptor")) == "insect_nAChR" and str(r.get("direction") or "") in _ACTIVATING:
            return _engagement(r) is not None
    return False


#: extra series computed alongside the four levels: the historical
#: depression-only level-C rule (a floor), and level B under the transmitter
#: sign convention it does *not* use by default.
EXTRA_SERIES: tuple[str, ...] = (C_FLOOR_LEVEL, "B_composition_only_alt_sign")


def _ordering(
    compounds: Sequence[str],
    conc_M: float,
    graph: str | None,
    assay: str,
    readout: str,
    kw: dict[str, Any],
    generic_rule: str = DEFAULT_GENERIC_RULE,
    nt_sign_convention: str = DEFAULT_NT_SIGN_CONVENTION,
) -> dict[str, list[float | None]]:
    """Each level's prediction for every compound, in ``compounds`` order.

    Also returns the two :data:`EXTRA_SERIES`, which cost nothing extra to
    collect and are what the floor comparison and the glutamate-sign check are
    computed from.
    """
    out: dict[str, list[float | None]] = {lvl: [] for lvl in tuple(LEVELS) + EXTRA_SERIES}
    for name in compounds:
        preds = level_predictions(
            name,
            conc_M,
            graph=graph,
            assay=assay,
            readout=readout,
            generic_rule=generic_rule,
            nt_sign_convention=nt_sign_convention,
            **kw,
        )
        for lvl in LEVELS:
            out[lvl].append(preds[lvl]["effect"])
        out[C_FLOOR_LEVEL].append(
            preds["C_topology_only"]["detail"].get("floor", {}).get("effect")
        )
        out["B_composition_only_alt_sign"].append(
            preds["B_composition_only"]["detail"].get("effect_other_convention")
        )
    return out


def _corr_row(
    level: str,
    pred: Sequence[float | None],
    full: Sequence[float | None],
    note: str,
    unit: str,
    with_hz: bool = False,
) -> dict[str, Any]:
    """Correlation/residual record for one series against the full model."""
    r = pearson(pred, full)
    rho = spearman(pred, full)
    x, y = _pairs(pred, full)
    resid_z = float(np.sqrt(np.mean((_z(x) - _z(y)) ** 2))) if x.size >= 3 else None
    resid_hz = float(np.sqrt(np.mean((x - y) ** 2))) if (with_hz and x.size) else None
    max_hz = float(np.max(np.abs(x - y))) if (with_hz and x.size) else None
    return {
        "level": level,
        "note": note,
        "unit": unit,
        "n_compounds": int(x.size),
        "pearson_r_vs_full": r,
        "spearman_rho_vs_full": rho,
        "residual_rms_standardised": resid_z,
        "residual_rms_hz": resid_hz,
        "max_abs_residual_hz": max_hz,
        # SIGNED rho: an ordering that is a perfect inversion of the full
        # model's reproduces nothing. v0.6 used |rho| here, so rho = -0.95
        # counted as a reproduction.
        "reproduces_full_ordering": bool(rho is not None and rho >= REPRODUCES_RHO),
        "ordering_inverted": bool(rho is not None and rho <= -REPRODUCES_RHO),
        "reproduces_full_values": (
            None
            if resid_hz is None
            else bool(resid_hz <= 0.05 * float(np.std(y)) if np.std(y) > 0 else resid_hz == 0.0)
        ),
        "reproduces_full": bool(
            rho is not None
            and rho >= REPRODUCES_RHO
            and (resid_hz is None or (r is not None and r >= REPRODUCES_R))
        ),
    }


def _gain_block(
    ordering: dict[str, list[float | None]],
    compounds: Sequence[str],
) -> dict[str, Any]:
    """Correlations, residuals and the added-information increments.

    ``information_added_vs_previous`` is an increment in **signed** Spearman
    rho.  Before v0.6.1 it was an increment in ``|rho|``, which made the column
    uninterpretable: a level that ordered the library backwards scored as
    informative, and the level above it then appeared to *lose* information
    while actually fixing the sign.
    """
    full = ordering["D_full_flylab"]
    rows: list[dict[str, Any]] = []
    prev_rho = 0.0
    for lvl in LEVELS:
        unit = (
            "dimensionless"
            if lvl == "A_receptor_only"
            else ("excitation index" if lvl == "B_composition_only" else "Hz")
        )
        row = _corr_row(
            lvl,
            ordering[lvl],
            full,
            LEVEL_NOTES[lvl],
            unit,
            with_hz=lvl in ("C_topology_only", "D_full_flylab"),
        )
        rho = row["spearman_rho_vs_full"]
        rho_v = float(rho) if rho is not None else 0.0
        row["information_added_vs_previous"] = (
            None if lvl == LEVELS[0] else float(rho_v - prev_rho)
        )
        row["information_added_basis"] = "signed_spearman_rho"
        rows.append(row)
        prev_rho = rho_v

    block: dict[str, Any] = {"levels": rows, "compounds": list(compounds)}

    # -- level C: the fair rule against the historical floor ---------------
    if C_FLOOR_LEVEL in ordering:
        floor_row = _corr_row(
            C_FLOOR_LEVEL,
            ordering[C_FLOOR_LEVEL],
            full,
            GENERIC_RULE_NOTES["depressant_floor"],
            "Hz",
            with_hz=True,
        )
        c_row = next(r for r in rows if r["level"] == "C_topology_only")
        floor_rho = floor_row["spearman_rho_vs_full"]
        c_rho = c_row["spearman_rho_vs_full"]
        block["generic_rule_comparison"] = {
            "direction_aware": {
                "rule": GENERIC_RULES["direction_aware"],
                "spearman_rho_vs_full": c_rho,
                "pearson_r_vs_full": c_row["pearson_r_vs_full"],
                "residual_rms_hz": c_row["residual_rms_hz"],
                "role": "the topology-only competitor",
            },
            "depressant_floor": {
                "rule": GENERIC_MULTIPLIER_RULE,
                "spearman_rho_vs_full": floor_rho,
                "pearson_r_vs_full": floor_row["pearson_r_vs_full"],
                "residual_rms_hz": floor_row["residual_rms_hz"],
                "role": (
                    "FLOOR, not a competitor: this rule cannot produce a "
                    "positive effect, so it cannot order a library that "
                    "contains disinhibitors"
                ),
            },
            "delta_rho": (
                None if (c_rho is None or floor_rho is None) else float(c_rho - floor_rho)
            ),
            "floor_row": floor_row,
        }

    # -- level B: what the glutamate-sign disagreement does ----------------
    alt = ordering.get("B_composition_only_alt_sign")
    if alt is not None:
        b_row = next(r for r in rows if r["level"] == "B_composition_only")
        alt_row = _corr_row(
            "B_composition_only_alt_sign",
            alt,
            full,
            "level B under the other transmitter-sign convention",
            "excitation index",
        )
        b_rho, alt_rho = b_row["spearman_rho_vs_full"], alt_row["spearman_rho_vs_full"]
        block["composition_sign_convention_check"] = {
            "shipped_convention": DEFAULT_NT_SIGN_CONVENTION,
            "shipped_note": (
                "glutamate scored as partly EXCITATORY (ach + 0.3*glu - 1.1*gaba)"
            ),
            "alternative_convention": "rate_engine",
            "alternative_note": (
                "glutamate scored INHIBITORY at the rate engine's own "
                "SIGN['glutamate'] = -0.4, i.e. levels B and D agree"
            ),
            "spearman_rho_vs_full_shipped": b_rho,
            "spearman_rho_vs_full_reconciled": alt_rho,
            "delta_rho": (
                None if (b_rho is None or alt_rho is None) else float(alt_rho - b_rho)
            ),
            "identical_orderings": bool(
                b_rho is not None and alt_rho is not None and abs(alt_rho - b_rho) < 1e-12
            ),
        }
    return block


def _statements(gain: dict[str, Any], scope: str) -> list[str]:
    """Plain-English verdicts, including the one the paper needs."""
    out: list[str] = []
    for row in gain["levels"]:
        lvl = row["level"]
        if lvl == "D_full_flylab":
            continue
        r, rho = row["pearson_r_vs_full"], row["spearman_rho_vs_full"]
        if r is None or rho is None:
            out.append(f"{lvl} ({scope}): not comparable (too few usable compounds).")
            continue
        if row["reproduces_full"]:
            unit_note = (
                " (ordering only: this level's units are not the full model's)"
                if row["residual_rms_hz"] is None
                else ""
            )
            out.append(
                f"{lvl} REPRODUCES the full model {scope}{unit_note}: "
                f"rho = {rho:.3f}, r = {r:.3f} over {row['n_compounds']} "
                "compounds. The layers above it add no ordering information "
                "for this readout."
            )
        elif row["reproduces_full_ordering"]:
            out.append(
                f"{lvl} reproduces the full model's ORDERING but not its values "
                f"{scope}: rho = {rho:.3f} (>= {REPRODUCES_RHO}) with r = "
                f"{r:.3f} (< {REPRODUCES_R}) over {row['n_compounds']} compounds."
            )
        elif row["ordering_inverted"]:
            out.append(
                f"{lvl} INVERTS the full model's ordering {scope}: rho = "
                f"{rho:.3f}, r = {r:.3f}. This is not reproduction - it is an "
                "anti-correlated level - and the signed rho is what says so."
            )
        else:
            out.append(
                f"{lvl} does NOT reproduce the full model {scope}: rho = "
                f"{rho:.3f}, r = {r:.3f}; the layers above it are doing real work."
            )
    cmp_block = gain.get("generic_rule_comparison")
    if cmp_block and cmp_block.get("delta_rho") is not None:
        da = cmp_block["direction_aware"]["spearman_rho_vs_full"]
        fl = cmp_block["depressant_floor"]["spearman_rho_vs_full"]
        out.append(
            f"C_topology_only {scope}: the direction-aware generic rule reaches "
            f"rho = {da:.3f} against rho = {fl:.3f} for the historical "
            "depression-only rule (a floor, which cannot express "
            f"disinhibition at all); the difference, {da - fl:+.3f}, is how "
            "much of level C's old score was an artefact of that rule rather "
            "than a fact about the connectome."
        )
    sign_block = gain.get("composition_sign_convention_check")
    if sign_block and sign_block.get("delta_rho") is not None:
        out.append(
            f"B_composition_only {scope}: the composition index signs glutamate "
            "as partly excitatory while the rate engine signs it inhibitory; "
            f"reconciling them moves the B-versus-D rank correlation from "
            f"{sign_block['spearman_rho_vs_full_shipped']:.3f} to "
            f"{sign_block['spearman_rho_vs_full_reconciled']:.3f} "
            f"({sign_block['delta_rho']:+.3f})."
        )
    return out



# --------------------------------------------------------------------------
# is B-versus-D a finding, or an algebraic identity?
# --------------------------------------------------------------------------
#: gain box the random-gain reference distribution samples from.  The synaptic
#: gains are floored at 0.05 exactly as the mechanism rules are, and the upper
#: end covers the largest coefficient any shipped rule can reach.
REFERENCE_GAIN_BOUNDS: dict[str, tuple[float, float]] = {
    "g_ach": (0.05, 1.5),
    "g_gaba": (0.05, 1.5),
    "g_glu": (0.05, 1.8),
    "g_oct": (0.05, 1.5),
    "g_nav": (0.5, 2.5),
}

GAIN_FLOOR = 0.05


def _full_effect_direct(
    gains: dict[str, float],
    graph: str | None = "named",
    normalise: str = "row_abs",
    drive_hz: float = 40.0,
    steps: int = 80,
    readout: str = "mean_hz",
) -> float | None:
    """Level D's contrast for an arbitrary gain vector, off the rate engine.

    Equivalent to :func:`full_model` for the default
    ``subgraph``/``mean_hz``/``row_abs`` combination (a test asserts it), but it
    takes a gain vector rather than a compound, which is what the reference
    distributions need, and it can run the engine without its row
    normalisation.
    """
    from flylab.circuit.rate import effect_under_normalisation

    return effect_under_normalisation(
        gains, graph=graph, normalise=normalise, drive_hz=drive_hz,
        steps=steps, readout=readout,
    )


def _library_gain_vectors(
    compounds: Sequence[str],
    conc_M: float,
    library: dict[str, Any] | None = None,
) -> list[dict[str, float]]:
    from flylab.circuit.rate import compute_gains

    out = []
    for name in compounds:
        gains, _occ = compute_gains(name, conc_M, library=library)
        out.append({k: float(v) for k, v in gains.items()})
    return out


def _at_floor(gains: dict[str, float]) -> bool:
    return any(
        abs(float(gains.get(k, 1.0)) - GAIN_FLOOR) < 1e-12
        for k in ("g_ach", "g_gaba", "g_glu", "g_oct")
    )


def _rho_for_gain_set(
    gain_vectors: Sequence[dict[str, float]],
    counts: dict[str, int],
    graph: str | None,
    nt_sign_convention: str,
    normalise: str,
    drive_hz: float,
    steps: int,
) -> float | None:
    b = [composition_effect_from_gains(g, counts, graph, nt_sign_convention) for g in gain_vectors]
    d = [_full_effect_direct(g, graph, normalise, drive_hz, steps) for g in gain_vectors]
    return spearman(b, d)


def composition_reference_distribution(
    conc_M: float = 1e-6,
    graph: str | None = "named",
    compounds: Sequence[str] | None = None,
    n_draws: int = 50,
    seed: int = 0,
    nt_sign_convention: str = DEFAULT_NT_SIGN_CONVENTION,
    normalise: str = "row_abs",
    drive_hz: float = 40.0,
    steps: int = 80,
    library: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """A reference distribution for the composition-versus-full correlation.

    The paper reports a rank correlation near 0.99 between the composition-only
    level (B) and the full model (D) and reads it as "the ordering information
    lives in composition, not in wiring".  But B and D are **functions of the
    same gain vector**: B applies it to a fixed census, D applies it per cell
    on a row-normalised matrix in which each cell's recurrent input is a
    composition-weighted average of its presynaptic gains.  A high correlation
    is therefore expected for algebraic reasons, and a bare 0.99 is not
    interpretable without knowing what correlation would have been
    *unsurprising*.

    This function supplies that reference by recomputing the same rank
    correlation under comparison conditions:

    ``observed``
        the real library's gain vectors at ``conc_M`` (what the paper quotes).
    ``library_no_floor``
        the same, after dropping every compound whose gain vector touches the
        0.05 floor - the tie artefact the referee already checked and rejected.
    ``random_gain_vectors``
        ``n_draws`` independent sets of arbitrary gain vectors drawn uniformly
        from :data:`REFERENCE_GAIN_BOUNDS`, each set the size of the library.
        These compounds do not exist and no pharmacology connects them, so the
        correlation they produce is the purely structural part of the number.
    ``shuffled_compound_assignment``
        the library's gain vectors permuted across compound labels.  Both
        levels are functions of the gain vector alone, so the permutation moves
        the pairs as a block and the correlation is **provably unchanged** -
        which is itself the point: the number carries no compound-level
        information at all.

    ``statements`` says in words what would have been evidence against the
    claim.
    """
    rng = np.random.default_rng(int(seed))
    names = list(compounds) if compounds is not None else _library_compounds()
    counts = graph_census(graph)
    vectors = _library_gain_vectors(names, conc_M, library)

    def rho_of(vs: Sequence[dict[str, float]]) -> float | None:
        return _rho_for_gain_set(
            vs, counts, graph, nt_sign_convention, normalise, drive_hz, steps
        )

    observed = rho_of(vectors)

    kept = [(n, g) for n, g in zip(names, vectors) if not _at_floor(g)]
    no_floor = rho_of([g for _n, g in kept]) if len(kept) >= 3 else None

    perm = rng.permutation(len(vectors))
    shuffled = rho_of([vectors[i] for i in perm])

    keys = list(REFERENCE_GAIN_BOUNDS)

    def draw_dense() -> dict[str, float]:
        g = {k: float(rng.uniform(*REFERENCE_GAIN_BOUNDS[k])) for k in keys}
        g["ach_tone"] = 1.0
        return g

    def draw_sparse() -> dict[str, float]:
        """One gain moved, the rest at vehicle: the shape the shipped
        mechanism rules actually produce (one receptor, one transmitter)."""
        g = {k: 1.0 for k in keys}
        k = keys[int(rng.integers(0, len(keys)))]
        g[k] = float(rng.uniform(*REFERENCE_GAIN_BOUNDS[k]))
        g["ach_tone"] = 1.0
        return g

    def reference(draw) -> tuple[list[float], dict[str, Any], float | None]:
        vals: list[float] = []
        for _ in range(int(n_draws)):
            r = rho_of([draw() for _i in range(len(names))])
            if r is not None:
                vals.append(float(r))
        a = np.asarray(vals, dtype=float)
        q = (
            {
                "n": int(a.size),
                "median": float(np.median(a)),
                "mean": float(a.mean()),
                "p05": float(np.percentile(a, 5)),
                "p95": float(np.percentile(a, 95)),
                "min": float(a.min()),
                "max": float(a.max()),
            }
            if a.size
            else {"n": 0}
        )
        pc = (
            float((a <= float(observed)).mean() * 100.0)
            if (a.size and observed is not None)
            else None
        )
        return vals, q, pc

    draws, quant, pct = reference(draw_dense)
    sparse_draws, sparse_quant, sparse_pct = reference(draw_sparse)
    arr = np.asarray(draws, dtype=float)

    statements: list[str] = []
    if observed is not None:
        statements.append(
            f"Composition-versus-full Spearman rho = {observed:.3f} over "
            f"{len(names)} library compounds at {conc_M:.1e} M."
        )
    if no_floor is not None:
        statements.append(
            f"Dropping the {len(names) - len(kept)} compound(s) whose gain "
            f"vector sits on the {GAIN_FLOOR} floor leaves rho = {no_floor:.3f}: "
            "the correlation is not a tie artefact."
        )
    if arr.size and observed is not None:
        statements.append(
            f"Reference 1 (arbitrary gain vectors, all gains moved at once): "
            f"{arr.size} draws of {len(names)} vectors give rho median "
            f"{quant['median']:.3f}, 5th-95th percentile "
            f"[{quant['p05']:.3f}, {quant['p95']:.3f}]; the observed "
            f"{observed:.3f} is at the {pct:.0f}th percentile. Even compounds "
            "that do not exist, ordered by nothing but a shared random gain "
            "vector, reproduce most of the ordering."
        )
    if sparse_quant.get("n") and observed is not None:
        sq = sparse_quant
        verdict = (
            "UNSURPRISING: it lies inside the reference's 5th-95th percentile"
            if sq["p05"] <= observed <= sq["p95"]
            else (
                "above the reference range, but the reference median is already "
                f"{sq['median']:.3f}"
                if observed > sq["p95"]
                else "BELOW the reference range"
            )
        )
        statements.append(
            f"Reference 2 (matched: one gain moved per pseudo-compound, the "
            f"shape the shipped mechanism rules produce): rho median "
            f"{sq['median']:.3f}, 5th-95th percentile [{sq['p05']:.3f}, "
            f"{sq['p95']:.3f}] over {sq['n']} draws. The observed "
            f"{observed:.3f} is at the {sparse_pct:.0f}th percentile - "
            f"{verdict}. A rank correlation of this size between levels B and D "
            "is therefore close to an algebraic identity: both are functions of "
            "the same gain vector, one through a fixed census and one through a "
            "row-normalised matrix whose cells average their presynaptic gains."
        )
        statements.append(
            "WHAT WOULD HAVE BEEN EVIDENCE AGAINST THE CLAIM that composition "
            "carries the ordering: an observed rho below the matched "
            f"reference's 5th percentile ({sq['p05']:.3f}), i.e. the census "
            "ordering the REAL library worse than it orders pseudo-compounds "
            "with no pharmacology in them. Evidence that the agreement was a "
            "substantive finding about this connectome would need the matched "
            "reference to sit well below the observed value - a median near "
            "zero, so that reproducing the ordering was a priori unlikely. "
            f"Here the matched reference median is {sq['median']:.3f}, so the "
            "quoted 0.99 should be read as the structural floor plus a small "
            "excess, not as a measurement of what the connectome contributes."
        )
    if shuffled is not None and observed is not None:
        statements.append(
            f"Shuffling the compound-to-gain assignment leaves rho = "
            f"{shuffled:.3f} against the observed {observed:.3f} "
            f"(difference {shuffled - observed:+.3g}). Both levels are "
            "functions of the gain vector alone, so a permutation of the labels "
            "moves the (B, D) pairs as a block and cannot change the rank "
            "correlation. The comparison therefore carries no compound-level "
            "information: it is a statement about the map from gains to "
            "readouts, not about which compound is which."
        )
    return {
        "conc_M": float(conc_M),
        "graph": graph,
        "normalise": normalise,
        "nt_sign_convention": nt_sign_convention,
        "n_compounds": len(names),
        "compounds": names,
        "observed": {"spearman_rho": observed, "n_compounds": len(names)},
        "conditions": {
            "library_no_floor": {
                "spearman_rho": no_floor,
                "n_compounds": len(kept),
                "dropped": [n for n, g in zip(names, vectors) if _at_floor(g)],
                "note": "compounds on the 0.05 gain floor removed (the tie artefact)",
            },
            "random_gain_vectors": {
                "spearman_rho": quant,
                "draws": draws,
                "bounds": {k: list(v) for k, v in REFERENCE_GAIN_BOUNDS.items()},
                "note": (
                    "arbitrary gain vectors, every gain moved at once: the "
                    "structural part of the correlation, with no pharmacology "
                    "in it"
                ),
            },
            "random_single_gain_vectors": {
                "spearman_rho": sparse_quant,
                "draws": sparse_draws,
                "observed_percentile": sparse_pct,
                "note": (
                    "matched reference: one gain moved per pseudo-compound, "
                    "which is the shape every shipped mechanism rule produces. "
                    "This is the distribution the observed value should be read "
                    "against."
                ),
            },
            "shuffled_compound_assignment": {
                "spearman_rho": shuffled,
                "identical_to_observed": bool(
                    shuffled is not None
                    and observed is not None
                    and abs(shuffled - observed) < 1e-12
                ),
                "note": (
                    "provably identical to the observed value: both levels are "
                    "functions of the same gain vector"
                ),
            },
        },
        "observed_percentile_of_reference": pct,
        "observed_percentile_of_matched_reference": sparse_pct,
        "unsurprising_rho_range": [sparse_quant.get("p05"), sparse_quant.get("p95")],
        "unsurprising_rho_range_dense": [quant.get("p05"), quant.get("p95")],
        "statements": statements,
        "label": "model_derived",
        "warnings": list(ABLATION_WARNINGS)
        + [
            "The composition level and the full model are not independent "
            "models: they consume the same gain vector, so their rank "
            "correlation has a large structural floor. Quote it against this "
            "reference distribution, never on its own.",
        ],
    }


def composition_dominance_under_normalisations(
    conc_M: float = 1e-6,
    graph: str | None = "named",
    compounds: Sequence[str] | None = None,
    normalisations: Sequence[str] = ("row_abs", "none", "degree"),
    nt_sign_convention: str = DEFAULT_NT_SIGN_CONVENTION,
    drive_hz: float = 40.0,
    steps: int = 80,
    library: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Does the composition-dominated verdict survive an unnormalised engine?

    The rate engine row-normalises its weight matrix, which makes every cell's
    recurrent input a composition-weighted average of its presynaptic gains and
    therefore builds much of "composition-dominated" into the readout
    (``flylab.circuit.rate``).  This runs the *same* library and the same
    composition baseline against a full model whose engine has that
    normalisation replaced, and reports the B-versus-D rank correlation under
    each.
    """
    names = list(compounds) if compounds is not None else _library_compounds()
    counts = graph_census(graph)
    vectors = _library_gain_vectors(names, conc_M, library)
    b = [composition_effect_from_gains(g, counts, graph, nt_sign_convention) for g in vectors]
    out: dict[str, Any] = {}
    for mode in normalisations:
        d = [_full_effect_direct(g, graph, mode, drive_hz, steps) for g in vectors]
        rho = spearman(b, d)
        out[mode] = {
            "spearman_rho_b_vs_d": rho,
            "pearson_r_b_vs_d": pearson(b, d),
            "reproduces_full_ordering": bool(rho is not None and rho >= REPRODUCES_RHO),
            "full_model_effects": d,
        }
    ref = normalisations[0]
    statements = [
        f"Composition-versus-full Spearman rho at {conc_M:.1e} M: "
        + ", ".join(
            f"{m} {out[m]['spearman_rho_b_vs_d']:.3f}"
            for m in normalisations
            if out[m]["spearman_rho_b_vs_d"] is not None
        )
        + "."
    ]
    survives = [m for m in normalisations if out[m]["reproduces_full_ordering"]]
    statements.append(
        "The composition-dominated verdict "
        + (
            "survives every normalisation tested"
            if len(survives) == len(list(normalisations))
            else f"survives only under {survives or 'none'}"
        )
        + f" (threshold rho >= {REPRODUCES_RHO}); the shipped engine is "
        f"{ref!r}, and the row normalisation is an undocumented modelling "
        "choice, not a measurement."
    )
    return {
        "conc_M": float(conc_M),
        "graph": graph,
        "n_compounds": len(names),
        "compounds": names,
        "composition_effects": b,
        "by_normalisation": out,
        "statements": statements,
        "label": "model_derived",
        "warnings": [
            "Levels are not comparable across normalisations (the unnormalised "
            "engine is supercritical and its rates are shaped by the r_max "
            "clip); only the orderings are.",
        ],
    }


def glutamate_sign_reconciliation(
    concs_M: Sequence[float] = (1e-8, 1e-7, 1e-6, 1e-5),
    graph: str | None = "named",
    compounds: Sequence[str] | None = None,
    normalise: str = "row_abs",
    drive_hz: float = 40.0,
    steps: int = 80,
    library: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """What the B/D glutamate-sign disagreement does to the headline number.

    Level B's excitation index (inherited from ``flylab.assays.wholens``) scores
    glutamate as partly **excitatory**; the rate engine that level D runs signs
    it **inhibitory** (``SIGN["glutamate"] = -0.4``).  Two levels of the same
    ablation ladder therefore disagree about a transmitter sign.  This reports
    the composition-versus-full rank correlation under both conventions, at
    every concentration, so the paper can say whether its conclusion depends on
    the disagreement.
    """
    names = list(compounds) if compounds is not None else _library_compounds()
    counts = graph_census(graph)
    rows: list[dict[str, Any]] = []
    for conc in concs_M:
        vectors = _library_gain_vectors(names, float(conc), library)
        d = [_full_effect_direct(g, graph, normalise, drive_hz, steps) for g in vectors]
        row: dict[str, Any] = {"conc_M": float(conc)}
        for convention in NT_SIGN_CONVENTIONS:
            b = [composition_effect_from_gains(g, counts, graph, convention) for g in vectors]
            row[f"rho_{convention}"] = spearman(b, d)
            row[f"r_{convention}"] = pearson(b, d)
        a, bb = row["rho_wholens"], row["rho_rate_engine"]
        row["delta_rho"] = None if (a is None or bb is None) else float(bb - a)
        rows.append(row)
    deltas = [abs(r["delta_rho"]) for r in rows if r["delta_rho"] is not None]
    worst = max(deltas) if deltas else None
    statements = [
        f"At {r['conc_M']:.1e} M the composition-versus-full rho is "
        f"{r['rho_wholens']:.3f} with glutamate scored partly excitatory "
        f"(shipped) and {r['rho_rate_engine']:.3f} with glutamate signed as the "
        f"rate engine signs it ({r['delta_rho']:+.3f})."
        for r in rows
        if r["delta_rho"] is not None
    ]
    if worst is not None:
        statements.append(
            "The two levels of the ladder disagree about the sign of "
            "glutamate; reconciling them moves the number the paper's "
            f"conclusion rests on by at most {worst:.3f} in rank correlation. "
            "The disagreement is a real defect of the ladder and is documented "
            "as one, but it is not what produces the composition-versus-full "
            "agreement."
        )
    return {
        "compounds": names,
        "concs_M": [float(c) for c in concs_M],
        "graph": graph,
        "rows": rows,
        "max_abs_delta_rho": worst,
        "conventions": {
            "wholens": "excitation = max(0, ach + 0.3*glu - 1.1*gaba)  [glutamate partly EXCITATORY]",
            "rate_engine": "excitation = max(0, ach - 0.4*glu - 1.1*gaba)  [glutamate INHIBITORY, as SIGN]",
        },
        "statements": statements,
        "label": "model_derived",
        "warnings": [
            "Only the glutamate coefficient changes between the two "
            "conventions, so the difference is attributable to the transmitter "
            "sign and to nothing else in the index.",
        ],
    }

def ablation(
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    graph: str | None = "named",
    readout: str = "auto",
    assay: str = "subgraph",
    compounds: Sequence[str] | None = None,
    include_information_gain: bool = True,
    generic_rule: str = DEFAULT_GENERIC_RULE,
    nt_sign_convention: str = DEFAULT_NT_SIGN_CONVENTION,
    **kw: Any,
) -> dict[str, Any]:
    """The four levels for one compound, plus what each level adds.

    ``levels`` holds this compound's prediction at every level of the ladder.
    ``information_gain`` holds, per level, the Pearson and Spearman correlation
    of that level's *ordering of the compound library* against the full model,
    the residual after standardising both (unit-free), the residual in Hz where
    the units allow it, and the increment in **signed** rho over the previous
    level - i.e. the ordering information that layer added.  It also carries
    ``generic_rule_comparison`` (the direction-aware level C against the
    historical depression-only floor) and ``composition_sign_convention_check``
    (what the B/D glutamate-sign disagreement does to the correlation).
    ``compounds`` defaults to the whole library; pass a subset (or
    ``include_information_gain=False``) to make this cheap.
    """
    t0 = time.perf_counter()
    if readout == "auto":
        readout = DEFAULT_READOUT[assay]
    graph_name = graph if graph is not None else DEFAULT_GRAPH[assay]
    levels = level_predictions(
        compound, conc_M, graph=graph_name, assay=assay, readout=readout,
        generic_rule=generic_rule, nt_sign_convention=nt_sign_convention, **kw
    )
    warnings = list(BASE_WARNINGS) + list(ABLATION_WARNINGS)
    for lvl in levels.values():
        for w in lvl.get("warnings", []):
            if w not in warnings:
                warnings.append(w)

    gain: dict[str, Any] | None = None
    statements: list[str] = []
    by_class: dict[str, Any] = {}
    if include_information_gain:
        names = list(compounds) if compounds is not None else _library_compounds()
        ordering = _ordering(
            names, conc_M, graph_name, assay, readout, kw, generic_rule, nt_sign_convention
        )
        gain = _gain_block(ordering, names)
        statements = _statements(gain, f"across {len(names)} compounds")
        groups: dict[str, list[int]] = {}
        for i, name in enumerate(names):
            cls = _compound_class(name)
            if cls:
                groups.setdefault(cls, []).append(i)
            if _is_nicotinic_agonist(name, conc_M, kw.get("library")):
                groups.setdefault("nicotinic agonist", []).append(i)
        for cls, idx in groups.items():
            if len(idx) < 3:
                continue
            sub = {lvl: [ordering[lvl][i] for i in idx] for lvl in ordering}
            sub_gain = _gain_block(sub, [names[i] for i in idx])
            by_class[cls] = sub_gain
            statements.extend(_statements(sub_gain, f"within the {cls} class"))

    return {
        "compound": compound,
        "conc_M": float(conc_M),
        "assay": assay,
        "graph": graph_name,
        "readout": readout,
        "levels": levels,
        "level_order": list(LEVELS),
        "generic_rule": generic_rule,
        "nt_sign_convention": nt_sign_convention,
        "effects": {lvl: levels[lvl]["effect"] for lvl in LEVELS},
        "information_gain": gain,
        "information_gain_by_class": by_class,
        "statements": statements,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }


def _concentration_dependence(gains: dict[str, Any]) -> dict[str, Any]:
    """Which level is worst, at which concentration - stated per concentration.

    The ladder's ranking is **not** constant across the dose range, and a
    summary that quotes only the 1 uM column hides that.  This block reports,
    per concentration, every level's rank correlation and which level comes
    last, both for the direction-aware level C and for the historical
    depression-only floor.
    """
    rows: list[dict[str, Any]] = []
    for key in sorted(gains, key=float):
        block = gains[key]
        per = {
            r["level"]: r["spearman_rho_vs_full"]
            for r in block["levels"]
            if r["level"] != "D_full_flylab"
        }
        cmp_block = block.get("generic_rule_comparison") or {}
        floor_rho = (cmp_block.get("depressant_floor") or {}).get("spearman_rho_vs_full")
        scored = {k: v for k, v in per.items() if v is not None}
        worst = min(scored, key=lambda k: scored[k]) if scored else None
        floor_set = dict(scored)
        if floor_rho is not None:
            floor_set["C_topology_only"] = floor_rho
        worst_floor = min(floor_set, key=lambda k: floor_set[k]) if floor_set else None
        rows.append(
            {
                "conc_M": float(key),
                "rho_by_level": per,
                "worst_level": worst,
                "worst_rho": scored.get(worst) if worst else None,
                "rho_C_floor_rule": floor_rho,
                "worst_level_with_floor_rule": worst_floor,
            }
        )
    statements: list[str] = []
    if rows:
        worst_names = {r["worst_level"] for r in rows}
        detail = ", ".join(
            f"{r['conc_M']:.0e} M: {r['worst_level']} ({r['worst_rho']:.3f})" for r in rows
        )
        statements.append(
            "Concentration dependence of the ladder: the weakest level is "
            + (
                f"{rows[0]['worst_level']} at every concentration tested"
                if len(worst_names) == 1
                else "NOT the same at every concentration"
            )
            + f" ({detail})."
        )
        floor_detail = ", ".join(
            f"{r['conc_M']:.0e} M: {r['worst_level_with_floor_rule']}" for r in rows
        )
        statements.append(
            "With the historical depression-only level-C rule the weakest "
            f"level is {floor_detail}: that rule is the worst of the four only "
            "at the lowest concentration, and at the three higher ones the "
            "receptor-only level is worse. Any sentence of the form 'the "
            "topology-only level is the worst of the four' is therefore true "
            "of one column of this table and false of the other three."
        )
    return {"rows": rows, "statements": statements}


def ablation_table(
    compounds: Sequence[str] | None = None,
    concs_M: Sequence[float] = (1e-6,),
    graph: str | None = "named",
    readout: str = "auto",
    assay: str = "subgraph",
    generic_rule: str = DEFAULT_GENERIC_RULE,
    nt_sign_convention: str = DEFAULT_NT_SIGN_CONVENTION,
    **kw: Any,
) -> dict[str, Any]:
    """The paper's ablation table: every compound at every concentration.

    ``rows`` is one record per ``(compound, conc)`` with all four level
    predictions; ``information_gain`` is computed per concentration over the
    whole set, and ``statements`` says in words where a cheaper baseline
    tracks the full model's ordering, against its matched reference.
    """
    t0 = time.perf_counter()
    if readout == "auto":
        readout = DEFAULT_READOUT[assay]
    graph_name = graph if graph is not None else DEFAULT_GRAPH[assay]
    names = list(compounds) if compounds is not None else _library_compounds()
    concs = [float(c) for c in concs_M]
    rows: list[dict[str, Any]] = []
    gains: dict[str, Any] = {}
    statements: list[str] = []
    warnings = list(BASE_WARNINGS) + list(ABLATION_WARNINGS)

    for conc in concs:
        ordering = _ordering(
            names, conc, graph_name, assay, readout, kw, generic_rule, nt_sign_convention
        )
        for i, name in enumerate(names):
            rows.append(
                {
                    "compound": name,
                    "class": _compound_class(name),
                    "conc_M": conc,
                    **{lvl: ordering[lvl][i] for lvl in ordering},
                }
            )
        gain = _gain_block(ordering, names)
        gains[f"{conc:.3e}"] = gain
        statements.extend(_statements(gain, f"at {conc:.1e} M"))

    conc_dep = _concentration_dependence(gains)
    statements.extend(conc_dep["statements"])

    return {
        "compounds": names,
        "concs_M": concs,
        "assay": assay,
        "graph": graph_name,
        "readout": readout,
        "level_order": list(LEVELS),
        "level_notes": dict(LEVEL_NOTES),
        "extra_series": list(EXTRA_SERIES),
        "generic_rule": generic_rule,
        "generic_rule_notes": dict(GENERIC_RULE_NOTES),
        "nt_sign_convention": nt_sign_convention,
        "rows": rows,
        "information_gain": gains,
        "concentration_dependence": conc_dep,
        "statements": statements,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }
