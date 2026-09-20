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
                         (``max(0.05, 1 - theta_max)`` applied to every
                         transmitter) instead of mechanism-specific gains.
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

The expected - and publishable - result is that for the nicotinic agonists the
composition-only baseline already reproduces the full model: their mean-rate
effect is an E/I-balance effect, which is exactly what the null models say
independently.  The output says so in words when it happens, and says the
opposite when a level fails.

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
    "REPRODUCES_R",
    "REPRODUCES_RHO",
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
        "topology only: the real cut and the real rate model with one generic "
        "multiplier per compound instead of mechanism-specific gains. Unit: Hz."
    ),
    "D_full_flylab": (
        "full FlyLab: mechanism-specific per-transmitter gains on the real cut. "
        "Unit: Hz."
    ),
}

#: Baseline C's stand-in for pharmacology: "the drug removes a fraction
#: theta_max of synaptic transmission, whatever the transmitter".
GENERIC_MULTIPLIER_RULE = "g_all = max(0.05, 1 - theta_max)"

#: a cheaper level "reproduces" the full model at or above these
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
    carried no number at all."""
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
        best = e if best is None else max(best, e)
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
def _settle(counts: dict[str, float], gains: dict[str, float], g_nav: float = 1.0) -> dict[str, float]:
    """Whole-CNS style excitation index, on whatever census it is handed.

    This is the ``flylab.assays.wholens`` index (imported when available so the
    two cannot drift), applied to the *graph's own* transmitter census instead
    of the 165k-cell MaleCNS one.
    """
    try:  # keep one implementation when the assay layer is importable
        from flylab.assays.wholens import _settle as _wholens_settle

        return _wholens_settle(counts, gains, g_nav)
    except Exception:  # pragma: no cover - fallback copy of the same formula
        tot = max(sum(counts.values()), 1)
        drive = {k: 40.0 * g_nav * (v / tot) for k, v in counts.items()}
        ach = drive.get("acetylcholine", 0.0) * gains["acetylcholine"]
        gaba = drive.get("gaba", 0.0) * gains["gaba"]
        glu = drive.get("glutamate", 0.0) * gains.get("glutamate", 1.0)
        excitation = max(0.0, ach + 0.3 * glu - 1.1 * gaba)
        inhibition = max(0.0, gaba + 0.7 * glu)
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


def baseline_composition_only(
    compound: str | None,
    conc_M: float = 1e-6,
    graph: str | None = "named",
    library: dict[str, Any] | None = None,
    rule_overrides: dict[str, float] | None = None,
    **_: Any,
) -> dict[str, Any]:
    """Baseline B: mechanism gains x the graph's global transmitter proportions.

    No edges: the cut is reduced to "this fraction of the cells is
    cholinergic, this fraction GABAergic", and the drug's gains are applied to
    those proportions.  ``effect`` is the change in the excitation index.
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
    treated = _settle(counts, nt_gains, gains["g_nav"])
    vehicle = _settle(counts, {k: 1.0 for k in nt_gains}, 1.0)
    total = max(sum(counts.values()), 1)
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
        },
        "warnings": [],
    }


# --------------------------------------------------------------------------
# level C: topology only
# --------------------------------------------------------------------------
def generic_multiplier(
    compound: str | None,
    conc_M: float = 1e-6,
    library: dict[str, Any] | None = None,
) -> tuple[float, dict[str, Any]]:
    """One number per compound: ``max(0.05, 1 - theta_max)``.

    ``theta_max`` is the largest modelled insect engagement with a direction.
    This is the "generic pharmacology" a connectome-only model would use: the
    drug perturbs synaptic transmission by that fraction, and the model has no
    idea which transmitter it acts on.
    """
    if not compound:
        return 1.0, {"theta_max": None, "n_not_modelled": 0}
    _occ, rows = _insect_rows(compound, conc_M, library)
    theta, not_modelled = _theta_max(rows)
    if theta is None:
        return 1.0, {"theta_max": None, "n_not_modelled": not_modelled}
    return max(0.05, 1.0 - float(theta)), {
        "theta_max": float(theta),
        "n_not_modelled": not_modelled,
    }


def baseline_topology_only(
    compound: str | None,
    conc_M: float = 1e-6,
    graph: str | None = "named",
    assay: str = "subgraph",
    readout: str = "auto",
    library: dict[str, Any] | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Baseline C: the real cut, the real rate model, one generic multiplier.

    The multiplier (:func:`generic_multiplier`, ``max(0.05, 1 - theta_max)``)
    is applied to **every** transmitter gain, so the model keeps all of the
    wiring and none of the mechanism.  ``effect`` is ``treated - vehicle`` in
    Hz, the same contrast and the same units as the full model.
    """
    from flylab.circuit.rate import DEFAULT_GAINS

    if readout == "auto":
        readout = DEFAULT_READOUT[assay]
    graph_name = graph if graph is not None else DEFAULT_GRAPH[assay]
    m, detail = generic_multiplier(compound, conc_M, library)
    gains = dict(DEFAULT_GAINS)
    for key in ("g_ach", "g_gaba", "g_glu", "g_oct"):
        gains[key] = float(m)
    effect = _effect_with_gains(assay, readout, graph_name, gains, kw)
    return {
        "level": "C_topology_only",
        "effect": effect,
        "unit": "Hz",
        "detail": {
            "rule": GENERIC_MULTIPLIER_RULE,
            "multiplier": float(m),
            "graph": graph_name,
            "readout": readout,
            **detail,
        },
        "warnings": (
            []
            if detail.get("theta_max") is not None
            else [
                f"no modelled insect engagement for {compound}: baseline C used "
                "a multiplier of 1.0, i.e. it predicts no effect."
            ]
        ),
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
    **kw: Any,
) -> dict[str, dict[str, Any]]:
    """All four levels for one compound at one concentration."""
    return {
        "A_receptor_only": baseline_receptor_only(compound, conc_M, **kw),
        "B_composition_only": baseline_composition_only(compound, conc_M, graph=graph, **kw),
        "C_topology_only": baseline_topology_only(
            compound, conc_M, graph=graph, assay=assay, readout=readout, **kw
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


def _ordering(
    compounds: Sequence[str],
    conc_M: float,
    graph: str | None,
    assay: str,
    readout: str,
    kw: dict[str, Any],
) -> dict[str, list[float | None]]:
    """Each level's prediction for every compound, in ``compounds`` order."""
    out: dict[str, list[float | None]] = {lvl: [] for lvl in LEVELS}
    for name in compounds:
        preds = level_predictions(name, conc_M, graph=graph, assay=assay, readout=readout, **kw)
        for lvl in LEVELS:
            out[lvl].append(preds[lvl]["effect"])
    return out


def _gain_block(
    ordering: dict[str, list[float | None]],
    compounds: Sequence[str],
) -> dict[str, Any]:
    """Correlations, residuals and the added-information increments."""
    full = ordering["D_full_flylab"]
    rows: list[dict[str, Any]] = []
    prev_abs_rho = 0.0
    for lvl in LEVELS:
        pred = ordering[lvl]
        r = pearson(pred, full)
        rho = spearman(pred, full)
        x, y = _pairs(pred, full)
        resid_z = float(np.sqrt(np.mean((_z(x) - _z(y)) ** 2))) if x.size >= 3 else None
        resid_hz = None
        max_hz = None
        if lvl in ("C_topology_only", "D_full_flylab") and x.size:
            resid_hz = float(np.sqrt(np.mean((x - y) ** 2)))
            max_hz = float(np.max(np.abs(x - y)))
        abs_rho = abs(rho) if rho is not None else 0.0
        rows.append(
            {
                "level": lvl,
                "note": LEVEL_NOTES[lvl],
                "unit": (
                    "dimensionless"
                    if lvl == "A_receptor_only"
                    else ("excitation index" if lvl == "B_composition_only" else "Hz")
                ),
                "n_compounds": int(x.size),
                "pearson_r_vs_full": r,
                "spearman_rho_vs_full": rho,
                "residual_rms_standardised": resid_z,
                "residual_rms_hz": resid_hz,
                "max_abs_residual_hz": max_hz,
                "information_added_vs_previous": (
                    None if lvl == LEVELS[0] else float(abs_rho - prev_abs_rho)
                ),
                # the levels do not share units, so "reproduces" is a statement
                # about the ORDERING; where the units do match (C vs D, both in
                # Hz) the values are checked too.
                "reproduces_full_ordering": bool(rho is not None and abs(rho) >= REPRODUCES_RHO),
                "reproduces_full_values": (
                    None
                    if resid_hz is None
                    else bool(resid_hz <= 0.05 * float(np.std(y)) if np.std(y) > 0 else resid_hz == 0.0)
                ),
                "reproduces_full": bool(
                    rho is not None
                    and abs(rho) >= REPRODUCES_RHO
                    and (
                        resid_hz is None
                        or (r is not None and abs(r) >= REPRODUCES_R)
                    )
                ),
            }
        )
        prev_abs_rho = abs_rho
    return {"levels": rows, "compounds": list(compounds)}


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
        else:
            out.append(
                f"{lvl} does NOT reproduce the full model {scope}: rho = "
                f"{rho:.3f}, r = {r:.3f}; the layers above it are doing real work."
            )
    return out


def ablation(
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    graph: str | None = "named",
    readout: str = "auto",
    assay: str = "subgraph",
    compounds: Sequence[str] | None = None,
    include_information_gain: bool = True,
    **kw: Any,
) -> dict[str, Any]:
    """The four levels for one compound, plus what each level adds.

    ``levels`` holds this compound's prediction at every level of the ladder.
    ``information_gain`` holds, per level, the Pearson and Spearman correlation
    of that level's *ordering of the compound library* against the full model,
    the residual after standardising both (unit-free), the residual in Hz where
    the units allow it, and the increment in ``|rho|`` over the previous level -
    i.e. the ordering information that layer added.  ``compounds`` defaults to the whole library;
    pass a subset (or ``include_information_gain=False``) to make this cheap.
    """
    t0 = time.perf_counter()
    if readout == "auto":
        readout = DEFAULT_READOUT[assay]
    graph_name = graph if graph is not None else DEFAULT_GRAPH[assay]
    levels = level_predictions(
        compound, conc_M, graph=graph_name, assay=assay, readout=readout, **kw
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
        ordering = _ordering(names, conc_M, graph_name, assay, readout, kw)
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
            sub = {lvl: [ordering[lvl][i] for i in idx] for lvl in LEVELS}
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
        "effects": {lvl: levels[lvl]["effect"] for lvl in LEVELS},
        "information_gain": gain,
        "information_gain_by_class": by_class,
        "statements": statements,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }


def ablation_table(
    compounds: Sequence[str] | None = None,
    concs_M: Sequence[float] = (1e-6,),
    graph: str | None = "named",
    readout: str = "auto",
    assay: str = "subgraph",
    **kw: Any,
) -> dict[str, Any]:
    """The paper's ablation table: every compound at every concentration.

    ``rows`` is one record per ``(compound, conc)`` with all four level
    predictions; ``information_gain`` is computed per concentration over the
    whole set, and ``statements`` says in words where a cheaper baseline
    already reproduces the full model.
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
        ordering = _ordering(names, conc, graph_name, assay, readout, kw)
        for i, name in enumerate(names):
            rows.append(
                {
                    "compound": name,
                    "class": _compound_class(name),
                    "conc_M": conc,
                    **{lvl: ordering[lvl][i] for lvl in LEVELS},
                }
            )
        gain = _gain_block(ordering, names)
        gains[f"{conc:.3e}"] = gain
        statements.extend(_statements(gain, f"at {conc:.1e} M"))

    return {
        "compounds": names,
        "concs_M": concs,
        "assay": assay,
        "graph": graph_name,
        "readout": readout,
        "level_order": list(LEVELS),
        "level_notes": dict(LEVEL_NOTES),
        "rows": rows,
        "information_gain": gains,
        "statements": statements,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }
