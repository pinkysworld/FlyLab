r"""At what network scale does a connectome-dependence verdict stabilise?

The question this module exists to answer
-----------------------------------------
FlyLab's dependence verdict reversed between its two committed cuts.  On the
``named`` cut (1 126 nodes, 1 360 edges -- an in-star in which 85 % of edges
land on the four seed cells and only 16 % of nodes have any input at all) most compounds come out *composition-dominated*.  On
``taste_motor`` (1 841 nodes, 19 066 edges, mean degree 10.4) not one landscape
cell is composition-dominated and everything with an effect is topology-
dependent.  The verdict is therefore a property of the substrate as much as of
the drug, and nothing in the method says where -- or whether -- it settles.

This module runs the *same* dependence profile across a ladder of nested cuts
of increasing size (:mod:`flylab.maps.extract`) and reports, per cut, the
verdict, the empirical permutation probability per null mode, the equivalence
gap, the necessary information level, and the cut's structural statistics, so a
reader can see **where** the verdict changes and **what changed structurally**
when it did.

It reuses :func:`flylab.analysis.dependence.dependence_profile` unchanged; no
statistic is reimplemented here.

Cost, honestly
--------------
A dependence profile costs, on one core,

.. math::

    T \approx n_\text{modes} \cdot n \cdot (a + b E)

with :math:`a = 1.09\times10^{-2}` s of per-shuffle overhead and
:math:`b = 6.86\times10^{-6}` s per edge, fitted on this machine over the
committed cuts and the scale ladder (:data:`PROFILE_COST_A`,
:data:`PROFILE_COST_B`; the fit reproduces the measured 1k-cut profile to
0.1 %).  Because :math:`E` grows faster than the node count along any realistic
cut ladder, a fixed ``n`` is not affordable at every rung, so
:func:`permutation_budget` scales ``n`` down with the cut and
:func:`estimate_scale_runtime` says what a study will cost **before** it is
started.  Scaling ``n`` down is not free: the permutation resolution is
:math:`1/(n+1)`, so below :math:`n = 1/\alpha - 1` no mode can be distinguished
at all and a "composition-dominated" verdict at the top of the ladder would be
uninterpretable.  :func:`permutation_budget` refuses to go below that floor and
every row records the ``n`` it was run at.
"""
from __future__ import annotations

import math
import time
from pathlib import Path
from typing import Any, Sequence

from flylab.analysis.dependence import (
    CONFIRMATORY_COMPOUNDS,
    DEFAULT_ALPHA,
    DEFAULT_MODES,
    PAPER_N,
    cut_census,
    dependence_profile,
)
from flylab.circuit.rate import DATA_DIRS, GRAPHS, load_graph
# the ladder's constants and naming live in a pandas-free module on purpose:
# this module is reachable from flylab.browser.bridge, which must import with
# pandas and pyarrow hidden (tests/test_browser_bridge.py enforces it).
# flylab.maps.extract, which needs pandas and pyarrow to *build* a cut, is
# never imported here.
from flylab.maps.ladder import LADDER_RECIPE, LADDER_SIZES, ladder_filename

# --------------------------------------------------------------------------
# the ladder of cuts
# --------------------------------------------------------------------------
#: Cut name -> filename, for every cut this module knows about: the two
#: committed neighbourhood cuts plus the scale ladder.  A rung that is not on
#: disk is skipped with a note rather than failing the study -- only the
#: smallest rung is small enough to commit; the rest are CI artifacts.
CUT_FILES: dict[str, str] = dict(GRAPHS)
CUT_FILES.update({f"scale_{k // 1000}k": ladder_filename(k) for k in LADDER_SIZES})

#: The default ladder, weakest substrate first.  ``named`` and ``taste_motor``
#: are included because they are the two cuts whose disagreement raised the
#: question; they are built by a *different* recipe from the ``scale_*`` rungs
#: and every report says so.
DEFAULT_CUTS: tuple[str, ...] = ("named", "taste_motor") + tuple(
    f"scale_{k // 1000}k" for k in LADDER_SIZES
)

#: Which recipe built each cut.  A verdict difference between two cuts built by
#: different recipes confounds scale with construction; within the ``scale_*``
#: rungs it does not, because they are nested views of one growth order.
CUT_RECIPE: dict[str, str] = {
    "named": "hops_neighborhood (1 hop from MN9/DNp01, min_weight 5, no closure)",
    "taste_motor": (
        "hops_neighborhood (1 hop from MN9/DNp01 + 6 labellar GRN types, "
        "min_weight 5, induced closure at weight 10)"
    ),
}
for _k in LADDER_SIZES:
    CUT_RECIPE[f"scale_{_k // 1000}k"] = "scale_ladder (ranked BFS, induced, min_weight 5)"

#: Fitted per-shuffle overhead, seconds (one core, this engine).
PROFILE_COST_A = 1.09e-2
#: Fitted marginal cost per edge per shuffle, seconds.
PROFILE_COST_B = 6.86e-6
#: Measured cost of one rate step per edge, seconds (sparse scatter), at
#: whole-CNS size.  It is *not* constant with size: the gather ``r[pre]``
#: misses cache once the rate vector leaves L2, so the same kernel costs
#: 3.4e-9 s/edge/step on a 5k cut and 1.48e-8 on the full connectome.
#: :data:`MEASURED_RUN_COSTS` has the measurements; :func:`run_cost_s`
#: interpolates between them instead of extrapolating one constant.
RUN_COST_PER_EDGE_STEP = 1.48e-8

#: ``(n_edges, seconds for one 80-step run)``, all measured on this machine
#: with the sparse engine.  These are measurements, not a model.
MEASURED_RUN_COSTS: tuple[tuple[int, float], ...] = (
    (1_360, 0.0020),
    (19_066, 0.0060),
    (22_857, 0.0110),
    (279_845, 0.0770),
    (621_600, 0.1710),
    (1_364_375, 0.4690),
    (2_216_881, 0.9670),
    (25_563_197, 30.2600),
)


def run_cost_s(n_edges: int, steps: int = 80) -> float:
    """Seconds for one rate run, log-interpolated over :data:`MEASURED_RUN_COSTS`."""
    pts = MEASURED_RUN_COSTS
    e = float(max(int(n_edges), 1))
    if e <= pts[0][0]:
        base = pts[0][1] * e / pts[0][0]
    elif e >= pts[-1][0]:
        base = pts[-1][1] * e / pts[-1][0]
    else:
        base = pts[-1][1]
        for (e0, t0), (e1, t1) in zip(pts, pts[1:]):
            if e0 <= e <= e1:
                f = (math.log(e) - math.log(e0)) / (math.log(e1) - math.log(e0))
                base = math.exp(math.log(t0) + f * (math.log(t1) - math.log(t0)))
                break
    return float(base * int(steps) / 80.0)

#: Below this ``n`` the permutation resolution 1/(n+1) is coarser than
#: alpha = 0.05 and NO mode can be distinguished, whatever the wiring does.
MIN_USEFUL_N = 20


def cut_path(name: str | Path) -> Path | None:
    """Resolve a cut name to a file on disk, or ``None`` if it is not there."""
    if isinstance(name, Path) or (isinstance(name, str) and name not in CUT_FILES):
        p = Path(name)
        return p if p.exists() else None
    fname = CUT_FILES[str(name)]
    for d in DATA_DIRS:
        p = Path(d) / fname
        if p.exists():
            return p
    return None


def available_cuts(cuts: Sequence[str | Path] = DEFAULT_CUTS) -> list[dict[str, Any]]:
    """The subset of ``cuts`` present on this machine, smallest first.

    Each entry carries ``name``, ``path``, ``n_nodes``, ``n_edges`` and
    ``recipe``; the missing ones are reported by :func:`missing_cuts`.
    """
    out: list[dict[str, Any]] = []
    for name in cuts:
        p = cut_path(name)
        if p is None:
            continue
        g = load_graph(p)
        out.append(
            {
                "name": str(name),
                "path": str(p),
                "n_nodes": int(g.get("n_nodes") or len(g["nodes"])),
                "n_edges": int(g.get("n_edges") or len(g["edges"])),
                "recipe": CUT_RECIPE.get(str(name), g.get("growth") or "unknown"),
            }
        )
    out.sort(key=lambda c: (c["n_nodes"], c["n_edges"]))
    return out


def missing_cuts(cuts: Sequence[str | Path] = DEFAULT_CUTS) -> list[str]:
    """Names in ``cuts`` with no file on disk (build them with the Action)."""
    return [str(c) for c in cuts if cut_path(c) is None]


# --------------------------------------------------------------------------
# cost
# --------------------------------------------------------------------------
def profile_cost_s(n_edges: int, n: int, n_modes: int = len(DEFAULT_MODES)) -> float:
    """Predicted serial runtime of one :func:`dependence_profile`, seconds."""
    return float(n_modes) * float(n) * (PROFILE_COST_A + PROFILE_COST_B * float(n_edges))


def permutation_budget(
    n_edges: int,
    seconds: float = 900.0,
    n_modes: int = len(DEFAULT_MODES),
    max_n: int = PAPER_N,
    min_n: int = MIN_USEFUL_N,
    alpha: float = DEFAULT_ALPHA,
) -> dict[str, Any]:
    """How many permutations fit in ``seconds`` on this cut, and what that costs.

    Returns the chosen ``n``, whether it hit the ceiling or the floor, the
    permutation resolution ``1/(n+1)`` it buys, and -- when the floor bound --
    a warning that the budget is not enough for a test at this ``alpha``.
    """
    per = PROFILE_COST_A + PROFILE_COST_B * float(n_edges)
    raw = float(seconds) / max(per * float(n_modes), 1e-12)
    n = int(max(min(math.floor(raw), int(max_n)), int(min_n)))
    res = 1.0 / (n + 1)
    warnings: list[str] = []
    if n >= int(max_n):
        note = "ceiling"
    elif raw < min_n:
        note = "floor"
        warnings.append(
            f"the {float(seconds):.0f} s budget buys only {raw:.1f} permutations on a "
            f"{int(n_edges)}-edge cut; n was raised to the floor {int(min_n)}, so this "
            f"profile will cost about {profile_cost_s(n_edges, n, n_modes):.0f} s, not "
            f"{float(seconds):.0f} s."
        )
    else:
        note = "budget"
    if res > float(alpha):
        warnings.append(
            f"n = {n} gives a permutation resolution of {res:.4f}, coarser than "
            f"alpha = {float(alpha)}: no mode can be distinguished from the real "
            "graph at this n, so a composition-dominated verdict here says "
            "nothing about the wiring."
        )
    if n < PAPER_N:
        warnings.append(
            f"n = {n} < {PAPER_N}: this rung is exploratory, not one of the "
            "prespecified confirmatory tests."
        )
    return {
        "n": n,
        "n_edges": int(n_edges),
        "bound_by": note,
        "seconds_budget": float(seconds),
        "predicted_s": profile_cost_s(n_edges, n, n_modes),
        "p_resolution": res,
        "resolution_coarser_than_alpha": bool(res > float(alpha)),
        "warnings": warnings,
    }


def estimate_scale_runtime(
    compounds: Sequence[str],
    cuts: Sequence[str | Path] = DEFAULT_CUTS,
    n: int | None = None,
    seconds_per_cut: float = 900.0,
    n_jobs: int = 1,
) -> dict[str, Any]:
    """What a :func:`dependence_vs_scale` call will cost, before starting it."""
    rows = []
    total = 0.0
    for c in available_cuts(cuts):
        budget = (
            {"n": int(n), "predicted_s": profile_cost_s(c["n_edges"], int(n)), "warnings": []}
            if n is not None
            else permutation_budget(c["n_edges"], seconds=seconds_per_cut)
        )
        per_cut = budget["predicted_s"] * len(compounds)
        total += per_cut
        rows.append({**c, "n": budget["n"], "predicted_s": per_cut})
    return {
        "compounds": list(compounds),
        "rows": rows,
        "n_jobs": int(n_jobs),
        "estimate_s_serial": total,
        "estimate_s": total / max(int(n_jobs), 1),
        "missing_cuts": missing_cuts(cuts),
        "note": (
            "fitted from measured profile runtimes on this engine "
            f"(a = {PROFILE_COST_A:.3g} s, b = {PROFILE_COST_B:.3g} s/edge); "
            "the parallel figure is optimistic because the weight-matched "
            "transmitter null runs serially."
        ),
    }


# --------------------------------------------------------------------------
# the study
# --------------------------------------------------------------------------
SCALE_WARNINGS: list[str] = [
    "A verdict compared across cuts is only a statement about scale when the "
    "cuts were built by the same recipe. 'named' and 'taste_motor' are hops-"
    "limited neighbourhoods and the 'scale_*' rungs are nested ranked-BFS "
    "induced cuts, so a difference between the two families confounds scale "
    "with construction; differences WITHIN the scale_* rungs do not.",
    "The permutation count is scaled down with cut size, so the rungs are not "
    "equally powered. A higher rung failing to reject where a lower one "
    "rejected can mean a weaker test rather than a weaker dependence; the "
    "p_resolution column is the check.",
    "The network-mean readout dilutes as the cut grows (the same seed cells "
    "drive ever more cells), so the raw effect size falls with scale by "
    "construction. The effect floor and the equivalence margin are both "
    "relative to the vehicle readout, which absorbs most of this, but the "
    "absolute Hz columns are not comparable across rungs.",
]


def _mode_rows(profile: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for row in profile.get("modes") or []:
        out.append(
            {
                "mode": row["mode"],
                "p_two_sided": row.get("p_two_sided"),
                "p_adjusted": row.get("p_adjusted"),
                "z": row.get("z"),
                "verdict": row.get("verdict"),
                "null_median": row.get("null_median"),
                "equivalence_gap": row.get("abs_gap_from_null_median"),
                "n_ok": row.get("n_ok"),
                "stabilised": row.get("stabilised"),
            }
        )
    return out


def dependence_vs_scale(
    compounds: Sequence[str] = CONFIRMATORY_COMPOUNDS,
    cuts: Sequence[str | Path] = DEFAULT_CUTS,
    conc_M: float = 1e-6,
    n: int | None = None,
    seconds_per_cut: float = 900.0,
    assay: str = "subgraph",
    readout: str = "auto",
    seed: int = 0,
    alpha: float = DEFAULT_ALPHA,
    n_jobs: int = 1,
    modes: Sequence[str] = DEFAULT_MODES,
    census: bool = True,
    progress: Any = None,
) -> dict[str, Any]:
    """The same dependence profile across a ladder of cuts of increasing size.

    Args:
        compounds: compounds to profile at every rung.
        cuts: cut names or paths; rungs absent from disk are skipped and named
            in ``missing_cuts``.
        conc_M: one concentration (the scaling question is about the graph, not
            the dose; use :func:`flylab.analysis.dependence.dependence_landscape`
            for the dose axis).
        n: fixed permutation count.  Left ``None``, each rung gets
            :func:`permutation_budget` for ``seconds_per_cut`` -- the honest
            default, because a fixed ``n`` is either unaffordable at the top or
            wasteful at the bottom.
        census: also compute each cut's
            :func:`~flylab.analysis.dependence.cut_census` (mean degree, share
            of edges onto seeds, share of nodes with input, in-star flag).

    Returns a block with one ``rows`` entry per (cut, compound) carrying the
    class, the per-mode verdicts and probabilities, the equivalence gaps, the
    necessary information level and the cut's structure, plus a
    ``stability`` block from :func:`verdict_stability`.
    """
    t0 = time.perf_counter()
    cut_list = available_cuts(cuts)
    if not cut_list:
        raise FileNotFoundError(
            "none of the requested cuts is on disk: "
            + ", ".join(str(c) for c in cuts)
            + ". Build the ladder with flylab.maps.extract.extract_ladder, or "
            "download the extract-malecns-subgraph Action's artifact."
        )
    rows: list[dict[str, Any]] = []
    warnings = list(SCALE_WARNINGS)
    census_by_cut: dict[str, Any] = {}
    for c in cut_list:
        if census:
            try:
                census_by_cut[c["name"]] = cut_census(c["path"])
            except Exception as exc:  # pragma: no cover - defensive
                warnings.append(f"cut_census failed for {c['name']}: {exc}")
        budget = (
            {"n": int(n), "p_resolution": 1.0 / (int(n) + 1), "bound_by": "fixed",
             "predicted_s": profile_cost_s(c["n_edges"], int(n)), "warnings": []}
            if n is not None
            else permutation_budget(c["n_edges"], seconds=seconds_per_cut, alpha=alpha)
        )
        for w in budget["warnings"]:
            tagged = f"{c['name']}: {w}"
            if tagged not in warnings:
                warnings.append(tagged)
        for compound in compounds:
            if progress is not None:
                progress(c["name"], compound, budget["n"])
            t = time.perf_counter()
            prof = dependence_profile(
                compound,
                conc_M,
                assay=assay,
                readout=readout,
                graph=c["path"],
                modes=tuple(modes),
                n=int(budget["n"]),
                seed=seed,
                alpha=alpha,
                n_jobs=n_jobs,
            )
            cen = census_by_cut.get(c["name"]) or {}
            lvl = prof.get("necessary_information_level") or {}
            rows.append(
                {
                    "cut": c["name"],
                    "path": c["path"],
                    "recipe": c["recipe"],
                    "compound": compound,
                    "conc_M": float(conc_M),
                    "n_nodes": c["n_nodes"],
                    "n_edges": c["n_edges"],
                    "mean_degree": cen.get("mean_degree"),
                    "share_onto_seeds": cen.get("share_onto_seeds"),
                    "share_nodes_with_in_degree": cen.get("share_nodes_with_in_degree"),
                    "share_edges_from_nodes_with_input": cen.get(
                        "share_edges_from_nodes_with_input"
                    ),
                    "in_star": cen.get("in_star"),
                    "n": int(budget["n"]),
                    "p_resolution": prof.get("p_resolution"),
                    "resolution_coarser_than_alpha": prof.get(
                        "resolution_coarser_than_alpha"
                    ),
                    "confirmatory": prof.get("confirmatory"),
                    "class": prof.get("class"),
                    "real_effect": prof.get("real_effect"),
                    "real_vehicle": prof.get("real_vehicle"),
                    "relative_effect": (
                        None
                        if not prof.get("real_vehicle") or prof.get("real_effect") is None
                        else float(prof["real_effect"]) / float(prof["real_vehicle"])
                    ),
                    "delta": prof.get("delta"),
                    "p": dict(prof.get("p") or {}),
                    "verdicts": dict(prof.get("verdicts") or {}),
                    "modes": _mode_rows(prof),
                    "necessary_information_level": lvl.get("level"),
                    "necessary_information_rank": lvl.get("rank"),
                    "necessary_information_verdict": lvl.get("verdict"),
                    "runtime_s": float(time.perf_counter() - t),
                }
            )
    result = {
        "compounds": list(compounds),
        "conc_M": float(conc_M),
        "assay": assay,
        "readout": readout,
        "alpha": float(alpha),
        "seed": int(seed),
        "modes": list(modes),
        "cuts": cut_list,
        "missing_cuts": missing_cuts(cuts),
        "census": census_by_cut,
        "rows": rows,
        "ladder_recipe": LADDER_RECIPE,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }
    result["stability"] = verdict_stability(result)
    return result


# --------------------------------------------------------------------------
# has the verdict settled?
# --------------------------------------------------------------------------
def verdict_stability(
    result: dict[str, Any] | Sequence[dict[str, Any]],
    key: str = "class",
    recipe: str | None = None,
) -> dict[str, Any]:
    """Has each compound's verdict settled by the largest affordable cut?

    A compound is ``settled`` only if its verdict is the same on the **two
    largest** cuts that were run and has not changed anywhere above
    ``settled_from``.  One rung of agreement is not stability; it is a
    coincidence with n = 1 replicate.

    Reports, per compound: the verdict sequence in increasing cut size, the
    number of changes, the last cut at which it changed, and -- when it never
    settles -- an explicit statement that the method's conclusion is
    scale-dependent with no stable point inside the range tested.  That finding
    is a legitimate result, and this function is written so that it cannot be
    quietly rounded into a reassuring one.

    ``recipe`` restricts the sequence to cuts whose ``recipe`` starts with that
    string.  Use ``recipe="scale_ladder"`` for the **clean** scale axis: the
    committed ``named`` and ``taste_motor`` cuts sit inside the ladder's size
    range but were built by a different rule, so a sequence that interleaves
    them mixes a size effect with a construction effect.
    """
    rows = list(result["rows"]) if isinstance(result, dict) else list(result)
    if recipe is not None:
        rows = [r for r in rows if str(r.get("recipe", "")).startswith(recipe)]
    by_compound: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        by_compound.setdefault(r["compound"], []).append(r)
    out: dict[str, Any] = {}
    any_unsettled = False
    for compound, rs in by_compound.items():
        rs = sorted(rs, key=lambda r: (r["n_nodes"], r["n_edges"]))
        seq = [r.get(key) for r in rs]
        cutnames = [r["cut"] for r in rs]
        changes = [
            {"from_cut": cutnames[i - 1], "to_cut": cutnames[i],
             "from": seq[i - 1], "to": seq[i]}
            for i in range(1, len(seq))
            if seq[i] != seq[i - 1]
        ]
        final = seq[-1] if seq else None
        # longest suffix equal to the final verdict
        i = len(seq)
        while i > 0 and seq[i - 1] == final:
            i -= 1
        settled_from = cutnames[i] if i < len(seq) else None
        n_agreeing = len(seq) - i
        settled = bool(len(seq) >= 2 and n_agreeing >= 2)
        under_powered = [
            r["cut"] for r in rs if r.get("resolution_coarser_than_alpha")
        ]
        if not settled:
            any_unsettled = True
        if len(seq) < 2:
            statement = (
                f"{compound}: only {len(seq)} cut was run, so stability cannot "
                f"be assessed at all; the verdict there is '{final}'."
            )
        elif settled:
            statement = (
                f"{compound}: the verdict is '{final}' on the largest "
                f"{n_agreeing} of {len(seq)} cuts tested, first reached at "
                f"'{settled_from}'."
            )
        else:
            statement = (
                f"{compound}: THE VERDICT HAS NOT SETTLED. It is still "
                f"changing at the largest affordable cut ('{cutnames[-2]}' -> "
                f"'{cutnames[-1]}': {seq[-2]} -> {final}). Within the range "
                "tested this method has no stable point, and its conclusions "
                "must be reported as scale-dependent."
            )
        out[compound] = {
            "compound": compound,
            "key": key,
            "cuts": cutnames,
            "n_nodes": [r["n_nodes"] for r in rs],
            "n_edges": [r["n_edges"] for r in rs],
            "n_permutations": [r["n"] for r in rs],
            "sequence": seq,
            "final": final,
            "n_changes": len(changes),
            "changes": changes,
            "settled": settled,
            "settled_from": settled_from if settled else None,
            "n_cuts_agreeing_at_top": n_agreeing,
            "under_powered_cuts": under_powered,
            "statement": statement,
        }
    return {
        "by_compound": out,
        "all_settled": bool(out) and not any_unsettled,
        "key": key,
        "recipe_filter": recipe,
        "criterion": (
            "settled = the same verdict on the two largest cuts run, with no "
            "change above that point. One rung of agreement is not stability."
        ),
        "caveat": (
            "Stability inside the range tested is not stability. The largest "
            "rung here is far short of the whole CNS, and the permutation "
            "count falls with scale, so 'settled' means 'did not move over the "
            "cuts this study could afford'."
        ),
    }


# --------------------------------------------------------------------------
# what does the verdict actually track?
# --------------------------------------------------------------------------
#: Structural statistics of a cut that a verdict might plausibly track, as
#: :func:`~flylab.analysis.dependence.cut_census` names them.  ``n_nodes`` and
#: ``n_edges`` are included precisely so that "it is just size" stays a
#: falsifiable alternative rather than an unexamined assumption.
STRUCTURE_KEYS: tuple[str, ...] = (
    "n_nodes",
    "n_edges",
    "mean_degree",
    "share_onto_seeds",
    "share_nodes_with_in_degree",
    "share_edges_from_nodes_with_input",
)


def verdict_vs_structure(
    result: dict[str, Any] | Sequence[dict[str, Any]],
    key: str = "class",
    statistics: Sequence[str] = STRUCTURE_KEYS,
) -> dict[str, Any]:
    """Which structural statistic, if any, separates the verdicts?

    For every statistic in ``statistics`` this reports the range taken by each
    verdict class and whether those ranges are **disjoint**, i.e. whether a
    single threshold on that statistic reproduces the verdict over the cuts
    that were run.  It also cross-tabulates the boolean ``in_star`` flag that
    :func:`~flylab.analysis.dependence.cut_census` already computes.

    This is an association over a handful of cuts, not a causal claim, and the
    returned ``strength`` says so: with a minority class of one or two cells
    every statistic that happens to be extreme on those cuts will separate,
    and the honest way to break the tie is a controlled cut -- the same node
    set with its edge structure changed -- not more statistics.
    """
    rows = list(result["rows"]) if isinstance(result, dict) else list(result)
    classes: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        classes.setdefault(str(r.get(key)), []).append(r)
    minority = min((len(v) for v in classes.values()), default=0)
    out: dict[str, Any] = {}
    for stat in statistics:
        per_class = {
            k: [r[stat] for r in v if r.get(stat) is not None]
            for k, v in classes.items()
        }
        per_class = {k: v for k, v in per_class.items() if v}
        if len(per_class) < 2:
            out[stat] = {"separates": None, "by_class": per_class}
            continue
        lo = {k: min(v) for k, v in per_class.items()}
        hi = {k: max(v) for k, v in per_class.items()}
        order = sorted(per_class, key=lambda k: lo[k])
        disjoint = all(
            hi[order[i]] < lo[order[i + 1]] for i in range(len(order) - 1)
        )
        out[stat] = {
            "separates": bool(disjoint),
            "by_class": {k: {"min": lo[k], "max": hi[k], "n": len(per_class[k])} for k in per_class},
            "threshold": (
                (hi[order[0]] + lo[order[1]]) / 2.0 if disjoint and len(order) == 2 else None
            ),
            "ascending_classes": order,
        }
    in_star = {}
    for k, v in classes.items():
        flags = [r.get("in_star") for r in v if r.get("in_star") is not None]
        in_star[k] = {"n": len(flags), "n_in_star": sum(1 for f in flags if f)}
    separating = [s for s, b in out.items() if b.get("separates")]
    if len(classes) < 2:
        statement = (
            f"every cut gave the same verdict ({next(iter(classes), None)!r}), so "
            "nothing here separates: the question of what the verdict tracks "
            "does not arise over these cuts."
        )
    elif minority <= 2:
        statement = (
            f"the minority verdict class has only {minority} cell(s), so the "
            + (
                "statistics " + ", ".join(separating)
                if separating
                else "statistics tested"
            )
            + " separate it trivially and cannot be told apart from one another. "
            "Read this as 'the cuts that differ are extreme on all of these at "
            "once', and settle it with a controlled cut (the same node set, a "
            "different edge structure) rather than with more statistics."
        )
    else:
        statement = (
            ("separated by " + ", ".join(separating))
            if separating
            else "no single statistic tested separates the verdict classes"
        ) + f" over {len(rows)} cut x compound cells."
    return {
        "key": key,
        "n_rows": len(rows),
        "classes": {k: len(v) for k, v in classes.items()},
        "minority_class_size": minority,
        "by_statistic": out,
        "separating_statistics": separating,
        "in_star_by_class": in_star,
        "strength": (
            "trivial (minority class <= 2)" if minority <= 2 and len(classes) > 1
            else "degenerate (one class)" if len(classes) < 2
            else "observational"
        ),
        "statement": statement,
        "label": "model_derived",
    }


# --------------------------------------------------------------------------
# the paper-ready table
# --------------------------------------------------------------------------
#: Columns of :func:`scale_report`, in order.
REPORT_COLUMNS: tuple[str, ...] = (
    "cut",
    "n_nodes",
    "n_edges",
    "mean_degree",
    "share_onto_seeds",
    "in_star",
    "compound",
    "n",
    "p_resolution",
    "real_effect",
    "relative_effect",
    "p_ER",
    "p_degree",
    "p_weight",
    "p_sign",
    "p_sign_wm",
    "necessary_information_level",
    "class",
)


def _fmt(v: Any) -> str:
    if v is None:
        return "n/a"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e5):
            return f"{v:.2e}"
        return f"{v:.4g}"
    return str(v)


def scale_report(
    result: dict[str, Any],
    columns: Sequence[str] = REPORT_COLUMNS,
) -> dict[str, Any]:
    """A paper-ready table of verdict versus scale, plus the stability verdict.

    Returns ``rows`` (dicts), ``columns``, ``markdown`` (a GitHub table) and
    ``statements`` -- one plain sentence per compound saying whether its
    verdict settled and where.
    """
    rows = []
    for r in sorted(result["rows"], key=lambda r: (r["n_nodes"], r["n_edges"], r["compound"])):
        flat = dict(r)
        flat.update(r.get("p") or {})
        rows.append({c: flat.get(c) for c in columns})
    head = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = ["| " + " | ".join(_fmt(row[c]) for c in columns) + " |" for row in rows]
    stability = result.get("stability") or verdict_stability(result)
    statements = [b["statement"] for b in stability["by_compound"].values()]
    return {
        "columns": list(columns),
        "rows": rows,
        "markdown": "\n".join([head, sep] + body),
        "statements": statements,
        "all_settled": stability["all_settled"],
        "stability": stability,
        "warnings": list(result.get("warnings") or []),
        "label": "model_derived",
    }


# --------------------------------------------------------------------------
# what is affordable at each scale
# --------------------------------------------------------------------------
#: Wall-clock beyond which an analysis is called infeasible on one core.
#: A day is generous for a single result and hopeless for a study that must be
#: repeated over 21 compounds and 4 concentrations.
FEASIBLE_SECONDS = 24 * 3600.0

#: LIF builds a dense ``N x N`` float32 matrix
#: (:class:`flylab.circuit.lif.LIFNetwork`), even when its sparse kernel is
#: selected, so its ceiling is memory, not time.
LIF_DENSE_BYTES_PER_NODE2 = 4


def feasibility_frontier(
    scales: Sequence[tuple[str, int, int]] | None = None,
    n_paper: int = PAPER_N,
    n_landscape: int = 100,
    landscape_cells: int = 84,
    memory_gb: float = 16.0,
    feasible_seconds: float = FEASIBLE_SECONDS,
) -> dict[str, Any]:
    """What is affordable at each scale: single run, profile, landscape, LIF.

    ``scales`` is ``[(name, n_nodes, n_edges), ...]``; the default is the two
    committed cuts, the scale ladder and the whole traced CNS at two synapse
    floors.  Costs come from the same measured constants the estimators use.
    The row where ``dependence_profile_feasible`` turns False is the point
    beyond which dependence testing stops being possible.
    """
    if scales is None:
        scales = DEFAULT_FRONTIER_SCALES
    rows = []
    for name, n_nodes, n_edges in scales:
        run_s = run_cost_s(n_edges)
        prof_s = profile_cost_s(n_edges, n_paper)
        land_s = profile_cost_s(n_edges, n_landscape) * landscape_cells
        lif_gb = LIF_DENSE_BYTES_PER_NODE2 * float(n_nodes) ** 2 / 1e9
        rows.append(
            {
                "scale": name,
                "n_nodes": int(n_nodes),
                "n_edges": int(n_edges),
                "single_run_s": run_s,
                "single_run_feasible": bool(run_s <= feasible_seconds),
                "dependence_profile_s": prof_s,
                "dependence_profile_h": prof_s / 3600.0,
                "dependence_profile_feasible": bool(prof_s <= feasible_seconds),
                "landscape_s": land_s,
                "landscape_days": land_s / 86400.0,
                "landscape_feasible": bool(land_s <= feasible_seconds),
                "lif_dense_gb": lif_gb,
                "lif_feasible": bool(lif_gb <= float(memory_gb)),
            }
        )
    feasible_profile = [r for r in rows if r["dependence_profile_feasible"]]
    feasible_lif = [r for r in rows if r["lif_feasible"]]
    return {
        "rows": rows,
        "n_paper": int(n_paper),
        "n_landscape": int(n_landscape),
        "landscape_cells": int(landscape_cells),
        "feasible_seconds": float(feasible_seconds),
        "memory_gb": float(memory_gb),
        "largest_cut_with_feasible_profile": (
            feasible_profile[-1]["scale"] if feasible_profile else None
        ),
        "largest_cut_with_feasible_lif": (
            feasible_lif[-1]["scale"] if feasible_lif else None
        ),
        "constants": {
            "profile_cost_a_s": PROFILE_COST_A,
            "profile_cost_b_s_per_edge": PROFILE_COST_B,
            "run_cost_s_per_edge_step": RUN_COST_PER_EDGE_STEP,
            "lif_bytes_per_node_squared": LIF_DENSE_BYTES_PER_NODE2,
        },
        "future_work": FUTURE_WORK,
        "label": "model_derived",
        "warnings": [
            "One core. Divide the profile and landscape columns by n_jobs for "
            "the parallel path, but not by more than ~3 on four cores: the "
            "weight-matched transmitter null runs serially.",
            "The LIF column is a memory bound, not a time bound: "
            "flylab.circuit.lif.LIFNetwork materialises a dense N x N float32 "
            "matrix even when its sparse kernel is chosen, so it fails on "
            "allocation rather than running slowly.",
            "Feasible here means 'one result in under a day on one core'. A "
            "landscape that is 'feasible' by that rule is still a day of "
            "compute per specification, and the robustness family has 25.",
        ],
    }


#: Default frontier rows: the committed cuts, the scale ladder, and the whole
#: traced CNS at two synapse floors.  Every node and edge count here is
#: **measured** on MaleCNS v1.0: 6 235 682 traced-to-traced edges at weight 5
#: and 25 563 197 at weight 1 (the file itself has 151 856 684 rows, but two
#: thirds of them touch a body that is not traced).
DEFAULT_FRONTIER_SCALES: tuple[tuple[str, int, int], ...] = (
    ("named", 1126, 1360),
    ("taste_motor", 1841, 19066),
    ("scale_1k", 1000, 22857),
    ("scale_5k", 5000, 279845),
    ("scale_10k", 10000, 621600),
    ("scale_25k", 25000, 1364375),
    ("scale_50k", 50000, 2216881),
    ("whole_cns_w5", 165122, 6235682),
    ("whole_cns_w1", 165122, 25563197),
)

FUTURE_WORK: list[dict[str, Any]] = [
    {
        "idea": "reuse one decomposition across the null ensemble",
        "what": (
            "Every shuffle currently rebuilds the signed weights, the row "
            "normaliser and the 80-step iteration from scratch. The weight-"
            "permuting and transmitter-permuting nulls change only a per-edge "
            "or per-node scalar, not the sparsity pattern, so the fixed point "
            "could be tracked from the unshuffled solution with a few "
            "Neumann/Krylov corrections instead of 80 full sweeps."
        ),
        "estimated_speedup": "5-15x on the weight and transmitter nulls",
        "would_move_frontier_to": "whole CNS at weight 5 for a 1000-permutation profile",
        "status": "not implemented",
    },
    {
        "idea": "analytic permutation distribution for the rate readout",
        "what": (
            "With row_abs normalisation each cell's input is a composition-"
            "weighted mean of presynaptic gains, and the network mean is close "
            "to a smooth function of the gain vector and the cut's weighted "
            "transmitter shares. Under the transmitter nulls those shares are "
            "a hypergeometric-type statistic whose first two moments are "
            "available in closed form, so the null distribution of the effect "
            "could be approximated without permuting at all, and the "
            "permutation test kept only as a spot check."
        ),
        "estimated_speedup": "O(1) instead of O(n) shuffles for the rank-3 rung",
        "would_move_frontier_to": "unbounded for the transmitter nulls; the "
        "wiring nulls still need permutation",
        "status": "not implemented; would need validating against the exact "
        "test on the committed cuts before any claim rests on it",
    },
    {
        "idea": "subsample the shuffle, not the graph",
        "what": (
            "A degree-preserving rewire currently does O(E) double edge swaps "
            "per shuffle. A partial rewire (a fixed fraction of swaps) is a "
            "weaker null but a continuous one, and the dependence verdict "
            "could be reported as a function of rewire depth rather than at "
            "one fully-rewired point."
        ),
        "estimated_speedup": "linear in the swap fraction",
        "would_move_frontier_to": "50k-100k cells at n = 1000",
        "status": "not implemented; changes the null's meaning, so it is a "
        "different question rather than a cheaper answer to this one",
    },
]


__all__ = [
    "CUT_FILES",
    "CUT_RECIPE",
    "DEFAULT_CUTS",
    "DEFAULT_FRONTIER_SCALES",
    "FEASIBLE_SECONDS",
    "FUTURE_WORK",
    "MIN_USEFUL_N",
    "PROFILE_COST_A",
    "PROFILE_COST_B",
    "REPORT_COLUMNS",
    "RUN_COST_PER_EDGE_STEP",
    "MEASURED_RUN_COSTS",
    "run_cost_s",
    "SCALE_WARNINGS",
    "STRUCTURE_KEYS",
    "verdict_vs_structure",
    "available_cuts",
    "cut_path",
    "dependence_vs_scale",
    "estimate_scale_runtime",
    "feasibility_frontier",
    "missing_cuts",
    "permutation_budget",
    "profile_cost_s",
    "scale_report",
    "verdict_stability",
]
