"""The scaling study: at what network scale does a dependence verdict settle?

What is being pinned down here:

* the cut **ladder** is deterministic, nested and built by one documented rule
  (:mod:`flylab.maps.extract`), so a verdict difference between two rungs is a
  statement about size and not about the recipe;
* the cost model that decides how many permutations each rung gets is the one
  the estimators publish, the permutation floor is respected, and a rung that
  cannot resolve ``alpha`` says so;
* :func:`verdict_stability` calls an unsettled verdict unsettled -- the
  uncomfortable answer must survive the reporting layer;
* the feasibility frontier names the point beyond which dependence testing
  stops being possible.

Everything here runs without the 1.0 GB MaleCNS weight matrix.  The rungs above
1k are CI artifacts, so the study itself is exercised on the committed cuts and
on a synthetic cut; the full-ladder run is marked ``slow``.
"""
from __future__ import annotations

import builtins
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from flylab.analysis.dependence import synthetic_cut
from flylab.analysis.scale import (
    CUT_FILES,
    STRUCTURE_KEYS,
    verdict_vs_structure,
    MEASURED_RUN_COSTS,
    run_cost_s,
    CUT_RECIPE,
    DEFAULT_CUTS,
    DEFAULT_FRONTIER_SCALES,
    FUTURE_WORK,
    MIN_USEFUL_N,
    PROFILE_COST_A,
    PROFILE_COST_B,
    REPORT_COLUMNS,
    available_cuts,
    cut_path,
    dependence_vs_scale,
    estimate_scale_runtime,
    feasibility_frontier,
    missing_cuts,
    permutation_budget,
    profile_cost_s,
    scale_report,
    verdict_stability,
)
from flylab.maps.extract import (
    COMMIT_BYTE_BUDGET,
    LADDER_MIN_WEIGHT,
    LADDER_RECIPE,
    LADDER_SEED_TYPES,
    LADDER_SIZES,
    _composition,
    _growth_order,
    ladder_filename,
    ladder_summary,
)

ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------
# the ladder rule: deterministic, nested, documented
# --------------------------------------------------------------------------
def _toy_graph(seed: int = 0):
    """A small weighted digraph with a dense core and a sparse rim."""
    rng = np.random.default_rng(seed)
    bodies = np.arange(100, 160, dtype=np.int64)
    pre = rng.choice(bodies, 400)
    post = rng.choice(bodies, 400)
    keep = pre != post
    w = rng.integers(5, 50, 400).astype(float)
    return pre[keep], post[keep], w[keep], bodies


def test_growth_order_is_deterministic():
    pre, post, w, bodies = _toy_graph()
    seeds = bodies[:3]
    a, ranks_a = _growth_order(pre, post, w, bodies, seeds, 30)
    b, ranks_b = _growth_order(pre, post, w, bodies, seeds, 30)
    assert list(a) == list(b)
    assert ranks_a == ranks_b
    assert sum(ranks_a) == len(a) == 30


def test_the_ladder_is_nested():
    """A smaller rung's node set is a prefix of a larger rung's."""
    pre, post, w, bodies = _toy_graph()
    seeds = bodies[:3]
    small, _ = _growth_order(pre, post, w, bodies, seeds, 12)
    large, _ = _growth_order(pre, post, w, bodies, seeds, 40)
    assert list(large[: len(small)]) == list(small)
    assert set(small) <= set(large)


def test_growth_starts_from_the_seeds_and_ranks_by_connection_strength():
    pre, post, w, bodies = _toy_graph()
    seeds = bodies[:3]
    order, ranks = _growth_order(pre, post, w, bodies, seeds, 20)
    assert ranks[0] == 3
    assert sorted(bodies[order[:3]]) == sorted(seeds.tolist())
    # the first cell taken after the seeds is the strongest-connected one
    sel = np.isin(pre, seeds) | np.isin(post, seeds)
    score = {}
    for a, b, v in zip(pre[sel], post[sel], w[sel]):
        for node in (a, b):
            if node not in seeds:
                score[int(node)] = score.get(int(node), 0.0) + float(v)
    best = max(sorted(score), key=lambda k: (score[k], -k))
    assert int(bodies[order[3]]) == best


def test_growth_stops_when_the_component_is_exhausted():
    pre = np.array([10, 11], dtype=np.int64)
    post = np.array([11, 12], dtype=np.int64)
    w = np.array([9.0, 9.0])
    bodies = np.array([10, 11, 12, 99], dtype=np.int64)  # 99 is isolated
    order, _ = _growth_order(pre, post, w, bodies, np.array([10], dtype=np.int64), 99)
    assert sorted(bodies[order].tolist()) == [10, 11, 12]


def test_the_recipe_is_recorded_and_says_why_not_the_alternatives():
    for phrase in ("More hops", "lower weight floor", "whole neuropil", "Determinism"):
        assert phrase.lower() in LADDER_RECIPE.lower(), phrase
    assert "induced" in LADDER_RECIPE


def test_ladder_constants_are_coherent():
    assert LADDER_SIZES == (1000, 5000, 10000, 25000, 50000)
    assert "MN9" in LADDER_SEED_TYPES and "DNp01" in LADDER_SEED_TYPES
    assert LADDER_MIN_WEIGHT == 5
    assert ladder_filename(5000) == "malecns_scale_5k.json"
    assert ladder_filename(1000).endswith("_1k.json")


def test_composition_reports_out_weight_not_just_cell_counts():
    """Weighted shares, because the gain patch acts through outgoing weight."""
    nodes = [
        {"bodyId": 1, "consensus_nt": "acetylcholine"},
        {"bodyId": 2, "consensus_nt": "gaba"},
        {"bodyId": 3, "consensus_nt": None},
    ]
    edges = [
        {"pre": 1, "post": 2, "weight": 90.0},
        {"pre": 2, "post": 3, "weight": 10.0},
    ]
    comp = _composition(nodes, edges)
    assert comp["cell_counts"]["acetylcholine"] == 1
    assert comp["cell_counts"]["unclear"] == 1
    assert comp["out_weight_shares"]["acetylcholine"] == pytest.approx(0.9)
    assert comp["out_weight_shares"]["gaba"] == pytest.approx(0.1)
    assert comp["total_out_weight"] == pytest.approx(100.0)


def test_ladder_summary_drops_the_bulk():
    payload = {"n_nodes": 3, "nodes": [1, 2, 3], "edges": [], "recipe": "x", "bytes": 9}
    s = ladder_summary(payload)
    assert "nodes" not in s and "edges" not in s and "recipe" not in s
    assert s["n_nodes"] == 3


# --------------------------------------------------------------------------
# which cuts exist, and what is committed
# --------------------------------------------------------------------------
def test_the_committed_ladder_rung_is_small_enough_to_be_committed():
    p = cut_path("scale_1k")
    if p is None:
        pytest.skip("scale_1k not built in this checkout")
    assert p.stat().st_size <= COMMIT_BYTE_BUDGET
    g = json.loads(p.read_text())
    assert g["growth"] == "ranked_bfs_induced"
    assert g["min_weight"] == LADDER_MIN_WEIGHT
    assert g["n_nodes"] == 1000
    assert g["census"]["mean_degree"] > 1.0
    assert g["composition"]["out_weight_shares"]["acetylcholine"] > 0.3
    assert g["seeds"]["MN9"] and g["seeds"]["DNp01"]


def test_the_default_ladder_covers_the_two_committed_cuts_and_the_rungs():
    assert DEFAULT_CUTS[:2] == ("named", "taste_motor")
    assert set(CUT_FILES) >= {"named", "taste_motor", "scale_1k", "scale_50k"}
    for name in DEFAULT_CUTS:
        assert CUT_RECIPE[name]
    # the two families are labelled differently, so the confound is visible
    assert CUT_RECIPE["named"] != CUT_RECIPE["scale_1k"]


def test_available_cuts_is_sorted_by_size_and_missing_ones_are_named():
    cuts = available_cuts()
    assert cuts, "no cut on disk at all"
    assert [c["n_nodes"] for c in cuts] == sorted(c["n_nodes"] for c in cuts)
    assert set(missing_cuts()) == set(DEFAULT_CUTS) - {c["name"] for c in cuts}


def test_an_explicit_path_is_accepted_and_a_missing_one_is_none(tmp_path):
    assert cut_path(tmp_path / "nope.json") is None
    p = tmp_path / "cut.json"
    p.write_text("{}")
    assert cut_path(p) == p


# --------------------------------------------------------------------------
# cost: be honest before spending the afternoon
# --------------------------------------------------------------------------
def test_profile_cost_is_the_published_linear_model():
    assert profile_cost_s(1000, 10, n_modes=2) == pytest.approx(
        2 * 10 * (PROFILE_COST_A + PROFILE_COST_B * 1000)
    )
    assert profile_cost_s(10_000, 100) > profile_cost_s(1_000, 100)


def test_permutation_budget_scales_n_down_with_the_cut():
    small = permutation_budget(1_360, seconds=600)
    big = permutation_budget(2_216_881, seconds=600)
    assert small["n"] > big["n"]
    assert big["n"] >= MIN_USEFUL_N


def test_permutation_budget_never_goes_below_the_resolution_floor():
    b = permutation_budget(50_000_000, seconds=1.0)
    assert b["n"] == MIN_USEFUL_N
    assert b["bound_by"] == "floor"
    assert any("raised to the floor" in w for w in b["warnings"])
    assert b["predicted_s"] > b["seconds_budget"]


def test_permutation_budget_caps_at_the_paper_n():
    b = permutation_budget(10, seconds=1e9)
    assert b["n"] == 1000
    assert b["bound_by"] == "ceiling"
    assert b["p_resolution"] == pytest.approx(1 / 1001)


def test_a_coarse_resolution_is_flagged_as_unable_to_reject():
    b = permutation_budget(1000, seconds=600, min_n=5, alpha=0.05)
    if b["n"] < 19:
        assert b["resolution_coarser_than_alpha"]
        assert any("no mode can be distinguished" in w for w in b["warnings"])


def test_run_cost_reproduces_its_own_measurements():
    """The run-cost curve interpolates measurements; it must hit them."""
    for n_edges, seconds in MEASURED_RUN_COSTS:
        assert run_cost_s(n_edges) == pytest.approx(seconds, rel=1e-6)
    assert run_cost_s(1_000_000) > run_cost_s(100_000)
    assert run_cost_s(19_066, steps=160) == pytest.approx(2 * run_cost_s(19_066))


def test_run_cost_per_edge_is_not_constant_with_size():
    """Cache behaviour, not arithmetic: a single constant would mislead."""
    small = run_cost_s(279_845) / 279_845
    large = run_cost_s(25_563_197) / 25_563_197
    assert large > 3 * small


def test_the_measured_whole_cns_run_is_in_the_table():
    edges = dict(MEASURED_RUN_COSTS)
    assert 25_563_197 in edges
    assert 25 < edges[25_563_197] < 40  # measured 30.3 s for 80 steps


def test_missing_cuts_are_reported_in_the_runtime_estimate():
    est = estimate_scale_runtime(["imidacloprid"], cuts=("named", "does_not_exist.json"))
    assert "does_not_exist.json" in est["missing_cuts"]
    assert all(r["name"] == "named" for r in est["rows"])


def test_estimate_scale_runtime_totals_over_compounds_and_cuts():
    est = estimate_scale_runtime(["imidacloprid", "fipronil"], seconds_per_cut=600)
    assert est["rows"]
    assert est["estimate_s_serial"] == pytest.approx(sum(r["predicted_s"] for r in est["rows"]))
    assert est["estimate_s"] == est["estimate_s_serial"]
    par = estimate_scale_runtime(["imidacloprid"], seconds_per_cut=600, n_jobs=4)
    assert par["estimate_s"] < par["estimate_s_serial"]


# --------------------------------------------------------------------------
# the study itself
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def tiny_cuts(tmp_path_factory):
    """Two synthetic cuts of different size, written where cut_path finds them."""
    d = tmp_path_factory.mktemp("cuts")
    out = []
    for n_nodes, n_edges in ((120, 400), (260, 1600)):
        g = synthetic_cut(n_nodes=n_nodes, n_edges=n_edges, loop_strength=1.0, seed=3)
        g.update(
            {
                "map": "male-cns:v1.0",
                "citation": "synthetic",
                "n_nodes": len(g["nodes"]),
                "n_edges": len(g["edges"]),
            }
        )
        p = d / f"synthetic_{n_nodes}.json"
        p.write_text(json.dumps(g))
        out.append(p)
    return out


def test_dependence_vs_scale_reports_one_row_per_cut_and_compound(tiny_cuts):
    res = dependence_vs_scale(
        compounds=("imidacloprid",), cuts=tiny_cuts, n=MIN_USEFUL_N, conc_M=1e-6
    )
    assert len(res["rows"]) == len(tiny_cuts)
    assert [r["n_nodes"] for r in res["rows"]] == sorted(r["n_nodes"] for r in res["rows"])
    for row in res["rows"]:
        assert row["class"] is not None
        assert row["n"] == MIN_USEFUL_N
        # the p per mode, the equivalence gap per mode, the level, the structure
        assert set(row["p"])
        assert row["necessary_information_level"]
        assert row["mean_degree"] is not None
        assert row["share_onto_seeds"] is not None
        assert row["in_star"] in (True, False)
        assert any(m["equivalence_gap"] is not None for m in row["modes"])
    json.dumps(res)  # serialisable


def test_the_study_warns_that_two_recipes_are_not_one_scale_axis(tiny_cuts):
    res = dependence_vs_scale(compounds=("imidacloprid",), cuts=tiny_cuts, n=MIN_USEFUL_N)
    assert any("confounds scale with construction" in w for w in res["warnings"])
    assert any("dilutes as the cut grows" in w for w in res["warnings"])
    assert any("not equally powered" in w for w in res["warnings"])


def test_a_ladder_with_nothing_on_disk_is_an_explicit_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="extract_ladder"):
        dependence_vs_scale(compounds=("imidacloprid",), cuts=[tmp_path / "absent.json"])


# --------------------------------------------------------------------------
# has it settled?  the uncomfortable answer must survive
# --------------------------------------------------------------------------
def _rows(classes):
    return [
        {
            "compound": "x",
            "cut": f"c{i}",
            "n_nodes": 100 * (i + 1),
            "n_edges": 1000 * (i + 1),
            "n": 100,
            "class": c,
        }
        for i, c in enumerate(classes)
    ]


def test_a_verdict_that_agrees_on_the_top_two_cuts_is_settled():
    s = verdict_stability({"rows": _rows(["composition-dominated"] * 2 + ["topology-dependent"] * 3)})
    b = s["by_compound"]["x"]
    assert b["settled"] is True
    assert b["settled_from"] == "c2"
    assert b["n_changes"] == 1
    assert s["all_settled"] is True


def test_a_verdict_still_moving_at_the_top_is_reported_as_unsettled():
    s = verdict_stability({"rows": _rows(["composition-dominated", "mixed", "topology-dependent"])})
    b = s["by_compound"]["x"]
    assert b["settled"] is False
    assert b["settled_from"] is None
    assert "HAS NOT SETTLED" in b["statement"]
    assert "scale-dependent" in b["statement"]
    assert s["all_settled"] is False


def test_one_rung_of_agreement_is_not_stability():
    """The last two cuts must agree; a single final rung proves nothing."""
    s = verdict_stability({"rows": _rows(["a", "a", "a", "b"])})
    assert s["by_compound"]["x"]["settled"] is False
    assert s["by_compound"]["x"]["n_cuts_agreeing_at_top"] == 1


def test_a_single_cut_cannot_be_called_stable():
    s = verdict_stability({"rows": _rows(["composition-dominated"])})
    b = s["by_compound"]["x"]
    assert b["settled"] is False
    assert "stability cannot be assessed" in b["statement"]
    assert s["all_settled"] is False


def test_the_sequence_can_be_restricted_to_one_recipe():
    """named/taste_motor sit inside the ladder's size range but are built by a
    different rule, so the clean scale axis is the scale_ladder rows alone."""
    rows = _rows(["a", "b", "a"])
    for r, rec in zip(rows, ["scale_ladder (x)", "hops_neighborhood (y)", "scale_ladder (x)"]):
        r["recipe"] = rec
    mixed = verdict_stability({"rows": rows})["by_compound"]["x"]
    clean = verdict_stability({"rows": rows}, recipe="scale_ladder")["by_compound"]["x"]
    assert mixed["sequence"] == ["a", "b", "a"] and mixed["n_changes"] == 2
    assert clean["sequence"] == ["a", "a"] and clean["settled"] is True
    assert clean["cuts"] == ["c0", "c2"]


def test_stability_carries_the_range_caveat():
    s = verdict_stability({"rows": _rows(["a", "a"])})
    assert "Stability inside the range tested is not stability" in s["caveat"]
    assert "two largest" in s["criterion"]


def test_under_powered_rungs_are_named():
    rows = _rows(["a", "a"])
    rows[1]["resolution_coarser_than_alpha"] = True
    s = verdict_stability({"rows": rows})
    assert s["by_compound"]["x"]["under_powered_cuts"] == ["c1"]


# --------------------------------------------------------------------------
# what does the verdict track?
# --------------------------------------------------------------------------
def _structured(specs):
    """rows as (class, mean_degree, in_star) triples."""
    out = []
    for i, (cls, deg, star) in enumerate(specs):
        out.append(
            {
                "compound": "x",
                "cut": f"c{i}",
                "n_nodes": 1000 + i,
                "n_edges": 1000 * (i + 1),
                "n": 100,
                "class": cls,
                "mean_degree": deg,
                "share_onto_seeds": 1.0 / deg,
                "share_nodes_with_in_degree": min(1.0, deg / 25.0),
                "share_edges_from_nodes_with_input": min(1.0, deg / 25.0),
                "in_star": star,
            }
        )
    return out


def test_a_statistic_that_separates_the_verdicts_is_named():
    rows = _structured(
        [("composition-dominated", 1.2, True), ("composition-dominated", 1.5, True),
         ("composition-dominated", 2.0, True), ("topology-dependent", 20.0, False),
         ("topology-dependent", 40.0, False), ("topology-dependent", 60.0, False)]
    )
    v = verdict_vs_structure({"rows": rows})
    assert "mean_degree" in v["separating_statistics"]
    assert v["by_statistic"]["mean_degree"]["separates"] is True
    assert 2.0 < v["by_statistic"]["mean_degree"]["threshold"] < 20.0
    assert v["strength"] == "observational"
    assert "separated by" in v["statement"]


def test_size_is_kept_as_a_falsifiable_alternative():
    """'It is just size' must be testable, so n_nodes is one of the statistics."""
    assert "n_nodes" in STRUCTURE_KEYS and "n_edges" in STRUCTURE_KEYS
    rows = _structured(
        [("composition-dominated", 1.2, True), ("composition-dominated", 1.5, True),
         ("composition-dominated", 2.0, True), ("topology-dependent", 20.0, False),
         ("topology-dependent", 40.0, False), ("topology-dependent", 60.0, False)]
    )
    v = verdict_vs_structure({"rows": rows})
    # n_nodes increases with the row index here, so it separates too, and the
    # report must not pretend one of the two explanations has been ruled out
    assert v["by_statistic"]["n_nodes"]["separates"] is True


def test_a_minority_class_of_one_is_called_trivial_not_a_finding():
    rows = _structured(
        [("composition-dominated", 1.2, True)] + [("topology-dependent", 20.0 + i, False) for i in range(5)]
    )
    v = verdict_vs_structure({"rows": rows})
    assert v["minority_class_size"] == 1
    assert v["strength"] == "trivial (minority class <= 2)"
    assert "cannot be told apart" in v["statement"]
    assert "controlled cut" in v["statement"]


def test_overlapping_ranges_do_not_separate():
    """An in-star cut that gives BOTH verdicts means the flag is not sufficient."""
    rows = _structured(
        [("composition-dominated", 1.2, True), ("topology-dependent", 1.2, True),
         ("topology-dependent", 20.0, False), ("topology-dependent", 40.0, False)]
    )
    v = verdict_vs_structure({"rows": rows})
    assert "mean_degree" not in v["separating_statistics"]
    assert v["by_statistic"]["mean_degree"]["separates"] is False
    assert v["by_statistic"]["mean_degree"]["threshold"] is None
    assert v["by_statistic"]["share_onto_seeds"]["separates"] is False
    # the in-star flag is necessary but not sufficient: one in-star cut lands
    # in each class, so the flag alone does not predict the verdict
    assert v["in_star_by_class"]["topology-dependent"]["n_in_star"] == 1
    assert v["in_star_by_class"]["composition-dominated"]["n_in_star"] == 1


def test_one_class_everywhere_is_reported_as_degenerate():
    v = verdict_vs_structure({"rows": _structured([("topology-dependent", 20.0, False)] * 3)})
    assert v["strength"] == "degenerate (one class)"
    assert "does not arise" in v["statement"]


# --------------------------------------------------------------------------
# the paper-ready table
# --------------------------------------------------------------------------
def test_scale_report_is_a_table_a_paper_can_paste(tiny_cuts):
    res = dependence_vs_scale(compounds=("imidacloprid",), cuts=tiny_cuts, n=MIN_USEFUL_N)
    rep = scale_report(res)
    assert rep["columns"] == list(REPORT_COLUMNS)
    lines = rep["markdown"].splitlines()
    assert lines[0].startswith("| cut |") and lines[0].rstrip().endswith("class |")
    assert len(lines) == 2 + len(rep["rows"])
    assert all(line.count("|") == len(REPORT_COLUMNS) + 1 for line in lines)
    assert rep["statements"] and all(isinstance(s, str) for s in rep["statements"])
    assert "class" in rep["rows"][0]


def test_the_report_never_hides_an_unsettled_verdict():
    res = {"rows": _rows(["composition-dominated", "topology-dependent"]), "warnings": []}
    for r in res["rows"]:
        r.update({"p": {}, "mean_degree": 1.0, "in_star": False, "compound": "x"})
    rep = scale_report(res)
    assert rep["all_settled"] is False
    assert any("HAS NOT SETTLED" in s for s in rep["statements"])


# --------------------------------------------------------------------------
# the feasibility frontier
# --------------------------------------------------------------------------
def test_the_frontier_names_where_dependence_testing_stops():
    f = feasibility_frontier()
    names = [r["scale"] for r in f["rows"]]
    assert names[:2] == ["named", "taste_motor"]
    assert "whole_cns_w1" in names
    assert f["largest_cut_with_feasible_profile"] is not None
    # a single run of the whole CNS is affordable; a 1000-permutation profile is not
    whole = [r for r in f["rows"] if r["scale"] == "whole_cns_w1"][0]
    assert whole["single_run_feasible"] is True
    assert whole["dependence_profile_feasible"] is False
    assert whole["landscape_feasible"] is False
    assert whole["lif_feasible"] is False


def test_the_frontier_is_monotone_in_edges():
    f = feasibility_frontier()
    by_edges = sorted(f["rows"], key=lambda r: r["n_edges"])
    costs = [r["dependence_profile_s"] for r in by_edges]
    assert costs == sorted(costs)
    # once infeasible, always infeasible further up
    flags = [r["dependence_profile_feasible"] for r in by_edges]
    assert flags == sorted(flags, reverse=True)


def test_lif_is_bounded_by_memory_not_time():
    f = feasibility_frontier(memory_gb=16.0)
    whole = [r for r in f["rows"] if r["scale"].startswith("whole_cns")][0]
    assert whole["lif_dense_gb"] > 100
    assert any("memory bound, not a time bound" in w for w in f["warnings"])


def test_the_frontier_scales_match_the_measured_cuts():
    named = [s for s in DEFAULT_FRONTIER_SCALES if s[0] == "named"][0]
    assert named[1:] == (1126, 1360)
    p = cut_path("scale_1k")
    if p is not None:
        rung = [s for s in DEFAULT_FRONTIER_SCALES if s[0] == "scale_1k"][0]
        g = json.loads(p.read_text())
        assert rung[1] == g["n_nodes"] and rung[2] == g["n_edges"]


def test_future_work_is_described_with_an_estimate_and_not_implemented():
    assert len(FUTURE_WORK) >= 2
    for item in FUTURE_WORK:
        assert item["estimated_speedup"]
        assert item["would_move_frontier_to"]
        assert "not implemented" in item["status"]


# --------------------------------------------------------------------------
# the real ladder (needs the CI artifacts)
# --------------------------------------------------------------------------
@pytest.mark.slow
def test_the_real_ladder_runs_end_to_end():
    cuts = [c["name"] for c in available_cuts()]
    if len(cuts) < 2:
        pytest.skip("fewer than two cuts on disk; run the extract Action")
    res = dependence_vs_scale(
        compounds=("imidacloprid", "fipronil"), cuts=cuts, seconds_per_cut=120, n_jobs=2
    )
    rep = scale_report(res)
    assert len(res["rows"]) == 2 * len(cuts)
    assert rep["markdown"].count("\n") >= len(res["rows"])


# --------------------------------------------------------------------------
# the browser guarantee: the core imports without the heavy stack
# --------------------------------------------------------------------------
#: what Pyodide does not have; the same list tests/test_browser_bridge.py uses
HEAVY = ("pydantic", "fastapi", "starlette", "typer", "uvicorn", "pandas", "pyarrow")


@pytest.fixture
def without_heavy_stack(monkeypatch):
    real_import = builtins.__import__

    def guarded(name, *args, **kw):
        if name.split(".")[0] in HEAVY:
            raise ImportError(f"{name} is not available in this environment")
        return real_import(name, *args, **kw)

    for mod in list(sys.modules):
        if mod.split(".")[0] in HEAVY or mod.startswith("flylab"):
            monkeypatch.delitem(sys.modules, mod, raising=False)
    monkeypatch.setattr(builtins, "__import__", guarded)
    yield
    for mod in list(sys.modules):
        if mod.startswith("flylab"):
            sys.modules.pop(mod, None)


@pytest.mark.parametrize(
    "module",
    [
        "flylab.maps.ladder",
        "flylab.analysis.scale",
        "flylab.assays.fullcns",
        "flylab.analysis",
        "flylab.assays",
    ],
)
def test_the_scaling_layer_imports_with_pandas_and_pyarrow_hidden(
    without_heavy_stack, module
):
    """``flylab.browser.bridge`` reaches all of these, and Pyodide has no pandas.

    The ladder's *constants* therefore live in ``flylab.maps.ladder`` (stdlib
    only) and ``flylab.maps.extract``, which needs pandas and pyarrow to build
    a cut, is never imported at module scope from the analysis or assay layer.
    """
    importlib.import_module(module)
    assert "pandas" not in sys.modules and "pyarrow" not in sys.modules
    assert "flylab.maps.extract" not in sys.modules


def test_the_ladder_constants_module_has_no_heavy_imports():
    """A regression guard that does not need the import machinery to run."""
    import ast

    tree = ast.parse((ROOT / "flylab" / "maps" / "ladder.py").read_text())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert not (imported & {"pandas", "pyarrow", "numpy"}), sorted(imported)


def test_extract_still_re_exports_the_ladder_constants():
    """Moving the constants must not break `from flylab.maps.extract import ...`."""
    from flylab.maps import extract, ladder

    for name in ladder.__all__:
        assert getattr(extract, name) is getattr(ladder, name), name
