"""Connectome null models: shuffles, null distributions and the summary score."""
from __future__ import annotations

import json
import time
from collections import Counter

import numpy as np
import pytest

from flylab.analysis.nullmodels import (
    MODES,
    connectome_information_score,
    null_distribution,
    null_panel,
    shuffle_graph,
)
from flylab.circuit.rate import load_graph


@pytest.fixture(scope="module")
def named():
    return load_graph("named")


def _pairs(graph):
    return [(e["pre"], e["post"]) for e in graph["edges"]]


@pytest.mark.parametrize("mode", MODES)
def test_shuffle_preserves_size_nodes_and_seeds(named, mode):
    s = shuffle_graph(named, mode, np.random.default_rng(0))
    assert s["n_nodes"] == named["n_nodes"] == len(s["nodes"])
    assert s["n_edges"] == named["n_edges"] == len(s["edges"])
    # seeds stay seeds: same block, and every seed body is still a node
    assert s["seeds"] == named["seeds"]
    bodies = {n["bodyId"] for n in s["nodes"]}
    for ids in s["seeds"].values():
        assert set(ids) <= bodies
    assert [n["bodyId"] for n in s["nodes"]] == [n["bodyId"] for n in named["nodes"]]
    # no self loops, no duplicated edges
    assert all(e["pre"] != e["post"] for e in s["edges"])
    assert len(set(_pairs(s))) == len(s["edges"])
    assert s["null_model"]["mode"] == mode
    # the cached source graph is never mutated
    assert named["n_edges"] == len(named["edges"])


def test_shuffle_is_deterministic_given_the_seed(named):
    for mode in MODES:
        a = shuffle_graph(named, mode, 7)
        b = shuffle_graph(named, mode, 7)
        c = shuffle_graph(named, mode, 8)
        assert json.dumps(a["edges"]) == json.dumps(b["edges"])
        assert json.dumps([n["consensus_nt"] for n in a["nodes"]]) == json.dumps(
            [n["consensus_nt"] for n in b["nodes"]]
        )
        differs = json.dumps(a["edges"]) != json.dumps(c["edges"]) or json.dumps(
            [n["consensus_nt"] for n in a["nodes"]]
        ) != json.dumps([n["consensus_nt"] for n in c["nodes"]])
        assert differs, f"{mode} gave the same graph for two different seeds"


def test_sign_permute_keeps_the_transmitter_histogram(named):
    s = shuffle_graph(named, "sign_permute", 0)
    assert Counter(n["consensus_nt"] for n in s["nodes"]) == Counter(
        n["consensus_nt"] for n in named["nodes"]
    )
    assert _pairs(s) == _pairs(named)  # topology untouched
    assert [n["consensus_nt"] for n in s["nodes"]] != [
        n["consensus_nt"] for n in named["nodes"]
    ]


def test_weight_permute_keeps_the_weight_histogram_and_topology(named):
    s = shuffle_graph(named, "weight_permute", 0)
    assert sorted(e["weight"] for e in s["edges"]) == sorted(
        e["weight"] for e in named["edges"]
    )
    assert _pairs(s) == _pairs(named)


def test_rewire_preserves_in_and_out_degree(named):
    s = shuffle_graph(named, "rewire_degree_preserving", 0)
    assert Counter(e["pre"] for e in s["edges"]) == Counter(e["pre"] for e in named["edges"])
    assert Counter(e["post"] for e in s["edges"]) == Counter(e["post"] for e in named["edges"])
    assert set(_pairs(s)) != set(_pairs(named))
    assert s["null_model"]["swaps_done"] >= 0.9 * 10 * named["n_edges"]


def test_erdos_renyi_draws_weights_from_the_empirical_list(named):
    s = shuffle_graph(named, "erdos_renyi", 0)
    empirical = set(float(e["weight"]) for e in named["edges"])
    assert all(float(e["weight"]) in empirical for e in s["edges"])
    assert set(_pairs(s)) != set(_pairs(named))


def test_null_distribution_shape_and_finite_z():
    t0 = time.perf_counter()
    res = null_distribution("subgraph", "imidacloprid", 1e-6, "sign_permute", n=10, seed=0)
    # Generous bound: this asserts the shuffles are not accidentally quadratic,
    # not a benchmark. CI runners and parallel test runs are slow and variable.
    assert time.perf_counter() - t0 < 60.0
    assert res["n"] == 10 and len(res["null_effects"]) == 10
    assert res["readout"] == "mean_hz"  # 'auto': MN9 is clamped on this assay
    assert res["real_effect"] is not None and np.isfinite(res["real_effect"])
    assert res["z"] is not None and np.isfinite(res["z"])
    assert 0.0 < res["p_two_sided"] <= 1.0
    assert res["p_two_sided"] >= 1.0 / (res["n_ok"] + 1)
    assert res["label"] == "model_derived" and res["warnings"]
    assert res["runtime_s"] > 0
    json.dumps(res)  # JSON-serialisable


def test_vehicle_and_treated_run_on_the_same_shuffled_graph():
    """The effect must be a within-graph contrast, so a null graph that kills
    the vehicle rate also kills the treated rate."""
    res = null_distribution("subgraph", "imidacloprid", 1e-6, "erdos_renyi", n=5, seed=0)
    assert all(v is not None for v in res["null_effects"])
    # an Erdos-Renyi graph has no MN9 neighborhood left to suppress
    assert abs(res["null_mean"]) < abs(res["real_effect"])


def test_imidacloprid_is_outside_the_sign_permute_null():
    """Finding, not a knob: the drug effect on the real wiring must sit off the
    centre of the NT-histogram-preserving null."""
    res = null_distribution("subgraph", "imidacloprid", 1e-6, "sign_permute", n=50, seed=0)
    assert abs(res["z"]) > 1.0, (
        f"imidacloprid mean_hz effect is inside the sign-permute null: "
        f"z={res['z']:.2f}, real={res['real_effect']:.3f}, "
        f"null={res['null_mean']:.3f}+/-{res['null_sd']:.3f}"
    )


def test_null_panel_and_information_score():
    panel = null_panel(
        "imidacloprid", 1e-6, modes=("sign_permute", "erdos_renyi"), n=5, n_taste=2, seed=0
    )
    assert {r["mode"] for r in panel["rows"]} == {"sign_permute", "erdos_renyi"}
    assert {r["assay"] for r in panel["rows"]} == {"subgraph", "taste_map"}
    assert all("z" in r and "p_two_sided" in r for r in panel["rows"])
    score = panel["summary"]
    assert score["score"] is None or 0.0 <= score["score"] <= 1.0
    assert score["n_rows"] == len(panel["rows"])
    assert panel["label"] == "model_derived" and panel["warnings"]
    json.dumps(panel)

    sub = connectome_information_score(panel=panel, assay="subgraph")
    assert sub["n_rows"] == 2
    assert sub["label"] == "model_derived"


def test_unknown_mode_and_assay_are_refused(named):
    with pytest.raises(ValueError):
        shuffle_graph(named, "not_a_mode", 0)
    with pytest.raises(ValueError):
        null_distribution("not_an_assay", "imidacloprid", 1e-6, "sign_permute", n=1)


# ==========================================================================
# v0.6: high-permutation runs, the empirical p, parallelism and convergence
# ==========================================================================
from flylab.analysis.nullmodels import (  # noqa: E402  (grouped with the v0.6 tests)
    DEFAULT_TOL,
    convergence_report,
    drug_effect,
    empirical_p,
    fast_drug_effect,
    p_resolution,
)


# --------------------------------------------------------------------------
# the shuffles keep their invariants at high n, and the array path used by the
# null loops is the same shuffle as the public dict API
# --------------------------------------------------------------------------
@pytest.mark.parametrize("mode", MODES)
def test_invariants_hold_over_many_shuffles(named, mode):
    """Many draws of every mode, not one: an invariant that only holds for
    seed 0 is not an invariant."""
    nts = Counter(n["consensus_nt"] for n in named["nodes"])
    weights = sorted(e["weight"] for e in named["edges"])
    pre_deg = Counter(e["pre"] for e in named["edges"])
    post_deg = Counter(e["post"] for e in named["edges"])
    # the degree-preserving rewire is ~20x the cost of the other modes
    reps = 60 if mode == "rewire_degree_preserving" else 200
    for i in range(reps):
        s = shuffle_graph(named, mode, np.random.default_rng([99, i]))
        assert s["n_nodes"] == named["n_nodes"] and s["n_edges"] == named["n_edges"]
        assert s["seeds"] == named["seeds"]
        assert [n["bodyId"] for n in s["nodes"]] == [n["bodyId"] for n in named["nodes"]]
        pairs = _pairs(s)
        assert len(set(pairs)) == len(pairs)
        assert all(a != b for a, b in pairs)
        assert Counter(n["consensus_nt"] for n in s["nodes"]) == nts
        if mode == "sign_permute":
            assert pairs == _pairs(named)
        elif mode == "weight_permute":
            assert pairs == _pairs(named)
            assert sorted(e["weight"] for e in s["edges"]) == weights
        elif mode == "rewire_degree_preserving":
            assert Counter(e["pre"] for e in s["edges"]) == pre_deg
            assert Counter(e["post"] for e in s["edges"]) == post_deg
        else:
            assert set(float(e["weight"]) for e in s["edges"]) <= set(weights)


def test_array_shuffles_match_the_dict_api(named):
    """The fast path must shuffle exactly the graph shuffle_graph documents."""
    from flylab.analysis.nullmodels import _GraphState, _shuffle_state

    base = _GraphState.from_graph(named)
    body = [n["bodyId"] for n in named["nodes"]]
    for mode in MODES:
        for i in (0, 1, 7):
            st = _shuffle_state(base, mode, np.random.default_rng([5, i]))
            sg = shuffle_graph(named, mode, np.random.default_rng([5, i]))
            assert [body[p] for p in st.pre.tolist()] == [e["pre"] for e in sg["edges"]]
            assert [body[q] for q in st.post.tolist()] == [e["post"] for e in sg["edges"]]
            assert st.w.tolist() == [float(e["weight"]) for e in sg["edges"]]
            assert list(st.nt) == [n["consensus_nt"] for n in sg["nodes"]]


# --------------------------------------------------------------------------
# the rate engine is the assay, only faster
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "assay,graph,readout,compound",
    [
        ("subgraph", "named", "mean_hz", "imidacloprid"),
        ("subgraph", "named", "mn9_hz", "fipronil"),
        ("taste_map", "taste_motor", "bitter_veto_ratio", "fipronil"),
        ("taste_map", "taste_motor", "mn9_hz", "fipronil"),
    ],
)
def test_fast_engine_reproduces_the_notebook_assay(assay, graph, readout, compound):
    slow = drug_effect(assay, compound, 1e-6, readout, graph=graph)
    fast = fast_drug_effect(assay, compound, 1e-6, readout, graph=graph)
    assert fast["effect"] == pytest.approx(slow["effect"], abs=1e-9)
    assert fast["treated"] == pytest.approx(slow["treated"], abs=1e-9)
    assert fast["vehicle"] == pytest.approx(slow["vehicle"], abs=1e-9)


def test_engine_choice_does_not_change_the_null(named):
    a = null_distribution("subgraph", "fipronil", 1e-6, "sign_permute", n=8, seed=11)
    b = null_distribution(
        "subgraph", "fipronil", 1e-6, "sign_permute", n=8, seed=11, engine="assay"
    )
    assert a["engine"] == "fast" and b["engine"] == "assay"
    assert a["real_effect"] == pytest.approx(b["real_effect"], abs=1e-9)
    for x, y in zip(a["null_effects"], b["null_effects"]):
        assert x == pytest.approx(y, abs=1e-9)
    assert a["p_two_sided"] == b["p_two_sided"]


def test_lif_assay_falls_back_to_the_notebook_path():
    res = null_distribution("spiking", "fipronil", 1e-6, "sign_permute", n=2, seed=0)
    assert res["engine"] == "assay"
    assert any("rate engine unavailable" in w for w in res["warnings"])


def test_a_thousand_shuffles_of_the_named_graph_are_practical():
    """The reviewer's requirement: the principal null analyses must run at
    n ~ 1000. Generous bound; this asserts the loop is not the old one."""
    t0 = time.perf_counter()
    res = null_distribution("subgraph", "imidacloprid", 1e-6, "sign_permute", n=1000, seed=0)
    elapsed = time.perf_counter() - t0
    assert res["n_ok"] == 1000
    assert elapsed < 90.0, f"1000 sign-permute shuffles took {elapsed:.1f} s"
    assert res["p_resolution"] == pytest.approx(1 / 1001)


# --------------------------------------------------------------------------
# empirical p, its resolution, and the convergence check
# --------------------------------------------------------------------------
def test_p_resolution_is_reported_and_correct():
    res = null_distribution("subgraph", "fipronil", 1e-6, "erdos_renyi", n=40, seed=0)
    assert res["p_resolution"] == pytest.approx(1.0 / (res["n_ok"] + 1))
    assert res["p_two_sided"] >= res["p_resolution"] - 1e-12
    # this arm is outside its null, so p sits exactly on the floor and says so
    assert res["p_two_sided"] == pytest.approx(res["p_resolution"])
    assert res["resolution_limited"] is True
    assert any("resolution floor" in w for w in res["warnings"])
    # and the panel-level statement is consistent with the row
    assert res["beats_null"] is (res["p_two_sided"] <= res["alpha"])


def test_empirical_p_counts_both_tails_and_includes_the_observation():
    nulls = [0.0] * 10
    assert empirical_p(5.0, nulls)["p"] == pytest.approx(1 / 11)
    assert empirical_p(0.0, nulls)["p"] == pytest.approx(1.0)
    # two-sided: the sign of the real effect must not matter
    sym = [-3.0, -1.0, 0.0, 1.0, 3.0]
    assert empirical_p(2.0, sym)["p"] == empirical_p(-2.0, sym)["p"]
    assert p_resolution(0) == 1.0 and p_resolution(999) == pytest.approx(1 / 1000)


def test_convergence_report_flags_a_moving_and_a_settled_p():
    settled = [0.0] * 400
    out = convergence_report(1.0, settled, tol=DEFAULT_TOL)
    assert [c["n"] for c in out["checkpoints"]] == [100, 200, 400]
    assert out["stabilised"] is True
    assert out["n_stabilised"] == 100
    assert out["verdict_stable"] is True

    moving = [0.0] * 100 + [5.0] * 100  # the null only widens in the second half
    out2 = convergence_report(1.0, moving, tol=0.001)
    assert out2["stabilised"] is False
    assert out2["p_max_shift"] > 0.001

    empty = convergence_report(None, [], tol=DEFAULT_TOL)
    assert empty["stabilised"] is False and empty["checkpoints"] == []


# --------------------------------------------------------------------------
# parallelism: identical, not merely similar
# --------------------------------------------------------------------------
@pytest.mark.parametrize("mode", ["sign_permute", "erdos_renyi"])
def test_parallel_is_identical_to_serial_for_the_same_seed(mode):
    serial = null_distribution("subgraph", "imidacloprid", 1e-6, mode, n=12, seed=5)
    parallel = null_distribution(
        "subgraph", "imidacloprid", 1e-6, mode, n=12, seed=5, n_jobs=3
    )
    assert serial["null_effects"] == parallel["null_effects"]
    assert serial["p_two_sided"] == parallel["p_two_sided"]
    assert serial["z"] == parallel["z"]
    assert serial["null_mean"] == parallel["null_mean"]
    assert parallel["n_jobs"] == 3
    assert any("identical to a serial run" in w for w in parallel["warnings"])


# --------------------------------------------------------------------------
# the panel now carries the draws F6 is drawn from
# --------------------------------------------------------------------------
def test_panel_rows_carry_their_null_effects_and_p_resolution():
    panel = null_panel(
        "fipronil", 1e-6, modes=("sign_permute", "erdos_renyi"), n=6, n_taste=2, seed=0
    )
    for row in panel["rows"]:
        assert len(row["null_effects"]) == row["n"]
        assert row["p_resolution"] == pytest.approx(1.0 / (row["n_ok"] + 1))
        assert "stabilised" in row and "engine" in row
        # the figure pipeline rebuilds the histogram from this alone
        draws = [v for v in row["null_effects"] if v is not None]
        assert len(draws) == row["n_ok"]
        if draws:
            assert row["null_mean"] == pytest.approx(float(np.mean(draws)))
    score = panel["summary"]
    assert score["score_p"] is None or 0.0 <= score["score_p"] <= 1.0
    assert score["alpha"] == 0.05
    json.dumps(panel)


def test_taste_map_arm_can_run_several_hundred_shuffles():
    """The taste-motor veto arm used ten shuffles in v0.5; it must now be able
    to run the hundreds the reviewer asked for."""
    t0 = time.perf_counter()
    res = null_distribution(
        "taste_map", "fipronil", 1e-6, "sign_permute", n=300, seed=0,
        readout="bitter_veto_ratio",
    )
    elapsed = time.perf_counter() - t0
    assert res["n"] == 300 and res["engine"] == "fast"
    # some shuffled taste graphs silence MN9, so the veto ratio has no
    # denominator there; those draws are dropped and the shortfall is reported.
    assert res["n_ok"] >= 200
    assert res["p_resolution"] == pytest.approx(1 / (res["n_ok"] + 1))
    if res["n_ok"] < 300:
        assert any("undefined readout" in w for w in res["warnings"])
    assert elapsed < 120.0, f"300 taste-map shuffles took {elapsed:.1f} s"


@pytest.mark.slow
def test_thousand_shuffle_panel_on_both_arms():
    panel = null_panel("fipronil", 1e-6, n=1000, n_taste=300, seed=0, n_jobs=1)
    assert all(row["n"] in (1000, 300) for row in panel["rows"])
    assert all(row["p_resolution"] <= 1 / 301 for row in panel["rows"])


def test_verdict_stability_is_reported_separately_from_p_stability():
    """A mid-range p can still be moving when the verdict at alpha is not:
    the two are reported separately so the paper can say which it means."""
    res = null_distribution("subgraph", "imidacloprid", 1e-6, "sign_permute", n=40, seed=0)
    assert res["stabilised"] is False  # p is still moving at n = 40
    assert res["verdict_stable"] is True  # but it was never near alpha
    assert any("the verdict at alpha never changed" in w for w in res["warnings"])
    assert res["convergence"]["verdict_stable"] is True
