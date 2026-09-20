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
    assert time.perf_counter() - t0 < 10.0
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
