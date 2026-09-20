"""Deterministic viewer layout."""
from __future__ import annotations

import json

import pytest

from flylab.analysis.layout import graph_for_viewer, graph_layout
from flylab.circuit.rate import compute_gains, load_graph


@pytest.mark.parametrize("name", ["named", "taste_motor"])
def test_every_node_gets_a_position(name):
    g = load_graph(name)
    pos = graph_layout(g)
    assert len(pos) == len(g["nodes"])
    for n in g["nodes"]:
        x, y = pos[n["bodyId"]]
        assert isinstance(x, float) and isinstance(y, float)


@pytest.mark.parametrize("name", ["named", "taste_motor"])
def test_layout_is_deterministic(name):
    g = load_graph(name)
    assert graph_layout(g, seed=0) == graph_layout(g, seed=11)


def test_named_graph_puts_seeds_in_the_middle():
    g = load_graph("named")
    pos = graph_layout(g)
    seeds = {b for ids in g["seeds"].values() for b in ids}
    assert all(pos[b][0] == 0.0 for b in seeds)
    xs = {p[0] for p in pos.values()}
    assert min(xs) < 0.0 < max(xs), "presynaptic left, postsynaptic right"


def test_taste_motor_puts_grns_left_and_motor_right():
    g = load_graph("taste_motor")
    pos = graph_layout(g)
    grns = [b for t in ("LB1a", "LB1b", "LB1c", "LB1d", "LB3b", "LB3c") for b in g["seeds"][t]]
    motors = list(g["seeds"]["MN9"]) + list(g["seeds"]["DNp01"])
    assert max(pos[b][0] for b in grns) < min(pos[b][0] for b in motors)


def test_viewer_payload_shape():
    g = load_graph("named")
    net_rates = {n["bodyId"]: 1.0 for n in g["nodes"]}
    gains, _ = compute_gains("imidacloprid", 1e-5)
    out = graph_for_viewer(g, rates=net_rates, gains=gains)
    assert out["n_nodes"] == len(g["nodes"])
    node = out["nodes"][0]
    assert set(node) >= {"id", "type", "superclass", "nt", "x", "y", "rate", "is_seed"}
    assert sum(1 for n in out["nodes"] if n["is_seed"]) == 4
    edge = out["edges"][0]
    assert set(edge) >= {"source", "target", "weight", "nt", "eff_weight"}
    assert all(e["weight"] >= 5 for e in out["edges"])
    ach = [e for e in out["edges"] if e["nt"] == "acetylcholine"]
    assert ach and all(abs(e["eff_weight"]) < e["weight"] for e in ach)
    json.dumps(out)


def test_viewer_payload_caps_the_dense_graph():
    out = graph_for_viewer(load_graph("taste_motor"), max_edges=2000)
    assert out["n_edges"] == 2000
    assert out["truncated_edges"] is True
