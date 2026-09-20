"""Per-edge / per-node / per-path drug impact."""
from __future__ import annotations

import pytest

from flylab.analysis.impact import (
    edge_impact,
    grn_to_mn9_paths,
    node_impact,
    nt_gain,
    path_impact,
    summarize_impact,
    totals_by_transmitter,
)
from flylab.circuit.rate import compute_gains, load_graph, rate_network


def test_nt_gain_applies_the_right_key():
    gains = {"g_ach": 0.2, "g_gaba": 0.5, "g_glu": 1.5, "ach_tone": 2.0}
    assert nt_gain("acetylcholine", gains) == pytest.approx(0.4)
    assert nt_gain("gaba", gains) == pytest.approx(0.5)
    assert nt_gain("glutamate", gains) == pytest.approx(1.5)
    assert nt_gain("unclear", gains) == pytest.approx(1.0)


def test_edge_impact_only_moves_ach_edges_for_a_neonicotinoid():
    g = load_graph("named")
    gains, _ = compute_gains("imidacloprid", 1e-5)
    rows = edge_impact(g, gains)
    assert len(rows) == len(g["edges"])
    assert abs(rows[0]["delta"]) >= abs(rows[-1]["delta"])
    for r in rows:
        if r["nt"] == "acetylcholine":
            assert r["eff_weight_after"] < r["eff_weight_before"]
        else:
            assert r["eff_weight_after"] == pytest.approx(r["eff_weight_before"])
    totals = totals_by_transmitter(rows)
    assert totals["acetylcholine"]["ratio"] == pytest.approx(gains["g_ach"], rel=1e-9)
    assert totals["gaba"]["ratio"] == pytest.approx(1.0)


def test_node_impact_is_sorted_and_labelled():
    g = load_graph("named")
    net = rate_network("named")
    gains, _ = compute_gains("imidacloprid", 1e-5)
    drive = net.drive_vector(net.seed_drive(40.0))
    veh = net.run(drive, None)
    tre = net.run(drive, gains)
    rows = node_impact(g, veh, tre, top=20)
    assert len(rows) == 20
    deltas = [abs(r["delta_hz"]) for r in rows]
    assert deltas == sorted(deltas, reverse=True)
    assert set(rows[0]) >= {"bodyId", "type", "superclass", "nt", "delta_hz"}


def test_path_impact_reaches_the_seeds():
    g = load_graph("named")
    gains, _ = compute_gains("imidacloprid", 1e-5)
    rows = path_impact(g, None, gains, max_len=2, top=10)
    assert 0 < len(rows) <= 10
    seeds = {b for ids in g["seeds"].values() for b in ids}
    for r in rows:
        assert r["path"][-1] in seeds
        assert 1 <= r["length"] <= 2
        assert r["min_weight"] >= g["min_weight"]
        assert set(r) >= {"product_before", "product_after", "types", "nts"}


def test_grn_to_mn9_paths_finds_sweet_and_bitter_routes():
    g = load_graph("taste_motor")
    gains, _ = compute_gains("imidacloprid", 1e-5)
    out = grn_to_mn9_paths(g, max_len=3, top=20, gains=gains)
    assert out["n_sweet"] > 100 and out["n_bitter"] > 100
    assert len(out["sweet"]) == 20 and len(out["bitter"]) == 20
    mn9 = set(g["seeds"]["MN9"])
    sweet_ids = {b for t in ("LB3b", "LB3c") for b in g["seeds"][t]}
    for r in out["sweet"]:
        assert r["path"][0] in sweet_ids
        assert r["path"][-1] in mn9
        assert r["modality"] == "sweet"
        assert len(r["path"]) - 1 <= 3
    assert "hypothesis" in out["note"].lower()
    # the ACh block must shrink ACh-fronted paths
    ach_sweet = [r for r in out["sweet"] if r["nts"][0] == "acetylcholine"]
    assert ach_sweet and all(abs(r["product_after"]) < abs(r["product_before"]) for r in ach_sweet)


def test_grn_paths_on_the_named_graph_explain_themselves():
    out = grn_to_mn9_paths(load_graph("named"))
    assert out["n_sweet"] == 0 and "taste_motor" in out["note"]


def test_summarize_impact_is_json_ready():
    import json

    out = summarize_impact("imidacloprid", 1e-6)
    assert out["compound"] == "imidacloprid"
    assert len(out["top_edges"]) == 25
    assert len(out["top_nodes"]) == 25
    assert len(out["top_paths"]) == 10
    assert set(out["totals_by_transmitter"]) >= {"acetylcholine", "gaba"}
    assert out["label"] == "model_derived"
    assert out["warnings"]
    json.dumps(out)
