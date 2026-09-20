"""Rate network: frozen regression values + the refactored API."""
from __future__ import annotations

import numpy as np
import pytest

from flylab.assays.subgraph import run_subgraph_assay
from flylab.circuit.rate import (
    DEFAULT_GAINS,
    GRAPHS,
    RateNetwork,
    compute_gains,
    fallback_gains,
    load_graph,
    rate_network,
    resolve_graph,
)

# Reference numbers captured from flylab 0.4.0 `assays/subgraph.py` BEFORE the
# extraction into circuit/rate.py. They are the contract: the refactor and any
# later gain plumbing must not move them.
REFERENCE = {
    None: {
        "mn9": [43.3870527190687, 33.39442109556258],
        "dnp01": [41.161738927609356, 42.07040399055634],
        "mean_hz": 6.631244285926711,
        "max_hz": 43.3870527190687,
        "g_ach": 1.0,
        "g_gaba": 1.0,
    },
    "imidacloprid": {
        "mn9": [39.59315224167709, 39.21820620512948],
        "dnp01": [39.83951798444125, 39.84941571135308],
        "mean_hz": 0.45846211979249263,
        "max_hz": 39.84941571135308,
        "g_ach": 0.05,
        "g_gaba": 1.0,
    },
    "nitenpyram": {
        "mn9": [39.59315224167709, 39.21820620512948],
        "dnp01": [39.83951798444125, 39.84941571135308],
        "mean_hz": 0.45846211979249263,
        "max_hz": 39.84941571135308,
        "g_ach": 0.05,
        "g_gaba": 1.0,
    },
    "nicotine": {
        "mn9": [38.88019224907244, 37.561053906809065],
        "dnp01": [39.54436053876982, 39.59235271376089],
        "mean_hz": 1.2272536152925329,
        "max_hz": 39.59235271376089,
        "g_ach": 0.17351561697271345,
        "g_gaba": 1.0,
    },
    "diazepam": {
        "mn9": [43.3870527190687, 33.39442109556258],
        "dnp01": [41.161738927609356, 42.07040399055634],
        "mean_hz": 6.631244285926711,
        "max_hz": 43.3870527190687,
        "g_ach": 1.0,
        "g_gaba": 1.0,
    },
    "fipronil": {
        "mn9": [55.37983063868171, 48.90009985416963],
        "dnp01": [44.78378759668054, 45.22314192475127],
        "mean_hz": 7.530851909481498,
        "max_hz": 55.37983063868171,
        "g_ach": 1.0,
        "g_gaba": 0.05,
    },
}


@pytest.mark.parametrize("compound", list(REFERENCE))
def test_rate_regression_against_v0_4_values(compound):
    ref = REFERENCE[compound]
    nb = run_subgraph_assay(compound, 1e-6 if compound else 0.0)
    r = nb["readouts"]
    assert [x["hz"] for x in r["named"]["MN9"]] == pytest.approx(ref["mn9"], rel=1e-12)
    assert [x["hz"] for x in r["named"]["DNp01"]] == pytest.approx(ref["dnp01"], rel=1e-12)
    assert r["mean_hz"] == pytest.approx(ref["mean_hz"], rel=1e-12)
    assert r["max_hz"] == pytest.approx(ref["max_hz"], rel=1e-12)
    assert r["g_ach"] == pytest.approx(ref["g_ach"], rel=1e-12)
    assert r["g_gaba"] == pytest.approx(ref["g_gaba"], rel=1e-12)


def test_notebook_carries_gains_superclass_and_top_changed():
    nb = run_subgraph_assay("imidacloprid", 1e-6)
    assert set(nb["gains"]) >= {"g_ach", "g_gaba", "g_glu", "g_oct", "g_nav", "ach_tone"}
    r = nb["readouts"]
    by_sc = r["by_superclass"]
    assert "cb_intrinsic" in by_sc and "descending_neuron" in by_sc
    assert all(isinstance(v, float) for v in by_sc.values())
    top = r["top_changed"]
    assert len(top) == 15
    assert set(top[0]) >= {"bodyId", "type", "superclass", "nt", "delta_hz"}
    deltas = [abs(t["delta_hz"]) for t in top]
    assert deltas == sorted(deltas, reverse=True)
    assert nb["warnings"], "hops-limited warning must survive the refactor"


def test_fipronil_moves_g_gaba_not_g_ach():
    fip = run_subgraph_assay("fipronil", 1e-6)
    imi = run_subgraph_assay("imidacloprid", 1e-6)
    assert fip["gains"]["g_gaba"] < 0.5
    assert fip["gains"]["g_ach"] == pytest.approx(1.0)
    assert imi["gains"]["g_ach"] < 0.5
    assert imi["gains"]["g_gaba"] == pytest.approx(1.0)


def test_imidacloprid_lowers_network_activity():
    veh = run_subgraph_assay(None, 0.0)["readouts"]
    imi = run_subgraph_assay("imidacloprid", 1e-5)["readouts"]
    assert imi["mean_hz"] < 0.2 * veh["mean_hz"]
    # MN9 itself is clamped by its own 40 Hz drive in this assay; that is the
    # documented limitation, not an accident.
    assert abs(imi["mn9_hz"] - veh["mn9_hz"]) < 5.0


def test_mn9_falls_when_mn9_is_not_driven_itself():
    """Drive MN9's presynaptic partners instead of MN9: the ACh block shows."""
    net = rate_network()
    g = load_graph()
    seeds = set(net.seed_ids)
    upstream = {e["pre"] for e in g["edges"] if e["post"] in seeds} - seeds
    drive = {b: 40.0 for b in upstream}
    gains, _ = compute_gains("imidacloprid", 1e-5)
    veh = net.run(drive, DEFAULT_GAINS)
    tre = net.run(drive, gains)
    mn9 = [net.node_index[b] for b in g["seeds"]["MN9"]]
    assert float(np.mean([tre[i] for i in mn9])) < float(np.mean([veh[i] for i in mn9]))
    assert tre.mean() < veh.mean()


def test_rate_network_exposes_contract_attributes():
    net = rate_network()
    assert net.W_signed.shape == (len(net), len(net))
    assert net.node_index[10331] < len(net)
    assert net.nt_of_node[10331] in {"acetylcholine", "gaba", "glutamate", "unclear", "histamine"}
    # rows are normalised to at most unit total absolute weight
    assert float(np.abs(net.W_signed).sum(axis=1).max()) <= 1.0 + 1e-9


def test_gain_keys_scale_the_right_transmitters():
    net = rate_network()
    idx_ach = [i for i, b in enumerate(net.body_ids) if net.nt_of_node[b] == "acetylcholine"][0]
    idx_gaba = [i for i, b in enumerate(net.body_ids) if net.nt_of_node[b] == "gaba"][0]
    s = net.scale_vector({"g_ach": 0.25, "g_gaba": 2.0, "ach_tone": 2.0})
    assert s[idx_ach] == pytest.approx(0.5)  # g_ach * ach_tone
    assert s[idx_gaba] == pytest.approx(2.0)


def test_g_nav_is_a_global_excitability_factor():
    net = rate_network()
    drive = net.seed_drive(40.0)
    base = net.run(drive, DEFAULT_GAINS)
    hot = net.run(drive, {**DEFAULT_GAINS, "g_nav": 1.5})
    assert hot.mean() > base.mean()


def test_fallback_rule_matches_the_frozen_coefficients():
    rows = [{"receptor": "insect_nAChR", "occupancy": 0.5, "direction": "agonist"}]
    assert fallback_gains(rows)["g_ach"] == pytest.approx(max(0.05, 1 + 0.4 * 0.5 - 1.6 * 0.25))
    rows = [{"receptor": "insect_RDL", "occupancy": 0.8, "direction": "antagonist"}]
    assert fallback_gains(rows)["g_gaba"] == pytest.approx(0.2)
    # the override hook only exists for sensitivity analysis
    rows = [{"receptor": "insect_nAChR", "occupancy": 0.5, "direction": "agonist"}]
    assert fallback_gains(rows, agonist_a=0.8, agonist_b=3.2)["g_ach"] == pytest.approx(0.6)


def test_graph_selection():
    assert set(GRAPHS) >= {"named", "taste_motor"}
    assert resolve_graph() == resolve_graph("named")
    assert resolve_graph("taste_motor").name == GRAPHS["taste_motor"]
    taste = rate_network("taste_motor")
    assert len(taste) > len(rate_network("named"))
    assert "LB3c" in taste.graph["seeds"]
    nb = run_subgraph_assay("imidacloprid", 1e-6, graph="taste_motor")
    assert nb["readouts"]["n_nodes"] == len(taste)


def test_explicit_graph_object_is_honoured():
    g = load_graph("named")
    thin = dict(g)
    thin["edges"] = [e for e in g["edges"] if e["weight"] >= 50]
    thin["n_edges"] = len(thin["edges"])
    net = RateNetwork(thin)
    assert int((np.abs(net.W_signed) > 0).sum()) < int((np.abs(rate_network().W_signed) > 0).sum())
