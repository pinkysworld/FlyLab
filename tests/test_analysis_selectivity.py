"""Receptor-level vs circuit-level selectivity."""
from __future__ import annotations

import json
import math
import time

import pytest

from flylab.analysis.selectivity import (
    DEFAULT_CONCS,
    circuit_selectivity_index,
    circuit_threshold_conc,
    graph_composition,
    receptor_selectivity_table,
    selectivity_landscape,
    vertebrate_threshold_conc,
)
from flylab.pharm.occupancy import list_compounds


def test_receptor_table_covers_every_compound():
    table = receptor_selectivity_table(1e-6)
    names = list_compounds()
    assert table["n_compounds"] == len(names)
    assert [r["compound"] for r in table["rows"]] == names
    for r in table["rows"]:
        assert 0.0 <= r["max_insect_occupancy"] <= 1.0
        assert 0.0 <= r["max_vertebrate_occupancy"] <= 1.0
        assert r["receptor_si_log10"] is not None
    assert table["label"] == "model_derived" and table["warnings"]
    json.dumps(table)


def test_imidacloprid_is_insect_selective_at_the_receptor_level():
    rows = {r["compound"]: r for r in receptor_selectivity_table(1e-6)["rows"]}
    imi = rows["imidacloprid"]
    assert imi["max_insect_occupancy"] >= 0.9
    assert imi["max_vertebrate_occupancy"] <= 0.2
    assert imi["receptor_si_log10"] > 1.0


def test_picrotoxin_receptor_si_is_about_zero():
    """The non-selective control: insect RDL and vertebrate GABA-A at the same
    order of magnitude."""
    rows = {r["compound"]: r for r in receptor_selectivity_table(1e-6)["rows"]}
    pic = rows["picrotoxin"]
    assert pic["best_pair"] == "GABA_A"
    assert abs(pic["receptor_si_log10"]) < 0.3
    assert not pic["placeholder"]


def test_vertebrate_threshold_conc_matches_the_hill_equation():
    v = vertebrate_threshold_conc("imidacloprid", 0.2)
    assert v["receptor"].startswith("vertebrate_")
    # theta = 0.2 -> C = EC50 * (0.2/0.8) = EC50/4 for n = 1
    assert v["conc_M"] == pytest.approx(v["ec50_M"] * 0.25, rel=1e-6)


def test_circuit_threshold_is_finite_for_imidacloprid_on_subgraph_mean_hz():
    res = circuit_threshold_conc("subgraph", "imidacloprid", "mean_hz", threshold_frac=0.5)
    assert res["conc_M"] is not None and math.isfinite(res["conc_M"])
    assert DEFAULT_CONCS[0] <= res["conc_M"] <= DEFAULT_CONCS[-1]
    assert len(res["curve"]) == len(DEFAULT_CONCS) == 13
    assert res["curve"][0]["rel_change"] < 0.5 <= res["curve"][-1]["rel_change"]
    assert res["label"] == "model_derived" and res["warnings"]
    json.dumps(res)


def test_circuit_threshold_returns_none_when_the_readout_never_moves_enough():
    res = circuit_threshold_conc("subgraph", "diazepam", "mean_hz", threshold_frac=0.5)
    assert res["conc_M"] is None
    assert any("never reaches" in w for w in res["warnings"])
    assert res["max_rel_change"] == pytest.approx(0.0, abs=1e-9)


def test_circuit_selectivity_index_reports_both_indices():
    si = circuit_selectivity_index("imidacloprid", assay="subgraph", readout="mean_hz")
    assert si["circuit_si_log10"] is not None
    assert si["receptor_si_log10"] is not None
    # circuit effect must start below the vertebrate occupancy threshold here
    assert si["circuit_threshold_M"] < si["vertebrate_threshold_M"]
    assert si["circuit_si_log10"] == pytest.approx(
        math.log10(si["vertebrate_threshold_M"]) - math.log10(si["circuit_threshold_M"])
    )
    assert si["si_gap_circuit_minus_receptor"] == pytest.approx(
        si["circuit_si_log10"] - si["receptor_si_log10"]
    )
    assert si["label"] == "model_derived"


def test_graph_composition_fractions_are_sane():
    for name in ("named", "taste_motor"):
        c = graph_composition(name)
        assert 0.0 < c["frac_ach_nodes"] < 1.0
        assert 0.0 < c["frac_ach_synapses"] < 1.0
        assert sum(c["node_fraction"].values()) == pytest.approx(1.0)
        assert sum(c["synapse_fraction"].values()) == pytest.approx(1.0)
        assert c["ei_ratio_synapses"] > 0


def test_selectivity_landscape_has_both_indices_and_composition():
    t0 = time.perf_counter()
    land = selectivity_landscape(
        assay="subgraph",
        graphs=("named", "taste_motor"),
        compounds=["imidacloprid", "picrotoxin", "diazepam", "deltamethrin"],
        readout="mean_hz",
    )
    assert time.perf_counter() - t0 < 30.0
    assert land["n_rows"] == 8
    by_key = {(r["compound"], r["graph"]): r for r in land["rows"]}
    for (compound, graph), r in by_key.items():
        assert r["receptor_si_log10"] is not None
        assert r["frac_ach_synapses"] > 0 and r["frac_gaba_synapses"] > 0
        if compound == "diazepam":
            # insect rows are all class placeholders -> skipped, not indexed
            assert r["skipped_reason"]
            assert r["circuit_si_log10"] is None
        else:
            assert r["skipped_reason"] is None
    imi = by_key[("imidacloprid", "named")]
    assert imi["circuit_si_log10"] is not None
    assert imi["circuit_beats_receptor"] is not None
    assert land["label"] == "model_derived" and land["warnings"]
    json.dumps(land)
