"""The whole-CNS rate assay: one run over all 165 122 traced cells.

What is being pinned down here:

* the assay refuses clearly, naming the download command, when the 1.0 GB
  weight matrix is absent, and refuses with numbers rather than thrashing when
  the edge arrays would exceed the memory ceiling;
* the streaming loader never goes through pandas for the matrix, and the
  sparse engine reproduces the dense one to numerical noise -- so "whole CNS"
  is the same model, not a different one;
* the notebook is a normal schema-0.3 notebook with the usual provenance, and
  its warnings state **in plain words** that no dependence testing, no null
  models and no specification robustness exist at this scale, that the gain
  patch is uniform over cells whose receptor expression is unmapped for about
  four fifths of them, and that the result cannot be interrogated;
* the measured runtime and peak memory are in the notebook, not in a comment.

Nothing here needs the weight matrix except the ``slow`` test at the end,
which skips when it is not on disk.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from flylab.assays.fullcns import (
    DEFAULT_MEMORY_CEILING_GB,
    DEFAULT_MIN_WEIGHT,
    FULLCNS_WARNINGS,
    RECEPTOR_COVERAGE,
    SECONDS_PER_EDGE_STEP,
    estimate_fullcns_cost,
    load_full_cns,
    run_fullcns_assay,
)
from flylab.circuit.rate import (
    DEFAULT_GAINS,
    SPARSE_CHUNK,
    RateNetwork,
    SparseRateNetwork,
    load_graph,
    sparse_memory_gb,
)
from flylab.maps.malecns import (
    TRACED_EDGES_BY_MIN_WEIGHT,
    edge_arrays,
    estimate_edge_memory_gb,
    weights_path,
)


# --------------------------------------------------------------------------
# the sparse engine is the same model as the dense one
# --------------------------------------------------------------------------
@pytest.mark.parametrize("graph", ["named", "taste_motor"])
@pytest.mark.parametrize("gains", [None, {"g_ach": 0.4, "g_gaba": 1.6}])
def test_sparse_matches_the_dense_engine(graph, gains):
    g = load_graph(graph)
    dense, sparse = RateNetwork(g), SparseRateNetwork.from_graph(g)
    d = dense.drive_vector(dense.seed_drive(40.0))
    s = sparse.drive_vector(sparse.seed_drive(40.0))
    assert np.array_equal(d, s)
    assert np.abs(dense.run(d, gains) - sparse.run(s, gains)).max() < 1e-9


@pytest.mark.parametrize("mode", ["row_abs", "none", "degree"])
def test_sparse_matches_the_dense_engine_under_every_normalisation(mode):
    g = load_graph("named")
    dense = RateNetwork(g, normalise=mode)
    sparse = SparseRateNetwork.from_graph(g, normalise=mode)
    assert np.allclose(dense.row_denominator, sparse.row_denominator)
    d = dense.drive_vector(dense.seed_drive(40.0))
    assert np.abs(dense.run(d, {"g_ach": 0.3}) - sparse.run(d, {"g_ach": 0.3})).max() < 1e-9


def test_sparse_matches_the_null_model_fast_path():
    """Three implementations of one model must agree, or one of them is wrong."""
    from flylab.analysis.nullmodels import _state_for

    state = _state_for("taste_motor")
    sparse = SparseRateNetwork.from_graph(load_graph("taste_motor"))
    d = np.zeros(state.n)
    for i in state.seed_indices(list(state.seeds)):
        d[i] = 40.0
    assert np.abs(state.run(d, {"g_ach": 0.4}) - sparse.run(d, {"g_ach": 0.4})).max() < 1e-9


def test_chunking_does_not_change_the_answer():
    g = load_graph("taste_motor")
    whole = SparseRateNetwork.from_graph(g)
    chunked = SparseRateNetwork.from_graph(g, chunk=37)
    d = whole.drive_vector(whole.seed_drive(40.0))
    assert np.abs(whole.run(d, {"g_ach": 0.4}) - chunked.run(d, {"g_ach": 0.4})).max() < 1e-9


def test_an_unknown_normalisation_is_rejected():
    with pytest.raises(ValueError, match="unknown normalise"):
        SparseRateNetwork.from_graph(load_graph("named"), normalise="sqrt")


def test_the_dense_engine_is_untouched():
    """The frozen regression: the shipped default must not have moved."""
    net = RateNetwork(load_graph("named"))
    r = net.run(net.drive_vector(net.seed_drive(40.0)), DEFAULT_GAINS, steps=80)
    assert float(r.mean()) == pytest.approx(6.6312, abs=1e-3)
    assert net.normalise == "row_abs"


# --------------------------------------------------------------------------
# refusals: no matrix, or not enough memory
# --------------------------------------------------------------------------
def test_a_missing_weight_matrix_names_the_download_command(tmp_path):
    with pytest.raises(FileNotFoundError) as exc:
        weights_path(tmp_path)
    assert "flylab download-malecns --full" in str(exc.value)
    assert "never committed" in str(exc.value)


def test_load_full_cns_fails_the_same_way(tmp_path):
    with pytest.raises(FileNotFoundError, match="download-malecns --full"):
        load_full_cns(dest=tmp_path)


def test_the_memory_ceiling_is_a_refusal_with_numbers(tmp_path):
    with pytest.raises(MemoryError) as exc:
        edge_arrays(tmp_path, min_weight=1, memory_ceiling_gb=0.001)
    msg = str(exc.value)
    assert "GB of edge arrays" in msg and "ceiling" in msg
    assert "min_weight" in msg


def test_memory_estimates_are_stated_and_fall_with_the_synapse_floor():
    assert estimate_edge_memory_gb(1, traced_only=True) > estimate_edge_memory_gb(5, traced_only=True)
    assert TRACED_EDGES_BY_MIN_WEIGHT[1] == 25_563_197
    assert TRACED_EDGES_BY_MIN_WEIGHT[1] > TRACED_EDGES_BY_MIN_WEIGHT[5]
    # the dense engine would need 218 GB for a 99.9 %-zero matrix, which is
    # the reason the sparse one exists; the sparse arrays need 0.3 GB
    assert 8 * 165_122 ** 2 / 1e9 == pytest.approx(218.1, rel=0.01)
    assert estimate_edge_memory_gb(1, traced_only=True) < 1.0


def test_the_cost_estimate_uses_the_traced_edge_count_and_the_measured_constant():
    est = estimate_fullcns_cost(1)
    assert est["n_edges_traced"] == 25_563_197
    assert est["predicted_runtime_s"] == pytest.approx(
        SECONDS_PER_EDGE_STEP * 25_563_197 * 80
    )
    # measured: 30.3 s for the real thing
    assert 20 < est["predicted_runtime_s"] < 45
    assert est["predicted_resident_gb"] == pytest.approx(
        sparse_memory_gb(165_122, 25_563_197, SPARSE_CHUNK)
    )


def test_the_default_floor_is_every_traced_edge():
    assert DEFAULT_MIN_WEIGHT == 1
    assert DEFAULT_MEMORY_CEILING_GB >= 4.0


# --------------------------------------------------------------------------
# the notebook and its warnings
# --------------------------------------------------------------------------
def _toy_cns(n=60, m=400, seed=0):
    """A small stand-in network with the fullcns interface (no 1 GB file)."""
    rng = np.random.default_rng(seed)
    nt = rng.choice(["acetylcholine", "gaba", "glutamate"], size=n, p=[0.6, 0.3, 0.1])
    pre = rng.integers(0, n, m)
    post = rng.integers(0, n, m)
    keep = pre != post
    body_ids = np.arange(9_000_000, 9_000_000 + n, dtype=np.int64)
    net = SparseRateNetwork(
        n=n,
        pre=pre[keep],
        post=post[keep],
        weight=rng.integers(1, 40, keep.sum()).astype(float),
        nt=nt,
        body_ids=body_ids,
        seeds={"MN9": [int(body_ids[0])], "DNp01": [int(body_ids[1])]},
    )
    info = {
        "types": ["MN9", "DNp01"] + [None] * (n - 2),
        "superclasses": ["vnc_motor", "descending_neuron"] + [None] * (n - 2),
        "seeds": {"MN9": 1, "DNp01": 1},
        "load_total_s": 0.1,
        "edge_arrays_gb": 1e-6,
        "predicted_resident_gb": net.memory_gb(),
    }
    return net, info


@pytest.fixture(scope="module")
def toy_notebook():
    net, info = _toy_cns()
    return run_fullcns_assay("imidacloprid", 1e-6, net=net, net_info=info)


def test_it_is_a_normal_notebook_with_the_usual_provenance(toy_notebook):
    nb = toy_notebook
    assert nb["flylab_notebook_version"] == "0.3"
    assert nb["assay"] == "malecns_full_cns"
    assert nb["map"]["name"] == "male-cns:v1.0"
    assert nb["map"]["version"] == "v1.0-full"
    assert nb["map"]["citation"]
    prov = nb["provenance"]
    assert prov["map"]["name"] == "male-cns:v1.0"
    assert "library_sha256" in prov and "platform" in prov and "flylab_version" in prov
    assert nb["compound"] == "imidacloprid"
    assert nb["gains"]["g_ach"] < 1.0  # the compound actually did something
    assert nb["label"] == "model_derived"
    json.dumps(nb)


def test_the_warnings_say_no_dependence_testing_is_available(toy_notebook):
    text = " ".join(toy_notebook["warnings"])
    assert "no dependence testing" in text.lower()
    assert "null model" in text.lower()
    assert "specification-robustness" in text.lower()
    assert "CANNOT BE INTERROGATED" in text


def test_the_warnings_say_the_gain_patch_is_uniform_and_expression_unmapped(toy_notebook):
    text = " ".join(toy_notebook["warnings"])
    assert "UNIFORMLY" in text
    assert "four fifths" in text
    assert "27.4" in text  # the measured nAChR class coverage
    assert RECEPTOR_COVERAGE["insect_nAChR"] < 0.3
    assert RECEPTOR_COVERAGE["insect_OctR"] == 0.0


def test_the_warnings_are_the_shipped_list(toy_notebook):
    for w in FULLCNS_WARNINGS:
        assert w in toy_notebook["warnings"]


def test_runtime_and_peak_memory_are_reported_not_asserted(toy_notebook):
    rt = toy_notebook["readouts"]["runtime"]
    mem = toy_notebook["readouts"]["memory"]
    assert rt["run_treated_s"] > 0 and rt["run_vehicle_s"] > 0
    assert rt["total_s"] >= rt["run_treated_s"]
    assert rt["seconds_per_edge_step"] > 0
    assert mem["peak_rss_gb"] > 0
    assert mem["memory_ceiling_gb"] == DEFAULT_MEMORY_CEILING_GB
    assert mem["chunk_edges"] == SPARSE_CHUNK


def test_the_readouts_carry_the_drug_contrast_and_the_structure(toy_notebook):
    r = toy_notebook["readouts"]
    assert r["n_nodes"] == 60
    assert r["mean_hz"] >= 0 and r["vehicle"]["mean_hz"] >= 0
    assert r["delta_mean_hz"] == pytest.approx(r["mean_hz"] - r["vehicle"]["mean_hz"])
    assert r["mn9_hz"] is not None and r["dnp01_hz"] is not None
    assert r["by_superclass"]["vnc_motor"] is not None
    assert "unknown" in r["by_superclass"]
    assert r["top_changed"] and "delta_hz" in r["top_changed"][0]
    assert r["receptor_expression_coverage"]["insect_nAChR"] < 0.3


def test_a_network_whose_activity_does_not_spread_says_so():
    """The whole-CNS run drives 76 cells out of 165 122; that must be visible."""
    net, info = _toy_cns(n=400, m=60, seed=1)  # far too sparse to propagate
    nb = run_fullcns_assay(None, 0.0, net=net, net_info=info)
    if nb["readouts"]["reach_vehicle"] < 0.05:
        assert any("ACTIVITY DOES NOT SPREAD" in w for w in nb["warnings"])
        assert any("NOT comparable with the same readout on a cut" in w for w in nb["warnings"])


def test_vehicle_only_runs_produce_no_effect():
    net, info = _toy_cns()
    nb = run_fullcns_assay(None, 0.0, net=net, net_info=info)
    assert nb["readouts"]["delta_mean_hz"] == pytest.approx(0.0, abs=1e-12)
    assert nb["occupancy"] == []


# --------------------------------------------------------------------------
# the real thing (needs the 1.0 GB matrix)
# --------------------------------------------------------------------------
def _matrix_present() -> bool:
    try:
        weights_path()
        return True
    except FileNotFoundError:
        return False


@pytest.mark.slow
def test_the_whole_cns_runs_within_the_stated_budget():
    if not _matrix_present():
        pytest.skip("MaleCNS weight matrix not downloaded")
    nb = run_fullcns_assay("imidacloprid", 1e-6)
    r = nb["readouts"]
    assert r["n_nodes"] == 165_122
    assert r["n_edges"] == 25_563_197
    assert r["memory"]["peak_rss_gb"] < DEFAULT_MEMORY_CEILING_GB
    assert r["runtime"]["total_s"] < 600
    assert r["mn9_hz"] is not None
