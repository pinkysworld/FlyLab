"""LIF network + spiking assay: determinism, runtime, and drug direction."""
from __future__ import annotations

import json
import time

import numpy as np
import pytest

from flylab.assays.spiking import ASSAY, run_spiking_assay
from flylab.circuit.lif import LIFNetwork, SpikeResult
from flylab.circuit.rate import DEFAULT_GAINS, compute_gains, load_graph

MN9 = (10331, 16949)


def _net_and_drive(t_hz=40.0, bg=26.0, seed=0):
    g = load_graph("named")
    net = LIFNetwork(g, seed=seed)
    drive = {n["bodyId"]: bg for n in g["nodes"]}
    for ids in g["seeds"].values():
        for b in ids:
            drive[b] = t_hz
    return net, drive


def test_lif_is_deterministic_for_a_seed():
    net, drive = _net_and_drive()
    a = net.run(drive, None, t_ms=300.0)
    b = net.run(drive, None, t_ms=300.0)
    assert isinstance(a, SpikeResult)
    assert a.spikes == b.spikes
    assert np.array_equal(a.rates_hz, b.rates_hz)
    net.seed = 7
    c = net.run(drive, None, t_ms=300.0)
    assert c.spikes != a.spikes, "a different seed must give a different realisation"


def test_500ms_run_is_fast():
    net, drive = _net_and_drive()
    t0 = time.perf_counter()
    res = net.run(drive, None, t_ms=500.0)
    elapsed = time.perf_counter() - t0
    assert res.t_ms == pytest.approx(500.0)
    assert elapsed < 3.0, f"500 ms of network took {elapsed:.2f}s"


def test_vehicle_mn9_is_in_a_plausible_band():
    nb = run_spiking_assay(None, 0.0, t_ms=500.0)
    mn9 = nb["readouts"]["mn9_hz"]
    assert 5.0 <= mn9 <= 60.0, f"vehicle MN9 calibration drifted: {mn9} Hz"
    assert 0.2 < nb["readouts"]["frac_active"] < 0.95


def test_imidacloprid_lowers_mn9_in_the_lif():
    veh = run_spiking_assay(None, 0.0, t_ms=500.0)["readouts"]
    imi = run_spiking_assay("imidacloprid", 1e-5, t_ms=500.0)["readouts"]
    assert imi["mn9_hz"] < veh["mn9_hz"]
    assert imi["mean_hz"] < veh["mean_hz"]
    assert imi["vehicle"]["mn9_hz"] == pytest.approx(veh["mn9_hz"])


def test_fipronil_moves_g_gaba_not_g_ach_in_the_spiking_notebook():
    fip = run_spiking_assay("fipronil", 1e-6, t_ms=200.0)
    imi = run_spiking_assay("imidacloprid", 1e-6, t_ms=200.0)
    assert fip["gains"]["g_gaba"] < 0.5 and fip["gains"]["g_ach"] == pytest.approx(1.0)
    assert imi["gains"]["g_ach"] < 0.5 and imi["gains"]["g_gaba"] == pytest.approx(1.0)


def test_spiking_notebook_shape_and_size():
    nb = run_spiking_assay("imidacloprid", 1e-6, t_ms=500.0, seed=3)
    assert nb["assay"] == ASSAY
    r = nb["readouts"]
    for key in ("named", "mn9_hz", "dnp01_hz", "by_superclass", "raster", "psth", "lif_params"):
        assert key in r
    seeds = {b for ids in load_graph("named")["seeds"].values() for b in ids}
    raster_ids = {row["bodyId"] for row in r["raster"]}
    assert seeds <= raster_ids, "named seed cells must always be in the raster"
    assert len(raster_ids) <= 110
    for row in r["raster"]:
        assert set(row) >= {"bodyId", "type", "nt", "spikes_ms"}
    psth = r["psth"]
    assert psth["bin_ms"] == 10.0 and psth["n_bins"] == 50
    assert len(psth["cells"]["MN9"][0]["counts"]) == 50
    size = len(json.dumps(nb).encode())
    assert size < 300_000, f"notebook JSON is {size} bytes"
    assert any("not" in w.lower() for w in nb["warnings"])
    assert any("fitted" in w.lower() for w in nb["warnings"])


def test_drive_dict_can_target_any_node():
    g = load_graph("named")
    other = next(n["bodyId"] for n in g["nodes"] if n["bodyId"] not in MN9)
    nb = run_spiking_assay(None, 0.0, t_ms=200.0, drive={other: 80.0})
    assert nb["readouts"]["n_spikes"] >= 1


def test_gains_change_effective_weights_only_for_their_transmitter():
    g = load_graph("named")
    net = LIFNetwork(g, seed=0)
    gains, _ = compute_gains("imidacloprid", 1e-5)
    base = net.effective_weights(DEFAULT_GAINS)
    drug = net.effective_weights(gains)
    ach = np.array([net.nt_of_node[b] == "acetylcholine" for b in net.body_ids])
    assert not np.allclose(base[:, ach], drug[:, ach])
    assert np.allclose(base[:, ~ach], drug[:, ~ach])
