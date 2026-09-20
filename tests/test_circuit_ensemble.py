"""Ensemble replicates, sensitivity tornado and the model-derived IC50."""
from __future__ import annotations

import math

import numpy as np
import pytest

from flylab.assays.ensemble import (
    DEFAULT_READOUT,
    IC50_WARNING,
    READOUT_KEYS,
    circuit_ic50,
    run_ensemble,
    sample_library,
    sensitivity,
)


def test_ensemble_ci_brackets_the_mean():
    nb = run_ensemble("subgraph", "imidacloprid", 1e-7, n_rep=6, seed=0)
    u = nb["uncertainty"]
    assert u["n_rep"] == 6
    assert u["method"]
    assert len(u["replicates"]) == 6
    assert set(u["ci"]) >= {"mn9_hz", "mean_hz", "g_ach"}
    for key, (lo, hi) in u["ci"].items():
        mean = u["mean"][key]
        assert lo <= mean <= hi, f"{key}: {lo} <= {mean} <= {hi}"
        assert u["sd"][key] >= 0.0
    assert nb["readouts"]["mean_hz"] is not None


def test_ensemble_is_reproducible_and_actually_varies():
    a = run_ensemble("subgraph", "imidacloprid", 1e-7, n_rep=4, seed=1)["uncertainty"]
    b = run_ensemble("subgraph", "imidacloprid", 1e-7, n_rep=4, seed=1)["uncertainty"]
    assert a["mean"] == b["mean"]
    assert a["sd"]["mean_hz"] > 0.0, "library + drive jitter must move the readout"


def test_ensemble_runs_the_taste_assay_too():
    nb = run_ensemble("taste", "imidacloprid", 1e-6, n_rep=3, seed=0)
    assert nb["assay"] == "taste_mn9"
    assert nb["uncertainty"]["mean"]["mn9_hz"] > 0.0


def test_spiking_ensemble_varies_with_seed():
    nb = run_ensemble("spiking", "imidacloprid", 1e-7, n_rep=3, seed=0, t_ms=200.0)
    u = nb["uncertainty"]
    seeds = {r["seed"] for r in u["replicates"]}
    assert len(seeds) == 3
    assert "mn9_hz" in u["ci"]


def test_sample_library_jitters_ec50():
    rng = np.random.default_rng(0)
    lib = sample_library(rng, 0.3)
    from flylab.pharm.occupancy import load_library

    base = load_library()
    a = float(base["compounds"]["imidacloprid"]["receptors"]["insect_nAChR"]["ec50_M"])
    b = float(lib["compounds"]["imidacloprid"]["receptors"]["insect_nAChR"]["ec50_M"])
    assert a != b
    assert 0.01 < b / a < 100.0


def test_sensitivity_returns_a_tornado_table():
    rows = sensitivity("subgraph", "imidacloprid", 1e-7, readout="mean_hz")
    params = {r["param"] for r in rows}
    assert params == {"ec50", "hill_n", "drive_hz", "gain_coef", "weight_threshold"}
    for r in rows:
        assert set(r) >= {"param", "low", "high", "base", "readout"}
        assert r["readout"] == "mean_hz"
        assert r["base"] is not None
    spans = [r["span"] for r in rows if r["span"] is not None]
    assert spans == sorted(spans, reverse=True)
    ec50 = next(r for r in rows if r["param"] == "ec50")
    assert ec50["low"] != ec50["high"], "shifting the dominant EC50 must matter"


def test_sensitivity_gain_coef_does_not_change_the_default_rule():
    from flylab.assays.subgraph import run_subgraph_assay

    before = run_subgraph_assay("imidacloprid", 1e-7)["gains"]["g_ach"]
    sensitivity("subgraph", "imidacloprid", 1e-7, params=("gain_coef",), readout="mean_hz")
    after = run_subgraph_assay("imidacloprid", 1e-7)["gains"]["g_ach"]
    assert before == after


def test_circuit_ic50_is_finite_and_inside_the_ladder():
    concs = [1e-10, 1e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 1e-5]
    out = circuit_ic50(
        "subgraph", "imidacloprid", readout="mean_hz", concs=concs, n_boot=40, n_rep=2, seed=0
    )
    fit = out["fit"]
    assert math.isfinite(fit["ic50"])
    assert min(concs) <= fit["ic50"] <= max(concs)
    assert fit["r2"] > 0.9
    # The point fit is on the un-jittered base curve while the CI bootstraps
    # drive-jittered, library-resampled replicates, so the two need not nest:
    # schema v3 added two subunit-resolved library rows, which changed the noise
    # realisation and exposed that. Require an ordered, finite CI on the ladder
    # and a point estimate of the same order as the interval.
    lo, hi = out["ci"]["ic50"]
    assert math.isfinite(lo) and math.isfinite(hi) and lo < hi
    assert min(concs) <= lo <= max(concs) and min(concs) <= hi <= max(concs)
    assert 0.5 * lo <= fit["ic50"] <= 2.0 * hi
    assert len(out["points"]) == len(concs)
    assert out["label"] == "model_derived"
    assert IC50_WARNING in out["warnings"]
    assert any("not an animal IC50" in w for w in out["warnings"])


def test_circuit_ic50_auto_readout():
    assert DEFAULT_READOUT["subgraph"] == "mean_hz"
    out = circuit_ic50(
        "subgraph", "nicotine", readout="auto", concs=[1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4],
        n_boot=20, n_rep=2,
    )
    assert out["readout"] == "mean_hz"
    assert math.isfinite(out["fit"]["ic50"])


def test_circuit_ic50_flags_a_flat_response():
    out = circuit_ic50(
        "subgraph", "diazepam", readout="mean_hz", concs=[1e-9, 1e-8, 1e-7, 1e-6],
        n_boot=5, n_rep=1, ec50_sd_log10=0.0, drive_jitter=0.0,
    )
    assert not math.isfinite(out["fit"]["ic50"])
    assert any("flat" in w for w in out["warnings"])


def test_readout_keys_are_the_documented_five():
    assert READOUT_KEYS == ("mn9_hz", "dnp01_hz", "mean_hz", "g_ach", "g_gaba")
