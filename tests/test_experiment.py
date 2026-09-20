"""Batch experiment designs."""
from __future__ import annotations

import csv
import io
import json
from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from flylab.assays.experiment import ExperimentDesign, design_from_yaml, run_experiment

DESIGN = {
    "assay": "subgraph",
    "compounds": ["imidacloprid", "fipronil"],
    "concs_M": [1e-9, 1e-8, 1e-7, 1e-6],
    "replicates": 3,
    "seed": 0,
    "readouts": ["mn9_hz", "dnp01_hz", "mean_hz", "g_ach", "g_gaba"],
    "jitter_log10": 0.2,
}


def test_row_count_is_compounds_x_concs_x_replicates():
    out = run_experiment(DESIGN)
    expected = len(DESIGN["compounds"]) * len(DESIGN["concs_M"]) * DESIGN["replicates"]
    assert len(out["rows"]) == expected == out["n_rows"]
    for r in out["rows"]:
        assert set(r) >= {"compound", "conc_M", "replicate", "mn9_hz", "g_ach"}
        assert 0 <= r["replicate"] < DESIGN["replicates"]


def test_vehicle_rows_are_kept_separately():
    out = run_experiment(DESIGN)
    assert len(out["vehicle_rows"]) == DESIGN["replicates"]
    assert all(r["compound"] is None and r["conc_M"] == 0.0 for r in out["vehicle_rows"])
    assert all(r["g_ach"] == pytest.approx(1.0) for r in out["vehicle_rows"])


def test_csv_is_valid_and_covers_every_row():
    out = run_experiment(DESIGN)
    rows = list(csv.DictReader(io.StringIO(out["csv"])))
    assert len(rows) == len(out["rows"]) + len(out["vehicle_rows"])
    assert list(rows[0]) == ["compound", "conc_M", "replicate"] + DESIGN["readouts"]
    assert float(rows[0]["mean_hz"]) >= 0.0
    assert out["csv"].endswith("\n")


def test_summary_is_per_compound_and_conc():
    out = run_experiment(DESIGN)
    n_cells = len(DESIGN["compounds"]) * len(DESIGN["concs_M"]) + 1  # + vehicle
    assert len(out["summary"]) == n_cells
    for e in out["summary"]:
        assert e["n"] == DESIGN["replicates"]
        assert "mean_hz_mean" in e and "mean_hz_sd" in e
    imi = [e for e in out["summary"] if e["compound"] == "imidacloprid"]
    veh = next(e for e in out["summary"] if e["compound"] is None)
    top = max(imi, key=lambda e: e["conc_M"])
    assert top["mean_hz_mean"] < veh["mean_hz_mean"]


def test_jitter_creates_replicate_spread_and_no_jitter_does_not():
    spread = run_experiment({**DESIGN, "replicates": 4, "concs_M": [1e-8]})
    vals = [r["mean_hz"] for r in spread["rows"] if r["compound"] == "imidacloprid"]
    assert len(set(vals)) > 1
    flat = run_experiment({**DESIGN, "replicates": 4, "concs_M": [1e-8], "jitter_log10": 0.0})
    vals = [r["mean_hz"] for r in flat["rows"] if r["compound"] == "imidacloprid"]
    assert len(set(vals)) == 1


def test_design_validation():
    with pytest.raises(ValidationError):
        ExperimentDesign(assay="not_an_assay")
    with pytest.raises(ValidationError):
        ExperimentDesign(compounds=[])
    with pytest.raises(ValidationError):
        ExperimentDesign(concs_M=[])
    with pytest.raises(ValidationError):
        ExperimentDesign(replicates=0)
    with pytest.raises(ValidationError):
        ExperimentDesign(concs_M=[-1.0])
    d = ExperimentDesign()
    assert d.assay == "subgraph" and d.replicates == 1


def test_design_from_yaml(tmp_path: Path):
    p = tmp_path / "design.yaml"
    p.write_text(yaml.safe_dump(DESIGN))
    d = design_from_yaml(p)
    assert d["assay"] == "subgraph"
    assert d["compounds"] == DESIGN["compounds"]
    out = run_experiment(d)
    assert len(out["rows"]) == 24
    # a wrapped {"design": {...}} file works too
    p2 = tmp_path / "wrapped.yaml"
    p2.write_text(yaml.safe_dump({"design": DESIGN}))
    assert design_from_yaml(p2)["replicates"] == 3


def test_result_is_json_serialisable_and_labelled():
    out = run_experiment({**DESIGN, "concs_M": [1e-7], "replicates": 1})
    assert out["label"] == "model_derived"
    assert out["warnings"]
    json.dumps(out)


def test_notebooks_are_optional():
    out = run_experiment({**DESIGN, "concs_M": [1e-7], "replicates": 1})
    assert out["notebooks"] == []
    kept = run_experiment({**DESIGN, "concs_M": [1e-7], "replicates": 1, "keep_notebooks": True})
    assert len(kept["notebooks"]) == 3  # 2 compounds + vehicle
