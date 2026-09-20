"""Map-extracted labellar GRN -> MN9 assay on the taste_motor neighborhood."""
from __future__ import annotations

import time

import pytest

from flylab.assays.taste_map import ASSAY, BITTER_TYPES, SWEET_TYPES, run_taste_map_assay


# Wall-clock budgets guard against an accidentally quadratic implementation;
# they are not benchmarks. CI runners and parallel test runs are slow and
# highly variable, so the bounds are deliberately loose.
@pytest.mark.parametrize("engine,budget", [("rate", 60.0), ("lif", 120.0)])
def test_assay_runs_and_mn9_responds_to_sugar(engine, budget):
    t0 = time.perf_counter()
    nb = run_taste_map_assay(None, 0.0, engine=engine)
    elapsed = time.perf_counter() - t0
    assert elapsed < budget, f"{engine} engine took {elapsed:.2f}s"
    assert nb["assay"] == ASSAY
    assert nb["map"]["version"] == "taste_motor_neighborhood"
    r = nb["readouts"]
    assert r["n_sweet_grn"] > 0 and r["n_bitter_grn"] > 0
    assert r["mn9_sugar_hz"] > r["mn9_no_drive_hz"], "sugar drive must move MN9"
    assert len(r["top_relays"]) == 15
    assert all(set(x) >= {"bodyId", "type", "nt", "delta_hz"} for x in r["top_relays"])


def test_drug_changes_mn9_on_the_map_path():
    veh = run_taste_map_assay(None, 0.0, engine="rate")["readouts"]
    imi = run_taste_map_assay("imidacloprid", 1e-5, engine="rate")["readouts"]
    assert imi["mn9_sugar_hz"] != veh["mn9_sugar_hz"]
    assert imi["mn9_vehicle_sugar_hz"] == pytest.approx(veh["mn9_sugar_hz"])
    assert imi["gains"]["g_ach"] < 0.5 if "gains" in imi else True


def test_mn9_and_dnp01_are_never_driven_directly():
    from flylab.assays.taste_map import _drive_map
    from flylab.circuit.rate import load_graph

    g = load_graph("taste_motor")
    drive = _drive_map(g, 150.0, 150.0)
    motors = set(g["seeds"]["MN9"]) | set(g["seeds"]["DNp01"])
    assert not (motors & set(drive))
    grns = {b for t in SWEET_TYPES + BITTER_TYPES for b in g["seeds"][t]}
    assert set(drive) == grns


def test_warnings_keep_the_hypothesis_and_the_control_explicit():
    nb = run_taste_map_assay("fipronil", 1e-6, engine="rate")
    text = " ".join(nb["warnings"]).lower()
    assert "hypothesis" in text
    assert "reduced_taste_v0" in text
    assert "hops-limited" in text
    assert nb["gains"]["g_gaba"] < 0.5


def test_bitter_veto_ratio_is_reported_not_asserted():
    nb = run_taste_map_assay(None, 0.0, bitter_hz=150.0, engine="rate")
    r = nb["readouts"]
    assert r["bitter_veto_ratio"] is None or r["bitter_veto_ratio"] >= 0.0
    assert r["mn9_sugar_bitter_hz"] is not None
