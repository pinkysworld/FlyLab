from pathlib import Path
import pytest
from flylab.maps.malecns import load_census
from flylab.assays.wholens import run_wholens_assay
ATLAS = Path("/tmp/malecns")
@pytest.mark.skipif(not (ATLAS / "body-annotations-male-cns-v1.0-minconf-0.5.feather").exists(), reason="no atlas")
def test_census_has_traced_ach_and_mn9():
    c = load_census(ATLAS)
    assert c["n_traced"] > 100_000
    assert c["neurotransmitter_counts"].get("acetylcholine", 0) > 50_000
    assert len(c["named_cells"]["MN9"]) >= 1
    assert len(c["named_cells"]["DNp01"]) >= 1
@pytest.mark.skipif(not (ATLAS / "body-annotations-male-cns-v1.0-minconf-0.5.feather").exists(), reason="no atlas")
def test_imidacloprid_lowers_excitation_index():
    veh = run_wholens_assay(compound=None, conc_M=0.0, dest=ATLAS)
    treat = run_wholens_assay(compound="imidacloprid", conc_M=1e-5, dest=ATLAS)
    assert treat["readouts"]["treated"]["cns_excitation_index"] < veh["readouts"]["vehicle"]["cns_excitation_index"]
    assert treat["map"]["name"] == "male-cns:v1.0"
