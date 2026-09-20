import json
from pathlib import Path

from flylab.maps.malecns import CENSUS_FALLBACK, GUSTATORY_SEEDS, load_census_fallback, load_gustatory_seeds
from flylab.assays.wholens import run_wholens_assay


def test_committed_census_is_small_and_complete():
    assert CENSUS_FALLBACK.exists()
    assert CENSUS_FALLBACK.stat().st_size < 200_000
    c = load_census_fallback()
    assert c["map"] == "male-cns:v1.0"
    assert c["n_traced"] > 150_000
    assert c["neurotransmitter_counts"]["acetylcholine"] > 50_000
    assert len(c["named_cells"]["MN9"]) == 2
    assert len(c["named_cells"]["DNp01"]) == 2


def test_gustatory_types_are_proven_on_malecns():
    """HANDOFF item 1: LB1a-d and LB3b/c exist as MaleCNS v1.0 `type` strings."""
    s = load_gustatory_seeds()
    for t in ("LB1a", "LB1b", "LB1c", "LB1d", "LB3b", "LB3c"):
        assert len(s["types"][t]) >= 5, t
    assert set(s["modality_hypothesis"]["bitter"]) == {"LB1a", "LB1b", "LB1c", "LB1d"}
    assert set(s["modality_hypothesis"]["sweet"]) == {"LB3b", "LB3c"}


def test_wholens_runs_without_local_atlas(tmp_path):
    nb = run_wholens_assay("imidacloprid", 1e-6, dest=tmp_path)
    assert nb["map"]["name"] == "male-cns:v1.0"
    assert nb["readouts"]["n_traced"] > 150_000
    assert nb["readouts"]["delta_excitation"] < 0
