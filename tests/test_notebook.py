"""Notebook schema 0.3: provenance, gains, backward compatibility."""
import json

from flylab.assays.taste import run_taste_assay
from flylab.notebook.schema import NOTEBOOK_VERSION, empty_notebook, git_sha, provenance, set_map

V02_KEYS = {
    "flylab_notebook_version", "assay", "map", "compound",
    "concentration_M", "occupancy", "readouts", "live_lab", "warnings",
}


def test_version_and_backward_compatible_keys():
    nb = empty_notebook("taste_mn9")
    assert NOTEBOOK_VERSION == "0.3"
    assert nb["flylab_notebook_version"] == "0.3"
    assert V02_KEYS <= set(nb)
    assert nb["map"] == {"name": None, "version": None, "citation": None}
    assert nb["live_lab"] is None


def test_new_03_blocks():
    nb = empty_notebook("taste_mn9", seed=11)
    assert nb["gains"] == {}
    assert nb["uncertainty"] is None
    assert nb["created_utc"].endswith("+00:00")
    prov = nb["provenance"]
    for key in ("flylab_version", "git_sha", "library_sha256", "map", "rng_seed", "platform"):
        assert key in prov
    assert prov["rng_seed"] == 11
    assert prov["library_sha256"] and len(prov["library_sha256"]) == 64
    assert prov["flylab_version"]
    assert "python-" in prov["platform"]


def test_git_sha_is_a_sha_or_none():
    sha = git_sha()
    assert sha is None or (len(sha) == 40 and all(c in "0123456789abcdef" for c in sha))
    assert provenance()["git_sha"] == sha


def test_set_map_updates_both_places():
    nb = set_map(empty_notebook("x"), "male-cns:v1.0", "v1.0", "citation")
    assert nb["map"]["name"] == "male-cns:v1.0"
    assert nb["provenance"]["map"] == nb["map"]


def test_notebook_is_json_serialisable():
    nb = run_taste_assay("imidacloprid", 1e-6)
    text = json.dumps(nb)
    assert json.loads(text)["flylab_notebook_version"] == "0.3"


def test_assay_notebook_carries_provenance_and_gains():
    nb = run_taste_assay("imidacloprid", 1e-6)
    assert nb["provenance"]["library_sha256"]
    assert nb["provenance"]["map"]["name"] == "reduced_taste_v0"
    assert set(nb["gains"]) == {"g_ach", "g_gaba", "g_glu", "g_oct", "g_nav", "ach_tone"}
    assert nb["uncertainty"] is None
    assert nb["live_lab"] is None
