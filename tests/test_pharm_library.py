"""Schema v2 validation of flylab/pharm/library.yaml."""
import math

import pytest

from flylab.pharm.occupancy import (
    ALLOWED_DIRECTIONS,
    ALLOWED_TIERS,
    PLACEHOLDER_ROW,
    compare_compound,
    hill_occupancy,
    library_sha256,
    list_compounds,
    load_library,
    occupancy_curve,
    receptor_spec,
    receptor_table,
    selectivity_pairs,
)

LIB = load_library()
V04_FROZEN = {
    "imidacloprid": {"insect_nAChR": (2.0e-08, 1.2, "agonist"),
                     "vertebrate_nAChR_a4b2": (1.0e-05, 1.0, "agonist")},
    "nitenpyram": {"insect_nAChR": (5.0e-08, 1.2, "agonist"),
                   "vertebrate_nAChR_a4b2": (2.0e-05, 1.0, "agonist")},
    "nicotine": {"insect_nAChR": (2.0e-07, 1.1, "agonist"),
                 "vertebrate_nAChR_a4b2": (1.0e-06, 1.2, "agonist")},
    "diazepam": {"vertebrate_GABA_A": (1.0e-08, 1.3, "positive_modulator")},
    "fipronil": {"insect_RDL": (3.0e-08, 1.2, "antagonist"),
                 "vertebrate_GABA_A": (1.0e-05, 1.0, "antagonist")},
}


def test_schema_and_version():
    assert LIB["schema_version"] == 2
    assert LIB["library_version"] == "0.5"
    assert "evidence_tiers" in LIB["notes"]
    assert 15 <= len(LIB["compounds"]) <= 25


def test_receptor_section_is_complete():
    receptors = receptor_table()
    assert len(receptors) >= 12
    for name, spec in receptors.items():
        assert spec["organism"] in {"insect", "vertebrate"}
        assert spec["family"]
        assert name.startswith(("insect_", "vertebrate_"))


@pytest.mark.parametrize("key", sorted(LIB["compounds"]))
def test_every_row_validates(key):
    entry = LIB["compounds"][key]
    assert entry["name"]
    assert entry["class"]
    rows = entry["receptors"]
    assert rows, f"{key} has no receptor rows"
    for receptor, spec in rows.items():
        assert receptor in LIB["receptors"], f"{key}:{receptor} not in receptors section"
        assert float(spec["ec50_M"]) > 0
        assert float(spec.get("n", 1.0)) > 0
        assert spec["direction"] in ALLOWED_DIRECTIONS
        assert spec["evidence_tier"] in ALLOWED_TIERS
        assert isinstance(spec["source"], str) and len(spec["source"].strip()) > 10
        if "efficacy" in spec:
            assert 0.0 <= float(spec["efficacy"]) <= 1.0
        if spec["evidence_tier"] == "class_placeholder":
            assert float(spec["ec50_M"]) == 0.01
            assert spec["direction"] == "none"
            assert "placeholder" in spec["source"].lower()


@pytest.mark.parametrize("key", sorted(LIB["compounds"]))
def test_each_compound_has_insect_and_vertebrate_rows(key):
    rows = LIB["compounds"][key]["receptors"]
    assert any(r.startswith("insect_") for r in rows)
    assert any(r.startswith("vertebrate_") for r in rows)


@pytest.mark.parametrize("key", sorted(LIB["compounds"]))
def test_missing_receptors_fall_back_to_placeholder(key):
    entry = LIB["compounds"][key]
    for receptor in receptor_table():
        spec = receptor_spec(entry, receptor)
        assert float(spec["ec50_M"]) > 0
        assert spec["direction"] in ALLOWED_DIRECTIONS
        if receptor not in entry["receptors"]:
            assert spec == PLACEHOLDER_ROW


def test_v04_numbers_are_frozen():
    for key, receptors in V04_FROZEN.items():
        for receptor, (ec50, n, direction) in receptors.items():
            spec = LIB["compounds"][key]["receptors"][receptor]
            assert float(spec["ec50_M"]) == ec50
            assert float(spec["n"]) == n
            assert spec["direction"] == direction


def test_named_sources_for_real_values():
    """Non-placeholder rows must name a publication, not just a class."""
    for key, entry in LIB["compounds"].items():
        for receptor, spec in entry["receptors"].items():
            if spec["evidence_tier"] == "class_placeholder":
                continue
            source = spec["source"]
            has_year = any(str(y) in source for y in range(1930, 2027))
            assert has_year, f"{key}:{receptor} source names no year: {source}"


def test_cas_numbers_look_like_cas():
    for key, entry in LIB["compounds"].items():
        cas = entry.get("cas")
        if cas is None:
            continue
        parts = cas.split("-")
        assert len(parts) == 3 and all(p.isdigit() for p in parts), f"{key}: {cas}"


def test_selectivity_pairs_reference_known_receptors():
    receptors = receptor_table()
    pairs = selectivity_pairs()
    assert len(pairs) == 5
    for pair in pairs.values():
        assert receptors[pair["insect"]]["organism"] == "insect"
        assert receptors[pair["vertebrate"]]["organism"] == "vertebrate"


def test_library_sha256_is_stable():
    assert library_sha256() == library_sha256()
    assert len(library_sha256()) == 64


def test_occupancy_curve_spans_the_grid():
    curve = occupancy_curve("imidacloprid")
    assert len(curve) == 17
    assert math.isclose(curve[0]["conc_M"], 1e-11, rel_tol=1e-9)
    assert math.isclose(curve[-1]["conc_M"], 1e-3, rel_tol=1e-9)
    occ = [p["receptors"]["insect_nAChR"] for p in curve]
    assert occ == sorted(occ)
    assert occ[0] < 0.01 < 0.99 < occ[-1]


def test_compare_compound_keeps_v04_keys():
    result = compare_compound("imidacloprid", 1e-6)
    for key in ("compound", "class", "concentration_M", "receptors", "disclaimer"):
        assert key in result
    row = result["receptors"][0]
    for key in ("receptor", "occupancy", "ec50_M", "direction", "source", "evidence_tier"):
        assert key in row


def test_list_compounds_matches_library():
    assert list_compounds() == sorted(LIB["compounds"])


def test_hill_sanity():
    assert hill_occupancy(0.0, 1e-6) == 0.0
    assert hill_occupancy(1e-6, 1e-6) == pytest.approx(0.5)
    assert hill_occupancy(1e-3, 1e-6) > 0.99
    steep = hill_occupancy(2e-6, 1e-6, n=3.0)
    shallow = hill_occupancy(2e-6, 1e-6, n=1.0)
    assert steep > shallow
    with pytest.raises(ValueError):
        hill_occupancy(-1.0, 1e-6)
    with pytest.raises(ValueError):
        hill_occupancy(1e-6, 0.0)
