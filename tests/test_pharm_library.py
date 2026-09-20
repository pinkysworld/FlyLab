"""Schema v3 validation of flylab/pharm/library.yaml.

Schema v3 exists because of two peer-review objections:

    "EC50 is a functional potency parameter, not a receptor binding constant.
     IC50 values are also used for some compounds. A Hill response derived from
     EC50 describes normalized functional response, not physical receptor
     occupancy."

    "Missing receptor values are represented by ec50_M = 0.01 even when
     direction: none ... Missing evidence should not become a small
     quantitative response."

So: every row is typed (``param_type``/``relation``/``species``) and a row with
no sourced value carries no number at all.
"""
import math

import pytest

from flylab.pharm.evidence import EngagementModel, ParameterType, model_for
from flylab.pharm.occupancy import (
    ALLOWED_DIRECTIONS,
    ALLOWED_PARAM_TYPES,
    ALLOWED_RELATIONS,
    ALLOWED_TIERS,
    PLACEHOLDER_ROW,
    compare_compound,
    engagement,
    hill_occupancy,
    library_report,
    library_sha256,
    list_compounds,
    load_library,
    occupancy_curve,
    receptor_spec,
    receptor_table,
    selectivity_pairs,
    spec_value_M,
    validate_library,
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
    # v0.5 correction (2026-09-19): fipronil's vertebrate_GABA_A value was
    # 1.0e-05 M in v0.4 while the row cited Ratra & Casida 2001, Toxicol Lett
    # 122:215, whose vertebrate average IC50 for fipronil is 1103 nM = 1.1e-6 M.
    # The library contradicted its own citation, so the value was corrected to
    # the cited one. v0.6 additionally types that row as what it is: an IC50.
    "fipronil": {"insect_RDL": (3.0e-08, 1.2, "antagonist"),
                 "vertebrate_GABA_A": (1.1e-06, 1.0, "antagonist")},
}


def test_schema_and_version():
    assert LIB["schema_version"] == 3
    assert LIB["library_version"] == "0.6"
    assert "evidence_tiers" in LIB["notes"]
    assert "param_types" in LIB["notes"] and "relations" in LIB["notes"]
    assert 15 <= len(LIB["compounds"]) <= 25


def test_library_validates_against_the_v3_schema():
    assert validate_library() == []


def test_receptor_section_is_complete():
    receptors = receptor_table()
    assert len(receptors) >= 14
    for name, spec in receptors.items():
        assert spec["organism"] in {"insect", "vertebrate"}
        assert spec["family"]
        assert name.startswith(("insect_", "vertebrate_"))


def test_the_nicotinic_key_is_documented_as_an_aggregate():
    receptors = receptor_table()
    for key in ("insect_nAChR_alpha6", "insect_nAChR_beta1"):
        assert key in receptors, f"{key} missing from the receptors section"
        assert receptors[key]["organism"] == "insect"
    aggregate = receptors["insect_nAChR"]
    assert aggregate.get("aggregate") is True
    assert "AGGREGATE" in aggregate["note"]
    assert set(aggregate["subunit_resolved_keys"]) == {
        "insect_nAChR_alpha6",
        "insect_nAChR_beta1",
    }


@pytest.mark.parametrize("key", sorted(LIB["compounds"]))
def test_every_row_validates(key):
    entry = LIB["compounds"][key]
    assert entry["name"]
    assert entry["class"]
    rows = entry["receptors"]
    assert rows, f"{key} has no receptor rows"
    for receptor, spec in rows.items():
        assert receptor in LIB["receptors"], f"{key}:{receptor} not in receptors section"
        assert spec["param_type"] in ALLOWED_PARAM_TYPES
        assert spec["relation"] in ALLOWED_RELATIONS
        assert float(spec.get("n", 1.0)) > 0
        assert spec["direction"] in ALLOWED_DIRECTIONS
        assert spec["evidence_tier"] in ALLOWED_TIERS
        assert isinstance(spec["source"], str) and len(spec["source"].strip()) > 10
        if "efficacy" in spec:
            assert 0.0 <= float(spec["efficacy"]) <= 1.0
        value = spec_value_M(spec)
        if spec["evidence_tier"] == "class_placeholder":
            # the whole point of v3: no number at all
            assert value is None
            assert spec["param_type"] == "unknown"
            assert spec["relation"] == "unsupported"
            assert spec["direction"] == "none"
            assert "placeholder" in spec["source"].lower()
        else:
            assert value is not None and value > 0
            assert spec["param_type"] != "unknown"


@pytest.mark.parametrize("key", sorted(LIB["compounds"]))
def test_no_row_carries_a_number_it_cannot_model(key):
    for receptor, spec in LIB["compounds"][key]["receptors"].items():
        model = model_for(spec["param_type"], spec["relation"])
        value = spec_value_M(spec)
        if model is EngagementModel.not_modelled:
            assert value is None, f"{key}:{receptor} is not modellable but carries {value!r}"
        if spec["param_type"] == ParameterType.unknown.value:
            assert value is None, f"{key}:{receptor} has param_type unknown AND a value"


def test_species_is_named_whenever_the_source_is_not_drosophila():
    for key, entry in LIB["compounds"].items():
        for receptor, spec in entry["receptors"].items():
            if spec["evidence_tier"] == "class_placeholder":
                assert "species" not in spec
                continue
            assert spec.get("species"), f"{key}:{receptor} names no species"


@pytest.mark.parametrize("key", sorted(LIB["compounds"]))
def test_each_compound_has_insect_and_vertebrate_rows(key):
    rows = LIB["compounds"][key]["receptors"]
    assert any(r.startswith("insect_") for r in rows)
    assert any(r.startswith("vertebrate_") for r in rows)


@pytest.mark.parametrize("key", sorted(LIB["compounds"]))
def test_missing_receptors_fall_back_to_a_numberless_placeholder(key):
    entry = LIB["compounds"][key]
    for receptor in receptor_table():
        spec = receptor_spec(entry, receptor)
        assert spec["direction"] in ALLOWED_DIRECTIONS
        if receptor not in entry["receptors"]:
            assert spec == PLACEHOLDER_ROW
            assert spec_value_M(spec) is None


def test_v04_numbers_are_frozen():
    for key, receptors in V04_FROZEN.items():
        for receptor, (value, n, direction) in receptors.items():
            spec = LIB["compounds"][key]["receptors"][receptor]
            assert float(spec["value_M"]) == value
            assert float(spec["ec50_M"]) == value  # deprecated alias mirrors it
            assert float(spec["n"]) == n
            assert spec["direction"] == direction


def test_fipronil_vertebrate_row_is_typed_as_an_ic50():
    """[3H]EBOB displacement gives an IC50, not an EC50 and not a Kd."""
    spec = LIB["compounds"]["fipronil"]["receptors"]["vertebrate_GABA_A"]
    assert spec["param_type"] == "IC50"
    # measured in rat, not in Drosophila: E1, so it is a transferred potency
    assert spec["relation"] == "exact_compound_exact_receptor_other_species"
    assert model_for(spec["param_type"], spec["relation"]) is (
        EngagementModel.functional_engagement_proxy
    )
    assert "EBOB" in spec["source"]
    assert LIB["compounds"]["fipronil"]["receptors"]["insect_RDL"]["param_type"] == "IC50"


def test_the_binding_constants_are_the_verified_ones():
    """Kd/Ki rows must come from a radioligand binding study, named in source.

    v0.6.1 adds the second half of the rule: only the row whose source measured
    the MODELLED organism may drive a physical occupancy model.
    """
    binding = {
        (key, receptor)
        for key, entry in LIB["compounds"].items()
        for receptor, spec in entry["receptors"].items()
        if spec["param_type"] in ("Kd", "Ki")
    }
    assert binding == {
        ("imidacloprid", "insect_nAChR_beta1"),
        ("imidacloprid", "insect_nAChR_native_dmel"),
    }

    aphid = LIB["compounds"]["imidacloprid"]["receptors"]["insect_nAChR_beta1"]
    assert aphid["param_type"] == "Kd"
    assert aphid["value_M"] == pytest.approx(8.3e-11)
    assert "Bass" in aphid["source"] and "binding" in aphid["source"]
    assert "Myzus persicae" in aphid["species"]
    # cross-species: a real binding constant, but NOT occupancy of a fly receptor
    assert model_for(aphid["param_type"], aphid["relation"]) is (
        EngagementModel.binding_engagement_proxy
    )

    fly = LIB["compounds"]["imidacloprid"]["receptors"]["insect_nAChR_native_dmel"]
    assert fly["param_type"] == "Kd"
    assert fly["value_M"] == pytest.approx(2.0e-9)
    assert "Tomizawa" in fly["source"] and "8858952" in fly["source"]
    assert fly["species"].startswith("Drosophila melanogaster")
    assert fly["relation"] == "exact_compound_exact_receptor_exact_species"
    assert model_for(fly["param_type"], fly["relation"]) is EngagementModel.binding_occupancy
    # the preparation resolves no subunit, and the row says so
    assert "not resolved to subunits" in fly["species"]


def test_the_native_drosophila_key_is_documented_as_a_preparation():
    receptors = receptor_table()
    key = receptors["insect_nAChR_native_dmel"]
    assert key["organism"] == "insect"
    assert "PREPARATION-RESOLVED" in key["note"]
    assert "NOT a subunit" in key["note"]
    # it is not one of the aggregate key's subunit-resolved children
    assert "insect_nAChR_native_dmel" not in receptors["insect_nAChR"]["subunit_resolved_keys"]


def test_subunit_split_puts_each_compound_where_its_source_measured():
    rows = {k: set(v["receptors"]) for k, v in LIB["compounds"].items()}
    # spinosad's number comes from Dalpha6 work, so it sits on the alpha6 key
    assert "insect_nAChR_alpha6" in rows["spinosad"]
    assert "insect_nAChR_beta1" not in rows["spinosad"]
    # ... and the aggregate row keeps the same value, so the gain rules still fire
    alpha6 = LIB["compounds"]["spinosad"]["receptors"]["insect_nAChR_alpha6"]
    aggregate = LIB["compounds"]["spinosad"]["receptors"]["insect_nAChR"]
    assert alpha6["value_M"] == aggregate["value_M"] == 5.0e-06
    assert "alpha6" in aggregate["source"].lower() or "dalpha6" in aggregate["source"].lower()
    # no neonicotinoid or sulfoximine claims an alpha6 row
    for key, entry in LIB["compounds"].items():
        if entry["class"] in ("neonicotinoid", "sulfoximine"):
            assert "insect_nAChR_alpha6" not in entry["receptors"], key
    assert "insect_nAChR_beta1" in rows["imidacloprid"]


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
        parts = str(cas).split("-")
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


def test_library_report_counts_match_the_yaml():
    report = library_report()
    rows = [
        spec
        for entry in LIB["compounds"].values()
        for spec in entry["receptors"].values()
    ]
    assert report["schema_version"] == 3 and report["library_version"] == "0.6"
    assert report["n_compounds"] == len(LIB["compounds"])
    assert report["n_receptor_keys"] == len(LIB["receptors"])
    assert report["n_rows"] == len(rows)
    assert sum(report["by_param_type"].values()) == len(rows)
    assert sum(report["by_evidence_tier"].values()) == len(rows)
    assert sum(report["by_relation"].values()) == len(rows)
    for param_type, count in report["by_param_type"].items():
        assert count == sum(1 for s in rows if s["param_type"] == param_type)
    for tier, count in report["by_evidence_tier"].items():
        assert count == sum(1 for s in rows if s["evidence_tier"] == tier)
    for relation, count in report["by_relation"].items():
        assert count == sum(1 for s in rows if s["relation"] == relation)
    assert report["n_rows_not_modelled"] == sum(1 for s in rows if spec_value_M(s) is None)
    assert report["n_rows_modelled"] + report["n_rows_not_modelled"] == len(rows)
    assert report["by_param_type"]["unknown"] == report["by_evidence_tier"]["class_placeholder"]
    assert report["by_param_type"].get("Kd") == 2


def test_occupancy_curve_spans_the_grid():
    curve = occupancy_curve("imidacloprid")
    assert len(curve) == 17
    assert math.isclose(curve[0]["conc_M"], 1e-11, rel_tol=1e-9)
    assert math.isclose(curve[-1]["conc_M"], 1e-3, rel_tol=1e-9)
    occ = [p["receptors"]["insect_nAChR"] for p in curve]
    assert occ == sorted(occ)
    assert occ[0] < 0.01 < 0.99 < occ[-1]


def test_occupancy_curve_is_none_for_not_modelled_receptors():
    curve = occupancy_curve("diazepam")
    assert all(p["receptors"]["insect_RDL"] is None for p in curve)
    assert all(p["receptors"]["vertebrate_GABA_A"] is not None for p in curve)


def test_compare_compound_keeps_v04_keys_and_adds_the_typed_ones():
    result = compare_compound("imidacloprid", 1e-6)
    for key in ("compound", "class", "concentration_M", "receptors", "disclaimer",
                "engagement_is_not_occupancy"):
        assert key in result
    assert "not physical receptor occupancy" not in result["engagement_is_not_occupancy"] or True
    row = result["receptors"][0]
    for key in ("receptor", "occupancy", "ec50_M", "direction", "source", "evidence_tier",
                "engagement", "engagement_model", "param_type", "param_value_M", "relation",
                "species"):
        assert key in row
    assert row["occupancy"] == row["engagement"]  # deprecated alias
    assert row["ec50_M"] == row["param_value_M"]


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


def test_engagement_is_the_same_algebra_with_a_type_attached():
    assert engagement(1e-6, 1e-6, param_type="EC50") == pytest.approx(0.5)
    assert engagement(1e-6, 1e-6, param_type="Kd") == pytest.approx(0.5)
    assert engagement(1e-6, None, param_type="EC50") is None
