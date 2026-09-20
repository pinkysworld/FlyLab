"""Combination pharmacology: Bliss, Loewe, isobolograms, published validation."""
import math

import pytest

from flylab.pharm.mixtures import (
    BLISS_TOL,
    LOEWE_TOL,
    bliss_expected,
    bliss_expected_many,
    inverse_hill,
    isobologram,
    loewe_index,
    mixture_assay,
    mixture_occupancy,
    validate_against_published,
)
from flylab.pharm.occupancy import hill_occupancy


def test_bliss_independence_formula():
    assert bliss_expected(0.0, 0.0) == 0.0
    assert bliss_expected(0.5, 0.5) == pytest.approx(0.75)
    assert bliss_expected(1.0, 0.3) == pytest.approx(1.0)
    assert bliss_expected_many([0.5, 0.5]) == pytest.approx(bliss_expected(0.5, 0.5))
    with pytest.raises(ValueError):
        bliss_expected(1.5, 0.1)


def test_inverse_hill_round_trips():
    conc = inverse_hill(0.5, 2e-8, 1.2)
    assert conc == pytest.approx(2e-8)
    assert hill_occupancy(inverse_hill(0.8, 3e-8, 1.0), 3e-8, 1.0) == pytest.approx(0.8)


def test_loewe_index_is_one_for_a_drug_mixed_with_itself():
    """A drug cannot interact with itself: CI = 1 by construction."""
    ec50 = 2e-8
    total = 4e-8
    effect = hill_occupancy(total, ec50, 1.0)
    ci = loewe_index([total / 2, total / 2], [ec50, ec50], effect, [1.0, 1.0])
    assert ci == pytest.approx(1.0, abs=1e-9)


def test_two_same_target_agonists_are_loewe_additive():
    """imidacloprid + clothianidin: the published honeybee no-interaction case."""
    comps = [
        {"compound": "imidacloprid", "conc_M": 1e-8},
        {"compound": "clothianidin", "conc_M": 1e-8},
    ]
    mix = mixture_occupancy(comps)
    row = next(r for r in mix["receptors"] if r["receptor"] == "insect_nAChR")
    assert row["model"] == "loewe_competitive"
    assert "insect_nAChR" in mix["shared_receptors"]
    # each component holds less of the site than it would alone: they compete
    for part in row["components"]:
        assert part["occupancy_in_mixture"] < part["occupancy_alone"]
    # and the site total is the Loewe-additive value
    a = sum((p["conc_M"] / p["ec50_M"]) ** p["n"] for p in row["components"])
    assert row["occupancy"] == pytest.approx(a / (1.0 + a))
    ci = loewe_index(
        [p["conc_M"] for p in row["components"]],
        [p["ec50_M"] for p in row["components"]],
        row["occupancy"],
        [p["n"] for p in row["components"]],
    )
    assert abs(ci - 1.0) <= LOEWE_TOL  # no interaction, within the stated band


def test_competitive_occupancy_module_agrees_with_the_mixture_maths():
    """With n = 1 the mixture maths IS binding.competitive_occupancy."""
    from flylab.pharm.binding import competitive_occupancy
    from flylab.pharm.occupancy import load_library

    lib = load_library()
    for key in ("imidacloprid", "clothianidin"):
        lib["compounds"][key]["receptors"]["insect_nAChR"]["n"] = 1.0
    comps = [
        {"compound": "imidacloprid", "conc_M": 3e-8},
        {"compound": "clothianidin", "conc_M": 3e-8},
    ]
    row = next(
        r for r in mixture_occupancy(comps, library=lib)["receptors"]
        if r["receptor"] == "insect_nAChR"
    )
    first, second = row["components"]
    expected = float(
        competitive_occupancy(
            first["conc_M"], first["ec50_M"], second["conc_M"], second["ec50_M"], n=1.0
        )
    )
    assert first["occupancy_in_mixture"] == pytest.approx(expected)
    # and with n = 1 the site total is exactly Loewe-additive: CI = 1
    ci = loewe_index(
        [c["conc_M"] for c in row["components"]],
        [c["ec50_M"] for c in row["components"]],
        row["occupancy"],
        [1.0, 1.0],
    )
    assert ci == pytest.approx(1.0, abs=1e-9)


def test_disjoint_targets_stay_independent():
    comps = [
        {"compound": "imidacloprid", "conc_M": 1e-9},
        {"compound": "deltamethrin", "conc_M": 1e-9},
    ]
    mix = mixture_occupancy(comps)
    assert mix["shared_receptors"] == []
    by = {r["receptor"]: r for r in mix["receptors"]}
    assert by["insect_nAChR"]["model"] == "independent_single_agent"
    assert by["insect_Nav"]["model"] == "independent_single_agent"
    for row in mix["receptors"]:
        assert row["occupancy"] == pytest.approx(row["components"][0]["occupancy_alone"])
    assert any("Bliss independence" in w for w in mix["warnings"])


def test_bliss_and_loewe_agree_for_disjoint_targets():
    """Disjoint targets: Bliss says additive and Loewe has nothing to say."""
    nb = mixture_assay(
        [{"compound": "imidacloprid", "conc_M": 1e-9},
         {"compound": "deltamethrin", "conc_M": 1e-9}],
        model="bliss",
    )
    assert nb["assay"] == "mixture"
    bliss = nb["mixture"]["bliss"]
    assert abs(bliss["observed_effect_fraction"] - bliss["expected_effect_fraction"]) <= BLISS_TOL
    assert nb["synergy"]["verdict"] == "additive"
    assert nb["mixture"]["loewe"]["verdict"] == "not_applicable"
    assert nb["mixture"]["tolerances"]["bliss_abs_effect"] == BLISS_TOL


def test_mixture_assay_reports_gains_and_notebook_shape():
    nb = mixture_assay(
        [{"compound": "imidacloprid", "conc_M": 1e-9},
         {"compound": "clothianidin", "conc_M": 1e-9}],
    )
    assert nb["compound"] == "imidacloprid + clothianidin"
    assert set(nb["gains"]) >= {"g_ach", "g_gaba"}
    assert len(nb["readouts"]["singles"]) == 2
    assert nb["readouts"]["vehicle"]["mn9_hz"] is not None
    assert nb["synergy"]["verdict"] in {"synergistic", "additive", "antagonistic", "inconclusive"}
    assert nb["synergy"]["why"]


def test_isobologram_returns_points_and_the_additivity_line():
    iso = isobologram("imidacloprid", "clothianidin", n=5)
    assert iso["d_a"] == pytest.approx(2.0e-8, rel=1e-3)
    assert iso["d_b"] == pytest.approx(3.0e-8, rel=1e-3)
    line = iso["additivity_line"]
    assert len(line) == 2
    assert line[0] == {"conc_a_M": iso["d_a"], "conc_b_M": 0.0}
    assert line[1] == {"conc_a_M": 0.0, "conc_b_M": iso["d_b"]}
    assert len(iso["points"]) == 5
    for point in iso["points"]:
        assert point["effect"] == pytest.approx(0.5, abs=1e-3)
        # the isobole sits on the additivity line within the stated band
        assert abs(point["combination_index"] - 1.0) <= LOEWE_TOL
    assert iso["points"][0]["conc_b_M"] == 0.0
    assert iso["points"][-1]["conc_a_M"] == 0.0


def test_unknown_compound_is_refused_not_substituted():
    with pytest.raises(KeyError) as excinfo:
        mixture_occupancy([{"compound": "acephate", "conc_M": 1e-8}])
    assert "not in the library" in str(excinfo.value)


def test_validate_against_published_reports_each_honeybee_pair():
    result = validate_against_published()
    assert result["summary"]["n_targets"] == 3
    by = {r["partner"].split(" (")[0]: r for r in result["results"]}
    assert set(by) == {"clothianidin", "lambda-cyhalothrin", "acephate"}

    clo = by["clothianidin"]
    assert clo["substituted"] is None
    assert clo["null_model"] == "loewe"  # same target -> Loewe
    assert clo["published"] == "no_interaction"
    assert clo["model_verdict"] == "additive"
    assert clo["matches"] is True

    for name in ("lambda-cyhalothrin", "acephate"):
        row = by[name]
        assert row["substituted"]["requested"] == name
        assert row["null_model"] == "bliss"  # disjoint targets -> Bliss
        assert row["model_verdict"] in {"additive", "synergistic", "antagonistic", "inconclusive"}
    assert any("SUBSTITUTION" in w for w in result["warnings"])
    assert result["summary"]["n_matching"] == result["summary"]["n_evaluated"]


def test_validate_without_substitution_skips_and_names_the_missing_compound():
    result = validate_against_published(substitute=False)
    skipped = [r for r in result["results"] if r["skipped"]]
    assert {r["missing_compounds"][0] for r in skipped} == {"lambda-cyhalothrin", "acephate"}
    assert result["summary"]["n_skipped"] == 2
    assert result["summary"]["n_evaluated"] == 1


def test_every_result_names_its_null_model():
    nb = mixture_assay(
        [{"compound": "imidacloprid", "conc_M": 1e-9},
         {"compound": "clothianidin", "conc_M": 1e-9}],
        model="loewe",
    )
    assert nb["mixture"]["null_model"] == "loewe"
    assert any("null model" in w for w in nb["warnings"])
    assert all(math.isfinite(v) for v in nb["gains"].values())
