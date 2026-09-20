"""Selectivity scorecard and the insect/vertebrate contrasts it must show.

Schema v3: a pair with a placeholder on either side has NO ratio ("Missing
evidence should not become a small quantitative response"), and the numbers are
engagements of typed parameters, not occupancies.
"""
import pytest

from flylab.pharm.occupancy import compare_compound


def by_receptor(result):
    return {row["receptor"]: row for row in result["receptors"]}


def test_imidacloprid_still_prefers_the_insect_receptor():
    result = compare_compound("imidacloprid", 1e-6)
    rows = by_receptor(result)
    assert rows["insect_nAChR"]["engagement"] > 0.9
    assert rows["vertebrate_nAChR_a4b2"]["engagement"] < 0.15
    # both are functional potencies, so neither is called an occupancy
    assert rows["insect_nAChR"]["engagement_model"] == "functional_engagement"
    assert rows["insect_nAChR"]["param_type"] == "EC50"
    # the subunit-resolved beta1 row is the one genuine binding constant
    beta1 = rows["insect_nAChR_beta1"]
    assert beta1["param_type"] == "Kd" and beta1["engagement_model"] == "binding_occupancy"
    assert beta1["engagement"] > 0.9
    pair = result["selectivity"]["nAChR"]
    assert pair["ec50_ratio_vert_over_insect"] == pytest.approx(500.0)
    assert pair["log10_ec50_ratio_vert_over_insect"] == pytest.approx(2.699, abs=1e-3)
    assert pair["occupancy_difference"] > 0.8


def test_selectivity_block_has_every_pair():
    result = compare_compound("imidacloprid", 1e-6)
    assert set(result["selectivity"]) == {"nAChR", "GABA_A", "GluCl", "AChE", "Nav"}
    for pair in result["selectivity"].values():
        for key in (
            "ec50_ratio_vert_over_insect",
            "log10_ec50_ratio_vert_over_insect",
            "occupancy_difference",
            "insect_occupancy",
            "vertebrate_occupancy",
            "evidence_tier",
            "ratio",
            "reason",
            "insect_param_type",
            "vertebrate_param_type",
            "comparable",
        ):
            assert key in pair


def test_picrotoxin_is_a_non_selective_control_on_its_sourced_pair():
    pair = compare_compound("picrotoxin", 1e-6)["selectivity"]["GABA_A"]
    assert pair["ratio"] == pytest.approx(1.0)
    assert pair["ec50_ratio_vert_over_insect"] == pytest.approx(1.0)
    assert abs(pair["log10_ec50_ratio_vert_over_insect"]) < 0.1
    assert abs(pair["occupancy_difference"]) < 0.05
    assert not pair["placeholder"]
    assert pair["reason"] is None
    # both sides are IC50s, so the ratio compares like with like
    assert pair["insect_param_type"] == pair["vertebrate_param_type"] == "IC50"
    assert pair["comparable"] is True


def test_selective_insecticides_beat_picrotoxin():
    pic = compare_compound("picrotoxin", 1e-6)["selectivity"]["GABA_A"]
    for key in ("imidacloprid", "deltamethrin", "ivermectin"):
        result = compare_compound(key, 1e-6)["selectivity"]
        # placeholder pairs carry no ratio and take no part in the comparison
        ratios = [
            p["log10_ec50_ratio_vert_over_insect"]
            for p in result.values()
            if p["log10_ec50_ratio_vert_over_insect"] is not None
        ]
        assert ratios, key
        assert max(ratios) > pic["log10_ec50_ratio_vert_over_insect"] + 0.5


def test_pairs_missing_from_a_compound_have_no_ratio_at_all():
    """The v0.5 bug: a placeholder pair used to come out as a ratio of 1.00."""
    pair = compare_compound("imidacloprid", 1e-6)["selectivity"]["Nav"]
    assert pair["placeholder"] is True
    assert pair["evidence_tier"] == "class_placeholder"
    assert pair["ratio"] is None and pair["reason"] == "placeholder"
    assert pair["ec50_ratio_vert_over_insect"] is None
    assert pair["log10_ec50_ratio_vert_over_insect"] is None
    assert pair["occupancy_difference"] is None
    assert pair["insect_occupancy"] is None and pair["vertebrate_occupancy"] is None


def test_a_pair_says_when_it_compares_two_different_kinds_of_parameter():
    for key in ("imidacloprid", "picrotoxin", "fipronil", "ivermectin"):
        for pair in compare_compound(key, 1e-6)["selectivity"].values():
            if pair["ratio"] is None:
                continue
            same = pair["insect_param_type"] == pair["vertebrate_param_type"]
            assert pair["comparable"] is same
            assert (pair["reason"] is None) is same


def test_evidence_tier_on_every_row():
    for row in compare_compound("ivermectin", 1e-7)["receptors"]:
        assert row["evidence_tier"] in {"literature_order", "class_placeholder", "measured_fit"}
        assert row["source"]


def test_the_result_explains_that_engagement_is_not_occupancy():
    result = compare_compound("fipronil", 1e-6)
    note = result["engagement_is_not_occupancy"]
    assert "occupancy" in note and "binding_occupancy" in note
    assert result["n_receptors_not_modelled"] == 2  # both nAChR rows
    assert result["schema_version"] == 3


def test_unknown_compound_raises_keyerror():
    with pytest.raises(KeyError):
        compare_compound("unobtainium", 1e-6)
