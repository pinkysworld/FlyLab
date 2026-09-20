"""Selectivity scorecard and the insect/vertebrate contrasts it must show."""
import pytest

from flylab.pharm.occupancy import compare_compound


def by_receptor(result):
    return {row["receptor"]: row for row in result["receptors"]}


def test_imidacloprid_still_prefers_the_insect_receptor():
    result = compare_compound("imidacloprid", 1e-6)
    rows = by_receptor(result)
    assert rows["insect_nAChR"]["occupancy"] > 0.9
    assert rows["vertebrate_nAChR_a4b2"]["occupancy"] < 0.2
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
        ):
            assert key in pair


def test_picrotoxin_is_a_non_selective_control():
    pair = compare_compound("picrotoxin", 1e-6)["selectivity"]["GABA_A"]
    assert pair["ec50_ratio_vert_over_insect"] == pytest.approx(1.0)
    assert abs(pair["log10_ec50_ratio_vert_over_insect"]) < 0.1
    assert abs(pair["occupancy_difference"]) < 0.05
    assert not pair["placeholder"]


def test_selective_insecticides_beat_picrotoxin():
    pic = compare_compound("picrotoxin", 1e-6)["selectivity"]["GABA_A"]
    for key in ("imidacloprid", "deltamethrin", "ivermectin"):
        result = compare_compound(key, 1e-6)["selectivity"]
        best = max(p["log10_ec50_ratio_vert_over_insect"] for p in result.values())
        assert best > pic["log10_ec50_ratio_vert_over_insect"] + 0.5


def test_pairs_missing_from_a_compound_are_flagged_placeholder():
    pair = compare_compound("imidacloprid", 1e-6)["selectivity"]["Nav"]
    assert pair["placeholder"] is True
    assert pair["evidence_tier"] == "class_placeholder"
    assert pair["ec50_ratio_vert_over_insect"] == pytest.approx(1.0)


def test_evidence_tier_on_every_row():
    for row in compare_compound("ivermectin", 1e-7)["receptors"]:
        assert row["evidence_tier"] in {"literature_order", "class_placeholder", "measured_fit"}
        assert row["source"]


def test_unknown_compound_raises_keyerror():
    with pytest.raises(KeyError):
        compare_compound("unobtainium", 1e-6)
