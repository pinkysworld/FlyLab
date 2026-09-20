"""Retrospective validation against published orderings, and the two known misses."""
import math

import pytest

from flylab.pharm.mechanisms import gains_from_occupancy
from flylab.pharm.occupancy import compare_compound, load_library
from flylab.validation import (
    KNOWN_DISCREPANCIES,
    kendall_tau,
    load_rank_orders,
    spearman_rho,
    validate_all,
    validate_entry,
    validation_report_markdown,
)


# --------------------------------------------------------------------------
# rank statistics
# --------------------------------------------------------------------------
def test_spearman_and_kendall_without_scipy():
    a = [1, 2, 3, 4, 5]
    assert spearman_rho(a, a) == pytest.approx(1.0)
    assert spearman_rho(a, list(reversed(a))) == pytest.approx(-1.0)
    assert kendall_tau(a, a) == pytest.approx(1.0)
    assert kendall_tau(a, list(reversed(a))) == pytest.approx(-1.0)
    # one swapped adjacent pair: rho = 1 - 6*2/(5*24) = 0.9
    assert spearman_rho(a, [2, 1, 3, 4, 5]) == pytest.approx(0.9)
    assert kendall_tau(a, [2, 1, 3, 4, 5]) == pytest.approx(0.8)
    # ties are handled, not crashed on
    assert math.isnan(spearman_rho([1, 1, 1], [1, 2, 3]))
    assert spearman_rho([1, 1, 2], [1, 2, 3]) == pytest.approx(0.8660254, abs=1e-6)


# --------------------------------------------------------------------------
# entries
# --------------------------------------------------------------------------
def test_dataset_loads():
    data = load_rank_orders()
    assert len(data["entries"]) >= 10
    assert data["not_sourced"]


def test_validate_entry_returns_a_finite_rho():
    result = validate_entry("neonic_potency_dmel_larval")
    assert result["skipped"] is False
    assert math.isfinite(result["spearman_rho"])
    assert math.isfinite(result["kendall_tau"])
    assert result["receptor"] == "insect_nAChR"
    assert result["pairwise"]
    assert result["n_concordant"] + result["n_discordant"] == len(result["pairwise"])


def test_alpha7_entry_uses_the_vertebrate_receptor_named_by_the_assay():
    result = validate_entry("neonic_alpha7_partial_agonism")
    if result.get("skipped"):
        # schema v3: thiamethoxam's alpha7 row is a placeholder (the source says
        # "no agonist effect", which is not a number), so the model cannot score
        # it. The entry is skipped by name instead of being ranked off a fake
        # ec50_M of 0.01.
        assert "vertebrate_nAChR_a7" in result["reason"]
        return
    assert result["receptor"] == "vertebrate_nAChR_a7"
    # thiamethoxam has no alpha7 agonist effect and must come last in both
    last = max(result["compounds"], key=lambda c: c["model_rank"])
    assert last["key"] == "thiamethoxam"
    # the EC50 ordering reproduces the publication exactly; the occupancy
    # ordering at 1e-6 M does not, and the result says why (Hill coefficients)
    assert result["ec50_ordering"]["spearman_rho_by_ec50"] == pytest.approx(1.0)
    if not result["ec50_ordering"]["agrees_with_occupancy_ordering"]:
        assert "Hill coefficient" in result["ec50_ordering"]["note"]


def test_entries_the_library_cannot_cover_are_skipped_with_names():
    result = validate_entry("benzodiazepine_inactive_at_insect_rdl")
    assert result["skipped"] is True
    assert "flunitrazepam" in result["missing_compounds"]

    unsourced = validate_entry("rdl_antagonist_order_drosophila")
    assert unsourced["skipped"] is True
    assert "no ranks or values" in unsourced["reason"]


def test_rank_order_only_entry_with_null_values_is_handled():
    result = validate_entry("neonic_insect_vs_vertebrate_nachr")
    assert result["skipped"] is False
    assert result["order_source"] == "direction_selectivity"
    assert math.isfinite(result["spearman_rho"])
    # nicotine is published as vertebrate-selective; the library makes it
    # (weakly) insect-selective, and that sign conflict is surfaced
    assert [c["compound"] for c in result["direction_conflicts"]] == ["nicotine"]


# --------------------------------------------------------------------------
# whole dataset
# --------------------------------------------------------------------------
def test_validate_all_summarises_both_assays():
    result = validate_all(assays=("occupancy", "subgraph"))
    summary = result["summary"]
    assert summary["n_entries"] == len(load_rank_orders()["entries"])
    assert summary["n_evaluated"] >= 1
    assert summary["n_evaluated"] + summary["n_skipped"] == summary["n_runs"]
    assert math.isfinite(summary["mean_rho"])
    occ = result["assays"]["occupancy"]
    assert occ["n_evaluated"] >= 4
    assert all(math.isfinite(v) for v in occ["rho_by_entry"].values())


def test_both_known_discrepancies_are_surfaced():
    block = validate_all(assays=("occupancy",))["known_discrepancies"]
    by_id = {d["id"]: d for d in block["discrepancies"]}
    assert set(by_id) == {d["id"] for d in KNOWN_DISCREPANCIES}

    nitenpyram = by_id["nitenpyram_rank_inverted"]
    assert nitenpyram["confirmed_present"] is True
    assert nitenpyram["model_rank_of_nitenpyram"] < nitenpyram["published_rank_of_nitenpyram"]
    assert "Perry" in nitenpyram["published_source"]

    clothianidin = by_id["clothianidin_vs_imidacloprid_inverted"]
    assert clothianidin["confirmed_present"] is True
    order = clothianidin["model_order_by_insect_nAChR_ec50"]
    assert order.index("imidacloprid") < order.index("clothianidin")

    # the explanation is stated, and the library is not tuned to hide it
    assert "DIFFERENT QUANTITIES" in block["explanation"]
    assert "not tuned" in block["policy"].lower()


def test_report_markdown_has_a_table_and_the_discrepancies():
    md = validation_report_markdown(assays=("occupancy",))
    assert "| entry |" in md
    assert "neonic_potency_dmel_larval" in md
    assert "nitenpyram_rank_inverted" in md
    assert "clothianidin_vs_imidacloprid_inverted" in md
    assert "mean Spearman rho" in md


# --------------------------------------------------------------------------
# the v0.5 library correction
# --------------------------------------------------------------------------
def test_fipronil_vertebrate_correction_is_in_compare_compound():
    """Ratra & Casida 2001 report 1103 nM, not the 1.0e-5 M v0.4 carried."""
    rows = {r["receptor"]: r for r in compare_compound("fipronil", 1e-6)["receptors"]}
    vert = rows["vertebrate_GABA_A"]
    assert vert["ec50_M"] == pytest.approx(1.1e-6)
    assert "1103" in vert["source"]
    assert "beta3" in vert["source"]  # the homopentamer caveat is stated
    assert vert["evidence_tier"] == "literature_order"
    # the insect row keeps its value and now records the 1500-fold spread
    insect = rows["insect_RDL"]
    assert insect["ec50_M"] == pytest.approx(3.0e-8)
    assert "Lees" in insect["source"] and "house-fly" in insect["source"]


def test_fipronil_selectivity_direction_and_gains_survive_the_correction():
    result = compare_compound("fipronil", 1e-6)
    pair = result["selectivity"]["GABA_A"]
    assert pair["ec50_ratio_vert_over_insect"] > 1.0, "fipronil must stay insect-selective"
    assert pair["insect_occupancy"] > pair["vertebrate_occupancy"]
    gains = gains_from_occupancy(result["receptors"])
    assert gains["g_gaba"] < 1.0
    assert gains["g_ach"] == pytest.approx(1.0)


def test_library_frozen_v04_values_are_unchanged_apart_from_that_row():
    lib = load_library()
    receptors = lib["compounds"]["imidacloprid"]["receptors"]
    assert float(receptors["insect_nAChR"]["ec50_M"]) == pytest.approx(2.0e-8)
    assert float(receptors["vertebrate_nAChR_a4b2"]["ec50_M"]) == pytest.approx(1.0e-5)
