"""The ablation ladder: four levels of model sophistication, and what each one
adds over the level below it."""
from __future__ import annotations

import json

import numpy as np
import pytest

from flylab.analysis.baselines import (
    C_FLOOR_LEVEL,
    DEFAULT_GENERIC_RULE,
    GENERIC_MULTIPLIER_RULE,
    GENERIC_RULES,
    LEVELS,
    NT_SIGN_CONVENTIONS,
    ablation,
    ablation_table,
    baseline_composition_only,
    baseline_receptor_only,
    baseline_topology_only,
    composition_dominance_under_normalisations,
    composition_effect_from_gains,
    composition_reference_distribution,
    full_model,
    generic_multiplier,
    glutamate_sign_reconciliation,
    graph_census,
    level_predictions,
    pearson,
    spearman,
)

PAPER_CONC = 1e-6
#: a small, mixed set: two nicotinic agonists, a GABA blocker, a pyrethroid and
#: a vertebrate-only ligand. Enough for a correlation, cheap enough for CI.
SET = ["imidacloprid", "nicotine", "fipronil", "deltamethrin", "diazepam", "picrotoxin"]


# --------------------------------------------------------------------------
# the four levels
# --------------------------------------------------------------------------
def test_ablation_returns_four_levels_in_order():
    out = ablation("imidacloprid", PAPER_CONC, compounds=SET)
    assert list(out["levels"]) == list(LEVELS) == out["level_order"]
    assert len(LEVELS) == 4
    for lvl in LEVELS:
        assert out["levels"][lvl]["level"] == lvl
        assert out["levels"][lvl]["unit"]
        assert out["effects"][lvl] is not None
        assert np.isfinite(float(out["effects"][lvl]))
    assert out["label"] == "model_derived" and out["warnings"]
    json.dumps(out)


def test_levels_are_what_they_claim_to_be():
    """Each baseline must be blind to the layer it ablates."""
    # A sees no graph at all: the same number on either cut.
    a_named = baseline_receptor_only("fipronil", PAPER_CONC)
    assert a_named["unit"] == "dimensionless"
    # B sees the census and no edges: it is a function of the node counts only.
    b = baseline_composition_only("imidacloprid", PAPER_CONC, graph="named")
    assert b["detail"]["census"] == graph_census("named")
    assert b["detail"]["census"] != graph_census("taste_motor")
    assert baseline_composition_only("imidacloprid", PAPER_CONC, graph="taste_motor")[
        "effect"
    ] != pytest.approx(b["effect"])
    # C sees the graph but only one generic multiplier.
    c = baseline_topology_only("imidacloprid", PAPER_CONC, graph="named")
    assert 0.05 <= c["detail"]["multiplier"] <= 2.0
    assert c["unit"] == "Hz"
    # D is the shipped model, and must equal the assay's own contrast.
    from flylab.analysis.nullmodels import drug_effect

    d = full_model("imidacloprid", PAPER_CONC, graph="named")
    ref = drug_effect("subgraph", "imidacloprid", PAPER_CONC, "mean_hz", graph="named")
    assert d["effect"] == pytest.approx(ref["effect"], abs=1e-9)


def test_vehicle_moves_nothing_at_any_level():
    preds = level_predictions("diazepam", PAPER_CONC)  # vertebrate-only ligand
    for lvl in ("B_composition_only", "C_topology_only", "D_full_flylab"):
        assert preds[lvl]["effect"] == pytest.approx(0.0, abs=1e-9)


def test_not_modelled_receptor_rows_are_skipped_not_read_as_zero():
    """The evidence layer reports None for a row with no sourced value; that
    means 'not modelled', and baseline A must never fold it in as a 0."""
    out = baseline_receptor_only("picrotoxin", PAPER_CONC)
    detail = out["detail"]
    assert detail["n_used"] + detail["n_not_modelled"] <= detail["n_rows"]
    assert all(v is None or np.isfinite(float(v)) for v in detail["per_receptor"].values())
    if detail["n_not_modelled"]:
        assert any("not modellable" in w or "no modellable" in w for w in out["warnings"])


# --------------------------------------------------------------------------
# information gain
# --------------------------------------------------------------------------
def test_information_gain_is_finite_and_ordered():
    out = ablation("imidacloprid", PAPER_CONC, compounds=SET)
    gain = out["information_gain"]
    rows = {r["level"]: r for r in gain["levels"]}
    assert set(rows) == set(LEVELS)
    for lvl, row in rows.items():
        for key in ("pearson_r_vs_full", "spearman_rho_vs_full"):
            val = row[key]
            assert val is None or (np.isfinite(val) and -1.0 <= val <= 1.0)
        assert row["residual_rms_standardised"] is None or np.isfinite(
            row["residual_rms_standardised"]
        )
    # the full model is its own reference
    assert rows["D_full_flylab"]["pearson_r_vs_full"] == pytest.approx(1.0)
    assert rows["D_full_flylab"]["spearman_rho_vs_full"] == pytest.approx(1.0)
    assert rows["D_full_flylab"]["residual_rms_hz"] == pytest.approx(0.0)
    # only the Hz levels get a residual in Hz
    assert rows["A_receptor_only"]["residual_rms_hz"] is None
    assert rows["C_topology_only"]["residual_rms_hz"] is not None
    assert out["statements"]


def test_composition_only_reproduces_the_full_model_for_nicotinic_agonists():
    """The publishable result: for the nicotinic agonists the mean-rate effect
    is an E/I-balance effect, so a census with no edges already orders them
    the way the full model does."""
    nicotinics = [
        "imidacloprid",
        "acetamiprid",
        "clothianidin",
        "nitenpyram",
        "nicotine",
        "thiamethoxam",
        "spinosad",
        "acetylcholine",
        "sulfoxaflor",
    ]
    out = ablation("imidacloprid", PAPER_CONC, compounds=nicotinics)
    rows = {r["level"]: r for r in out["information_gain"]["levels"]}
    b = rows["B_composition_only"]
    assert b["spearman_rho_vs_full"] >= 0.9
    assert b["reproduces_full_ordering"] is True
    assert any("B_composition_only REPRODUCES" in s for s in out["statements"])
    # and it is the composition layer specifically: the receptor scorecard on
    # its own gets the ordering wrong (it has no gain rule).
    assert rows["A_receptor_only"]["reproduces_full_ordering"] is False


def test_receptor_only_does_not_reproduce_the_full_model():
    """A scorecard with no gain rule cannot know that a saturating agonist
    silences its target, so its ordering is not the full model's."""
    out = ablation("imidacloprid", PAPER_CONC, compounds=SET)
    rows = {r["level"]: r for r in out["information_gain"]["levels"]}
    assert rows["A_receptor_only"]["reproduces_full"] is False
    assert any("A_receptor_only does NOT reproduce" in s for s in out["statements"])


def test_information_gain_can_be_switched_off():
    out = ablation("fipronil", PAPER_CONC, include_information_gain=False)
    assert out["information_gain"] is None
    assert out["effects"]["D_full_flylab"] is not None


# --------------------------------------------------------------------------
# correlation helpers
# --------------------------------------------------------------------------
def test_pearson_and_spearman_handle_ties_and_missing_values():
    assert pearson([1, 2, 3], [2, 4, 6]) == pytest.approx(1.0)
    assert spearman([1, 2, 3], [10, 20, 30]) == pytest.approx(1.0)
    assert spearman([1, 2, 3], [30, 20, 10]) == pytest.approx(-1.0)
    assert spearman([1, 1, 2, 3], [5, 5, 7, 9]) == pytest.approx(1.0)
    # pairwise-complete: the None pair is dropped, not treated as 0
    assert spearman([1, None, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert pearson([1, 1, 1], [1, 2, 3]) is None
    assert spearman([1, 2], [1, 2]) is None


# --------------------------------------------------------------------------
# the table
# --------------------------------------------------------------------------
def test_ablation_table_shape_and_statements():
    tab = ablation_table(compounds=SET, concs_M=(1e-7, 1e-6))
    assert len(tab["rows"]) == len(SET) * 2
    assert {r["compound"] for r in tab["rows"]} == set(SET)
    for row in tab["rows"]:
        for lvl in LEVELS:
            assert lvl in row
    assert set(tab["information_gain"]) == {"1.000e-07", "1.000e-06"}
    assert tab["statements"]
    assert tab["label"] == "model_derived" and tab["warnings"]
    json.dumps(tab)


@pytest.mark.slow
def test_ablation_table_over_the_whole_library():
    tab = ablation_table()
    assert len(tab["rows"]) == len(tab["compounds"])
    gain = list(tab["information_gain"].values())[0]
    rows = {r["level"]: r for r in gain["levels"]}
    assert rows["B_composition_only"]["spearman_rho_vs_full"] is not None


# --------------------------------------------------------------------------
# the topology-only baseline used to be sign-broken
# --------------------------------------------------------------------------
def test_generic_rule_is_direction_aware_by_default():
    """The old rule could only depress, so it could not express disinhibition;
    the new one takes the sign - and only the sign - from the mechanism table."""
    assert DEFAULT_GENERIC_RULE == "direction_aware"
    assert GENERIC_RULES["depressant_floor"] == GENERIC_MULTIPLIER_RULE

    # fipronil blocks RDL: the full model disinhibits, so a fair generic
    # baseline must be able to go up.
    m_new, det_new = generic_multiplier("fipronil", PAPER_CONC)
    m_old, det_old = generic_multiplier("fipronil", PAPER_CONC, rule="depressant_floor")
    assert m_new > 1.0 and det_new["direction"] == 1
    assert m_old <= 1.0 and det_old["direction"] == -1
    # same generic magnitude: only the sign is new
    assert det_new["theta_max"] == pytest.approx(det_old["theta_max"])
    assert m_new == pytest.approx(1.0 + det_new["theta_max"])
    assert m_old == pytest.approx(max(0.05, 1.0 - det_old["theta_max"]))

    # a nicotinic agonist at a silencing concentration still goes down
    m_imi, det_imi = generic_multiplier("imidacloprid", PAPER_CONC)
    assert m_imi < 1.0 and det_imi["direction"] == -1


def test_direction_aware_baseline_can_express_disinhibition():
    fip_new = baseline_topology_only("fipronil", PAPER_CONC, graph="named")
    fip_full = full_model("fipronil", PAPER_CONC, graph="named")
    assert fip_full["effect"] > 0
    assert fip_new["effect"] > 0, "a fair topology-only baseline must be able to go up"
    # the old rule is kept, reported alongside, and labelled a floor
    floor = fip_new["detail"]["floor"]
    assert floor["level"] == C_FLOOR_LEVEL
    assert floor["rule"] == GENERIC_MULTIPLIER_RULE
    assert floor["effect"] <= 0, "the historical rule cannot produce a positive effect"
    assert "FLOOR" in floor["rule_note"] or "floor" in floor["rule_note"].lower()
    old = baseline_topology_only(
        "fipronil", PAPER_CONC, graph="named", rule="depressant_floor"
    )
    assert old["effect"] == pytest.approx(floor["effect"])


def test_every_floor_rule_entry_is_non_positive_and_the_fix_moves_rho_up():
    tab = ablation_table(compounds=SET, concs_M=(PAPER_CONC,))
    for row in tab["rows"]:
        assert row[C_FLOOR_LEVEL] <= 1e-9, (
            "the historical generic rule is structurally incapable of a "
            "positive effect: that is why it is a floor, not a competitor"
        )
    cmp_block = tab["information_gain"][f"{PAPER_CONC:.3e}"]["generic_rule_comparison"]
    assert cmp_block["delta_rho"] > 0
    assert "FLOOR" in cmp_block["depressant_floor"]["role"]
    assert any("direction-aware generic rule" in s for s in tab["statements"])


def test_concentration_dependence_is_reported_per_concentration():
    tab = ablation_table(compounds=SET, concs_M=(1e-8, PAPER_CONC))
    dep = tab["concentration_dependence"]
    assert [r["conc_M"] for r in dep["rows"]] == [1e-8, PAPER_CONC]
    for row in dep["rows"]:
        assert row["worst_level"] in LEVELS
        assert row["worst_level_with_floor_rule"] in LEVELS
        assert row["rho_C_floor_rule"] is not None
    assert any("worst of the four only at the lowest" in s for s in dep["statements"])


# --------------------------------------------------------------------------
# signed rho (4a): an inverted ordering is not a reproduction
# --------------------------------------------------------------------------
def test_reproduction_flag_uses_signed_rho():
    from flylab.analysis.baselines import _corr_row

    forward = list(range(10))
    inverted = _corr_row("X", [-v for v in forward], forward, "n", "u")
    assert inverted["spearman_rho_vs_full"] == pytest.approx(-1.0)
    assert inverted["reproduces_full_ordering"] is False
    assert inverted["reproduces_full"] is False
    assert inverted["ordering_inverted"] is True
    same = _corr_row("X", forward, forward, "n", "u")
    assert same["reproduces_full_ordering"] is True
    assert same["ordering_inverted"] is False


def test_information_added_column_is_signed_and_says_so():
    out = ablation("imidacloprid", PAPER_CONC, compounds=SET)
    rows = {r["level"]: r for r in out["information_gain"]["levels"]}
    assert rows["A_receptor_only"]["information_added_vs_previous"] is None
    for lvl in LEVELS[1:]:
        assert rows[lvl]["information_added_basis"] == "signed_spearman_rho"
    # signed: the increment is exactly the difference of the signed rhos
    a = rows["A_receptor_only"]["spearman_rho_vs_full"]
    b = rows["B_composition_only"]["spearman_rho_vs_full"]
    assert rows["B_composition_only"]["information_added_vs_previous"] == pytest.approx(b - a)


# --------------------------------------------------------------------------
# B versus D: algebraic identity or finding? (3)
# --------------------------------------------------------------------------
def test_composition_reference_distribution_gives_an_unsurprising_range():
    out = composition_reference_distribution(PAPER_CONC, compounds=SET, n_draws=5, seed=1)
    obs = out["observed"]["spearman_rho"]
    assert obs is not None
    cond = out["conditions"]
    ref = cond["random_single_gain_vectors"]["spearman_rho"]
    assert ref["n"] == 5 and -1.0 <= ref["median"] <= 1.0
    assert out["unsurprising_rho_range"] == [ref["p05"], ref["p95"]]
    assert cond["random_gain_vectors"]["spearman_rho"]["n"] == 5
    # shuffling compound labels cannot move the correlation: both levels are
    # functions of the same gain vector. That is the point of the condition.
    assert cond["shuffled_compound_assignment"]["identical_to_observed"] is True
    assert any("EVIDENCE AGAINST" in s for s in out["statements"])
    assert any("floor" in s for s in out["statements"])
    json.dumps(out)


def test_composition_dominance_under_normalisations():
    out = composition_dominance_under_normalisations(PAPER_CONC, compounds=SET)
    assert set(out["by_normalisation"]) == {"row_abs", "none", "degree"}
    for mode, block in out["by_normalisation"].items():
        assert block["spearman_rho_b_vs_d"] is None or -1 <= block["spearman_rho_b_vs_d"] <= 1
    assert any("normalisation" in s for s in out["statements"])
    json.dumps(out)


def test_full_effect_direct_matches_the_full_model():
    from flylab.analysis.baselines import _full_effect_direct
    from flylab.circuit.rate import compute_gains

    for compound in ("imidacloprid", "fipronil"):
        gains, _ = compute_gains(compound, PAPER_CONC)
        assert _full_effect_direct(gains) == pytest.approx(
            full_model(compound, PAPER_CONC)["effect"], rel=1e-9
        )


# --------------------------------------------------------------------------
# the glutamate sign disagreement (4b)
# --------------------------------------------------------------------------
def test_levels_b_and_d_disagree_on_glutamate_and_the_gap_is_quantified():
    from flylab.circuit.rate import SIGN

    assert SIGN["glutamate"] < 0, "the rate engine signs glutamate inhibitory"
    assert set(NT_SIGN_CONVENTIONS) == {"wholens", "rate_engine"}
    # a gain on glutamate alone moves the index in opposite directions
    gains = {"g_ach": 1.0, "g_gaba": 1.0, "g_glu": 1.5, "g_oct": 1.0, "g_nav": 1.0, "ach_tone": 1.0}
    wholens = composition_effect_from_gains(gains, None, "named", "wholens")
    reconciled = composition_effect_from_gains(gains, None, "named", "rate_engine")
    assert wholens > 0 > reconciled

    b = baseline_composition_only("ivermectin", PAPER_CONC, graph="named")
    assert b["detail"]["nt_sign_convention"] == "wholens"
    assert b["detail"]["effect_other_convention"] is not None
    assert any("sign of glutamate" in w for w in b["warnings"])

    out = glutamate_sign_reconciliation(concs_M=(PAPER_CONC,), compounds=SET)
    row = out["rows"][0]
    assert row["rho_wholens"] is not None and row["rho_rate_engine"] is not None
    assert out["max_abs_delta_rho"] is not None
    assert any("glutamate" in s for s in out["statements"])
    json.dumps(out)


def test_information_gain_carries_the_sign_convention_check():
    out = ablation("imidacloprid", PAPER_CONC, compounds=SET)
    check = out["information_gain"]["composition_sign_convention_check"]
    assert check["shipped_convention"] == "wholens"
    assert check["alternative_convention"] == "rate_engine"
    assert check["delta_rho"] is not None
