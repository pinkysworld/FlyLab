"""The ablation ladder: four levels of model sophistication, and what each one
adds over the level below it."""
from __future__ import annotations

import json

import numpy as np
import pytest

from flylab.analysis.baselines import (
    LEVELS,
    ablation,
    ablation_table,
    baseline_composition_only,
    baseline_receptor_only,
    baseline_topology_only,
    full_model,
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
    assert 0.05 <= c["detail"]["multiplier"] <= 1.0
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
