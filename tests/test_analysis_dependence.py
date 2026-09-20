"""Connectome-Dependence Analysis: the profile D, the classification rule,
the information ladder and the landscape.

Every test here uses a small ``n`` so the default suite stays fast; the
high-permutation versions (n >= 500) are marked ``slow``.
"""
from __future__ import annotations

import json
import time

import pytest

from flylab.analysis.dependence import (
    CLASSES,
    CONFIRMATORY_COMPOUNDS,
    DEFAULT_DELTA_FRAC,
    BALANCED_TRANSMITTERS,
    DEFAULT_MODES,
    INFORMATION_LADDER,
    LADDER_FALLBACKS,
    MODE_INFORMATION,
    PAPER_N,
    STRUCTURAL_MODES,
    VERDICTS,
    _shuffle_state_any,
    balance_report,
    benjamini_hochberg,
    classify,
    cut_census,
    dependence_landscape,
    dependence_profile,
    effective_effect_floor,
    equivalence_margin,
    estimate_landscape_runtime,
    ladder_power,
    ladder_recovery,
    mode_verdict,
    necessary_information_level,
    synthetic_cut,
)
from flylab.analysis.nullmodels import MODES

PAPER_CONC = 1e-6


@pytest.fixture(scope="module")
def imidacloprid_profile():
    return dependence_profile("imidacloprid", PAPER_CONC, n=60, seed=0)


@pytest.fixture(scope="module")
def fipronil_profile():
    return dependence_profile("fipronil", PAPER_CONC, n=60, seed=0)


# --------------------------------------------------------------------------
# shape and serialisability
# --------------------------------------------------------------------------
def test_profile_reports_D_p_and_resolution(imidacloprid_profile):
    p = imidacloprid_profile
    # the four nullmodels degradations are always keyed, and the weight-matched
    # transmitter null this module adds is keyed beside them
    assert {"z_sign", "z_weight", "z_degree", "z_ER"} <= set(p["D"])
    assert {"p_sign", "p_weight", "p_degree", "p_ER"} <= set(p["p"])
    assert "z_sign_wm" in p["D"] and "p_sign_wm" in p["p"]
    assert p["D_keys"] == ["z_sign", "z_weight", "z_degree", "z_ER"]
    assert len(p["modes"]) == len(DEFAULT_MODES) == len(MODES) + 1
    for row in p["modes"]:
        # the empirical p is the headline statistic and knows its own resolution
        assert 0.0 < row["p_two_sided"] <= 1.0
        assert row["p_resolution"] == pytest.approx(1.0 / (row["n_ok"] + 1))
        assert row["p_two_sided"] >= row["p_resolution"] - 1e-12
        assert row["information_kept"] == MODE_INFORMATION[row["mode"]]["keeps"]
    assert p["label"] == "model_derived" and p["warnings"]
    json.dumps(p)


def test_profile_p_ordering_matches_the_null_distribution_arm():
    """A profile row must be the same number as the null_distribution it wraps."""
    from flylab.analysis.nullmodels import null_distribution

    prof = dependence_profile(
        "fipronil", PAPER_CONC, modes=("weight_permute",), n=25, seed=4
    )
    direct = null_distribution(
        "subgraph", "fipronil", PAPER_CONC, "weight_permute", n=25, seed=4, readout="mean_hz"
    )
    row = prof["modes"][0]
    assert row["p_two_sided"] == direct["p_two_sided"]
    assert row["z"] == direct["z"]
    assert prof["real_effect"] == direct["real_effect"]


# --------------------------------------------------------------------------
# the two headline findings
# --------------------------------------------------------------------------
def test_imidacloprid_mean_rate_is_composition_dominated(imidacloprid_profile):
    """v0.5 finding, now with a permutation p: the nicotinic agonist effect
    beats the Erdos-Renyi null and no structure-preserving null."""
    p = imidacloprid_profile
    assert p["real_effect"] < 0  # the network mean rate falls
    assert p["class"] == "composition-dominated", p["classification"]["reason"]
    assert p["p"]["p_ER"] <= 0.05
    assert p["p"]["p_weight"] > 0.05 and p["p"]["p_degree"] > 0.05
    level = p["necessary_information_level"]
    assert level["mode"] in ("rewire_degree_preserving", "weight_permute")
    assert level["rank"] < 4


def test_fipronil_is_topology_sensitive_and_classified_differently(
    fipronil_profile, imidacloprid_profile
):
    """The contrast the paper rests on: same cut, same readout, other verdict."""
    f = fipronil_profile
    assert f["real_effect"] > 0  # disinhibition raises the mean rate
    assert f["class"] == "topology-dependent", f["classification"]["reason"]
    assert f["classification"]["structural_distinguishable"]
    assert f["class"] != imidacloprid_profile["class"]
    # and the information level is strictly higher up the ladder
    assert (
        f["necessary_information_level"]["rank"]
        > imidacloprid_profile["necessary_information_level"]["rank"]
    )


# --------------------------------------------------------------------------
# the classification rule itself (pure function, no circuit runs)
# --------------------------------------------------------------------------
def _rows(**p):
    return {mode: {"p_two_sided": p.get(mode)} for mode in MODES}


def test_classification_rule_is_the_documented_one():
    beats_er_only = _rows(
        sign_permute=0.4, weight_permute=0.6, rewire_degree_preserving=0.3, erdos_renyi=0.001
    )
    assert classify(-6.0, beats_er_only)["class"] == "composition-dominated"
    assert classify(-6.0, beats_er_only)["er_distinguishable"] is True

    beats_structural = _rows(
        sign_permute=0.2, weight_permute=0.01, rewire_degree_preserving=0.3, erdos_renyi=0.001
    )
    assert classify(0.9, beats_structural)["class"] == "topology-dependent"

    sign_only = _rows(
        sign_permute=0.01, weight_permute=0.4, rewire_degree_preserving=0.5, erdos_renyi=0.2
    )
    assert classify(0.9, sign_only)["class"] == "mixed"

    nothing = _rows(
        sign_permute=0.9, weight_permute=0.9, rewire_degree_preserving=0.9, erdos_renyi=0.9
    )
    out = classify(0.9, nothing)
    assert out["class"] == "composition-dominated" and out["network_insensitive"] is True

    assert classify(0.0, beats_structural)["class"] == "no-effect"
    assert classify(None, beats_structural)["class"] == "undefined"
    assert set(CLASSES) >= {c["class"] for c in [classify(None, nothing), classify(0.0, nothing)]}


def test_necessary_information_level_walks_the_ladder():
    rows = [
        {"mode": "erdos_renyi", "p_two_sided": 0.001, "n": 50, "stabilised": True},
        {"mode": "rewire_degree_preserving", "p_two_sided": 0.001, "n": 50, "stabilised": True},
        {"mode": "weight_permute", "p_two_sided": 0.001, "n": 50, "stabilised": True},
        {"mode": "sign_permute", "p_two_sided": 0.4, "n": 50, "stabilised": True},
    ]
    profile = {"modes": rows, "real_effect": 1.0, "alpha": 0.05, "effect_floor": 1e-6}
    lvl = necessary_information_level(profile)
    assert lvl["mode"] == "sign_permute" and lvl["rank"] == 3

    for row in rows:
        row["p_two_sided"] = 0.001
    assert necessary_information_level(profile)["level"] == "real_connectome"

    for row in rows:
        row["p_two_sided"] = 0.9
    assert necessary_information_level(profile)["mode"] == "erdos_renyi"

    profile["real_effect"] = 0.0
    assert necessary_information_level(profile)["level"] == "none_required"
    assert INFORMATION_LADDER[-1] == "real_connectome"


# --------------------------------------------------------------------------
# convergence / stabilisation
# --------------------------------------------------------------------------
def test_convergence_flag_and_checkpoints_behave(fipronil_profile):
    for row in fipronil_profile["modes"]:
        assert isinstance(row["stabilised"], bool)
        assert row["n_stabilised"] is None or row["n_stabilised"] < row["n"]
    # The erdos-renyi arm sits right outside its null, so its p is pinned at
    # the resolution floor: at n = 60 it is still falling with n, and the
    # convergence check must say so rather than call it stable.
    er = [r for r in fipronil_profile["modes"] if r["mode"] == "erdos_renyi"][0]
    assert er["p_two_sided"] <= 0.05
    assert er["p_two_sided"] == pytest.approx(er["p_resolution"])
    assert er["resolution_limited"] is True
    assert er["stabilised"] is False


def test_permutations_needed_before_the_conclusion_stops_moving():
    """The advisor's question, answered by the flag itself: a sign-permute arm
    that has not settled at n = 40 has settled by n = 600, and the result says
    from which n it was stable."""
    small = dependence_profile(
        "fipronil", PAPER_CONC, modes=("sign_permute",), n=40, seed=0
    )["modes"][0]
    big = dependence_profile(
        "fipronil", PAPER_CONC, modes=("sign_permute",), n=600, seed=0
    )["modes"][0]
    assert small["stabilised"] is False and small["n_stabilised"] is None
    assert big["stabilised"] is True
    assert big["n_stabilised"] is not None and big["n_stabilised"] < big["n"]
    # both agree on the verdict; it is the precision that changed
    assert (small["p_two_sided"] > 0.05) == (big["p_two_sided"] > 0.05)
    assert big["p_resolution"] < small["p_resolution"]


# --------------------------------------------------------------------------
# the landscape
# --------------------------------------------------------------------------
def test_landscape_shape_matches_its_inputs():
    compounds = ["imidacloprid", "fipronil", "diazepam"]
    concs = (1e-7, 1e-6)
    t0 = time.perf_counter()
    land = dependence_landscape(compounds=compounds, concs_M=concs, n=25, seed=0)
    # Wall-clock budget: this guards against an accidentally quadratic
    # implementation, it is not a benchmark. CI runners and parallel test
    # runs are slow and highly variable, so the bound is deliberately loose.
    assert time.perf_counter() - t0 < 300.0
    assert land["shape"] == [3, 2]
    assert land["n_cells"] == 6 == len(land["cells"]) == len(land["table"])
    assert [c["compound"] for c in land["cells"]][:2] == ["imidacloprid", "imidacloprid"]
    assert land["p_resolution"] == pytest.approx(1.0 / 26)
    counts = land["summary"]["class_counts"]
    assert sum(counts.values()) == 6
    assert set(counts) <= set(CLASSES)
    assert all(row["class"] in CLASSES for row in land["table"])
    assert land["label"] == "model_derived" and land["warnings"]
    json.dumps(land)


def test_landscape_cells_equal_standalone_profiles():
    """The shared shuffle stream must not change a cell's numbers."""
    land = dependence_landscape(
        compounds=["fipronil"], concs_M=(1e-6,), modes=("sign_permute", "erdos_renyi"),
        n=15, seed=2,
    )
    solo = dependence_profile(
        "fipronil", 1e-6, modes=("sign_permute", "erdos_renyi"), n=15, seed=2
    )
    cell = land["cells"][0]
    assert cell["real_effect"] == solo["real_effect"]
    for a, b in zip(cell["modes"], solo["modes"]):
        assert a["mode"] == b["mode"]
        assert a["p_two_sided"] == b["p_two_sided"]
        assert a["z"] == pytest.approx(b["z"])
    assert cell["class"] == solo["class"]


def test_landscape_separates_the_two_known_compounds():
    """The raw screen separates them; the correction is asserted separately."""
    land = dependence_landscape(
        compounds=["imidacloprid", "fipronil"], concs_M=(1e-6,), n=60, seed=0
    )
    by = {row["compound"]: row for row in land["table"]}
    assert by["imidacloprid"]["class_raw"] == "composition-dominated"
    assert by["fipronil"]["class_raw"] == "topology-dependent"
    assert "fipronil" in land["summary"]["topology_dependent_compounds_raw"]
    assert "imidacloprid" in land["summary"]["composition_dominated_compounds"]


def test_too_few_shuffles_cannot_classify_and_says_so():
    """The reviewer's objection in one test: at n = 10 the resolution floor is
    1/11 > alpha, so no cell can be topology-dependent whatever the wiring does,
    and the result has to say that rather than quietly mislabel it."""
    land = dependence_landscape(compounds=["fipronil"], concs_M=(1e-6,), n=10)
    assert land["p_resolution"] > land["alpha"]
    assert land["resolution_coarser_than_alpha"] is True
    assert land["table"][0]["class"] != "topology-dependent"
    assert any("no cell in this landscape" in w for w in land["warnings"])
    prof = dependence_profile("fipronil", 1e-6, n=10)
    assert prof["resolution_coarser_than_alpha"] is True
    assert any("NO mode can be distinguished" in w for w in prof["warnings"])
    # ... and with enough shuffles the same cell is classified topology-dependent
    # on the raw screen (the FDR-corrected verdict needs more shuffles still,
    # which is the point of test_landscape_fdr_guards_an_unresolvable_grid)
    ok = dependence_landscape(compounds=["fipronil"], concs_M=(1e-6,), n=60)
    assert ok["resolution_coarser_than_alpha"] is False
    assert ok["table"][0]["class_raw"] == "topology-dependent"


def test_runtime_estimate_scales_with_n_and_cells():
    a = estimate_landscape_runtime(n_cells=10, n=100)
    b = estimate_landscape_runtime(n_cells=10, n=200)
    c = estimate_landscape_runtime(n_cells=20, n=100)
    assert b["estimate_s"] == pytest.approx(2 * a["estimate_s"])
    assert c["estimate_s"] > a["estimate_s"]
    json.dumps(a)


@pytest.mark.slow
def test_high_permutation_profile_is_stable():
    """The reviewer's ask: 500+ permutations, led by the empirical p.

    At n = 500 the verdict at alpha is settled for every mode from an earlier
    checkpoint; the *value* of a mid-range p (the sign-permute arm sits near
    0.27) is still moving, and the two flags say which is which.
    """
    p = dependence_profile("imidacloprid", PAPER_CONC, n=500, seed=0)
    assert p["p_resolution"] <= 1 / 500
    assert p["resolution_coarser_than_alpha"] is False
    assert p["class"] == "composition-dominated"
    assert p["verdict_stable"] is True
    settled = [r for r in p["modes"] if r["stabilised"]]
    assert len(settled) >= 3
    assert all(r["n_stabilised"] <= 250 for r in settled)


def test_the_n_ladder_is_usable_at_every_rung():
    """FAST_N must still be able to resolve alpha = 0.05 (n = 19 is the bare
    minimum, 1/(19+1) = 0.05; FAST_N sits just above it)."""
    from flylab.analysis.dependence import DEFAULT_N, FAST_N, PAPER_N
    from flylab.analysis.nullmodels import p_resolution

    assert p_resolution(FAST_N) <= 0.05
    assert p_resolution(15) > 0.05  # anything much smaller cannot resolve alpha
    assert FAST_N <= 25
    assert FAST_N < DEFAULT_N < PAPER_N
    land = dependence_landscape(
        compounds=["fipronil"], concs_M=(1e-6,), modes=("erdos_renyi",), n=FAST_N
    )
    assert land["resolution_coarser_than_alpha"] is False


def test_non_monotone_ladder_is_flagged_not_hidden():
    """weight_permute and sign_permute are not nested, so a degree-preserving
    rewire can reproduce an effect that a weight permutation does not. The
    level is still the cheapest model that works, and it says it is odd."""
    rows = [
        {"mode": "erdos_renyi", "p_two_sided": 0.001, "n": 500, "stabilised": True},
        {"mode": "rewire_degree_preserving", "p_two_sided": 0.40, "n": 500, "stabilised": True},
        {"mode": "weight_permute", "p_two_sided": 0.01, "n": 500, "stabilised": True},
        {"mode": "sign_permute", "p_two_sided": 0.30, "n": 500, "stabilised": True},
    ]
    profile = {"modes": rows, "real_effect": 1.0, "alpha": 0.05, "effect_floor": 1e-6}
    lvl = necessary_information_level(profile)
    assert lvl["mode"] == "rewire_degree_preserving"
    assert lvl["non_monotone"] is True
    assert lvl["richer_models_distinguishable"] == ["weight_permute"]
    assert any("not monotone" in w for w in lvl["warnings"])
    # and the cell is still classified topology-dependent, because a
    # structure-preserving null was beaten
    from flylab.analysis.dependence import classify

    assert classify(1.0, {r["mode"]: r for r in rows})["class"] == "topology-dependent"


# --------------------------------------------------------------------------
# Defect 1: the three-way verdict, and the margin it rests on
# --------------------------------------------------------------------------
def test_every_mode_returns_a_three_way_verdict_with_its_delta(imidacloprid_profile):
    """A mode is distinguishable, equivalent within tolerance, or
    indeterminate. Nothing else, and never "reproduces the effect"."""
    p = imidacloprid_profile
    assert p["delta"] is not None and p["delta"] > 0
    assert p["delta_frac"] == pytest.approx(DEFAULT_DELTA_FRAC)
    # the margin is an absolute value in readout units, and it is recorded
    assert p["delta"] == pytest.approx(DEFAULT_DELTA_FRAC * abs(p["real_vehicle"]))
    assert str(round(p["delta"], 6)) or p["delta_scale"]
    assert "vehicle readout" in p["delta_scale"]

    for row in p["modes"]:
        assert row["verdict"] in VERDICTS
        assert row["delta"] == pytest.approx(p["delta"])
        assert row["abs_gap_from_null_median"] == pytest.approx(
            abs(row["real_effect"] - row["null_median"])
        )
        # the verdict is exactly the documented rule, recomputed here
        if row["p_two_sided"] <= p["alpha"]:
            assert row["verdict"] == "distinguishable"
        elif row["abs_gap_from_null_median"] < p["delta"]:
            assert row["verdict"] == "equivalent_within_tolerance"
        else:
            assert row["verdict"] == "indeterminate"
        assert row["verdict_reason"]
    assert set(p["verdicts"]) == {r["mode"] for r in p["modes"]}


def test_non_significance_is_never_reported_as_equivalence():
    """The blocking defect: a mode the real effect does not beat must not be
    called equivalent unless the gap is actually inside the margin."""
    # far from the null median, but nowhere near significant
    v = mode_verdict(real_effect=-6.0, null_median=-3.0, p=0.275, alpha=0.05, delta=0.33)
    assert v["verdict"] == "indeterminate"
    assert v["within_tolerance"] is False
    assert "not evidence of equivalence" in v["verdict_reason"]

    # inside the margin: only now may equivalence be claimed
    v = mode_verdict(real_effect=-6.0, null_median=-6.1, p=0.275, alpha=0.05, delta=0.33)
    assert v["verdict"] == "equivalent_within_tolerance"

    # rejected, and small: the test decides, but the tension is visible
    v = mode_verdict(real_effect=0.9, null_median=0.65, p=0.006, alpha=0.05, delta=0.33)
    assert v["verdict"] == "distinguishable" and v["within_tolerance"] is True

    # no margin at all: equivalence is unavailable, not assumed
    v = mode_verdict(real_effect=-6.0, null_median=-6.1, p=0.275, alpha=0.05, delta=None)
    assert v["verdict"] == "indeterminate"


def test_the_equivalence_margin_is_a_fraction_of_the_vehicle_readout():
    m = equivalence_margin(6.631244285926712)
    assert m["delta"] == pytest.approx(0.05 * 6.631244285926712)
    assert m["prespecified"] is True
    # an absolute override is honoured and recorded as such
    m2 = equivalence_margin(6.63, delta=0.5)
    assert m2["delta"] == 0.5 and "absolute" in m2["delta_scale"]
    # a silenced or missing baseline yields no margin rather than a fake one
    for bad in (None, 0.0):
        assert equivalence_margin(bad)["delta"] is None


def test_no_returned_string_says_the_shuffle_reproduces_the_effect(
    imidacloprid_profile, fipronil_profile
):
    """The exact wording the reviewer objected to must be gone from output."""
    for profile in (imidacloprid_profile, fipronil_profile):
        blob = json.dumps(profile).lower()
        assert "reproduces" not in blob
        assert "already reproduce" not in blob
    lvl = imidacloprid_profile["necessary_information_level"]
    assert lvl["verdict"] in VERDICTS
    assert "equivalent_within_tolerance" in VERDICTS
    assert isinstance(lvl["equivalent_within_tolerance"], bool)


def test_necessary_level_separates_equivalence_from_indeterminacy():
    """The ladder answer must say which of the two it is."""
    rows = [
        {"mode": "erdos_renyi", "p_two_sided": 0.001, "n": 500, "stabilised": True,
         "real_effect": -6.0, "null_median": 0.0},
        {"mode": "rewire_degree_preserving", "p_two_sided": 0.47, "n": 500,
         "stabilised": True, "real_effect": -6.0, "null_median": -5.95},
        {"mode": "weight_permute", "p_two_sided": 0.59, "n": 500, "stabilised": True,
         "real_effect": -6.0, "null_median": -5.9},
        {"mode": "sign_permute", "p_two_sided": 0.27, "n": 500, "stabilised": True,
         "real_effect": -6.0, "null_median": -3.5},
    ]
    profile = {
        "modes": rows, "real_effect": -6.0, "alpha": 0.05,
        "effect_floor": 1e-6, "delta": 0.33,
    }
    lvl = necessary_information_level(profile)
    assert lvl["mode"] == "rewire_degree_preserving"
    assert lvl["verdict"] == "equivalent_within_tolerance"
    assert lvl["equivalent_within_tolerance"] is True
    assert "equivalent to the real cut" in lvl["description"]

    # widen the gap past the margin: the same non-rejection is now indeterminate
    rows[1]["null_median"] = -5.0
    lvl = necessary_information_level(profile)
    assert lvl["verdict"] == "indeterminate"
    assert lvl["equivalent_within_tolerance"] is False
    assert "not a demonstration that the level suffices" in lvl["description"]
    assert any("INDETERMINATE, not equivalent" in w for w in lvl["warnings"])


# --------------------------------------------------------------------------
# Defect 2: multiplicity across the landscape
# --------------------------------------------------------------------------
def test_benjamini_hochberg_matches_the_textbook_procedure():
    p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
    out = benjamini_hochberg(p, alpha=0.05)
    # step-up by hand: the largest i with p_(i) <= i*alpha/m is i = 2
    # (0.008 <= 0.010, while 0.039 > 0.015), so exactly two are rejected
    assert out["m"] == 10
    assert out["n_rejected"] == 2
    assert out["rejected"][:2] == [True, True]
    assert out["rejected"][2] is False
    assert out["adjusted"][0] == pytest.approx(0.01)
    assert out["adjusted"][1] == pytest.approx(0.04)
    # adjusted values are monotone and never below the raw p
    adj = out["adjusted"]
    assert all(a >= b - 1e-12 for a, b in zip(adj, p))
    assert adj == sorted(adj)
    # None passes through and does not count toward m
    out2 = benjamini_hochberg([0.01, None, 0.9], alpha=0.05)
    assert out2["m"] == 2 and out2["adjusted"][1] is None


def test_benjamini_hochberg_knows_when_resolution_forbids_rejection():
    """With 168 structural tests at 100 shuffles nothing can be rejected."""
    out = benjamini_hochberg([1 / 101] * 168, alpha=0.05, resolution=1 / 101)
    assert out["min_rejections"] == 34
    assert out["n_rejected"] == 168  # all at the floor together, so all clear it
    assert out["resolution_supports_single_rejection"] is False
    # one isolated strong cell cannot survive the same correction
    lonely = [1 / 101] + [0.9] * 167
    out2 = benjamini_hochberg(lonely, alpha=0.05, resolution=1 / 101)
    assert out2["n_rejected"] == 0
    # and a resolution too coarse for the family says so outright
    out3 = benjamini_hochberg([0.5] * 10, alpha=0.05, resolution=0.5)
    assert out3["can_reject"] is False


def test_landscape_reports_both_raw_and_fdr_controlled_counts():
    land = dependence_landscape(
        compounds=["imidacloprid", "fipronil", "diazepam"],
        concs_M=(1e-7, 1e-6),
        n=60,
        seed=0,
    )
    summary = land["summary"]
    # both counts are present, and the corrected one can only be smaller
    assert summary["n_topology_dependent_raw"] >= summary["n_topology_dependent"]
    assert summary["class_counts_raw"]["topology-dependent"] == summary[
        "n_topology_dependent_raw"
    ]
    assert summary["class_counts"]["topology-dependent"] == summary["n_topology_dependent"]
    assert sum(summary["class_counts"].values()) == land["n_cells"]
    assert sum(summary["class_counts_raw"].values()) == land["n_cells"]
    # the family, its size and the threshold are all on the record
    assert land["fdr"]["m"] == land["n_cells"] * len(STRUCTURAL_MODES)
    assert land["fdr_alpha"] == land["alpha"]
    assert set(land["fdr_modes"]) == set(STRUCTURAL_MODES)
    # per cell: the raw p, the adjusted p and the q-value, all three
    for cell in land["cells"]:
        assert cell["class_raw"] in CLASSES and cell["class"] in CLASSES
        for row in cell["modes"]:
            if row["mode"] in STRUCTURAL_MODES:
                assert row["p_adjusted"] is not None
                assert row["q_value"] == row["p_adjusted"]
                assert row["q_value"] >= row["p_two_sided"] - 1e-12
                assert row["fdr_alpha"] == land["fdr_alpha"]
                assert row["p_basis"] == "fdr_adjusted"
            else:
                assert row.get("p_adjusted") is None
        assert cell["q_value"] is not None
    # the summary states what kind of run this was
    assert land["summary"]["design"] in ("exploratory", "confirmatory")
    assert "Benjamini-Hochberg" in land["summary"]["statement"]
    json.dumps(land)


def test_landscape_says_when_its_count_change_is_the_correction():
    land = dependence_landscape(
        compounds=["imidacloprid", "fipronil"], concs_M=(1e-6,), n=60, seed=0
    )
    assert land["summary"]["n_topology_dependent_raw"] == 1
    assert land["summary"]["n_topology_dependent"] == 0
    assert any("multiplicity changes the headline count" in w for w in land["warnings"])


def test_landscape_fdr_guards_an_unresolvable_grid():
    """At 10 shuffles nothing can be rejected; that must be said, not returned
    as a count of zero topology-dependent cells."""
    land = dependence_landscape(
        compounds=["fipronil", "imidacloprid"], concs_M=(1e-6,), n=10, seed=0
    )
    assert land["fdr"]["can_reject"] is False
    assert land["summary"]["fdr_can_reject"] is False
    assert land["summary"]["n_topology_dependent"] == 0
    assert any("FDR CANNOT REJECT ANYTHING" in w for w in land["warnings"])
    assert "uninformative rather than zero" in land["summary"]["statement"]


def test_a_landscape_is_exploratory_and_a_prespecified_profile_is_not():
    # a real screen: 6 cells x 2 structural modes, so the adjusted threshold
    # 0.05/12 sits below what 60 shuffles can resolve
    land = dependence_landscape(
        compounds=["fipronil", "imidacloprid", "diazepam"],
        concs_M=(1e-7, 1e-6),
        n=60,
    )
    assert land["confirmatory"] is False
    assert land["design"] == "exploratory"
    assert any("EXPLORATORY" in w for w in land["warnings"])

    small = dependence_profile("fipronil", 1e-6, modes=("erdos_renyi",), n=25)
    assert small["confirmatory"] is False
    assert any("exploratory" in w for w in small["warnings"])
    # the prespecified confirmatory tests are named, and only they qualify
    assert set(CONFIRMATORY_COMPOUNDS) == {"imidacloprid", "fipronil"}
    forced = dependence_profile(
        "fipronil", 1e-6, modes=("erdos_renyi",), n=25, confirmatory=True
    )
    assert forced["confirmatory"] is True and forced["design"] == "confirmatory (prespecified)"
    assert PAPER_N == 1000


@pytest.mark.slow
def test_fdr_keeps_fipronil_at_a_confirmatory_permutation_count():
    """The correction is not a blunt instrument: at 1000 shuffles the real
    contrast survives it."""
    land = dependence_landscape(
        compounds=["imidacloprid", "fipronil"], concs_M=(1e-6,), n=1000, seed=0, n_jobs=4
    )
    by = {row["compound"]: row for row in land["table"]}
    assert by["fipronil"]["class"] == "topology-dependent"
    assert by["imidacloprid"]["class"] == "composition-dominated"
    assert land["summary"]["n_topology_dependent"] == 1


# --------------------------------------------------------------------------
# the transmitter null has to preserve the WEIGHTED excitation/inhibition
# balance, not just the label histogram
# --------------------------------------------------------------------------
def test_plain_sign_permute_moves_the_weighted_balance_and_says_so():
    """The joint-null defect, measured. `sign_permute` permutes the label
    array, so the share of outgoing synaptic weight carried by each
    transmitter moves - and a gain patch acts through exactly that share."""
    rep = balance_report(graph="named", n=100, seed=0)
    plain = rep["modes"]["sign_permute"]["acetylcholine"]
    assert plain["real"] == pytest.approx(0.6188, abs=5e-3)
    # the real graph sits outside its own null, so the null is not centred
    assert plain["percentile_of_real"] > 95.0
    assert plain["null_mean"] < plain["real"] - 0.05
    assert plain["within_tol"] is False
    assert rep["preserves_weighted_balance"]["sign_permute"] is False
    # and the mode says so in its own metadata rather than only here
    assert MODE_INFORMATION["sign_permute"]["joint"] is True
    assert "weighted" in MODE_INFORMATION["sign_permute"]["destroys"]
    assert "JOINT null" in MODE_INFORMATION["sign_permute"]["caveat"]


def test_the_weight_matched_null_holds_the_balance_and_still_shuffles():
    rep = balance_report(graph="named", n=100, seed=0)
    assert rep["preserves_weighted_balance"]["sign_permute_weight_matched"] is True
    for nt in BALANCED_TRANSMITTERS:
        row = rep["modes"]["sign_permute_weight_matched"][nt]
        assert row["max_abs_deviation"] <= rep["tol"] + 1e-12
    # it must still be a null: a useful fraction of the labels actually move,
    # and the label histogram is preserved exactly
    from collections import Counter

    from flylab.analysis.nullmodels import _state_for

    base = _state_for("named")
    shuffled = _shuffle_state_any(base, "sign_permute_weight_matched", 0)
    assert Counter(base.nt.tolist()) == Counter(shuffled.nt.tolist())
    moved = sum(1 for a, b in zip(base.nt.tolist(), shuffled.nt.tolist()) if a != b)
    assert moved > 0.40 * base.n, "the null must relabel a substantial share of nodes"
    assert shuffled.meta["mode"] == "sign_permute_weight_matched"
    assert shuffled.meta["swaps_accepted"] > 0


def test_the_ladders_rank_three_rung_is_the_weight_matched_null():
    assert INFORMATION_LADDER[3] == "sign_permute_weight_matched"
    assert "sign_permute" not in INFORMATION_LADDER
    assert "sign_permute" in DEFAULT_MODES  # kept: the contrast is informative
    assert LADDER_FALLBACKS["sign_permute_weight_matched"] == ("sign_permute",)


def test_a_ladder_without_the_matched_rung_falls_back_and_names_the_caveat():
    prof = dependence_profile("fipronil", PAPER_CONC, n=25, modes=MODES, seed=0)
    lvl = prof["necessary_information_level"]
    if lvl["mode"] == "sign_permute":
        assert lvl["ladder_substitutions"] == ["sign_permute_weight_matched -> sign_permute"]
        assert any("rung was not run" in w for w in prof["warnings"])
        assert any("JOINT null" in w for w in prof["warnings"])


# --------------------------------------------------------------------------
# what the cut under test looks like, and whether a negative can mean anything
# --------------------------------------------------------------------------
def test_the_named_cut_is_reported_as_an_in_star():
    c = cut_census("named")
    assert c["n_nodes"] == 1126 and c["n_edges"] == 1360
    assert c["edges_onto_seeds"] == 1154
    assert c["share_onto_seeds"] == pytest.approx(0.8485, abs=1e-3)
    assert c["n_nodes_with_in_degree"] == 184
    assert c["share_edges_from_nodes_with_input"] == pytest.approx(0.2493, abs=1e-3)
    assert c["in_star"] is True
    hubs = [h["in_degree"] for h in c["top_in_degrees"][:4]]
    assert hubs == [556, 477, 94, 27]
    assert "in-star" in c["note"]


def test_the_taste_motor_cut_is_not_an_in_star():
    c = cut_census("taste_motor")
    assert c["n_nodes"] == 1841 and c["n_edges"] == 19066
    assert c["mean_degree"] > 10.0
    assert c["share_edges_from_nodes_with_input"] > 0.9
    assert c["in_star"] is False
    # the contrast is the whole point: one cut can barely be rewired
    assert cut_census("named")["mean_degree"] < 2.0


# --------------------------------------------------------------------------
# does the instrument work at all? ground truth, planted and recovered
# --------------------------------------------------------------------------
def test_the_ladder_recovers_a_planted_topology_dependent_effect():
    """Without this, a negative from the ladder is indistinguishable from an
    under-powered test."""
    out = ladder_recovery(strengths=(0.0, 1.0), n=60, seed=0)
    by = {r["loop_strength"]: r for r in out["rows"]}
    assert by[0.0]["topology_dependent"] is False, by[0.0]["p"]
    assert by[1.0]["topology_dependent"] is True, by[1.0]["p"]
    assert out["recovered"] is True
    assert out["false_positive_on_control"] is False
    assert "recovered 1 of 1" in out["statement"]
    json.dumps(out)


def test_the_synthetic_cut_plants_what_it_claims_to_plant():
    planted = synthetic_cut(loop_strength=1.0, seed=0)
    control = synthetic_cut(loop_strength=0.0, seed=0)
    assert planted["meta"]["loop_nodes"] and not control["meta"]["loop_nodes"]
    assert len(planted["edges"]) == len(control["edges"]) + len(
        planted["meta"]["loop_nodes"]
    )
    # the two graphs are otherwise the same background
    assert planted["nodes"][-1]["bodyId"] == control["nodes"][-1]["bodyId"]
    assert set(planted["seeds"]) == {"MN9", "DNp01"}


def test_ladder_power_rises_with_effect_size():
    out = ladder_power(strengths=(0.0, 1.5), ns=(60,), replicates=3, seed=0)
    by = {(g["loop_strength"], g["n"]): g for g in out["grid"]}
    assert by[(0.0, 60)]["detection_rate"] <= 0.34   # the false-positive arm
    assert by[(1.5, 60)]["detection_rate"] >= 0.66   # the planted arm
    assert out["false_positive_rate"] is not None
    json.dumps(out)


@pytest.mark.slow
def test_ladder_power_curve_over_strength_and_permutations():
    out = ladder_power(
        strengths=(0.0, 0.25, 0.5, 1.0, 2.0), ns=(50, 200), replicates=5, seed=0
    )
    rates = {(g["loop_strength"], g["n"]): g["detection_rate"] for g in out["grid"]}
    # power is monotone in effect size at the larger permutation count
    strong = [rates[(s, 200)] for s in (0.0, 0.25, 0.5, 1.0, 2.0)]
    assert strong[-1] >= strong[0]
    assert out["false_positive_rate"] <= 0.2


# --------------------------------------------------------------------------
# the floors and the ladder's own monotonicity, surfaced for the paper
# --------------------------------------------------------------------------
def test_a_relative_effect_floor_keeps_numerical_noise_out_of_the_counts():
    f = effective_effect_floor(6.631, effect_floor=1e-6, effect_floor_frac=0.01)
    assert f["effect_floor"] == pytest.approx(0.06631)
    assert f["binding"] == "relative"
    # an effect of 2.7e-05 Hz against a 6.6 Hz baseline is not a drug effect
    assert abs(2.7e-05) < f["effect_floor"]
    # with no baseline the absolute floor still applies
    assert effective_effect_floor(None)["binding"] == "absolute"

    land = dependence_landscape(
        compounds=["acetylcholine", "fipronil"], concs_M=(1e-8,), n=25, seed=0
    )
    s = land["summary"]
    assert s["class_counts_absolute_floor"] is not None
    assert sum(s["class_counts_absolute_floor"].values()) == land["n_cells"]
    assert s["n_below_relative_effect_floor"] >= 1
    assert any("relative effect floor" in w for w in land["warnings"])


def test_the_landscape_counts_its_non_monotone_ladders():
    land = dependence_landscape(
        compounds=["fipronil", "imidacloprid"], concs_M=(1e-6,), n=40, seed=0
    )
    s = land["summary"]
    assert s["n_non_monotone_ladders"] == len(s["non_monotone_cells"])
    assert s["n_non_monotone_ladders"] <= land["n_cells"]
    for row in s["non_monotone_cells"]:
        assert row["richer_models_distinguishable"]
