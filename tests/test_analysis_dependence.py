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
    INFORMATION_LADDER,
    MODE_INFORMATION,
    classify,
    dependence_landscape,
    dependence_profile,
    estimate_landscape_runtime,
    necessary_information_level,
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
    assert set(p["D"]) == {"z_sign", "z_weight", "z_degree", "z_ER"}
    assert set(p["p"]) == {"p_sign", "p_weight", "p_degree", "p_ER"}
    assert p["D_keys"] == ["z_sign", "z_weight", "z_degree", "z_ER"]
    assert len(p["modes"]) == len(MODES)
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
    assert f["classification"]["structural_beaten"]
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
    assert classify(-6.0, beats_er_only)["er_beaten"] is True

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
    assert time.perf_counter() - t0 < 60.0
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
    land = dependence_landscape(
        compounds=["imidacloprid", "fipronil"], concs_M=(1e-6,), n=60, seed=0
    )
    by = {row["compound"]: row for row in land["table"]}
    assert by["imidacloprid"]["class"] == "composition-dominated"
    assert by["fipronil"]["class"] == "topology-dependent"
    assert "fipronil" in land["summary"]["topology_dependent_compounds"]
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
    assert any("NO mode can be beaten" in w for w in prof["warnings"])
    # ... and with enough shuffles the same cell is classified topology-dependent
    ok = dependence_landscape(compounds=["fipronil"], concs_M=(1e-6,), n=60)
    assert ok["resolution_coarser_than_alpha"] is False
    assert ok["table"][0]["class"] == "topology-dependent"


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
    assert lvl["richer_models_beaten"] == ["weight_permute"]
    assert any("not monotone" in w for w in lvl["warnings"])
    # and the cell is still classified topology-dependent, because a
    # structure-preserving null was beaten
    from flylab.analysis.dependence import classify

    assert classify(1.0, {r["mode"]: r for r in rows})["class"] == "topology-dependent"
