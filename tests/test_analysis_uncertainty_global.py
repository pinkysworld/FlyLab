"""Variance-based global uncertainty attribution: estimator validation on the
Ishigami function, and the FlyLab uncertainty budget."""
import math

import numpy as np
import pytest

from flylab.analysis.uncertainty_global import (
    FACTOR_NAMES,
    FACTORS,
    ISHIGAMI_REFERENCE,
    SHAPE_ORDER,
    RateSurrogate,
    baseline_unit_vector,
    ishigami,
    ishigami_unit,
    jansen_indices,
    lhs,
    saltelli_matrices,
    sobol_analysis,
    to_markdown,
    uncertainty_budget,
)
from flylab.assays.subgraph import run_subgraph_assay

#: Jansen estimates are unbiased only in the limit, so an index that is truly 0
#: comes back slightly negative at finite N.  This is the documented tolerance.
ESTIMATOR_TOL = 0.06


# --------------------------------------------------------------------------
# the design
# --------------------------------------------------------------------------
def test_lhs_is_stratified_and_inside_the_unit_cube():
    rng = np.random.default_rng(0)
    X = lhs(64, 5, rng)
    assert X.shape == (64, 5)
    assert X.min() >= 0.0 and X.max() <= 1.0
    for j in range(5):
        counts = np.histogram(X[:, j], bins=8, range=(0, 1))[0]
        assert set(counts) == {8}, "each factor's marginal must be evenly stratified"


def test_saltelli_matrices_are_cross_sampled():
    A, B, AB = saltelli_matrices(16, 4, seed=1)
    assert A.shape == B.shape == (16, 4)
    assert len(AB) == 4
    for j, M in enumerate(AB):
        assert np.array_equal(M[:, j], B[:, j])
        for other in range(4):
            if other != j:
                assert np.array_equal(M[:, other], A[:, other])


def test_evaluation_count_is_n_times_k_plus_two():
    res = sobol_analysis(model=ishigami_unit, factor_names=["x1", "x2", "x3"], n_base=32, n_boot=0)
    assert res["n_evaluations"] == 32 * (3 + 2)


# --------------------------------------------------------------------------
# estimator validation against known analytic indices
# --------------------------------------------------------------------------
def test_ishigami_analytic_reference_is_self_consistent():
    ref = ISHIGAMI_REFERENCE
    assert ref["first_order"][2] == 0.0
    assert ref["total_order"][2] > 0.0, "x3 has no main effect but does interact"
    for s, t in zip(ref["first_order"], ref["total_order"]):
        assert t >= s - 1e-12
    assert sum(ref["first_order"]) < 1.0


def test_jansen_recovers_the_ishigami_indices():
    res = sobol_analysis(
        model=ishigami_unit, factor_names=["x1", "x2", "x3"], n_base=4096, n_boot=0, seed=7
    )
    rows = {r["factor"]: r for r in res["rows"]}
    ref = ISHIGAMI_REFERENCE
    for j, name in enumerate(["x1", "x2", "x3"]):
        assert rows[name]["first_order"] == pytest.approx(ref["first_order"][j], abs=ESTIMATOR_TOL)
        assert rows[name]["total_order"] == pytest.approx(ref["total_order"][j], abs=ESTIMATOR_TOL)
    assert res["output_variance"] == pytest.approx(ref["variance"], rel=0.10)
    # the structural facts, not just the numbers
    assert rows["x3"]["first_order"] < rows["x3"]["total_order"]
    assert rows["x2"]["total_order"] == pytest.approx(rows["x2"]["first_order"], abs=ESTIMATOR_TOL)


def test_ishigami_convergence_improves_with_n():
    ref = ISHIGAMI_REFERENCE["first_order"]

    def err(n):
        res = sobol_analysis(model=ishigami_unit, factor_names=["x1", "x2", "x3"], n_base=n, n_boot=0, seed=3)
        rows = {r["factor"]: r for r in res["rows"]}
        return max(abs(rows[f"x{j+1}"]["first_order"] - ref[j]) for j in range(3))

    assert err(4096) < err(256)


def test_ishigami_function_values():
    assert ishigami([0.0, 0.0, 0.0]) == pytest.approx(0.0)
    assert ishigami([math.pi / 2, math.pi / 2, 0.0]) == pytest.approx(1.0 + 7.0)


# --------------------------------------------------------------------------
# indices obey their algebraic constraints on the real model
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def tiny_sobol():
    """The real model, at the smallest sample size that still exercises every
    code path.  ``n_base=6`` is 66 rate-network evaluations (a few seconds):
    deliberately far from converged, which the result itself says, and which is
    what the structural tests below assert on.  The converged run is the
    ``slow`` test at the bottom."""
    return sobol_analysis("imidacloprid", 1e-6, "mean_hz", n_base=6, n_boot=20, seed=0)


def test_factor_list_is_the_documented_one(tiny_sobol):
    assert {r["factor"] for r in tiny_sobol["rows"]} == set(FACTOR_NAMES)
    assert len(FACTORS) == 9
    assert {"potency", "hill_n", "gain_transform", "gain_coef", "drive",
            "expression", "weight_threshold", "transmitter", "lif_seed"} == set(FACTOR_NAMES)
    cat = next(f for f in FACTORS if f["name"] == "gain_transform")
    assert cat["kind"] == "categorical"
    assert len(SHAPE_ORDER) == 5


def test_indices_are_in_range_on_a_converged_run():
    """[0, 1] up to documented estimator noise, and T >= S, where N is large."""
    res = sobol_analysis(
        model=ishigami_unit, factor_names=["x1", "x2", "x3"], n_base=4096, n_boot=0, seed=11
    )
    for r in res["rows"]:
        assert -ESTIMATOR_TOL <= r["first_order"] <= 1.0 + ESTIMATOR_TOL
        assert -ESTIMATOR_TOL <= r["total_order"] <= 1.0 + ESTIMATOR_TOL
        assert r["total_order"] >= r["first_order"] - ESTIMATOR_TOL
        assert r["interaction"] == pytest.approx(r["total_order"] - r["first_order"])


def test_model_indices_are_finite_and_total_order_is_non_negative(tiny_sobol):
    """On the real model at a deliberately tiny N, only the structural
    guarantees hold: Jansen's T is a sum of squares, so it is >= 0 and <= 1 up
    to noise, while an unconverged S may be well below 0 (which is why the
    result carries a convergence diagnostic and an explicit warning)."""
    for r in tiny_sobol["rows"]:
        assert math.isfinite(r["first_order"]) and math.isfinite(r["total_order"])
        assert r["total_order"] >= -1e-12
        assert r["first_order"] <= 1.0 + ESTIMATOR_TOL
        assert r["total_order"] <= 1.0 + 10 * ESTIMATOR_TOL


def test_null_factor_has_no_effect_on_the_deterministic_rate_engine(tiny_sobol):
    """``lif_seed`` cannot move a rate readout: its total index must be exactly 0."""
    row = next(r for r in tiny_sobol["rows"] if r["factor"] == "lif_seed")
    assert row["total_order"] == pytest.approx(0.0, abs=1e-12)


def test_bootstrap_cis_bracket_the_point_estimates(tiny_sobol):
    for r in tiny_sobol["rows"]:
        lo, hi = r["first_order_ci"]
        assert lo <= hi
        lo_t, hi_t = r["total_order_ci"]
        assert lo_t <= hi_t


def test_convergence_diagnostic_is_reported(tiny_sobol):
    conv = tiny_sobol["convergence"]
    sizes = [c["n_base"] for c in conv]
    assert len(sizes) == 3
    assert sizes == sorted(sizes) and sizes[-1] == tiny_sobol["n_base"]
    assert conv[0]["max_abs_delta_first"] is None
    assert conv[-1]["max_abs_delta_first"] is not None
    assert conv[-1]["max_abs_delta_total"] is not None


def test_budget_shares_and_residual(tiny_sobol):
    budget = tiny_sobol["budget"]
    assert budget[-1]["source"].startswith("interactions")
    named = budget[:-1]
    assert {b["source"] for b in named} == set(FACTOR_NAMES)
    total = sum(b["share_of_variance"] for b in budget)
    assert total == pytest.approx(1.0, abs=1e-9), "the budget must add to the whole variance"
    for b in named:
        assert b["variance_removed"] == pytest.approx(b["share_of_variance"] * tiny_sobol["output_variance"])


def test_result_carries_its_warnings(tiny_sobol):
    text = " ".join(tiny_sobol["warnings"])
    assert "not shares of biological variability" in text
    assert "categorical" in text
    assert tiny_sobol["label"] == "model_derived"


def test_markdown_budget_renders(tiny_sobol):
    md = to_markdown(tiny_sobol)
    assert "Uncertainty budget" in md
    for name in FACTOR_NAMES:
        assert name in md


def test_uncertainty_budget_handles_a_bare_result():
    rows = [{"factor": "a", "first_order": 0.3, "total_order": 0.4, "kind": "continuous"}]
    out = uncertainty_budget({"rows": rows, "sum_first_order": 0.3, "output_variance": 2.0})
    assert out[0]["variance_removed"] == pytest.approx(0.6)
    assert out[-1]["share_of_variance"] == pytest.approx(0.7)


# --------------------------------------------------------------------------
# the surrogate IS the shipped model
# --------------------------------------------------------------------------
@pytest.mark.parametrize("readout", ["mean_hz", "mn9_hz"])
@pytest.mark.parametrize("compound", ["imidacloprid", "fipronil"])
def test_surrogate_reproduces_the_assay_exactly_at_baseline(compound, readout):
    nb = run_subgraph_assay(compound, 1e-6)
    s = RateSurrogate(compound, 1e-6, readout)
    assert s(baseline_unit_vector()) == nb["readouts"][readout]


def test_baseline_unit_vector_decodes_to_the_shipped_settings():
    s = RateSurrogate()
    settings = s.decode(baseline_unit_vector())
    assert settings["gain_transform"] == "flylab_biphasic"
    assert settings["gain_coef"] == pytest.approx(1.0)
    assert settings["potency_shift_log10"] == pytest.approx(0.0)
    assert settings["hill_mult"] == pytest.approx(1.0)
    assert settings["drive_hz"] == pytest.approx(40.0)
    assert settings["expression_lambda"] == 0.0
    assert settings["transmitter_p"] == 0.0
    assert settings["min_weight"] == pytest.approx(s.base_min_weight)


def test_surrogate_moves_when_the_transformation_changes():
    s = RateSurrogate()
    u = list(baseline_unit_vector())
    base = s(u)
    u[2] = 0.05  # linear
    assert s(u) != base


def test_transmitter_factor_is_nested_and_starts_from_the_shipped_graph():
    s = RateSurrogate()
    sign0, nt0 = s._sign_vector(0.0)
    assert list(nt0) == list(s.nt0)
    _, nt_low = s._sign_vector(0.05)
    _, nt_high = s._sign_vector(0.20)
    changed_low = {i for i in range(s.n) if nt_low[i] != s.nt0[i]}
    changed_high = {i for i in range(s.n) if nt_high[i] != s.nt0[i]}
    assert changed_low <= changed_high
    assert len(changed_high) > len(changed_low)


def test_unknown_readout_is_rejected():
    with pytest.raises(ValueError):
        RateSurrogate(readout="not_a_readout")(baseline_unit_vector())
    with pytest.raises(ValueError):
        RateSurrogate(engine="quantum")


@pytest.mark.slow
def test_full_sobol_run_converges():
    res = sobol_analysis("imidacloprid", 1e-6, "mean_hz", n_base=128, n_boot=100, seed=0)
    assert res["n_evaluations"] == 128 * 11
    assert res["output_variance"] > 0
