"""Value of information: ranking, experiment mapping and the standing caveat."""
import pytest

from flylab.analysis.uncertainty_global import FACTOR_NAMES, sobol_analysis
from flylab.analysis.voi import (
    EXPERIMENTS,
    NEGATIVE_VOI_NOTE,
    NULL_CONTROL,
    RESOLVED,
    UNRESOLVED,
    VOI_WARNING,
    to_markdown,
    voi,
)


def _surrogate(u):
    """A cheap analytic stand-in with the same nine factors.

    The VOI layer is pure algebra on a sobol_analysis result, so it is tested
    against a function whose factor structure is known and which costs
    microseconds instead of the rate network's tens of milliseconds: factor 3
    (``gain_coef``) dominates, factor 2 (``gain_transform``) is next and
    interacts with it, factor 8 (``lif_seed``) is a null factor.  One real
    end-to-end run is kept as the ``slow`` test at the bottom.
    """
    return (
        4.0 * u[3]                 # gain_coef      - dominant
        + 2.0 * u[2]               # gain_transform - second
        + 3.0 * u[2] * u[3]        # and they interact
        + 1.5 * u[0]               # potency
        + 1.2 * u[4]               # drive
        + 0.9 * u[6]               # weight_threshold
        + 0.6 * u[7]               # transmitter
        + 0.35 * u[5]              # expression
        + 0.2 * u[1]               # hill_n
        + 0.0 * u[8]               # lif_seed       - exact null factor
    )


@pytest.fixture(scope="module")
def tiny_sobol():
    return sobol_analysis(
        model=_surrogate, factor_names=list(FACTOR_NAMES), n_base=2048, n_boot=60, seed=0
    )


@pytest.fixture(scope="module")
def tiny_voi(tiny_sobol):
    return voi(result=tiny_sobol)


def test_every_factor_has_a_concrete_experiment():
    assert set(EXPERIMENTS) == set(FACTOR_NAMES)
    for name, exp in EXPERIMENTS.items():
        assert exp["experiment"] and exp["design"] and exp["resolves"]
        assert exp["cost"] and exp["feasibility"]
        assert len(exp["design"]) > 40, f"{name}: the design must be specific enough to act on"


def test_named_experiments_cover_the_advisor_examples():
    text = " ".join(e["experiment"] + " " + e["design"] for e in EXPERIMENTS.values()).lower()
    assert "immunostain" in text  # transmitter identity
    assert "fish" in text or "scrna" in text  # expression
    assert "dose-response" in text  # potency
    assert "calibration" in text  # gain rule
    assert EXPERIMENTS["expression"]["blocking_gate"]


def test_voi_rows_are_finite_and_ranked(tiny_voi):
    rows = tiny_voi["rows"]
    assert len(rows) == len(FACTOR_NAMES)
    assert [r["rank"] for r in rows] == list(range(1, len(rows) + 1))
    values = [r["voi_var"] for r in rows]
    assert all(v == v and abs(v) < float("inf") for v in values)
    assert values == sorted(values, reverse=True)


def test_voi_is_the_first_order_index_times_the_variance(tiny_sobol, tiny_voi):
    """The raw quantity is S_j * Var(Y); the decision quantity clips it at 0."""
    var = tiny_sobol["output_variance"]
    by_name = {r["factor"]: r for r in tiny_sobol["rows"]}
    for r in tiny_voi["rows"]:
        s1 = by_name[r["factor"]]["first_order"]
        assert r["voi_var_raw"] == pytest.approx(s1 * var)
        assert r["voi_fraction_raw"] == pytest.approx(s1)
        assert r["voi_var"] == pytest.approx(max(0.0, s1) * var)
        assert r["voi_fraction"] == pytest.approx(max(0.0, s1))
        assert r["decision_voi_var"] == pytest.approx(r["voi_var"])


def test_voi_never_exceeds_its_total_order_upper_bound(tiny_voi):
    """Resolving a factor cannot remove more than its total-order share.

    The tolerance is the estimator's, not the identity's: for a factor whose
    true indices are both 0 the Jansen total order is exactly 0 while the
    first order is a small signed noise term, so S can sit a hair above T.
    """
    for r in tiny_voi["rows"]:
        assert r["voi_fraction"] <= r["voi_upper_bound_fraction"] + 0.02


def test_voi_ranking_agrees_with_the_total_order_ranking(tiny_sobol, tiny_voi):
    """The top factor by VOI must also be near the top by total order."""
    by_total = sorted(tiny_sobol["rows"], key=lambda r: -r["total_order"])
    top_total = [r["factor"] for r in by_total[:2]]
    assert tiny_voi["rows"][0]["factor"] in top_total

    # and the two orderings must be positively rank-correlated over the
    # factors that carry a material share (factors at the noise floor have no
    # meaningful rank in either ordering)
    material = [r["factor"] for r in by_total if r["total_order"] >= 0.01]
    voi_rank = {f: i for i, f in enumerate(r["factor"] for r in tiny_voi["rows"]) if f in material}
    tot_rank = {f: i for i, f in enumerate(r["factor"] for r in by_total) if f in material}
    voi_rank = {f: i for i, f in enumerate(sorted(voi_rank, key=voi_rank.get))}
    tot_rank = {f: i for i, f in enumerate(sorted(tot_rank, key=tot_rank.get))}
    n = len(voi_rank)
    d2 = sum((voi_rank[f] - tot_rank[f]) ** 2 for f in voi_rank)
    rho = 1.0 - 6.0 * d2 / (n * (n * n - 1))
    assert n >= 5
    assert rho > 0.8, f"VOI and total-order rankings disagree (rho={rho:.2f})"


def test_recommendation_names_an_experiment(tiny_voi):
    rec = tiny_voi["recommendation"]
    assert rec and rec.startswith("Next experiment:")
    assert any(e["experiment"] in rec for e in EXPERIMENTS.values())


def test_warning_states_this_is_model_variance(tiny_voi):
    joined = " ".join(tiny_voi["warnings"])
    assert VOI_WARNING in tiny_voi["warnings"]
    assert "MODEL variance, not biological variance" in joined
    assert "FIRST-ORDER" in joined
    assert tiny_voi["label"] == "model_derived"


def test_factor_subset_is_honoured(tiny_sobol):
    out = voi(factors=["potency", "gain_coef"], result=tiny_sobol)
    assert {r["factor"] for r in out["rows"]} == {"potency", "gain_coef"}


def test_sd_reduction_is_consistent(tiny_voi):
    sd = tiny_voi["output_sd"]
    for r in tiny_voi["rows"]:
        assert r["residual_sd"] <= sd + 1e-9
        assert r["sd_reduction"] == pytest.approx(sd - r["residual_sd"])


def test_markdown_renders(tiny_voi):
    md = to_markdown(tiny_voi)
    assert "Value of information" in md
    assert "Next experiment" in md
    for name in FACTOR_NAMES:
        assert name in md


@pytest.mark.slow
def test_voi_runs_end_to_end():
    out = voi(n_base=64, n_boot=50)
    assert out["rows"] and out["recommendation"]


# --------------------------------------------------------------------------
# Defect 4: a value of information is never negative
# --------------------------------------------------------------------------
def _sparse_surrogate(u):
    """Only two factors matter; the other seven are exactly zero, so at a small
    sample their Jansen estimates land on both sides of zero."""
    return 4.0 * u[3] + 2.0 * u[2] + 3.0 * u[2] * u[3]


def _noisy_sobol(seed=3):
    """A sample small enough that several true-zero factors come back negative."""
    return sobol_analysis(
        model=_sparse_surrogate,
        factor_names=list(FACTOR_NAMES),
        n_base=128,
        n_boot=60,
        seed=seed,
    )


def test_no_voi_value_is_ever_negative():
    """The blocking arithmetic: resolving an assumption cannot add variance."""
    for res in (_noisy_sobol(3), _noisy_sobol(4)):
        out = voi(result=res)
        raw = [r["voi_fraction_raw"] for r in out["rows"]]
        assert any(v < 0 for v in raw), "this fixture must contain negative estimates"
        for r in out["rows"]:
            assert r["voi_var"] >= 0.0
            assert r["voi_fraction"] >= 0.0
            assert r["decision_voi_var"] >= 0.0
            assert r["voi_upper_bound_var"] >= 0.0
            assert r["sd_reduction"] >= -1e-12
            assert r["residual_sd"] <= out["output_sd"] + 1e-9
        # and the markdown a reader sees carries no negative VOI either
        md = to_markdown(out)
        body = [ln for ln in md.splitlines() if ln.startswith("| ") and "factor" not in ln]
        for line in body:
            cells = [c.strip() for c in line.split("|")]
            assert not cells[3].startswith("-"), line   # decision VOI share
            assert not cells[4].startswith("-"), line   # decision VOI in Hz^2


def test_a_negative_estimate_is_kept_but_labelled_unresolved():
    res = _noisy_sobol(3)
    out = voi(result=res)
    negative = [r for r in out["rows"] if r["voi_fraction_raw"] < 0]
    assert negative
    for r in negative:
        # the raw estimate survives as a convergence diagnostic ...
        assert r["clipped_from_negative"] is True
        assert r["voi_var_raw"] < 0
        # ... and the factor is never presented as a ranked experiment
        assert r["status"] == UNRESOLVED
        assert r["factor"] in out["unresolved_factors"]
        assert r["factor"] not in out["resolved_factors"]
    assert any("estimator noise" in w for w in out["warnings"])
    assert any(UNRESOLVED in w for w in out["warnings"])


def test_a_factor_whose_ci_spans_zero_is_unresolved_not_ranked():
    res = _noisy_sobol(3)
    out = voi(result=res)
    for r in out["rows"]:
        ci = r["first_order_ci"]
        if ci is not None and ci[0] <= 0.0 <= ci[1]:
            assert r["status"] == UNRESOLVED
            assert r["ci_includes_zero"] is True
        elif ci is not None:
            assert r["status"] == RESOLVED
    # the recommendation names only resolved factors, and says the rest are not
    rec = out["recommendation"]
    if out["unresolved_factors"]:
        assert UNRESOLVED in rec
        for name in out["unresolved_factors"]:
            assert name in rec


def test_voi_states_are_the_three_documented_ones(tiny_voi):
    assert set(tiny_voi["factor_states"].values()) <= {UNRESOLVED, RESOLVED, NULL_CONTROL}
    assert tiny_voi["factor_states"]["lif_seed"] == NULL_CONTROL
    assert "lif_seed" in tiny_voi["null_control_factors"]
    assert NEGATIVE_VOI_NOTE in tiny_voi["warnings"]
