"""Monte-Carlo over the log10 of a sourced library value (schema v3).

Placeholder rows have no value, so they get no interval and no point estimate:
"Missing evidence should not become a small quantitative response."
"""
import numpy as np
import pytest

from flylab.pharm.occupancy import load_library
from flylab.pharm.uncertainty import occupancy_ci, sample_library


def test_ci_brackets_the_point_estimate():
    ci = occupancy_ci("imidacloprid", 1e-6, n=400, seed=0)
    modelled = 0
    for receptor, stats in ci["receptors"].items():
        if stats["point"] is None:
            # not modelled: no point, no interval, no number of any kind
            assert stats["p2_5"] is None and stats["p97_5"] is None and stats["mean"] is None
            assert stats["perturbed"] is False
            continue
        modelled += 1
        assert stats["p2_5"] <= stats["point"] <= stats["p97_5"], receptor
        assert 0.0 <= stats["p2_5"] <= stats["p97_5"] <= 1.0
    assert modelled >= 2


def test_wider_sd_gives_wider_interval():
    narrow = occupancy_ci("nicotine", 2e-7, sd_log10=0.1, n=400, seed=1)
    wide = occupancy_ci("nicotine", 2e-7, sd_log10=0.6, n=400, seed=1)
    key = "insect_nAChR"
    w_n = narrow["receptors"][key]["p97_5"] - narrow["receptors"][key]["p2_5"]
    w_w = wide["receptors"][key]["p97_5"] - wide["receptors"][key]["p2_5"]
    assert w_w > w_n > 0


def test_placeholder_rows_report_none_not_a_point_value():
    ci = occupancy_ci("imidacloprid", 1e-6, n=100, seed=0)
    placeholder = ci["receptors"]["insect_RDL"]
    assert placeholder["perturbed"] is False
    assert placeholder["point"] is None
    assert placeholder["p2_5"] is None and placeholder["p97_5"] is None
    assert placeholder["engagement_model"] == "not_modelled"


def test_diazepam_rdl_has_no_interval_at_all():
    """The peer-review example: diazepam/RDL must not be 1e-4 with a CI."""
    ci = occupancy_ci("diazepam", 1e-6, n=100, seed=0)["receptors"]["insect_RDL"]
    assert ci["point"] is None and ci["mean"] is None


def test_seed_is_reproducible():
    a = occupancy_ci("fipronil", 1e-7, n=200, seed=7)
    b = occupancy_ci("fipronil", 1e-7, n=200, seed=7)
    c = occupancy_ci("fipronil", 1e-7, n=200, seed=8)
    assert a["receptors"] == b["receptors"]
    assert a["receptors"]["insect_RDL"] != c["receptors"]["insect_RDL"]


def test_sample_library_perturbs_only_estimates():
    rng = np.random.default_rng(0)
    base = load_library()
    sampled = sample_library(rng, 0.3, base)
    real = sampled["compounds"]["fipronil"]["receptors"]["insect_RDL"]
    placeholder = sampled["compounds"]["fipronil"]["receptors"]["insect_nAChR"]
    assert real["value_M"] != base["compounds"]["fipronil"]["receptors"]["insect_RDL"]["value_M"]
    assert real["ec50_M"] == real["value_M"]  # deprecated alias kept in step
    # schema v3: a placeholder has no number to perturb and gains none
    assert placeholder["value_M"] is None
    assert "ec50_M" not in placeholder
    # the source library is untouched
    assert base["compounds"]["fipronil"]["receptors"]["insect_RDL"]["value_M"] == 3.0e-08


def test_sample_library_is_lognormal_around_the_point():
    rng = np.random.default_rng(3)
    base = load_library()
    draws = [
        sample_library(rng, 0.3, base)["compounds"]["imidacloprid"]["receptors"]["insect_nAChR"]["ec50_M"]
        for _ in range(300)
    ]
    mean_log = np.mean(np.log10(np.array(draws) / 2.0e-08))
    assert abs(mean_log) < 0.1
    assert np.std(np.log10(np.array(draws) / 2.0e-08)) == pytest.approx(0.3, abs=0.06)


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        occupancy_ci("imidacloprid", 1e-6, n=1)
    with pytest.raises(KeyError):
        occupancy_ci("unobtainium", 1e-6)
    with pytest.raises(ValueError):
        sample_library(np.random.default_rng(0), -0.1)
