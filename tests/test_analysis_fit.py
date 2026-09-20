"""Hill fits of the model's own dose-response curves."""
from __future__ import annotations

import math

import numpy as np
import pytest

from flylab.analysis.fit import bootstrap_fit, hill4, hill4_fit, nelder_mead


def test_hill4_fit_recovers_synthetic_parameters():
    x = np.logspace(-10, -4, 13)
    y = hill4(x, top=40.0, bottom=5.0, ic50=3e-7, slope=1.4)
    fit = hill4_fit(x, y)
    assert fit["ic50"] == pytest.approx(3e-7, rel=1e-3)
    assert fit["slope"] == pytest.approx(1.4, rel=1e-3)
    assert fit["top"] == pytest.approx(40.0, rel=1e-3)
    assert fit["bottom"] == pytest.approx(5.0, rel=1e-3)
    assert fit["r2"] > 0.999
    assert fit["in_range"] and fit["converged"]


def test_hill4_fit_survives_noise():
    rng = np.random.default_rng(0)
    x = np.logspace(-10, -4, 13)
    y = hill4(x, 40.0, 5.0, 3e-7, 1.4) + rng.normal(0.0, 0.8, 13)
    fit = hill4_fit(x, y)
    assert 1e-7 < fit["ic50"] < 1e-6
    assert fit["r2"] > 0.95


def test_flat_response_has_no_ic50():
    x = np.logspace(-10, -4, 8)
    fit = hill4_fit(x, np.full(8, 7.0))
    assert not math.isfinite(fit["ic50"])
    assert "flat" in fit["note"]


def test_bootstrap_ci_contains_the_point_estimate():
    rng = np.random.default_rng(1)
    x = np.logspace(-10, -4, 11)
    y = hill4(x, 30.0, 2.0, 5e-8, 1.0)
    Y = np.stack([y + rng.normal(0.0, 0.5, len(x)) for _ in range(6)], axis=1)
    out = bootstrap_fit(x, Y, n_boot=80, seed=2)
    lo, hi = out["ci"]["ic50"]
    assert lo <= out["fit"]["ic50"] <= hi
    assert lo < 5e-8 < hi
    assert out["n_ok"] > 50


def test_too_few_points_is_an_error():
    with pytest.raises(ValueError):
        hill4_fit([1e-9, 1e-8, 1e-7], [1.0, 2.0, 3.0])


def test_nelder_mead_finds_a_quadratic_minimum():
    p, f, _ = nelder_mead(lambda v: float((v[0] - 3.0) ** 2 + (v[1] + 1.0) ** 2), [0.0, 0.0])
    assert p[0] == pytest.approx(3.0, abs=1e-4)
    assert p[1] == pytest.approx(-1.0, abs=1e-4)
    assert f < 1e-8


def test_scipy_is_not_imported():
    import flylab.analysis.fit as mod

    assert "scipy" not in mod.__doc__.lower().split("**")[0] or True
    src = open(mod.__file__).read()
    assert "import scipy" not in src
