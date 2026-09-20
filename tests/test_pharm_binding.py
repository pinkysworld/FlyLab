"""Gaddum / Schild / Black-and-Leff helpers."""
import numpy as np
import pytest

from flylab.pharm.binding import (
    competitive_occupancy,
    effective_ach_gain,
    operational_response,
    schild_shift,
)
from flylab.pharm.occupancy import compare_compound, hill_occupancy


def test_no_competitor_reduces_to_hill():
    occ = competitive_occupancy(1e-6, 1e-6, 0.0, 1e-6)
    assert float(occ) == pytest.approx(hill_occupancy(1e-6, 1e-6))


def test_competitor_shifts_the_curve_right():
    plain = float(competitive_occupancy(1e-6, 1e-6))
    shifted = float(competitive_occupancy(1e-6, 1e-6, 1e-5, 1e-6))
    assert shifted < plain
    # occupancy at Kd*(1+[B]/Kb) is back to 0.5
    ratio = float(schild_shift(1e-5, 1e-6))
    assert float(competitive_occupancy(1e-6 * ratio, 1e-6, 1e-5, 1e-6)) == pytest.approx(0.5)


def test_schild_dose_ratio():
    assert float(schild_shift(0.0, 1e-6)) == 1.0
    assert float(schild_shift(1e-6, 1e-6)) == pytest.approx(2.0)
    b = np.array([1e-8, 1e-7, 1e-6])
    r = schild_shift(b, 1e-7)
    assert r.shape == b.shape
    slope = np.polyfit(np.log10(b), np.log10(r - 1.0), 1)[0]
    assert slope == pytest.approx(1.0, abs=0.05)  # Schild slope of 1


def test_operational_model_is_bounded_and_ordered():
    concs = np.logspace(-10, -3, 30)
    full = operational_response(concs, 1e-6, tau=10.0)
    partial = operational_response(concs, 1e-6, tau=0.5)
    assert np.all((full >= 0) & (full <= 1))
    assert np.all(np.diff(full) >= -1e-12)     # monotone in concentration
    assert full[-1] > partial[-1]              # higher efficacy, higher ceiling
    assert np.all(operational_response(concs, 1e-6, tau=0.0) == 0.0)


def test_receptor_reserve_moves_ec50_below_ka():
    concs = np.logspace(-10, -3, 400)
    resp = operational_response(concs, 1e-6, tau=10.0)
    ec50 = concs[int(np.argmin(np.abs(resp - resp.max() / 2)))]
    assert ec50 < 1e-6  # spare receptors: half-maximal response below Ka


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        competitive_occupancy(-1.0, 1e-6)
    with pytest.raises(ValueError):
        schild_shift(1e-6, 0.0)
    with pytest.raises(ValueError):
        operational_response(1e-6, 0.0, 1.0)


def test_effective_ach_gain_scales_with_tone():
    rows = compare_compound("imidacloprid", 1e-9)["receptors"]
    base = effective_ach_gain(rows)
    doubled = effective_ach_gain(rows, ach_tone=2.0)
    assert doubled == pytest.approx(2.0 * base)
    assert effective_ach_gain(rows, ach_tone=1e-6) == 0.05  # floor
