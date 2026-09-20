"""One-compartment teaching exposure model."""
import numpy as np
import pytest

from flylab.pharm.exposure import ROUTE_DEFAULTS, exposure_profile


def test_profile_shape_and_peak():
    p = exposure_profile("imidacloprid", dose=0.05, route="feeding", t_h=24, dt_h=0.05)
    assert len(p["t_h"]) == len(p["conc_M"]) == 481
    assert p["cmax"] > 0
    assert 0 < p["tmax"] < 24
    assert p["auc"] > 0
    assert p["conc_M"][0] == pytest.approx(0.0)


def test_decay_is_monotone_after_tmax():
    p = exposure_profile("imidacloprid", dose=0.05, route="feeding")
    t = np.array(p["t_h"])
    c = np.array(p["conc_M"])
    after = c[t >= p["tmax"]]
    assert np.all(np.diff(after) <= 1e-18)
    assert after[-1] < 0.5 * p["cmax"]
    rise = c[t <= p["tmax"]]
    assert np.all(np.diff(rise) >= -1e-18)


def test_occupancy_tracks_concentration():
    p = exposure_profile("imidacloprid", dose=0.05, route="feeding")
    occ = np.array(p["occupancy"]["insect_nAChR"])
    assert len(occ) == len(p["conc_M"])
    assert occ.max() > 0.5
    assert occ[0] == pytest.approx(0.0)
    assert np.argmax(occ) == np.argmax(p["conc_M"])
    vert = np.array(p["occupancy"]["vertebrate_nAChR_a4b2"])
    assert vert.max() < occ.max()


def test_routes_differ_in_the_documented_direction():
    peaks = {
        route: exposure_profile("imidacloprid", dose=0.05, route=route)
        for route in ROUTE_DEFAULTS
    }
    assert peaks["bath"]["cmax"] > peaks["feeding"]["cmax"] > peaks["topical"]["cmax"]
    assert peaks["bath"]["tmax"] < peaks["topical"]["tmax"]


def test_dose_scales_linearly():
    one = exposure_profile("imidacloprid", dose=0.01, route="feeding")
    ten = exposure_profile("imidacloprid", dose=0.1, route="feeding")
    assert ten["cmax"] == pytest.approx(10 * one["cmax"], rel=1e-9)
    assert ten["tmax"] == pytest.approx(one["tmax"])


def test_warnings_say_teaching_order():
    p = exposure_profile("fipronil", dose=0.02, route="topical")
    assert any("teaching-order" in w.lower() for w in p["warnings"])
    assert any("class_placeholder" in w for w in p["warnings"])
    assert any("haemolymph" in w for w in p["warnings"])


def test_custom_params_and_flip_flop_case():
    p = exposure_profile("imidacloprid", dose=0.05, route="feeding", params={"ka": 0.5, "ke": 0.5})
    assert p["cmax"] > 0 and p["tmax"] == pytest.approx(2.0, abs=0.2)
    assert p["params"]["ka_per_h"] == 0.5


def test_bad_inputs_raise():
    with pytest.raises(ValueError):
        exposure_profile("imidacloprid", 0.05, "injection")
    with pytest.raises(ValueError):
        exposure_profile("imidacloprid", -1.0, "feeding")
    with pytest.raises(KeyError):
        exposure_profile("unobtainium", 0.05, "feeding")
