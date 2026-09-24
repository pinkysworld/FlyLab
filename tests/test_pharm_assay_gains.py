"""Assays must take their gains from flylab.pharm.mechanisms."""
import pytest

from flylab.assays.taste import run_taste_assay
from flylab.assays.wholens import GABA_SOURCE_CHANGE, run_wholens_assay
from flylab.maps.malecns import CENSUS_FALLBACK
from flylab.pharm.mechanisms import gains_from_occupancy
from flylab.pharm.occupancy import compare_compound

census = pytest.mark.skipif(not CENSUS_FALLBACK.exists(), reason="no census fallback")


def test_taste_gains_match_mechanisms():
    nb = run_taste_assay("imidacloprid", 1e-6)
    expected = gains_from_occupancy(compare_compound("imidacloprid", 1e-6)["receptors"])
    assert nb["gains"] == expected
    assert nb["readouts"]["g_ach"] == expected["g_ach"]


def test_taste_vehicle_is_unity():
    nb = run_taste_assay(None, 0.0)
    assert nb["gains"] == {k: 1.0 for k in nb["gains"]}
    assert nb["readouts"]["mn9_sugar_hz"] == nb["readouts"]["mn9_vehicle_sugar_hz"]


def test_taste_injects_the_gaba_gain_so_an_rdl_blocker_lifts_the_bitter_veto():
    """fipronil used to leave the reduced plot identical to vehicle.

    ``run_taste_circuit`` always accepted ``g_gaba`` (it scales the bitter
    inhibitory unit) but the assay never passed it, so an RDL antagonist read as
    "does nothing to feeding motor output". Lower g_gaba weakens the inhibitory
    unit, so sugar+bitter MN9 must rise and the veto must lift - the direction
    the map-based assay already shows.
    """
    veh = run_taste_assay(None, 0.0)
    nb = run_taste_assay("fipronil", 1e-6)
    assert nb["gains"]["g_gaba"] < 1.0
    assert nb["readouts"]["g_gaba"] == nb["gains"]["g_gaba"]
    assert nb["readouts"]["mn9_sugar_bitter_hz"] > veh["readouts"]["mn9_sugar_bitter_hz"]
    assert nb["readouts"]["bitter_veto_ratio"] > veh["readouts"]["bitter_veto_ratio"]
    # the vehicle arm is frozen: same drive, both gains at 1.0
    assert veh["readouts"]["mn9_sugar_hz"] == 300.0
    assert veh["readouts"]["mn9_sugar_bitter_hz"] == 0.0
    assert nb["readouts"]["mn9_vehicle_sugar_hz"] == veh["readouts"]["mn9_sugar_hz"]
    # no warning may still claim the circuit ignores RDL
    assert not any("injects only g_ach" in w for w in nb["warnings"])


def test_taste_still_warns_about_gains_it_really_cannot_inject():
    nb = run_taste_assay("deltamethrin", 1e-6)
    assert nb["gains"]["g_nav"] > 1.0
    assert any("g_nav" in w and "not applied" in w for w in nb["warnings"])


def test_unsourced_ache_parameter_cannot_change_the_taste_circuit():
    """Withheld chlorpyrifos-oxon potency must yield no insect AChE gain."""
    low = run_taste_assay("chlorpyrifos_oxon", 3e-9, sugar_hz=20.0)
    high = run_taste_assay("chlorpyrifos_oxon", 1e-6, sugar_hz=20.0)
    assert low["gains"]["ach_tone"] == high["gains"]["ach_tone"] == 1.0
    assert low["readouts"]["mn9_sugar_hz"] == low["readouts"]["mn9_vehicle_sugar_hz"]
    assert high["readouts"]["mn9_sugar_hz"] == high["readouts"]["mn9_vehicle_sugar_hz"]


@census
def test_wholens_gains_come_from_insect_rdl():
    fip = run_wholens_assay("fipronil", 1e-6)
    assert fip["gains"]["g_gaba"] < 0.5
    assert fip["gains"]["g_ach"] == 1.0
    assert fip["readouts"]["gains"]["gaba"] == fip["gains"]["g_gaba"]
    assert any(GABA_SOURCE_CHANGE == w for w in fip["warnings"])


@census
def test_wholens_diazepam_no_longer_moves_the_fly_gaba_gain():
    nb = run_wholens_assay("diazepam", 1e-6)
    assert nb["gains"]["g_gaba"] == 1.0
    assert any("no insect RDL row" in w for w in nb["warnings"])


@census
def test_wholens_keeps_v04_readout_keys():
    nb = run_wholens_assay("imidacloprid", 1e-5)
    for key in ("n_traced", "neurotransmitter_counts", "superclass_counts", "named_cells",
                "gains", "vehicle", "treated", "delta_excitation", "weights_present"):
        assert key in nb["readouts"]
    assert set(nb["readouts"]["gains"]) >= {"acetylcholine", "gaba", "glutamate", "other"}
    assert nb["readouts"]["treated"]["cns_excitation_index"] < nb["readouts"]["vehicle"]["cns_excitation_index"]
    assert nb["map"]["name"] == "male-cns:v1.0"


@census
def test_wholens_nav_modulator_raises_excitation():
    nb = run_wholens_assay("deltamethrin", 1e-6)
    assert nb["gains"]["g_nav"] > 1.0
    assert nb["readouts"]["delta_excitation"] > 0
