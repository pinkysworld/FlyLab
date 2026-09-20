"""The gain patch rules, including the four legacy v0.4 formulas."""
import pytest

from flylab.pharm.mechanisms import (
    GAIN_KEYS,
    MECHANISM_TABLE,
    default_gains,
    gains_from_occupancy,
    mechanism_table_rows,
)
from flylab.pharm.occupancy import compare_compound

GRID = [0.0, 0.05, 0.2, 0.5, 0.75, 0.9, 1.0]


def _row(receptor, occupancy, direction):
    return {"receptor": receptor, "occupancy": occupancy, "direction": direction}


def test_defaults_are_one():
    g = default_gains()
    assert set(g) == set(GAIN_KEYS)
    assert all(v == 1.0 for v in g.values())
    assert gains_from_occupancy([]) == g


@pytest.mark.parametrize("th", GRID)
def test_legacy_nachr_agonist_formula(th):
    legacy = max(0.05, 1.0 + 0.4 * th - 1.6 * th * th)
    assert gains_from_occupancy([_row("insect_nAChR", th, "agonist")])["g_ach"] == pytest.approx(legacy)


@pytest.mark.parametrize("th", GRID)
def test_legacy_nachr_antagonist_formula(th):
    legacy = max(0.05, 1.0 - th)
    assert gains_from_occupancy([_row("insect_nAChR", th, "antagonist")])["g_ach"] == pytest.approx(legacy)


@pytest.mark.parametrize("th", GRID)
def test_legacy_rdl_antagonist_formula(th):
    legacy = max(0.05, 1.0 - th)
    assert gains_from_occupancy([_row("insect_RDL", th, "antagonist")])["g_gaba"] == pytest.approx(legacy)


@pytest.mark.parametrize("th", GRID)
def test_legacy_rdl_agonist_formula(th):
    legacy = max(0.05, 1.0 + 0.4 * th)
    assert gains_from_occupancy([_row("insect_RDL", th, "agonist")])["g_gaba"] == pytest.approx(legacy)


@pytest.mark.parametrize("th", GRID)
def test_new_rules_are_bounded_and_monotone(th):
    glu = gains_from_occupancy([_row("insect_GluCl", th, "agonist")])["g_glu"]
    nav = gains_from_occupancy([_row("insect_Nav", th, "positive_modulator")])["g_nav"]
    oct_ = gains_from_occupancy([_row("insect_OctR", th, "agonist")])["g_oct"]
    assert glu == pytest.approx(max(0.05, 1.0 + 0.8 * th))
    assert nav == pytest.approx(1.0 + 1.5 * th)
    assert oct_ == pytest.approx(1.0 + 0.5 * th)
    assert 0.05 <= glu <= 1.8 and 1.0 <= nav <= 2.5 and 1.0 <= oct_ <= 1.5


def test_ache_inhibitor_raises_tone_then_blocks():
    low = gains_from_occupancy([_row("insect_AChE", 0.25, "inhibitor")])
    high = gains_from_occupancy([_row("insect_AChE", 1.0, "inhibitor")])
    assert low["ach_tone"] == pytest.approx(1.5)
    assert high["ach_tone"] == pytest.approx(3.0)
    assert low["g_ach"] > 1.0           # excitation phase
    assert high["g_ach"] == pytest.approx(0.05)  # desensitisation / block phase
    assert low["g_gaba"] == 1.0


def test_vertebrate_rows_never_move_fly_gains():
    rows = [
        _row("vertebrate_nAChR_a4b2", 1.0, "agonist"),
        _row("vertebrate_GABA_A", 1.0, "antagonist"),
        _row("vertebrate_GlyR", 1.0, "positive_modulator"),
    ]
    assert gains_from_occupancy(rows) == default_gains()


def test_unknown_receptor_and_none_direction_are_ignored():
    rows = [_row("insect_nAChR", 1.0, "none"), _row("insect_TRPA1", 1.0, "agonist")]
    assert gains_from_occupancy(rows) == default_gains()


def test_gains_are_floored_not_zero():
    for direction, receptor in [("antagonist", "insect_nAChR"), ("antagonist", "insect_RDL")]:
        g = gains_from_occupancy([_row(receptor, 1.0, direction)])
        assert min(g.values()) >= 0.05


def test_imidacloprid_is_cholinergic_and_fipronil_is_gabaergic():
    imi = gains_from_occupancy(compare_compound("imidacloprid", 1e-6)["receptors"])
    fip = gains_from_occupancy(compare_compound("fipronil", 1e-6)["receptors"])
    assert imi["g_ach"] < 0.5 and imi["g_gaba"] == 1.0
    assert fip["g_gaba"] < 0.5 and fip["g_ach"] == 1.0


def test_library_compounds_reach_each_gain():
    reached = {
        "g_ach": gains_from_occupancy(compare_compound("imidacloprid", 1e-6)["receptors"])["g_ach"],
        "g_gaba": gains_from_occupancy(compare_compound("fipronil", 1e-6)["receptors"])["g_gaba"],
        "g_glu": gains_from_occupancy(compare_compound("ivermectin", 1e-6)["receptors"])["g_glu"],
        "g_nav": gains_from_occupancy(compare_compound("deltamethrin", 1e-6)["receptors"])["g_nav"],
        "g_oct": gains_from_occupancy(compare_compound("chlordimeform", 1e-4)["receptors"])["g_oct"],
        "ach_tone": gains_from_occupancy(compare_compound("chlorpyrifos_oxon", 1e-7)["receptors"])["ach_tone"],
    }
    for key, value in reached.items():
        assert value != 1.0, f"{key} never leaves 1.0"


def test_mechanism_table_is_renderable():
    rows = mechanism_table_rows()
    assert len(rows) == len(MECHANISM_TABLE) >= 12
    for row in rows:
        assert row["receptor"].startswith("insect_")
        assert row["gain"] in GAIN_KEYS
        assert row["formula"] and row["rationale"] and row["applies_to"]
    import json

    json.dumps(MECHANISM_TABLE)  # must stay JSON-serialisable for the UI
