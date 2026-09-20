"""Genotype toggle: per-compound EC50 shifts from resistance_alleles.yaml."""
import pytest

from flylab.pharm.genotype import (
    WILD_TYPE,
    apply_genotype,
    genotype_assay,
    genotype_occupancy,
    genotype_panel,
    list_genotypes,
)
from flylab.pharm.occupancy import compare_compound, load_library

BASE = load_library()


def _ec50(lib, compound, receptor):
    return float(lib["compounds"][compound]["receptors"][receptor]["ec50_M"])


def test_list_genotypes_starts_with_wild_type():
    rows = list_genotypes()
    assert rows[0]["id"] == WILD_TYPE
    ids = {r["id"] for r in rows}
    assert {"rdl_A301S_nlug", "nachr_beta1_R81T_mper", "para_M918T_superkdr"} <= ids
    for row in rows:
        assert set(row) >= {
            "id", "gene", "allele", "species", "target_receptor",
            "compounds_affected", "n_rows", "caveats",
        }


def test_wild_type_is_identity():
    lib = apply_genotype(None, WILD_TYPE)
    assert lib["compounds"] == BASE["compounds"]
    assert lib["genotype"]["id"] == WILD_TYPE
    assert lib["genotype"]["shifted"] == []
    assert compare_compound("imidacloprid", 1e-7, library=lib) == compare_compound("imidacloprid", 1e-7)


def test_shift_applies_only_to_sourced_compounds():
    """R81T has a published imidacloprid shift and no clothianidin number."""
    lib = apply_genotype(None, "nachr_beta1_R81T_mper")
    assert _ec50(lib, "imidacloprid", "insect_nAChR") == pytest.approx(
        _ec50(BASE, "imidacloprid", "insect_nAChR") * 50.0
    )
    for untouched in ("clothianidin", "nitenpyram", "nicotine", "acetamiprid"):
        assert _ec50(lib, untouched, "insect_nAChR") == _ec50(BASE, untouched, "insect_nAChR")
    unshifted = {row["compound"] for row in lib["genotype"]["unshifted"]}
    assert {"clothianidin", "nitenpyram", "nicotine"} <= unshifted
    assert lib["unshifted"] == lib["genotype"]["unshifted"]


def test_rdl_a301s_shifts_gaba_but_barely_moves_fipronil():
    """Garrood 2017: A301S shifts GABA 33.9x and fipronil only 1.1x."""
    lib = apply_genotype(None, "rdl_A301S_nlug")
    assert _ec50(lib, "gaba", "insect_RDL") == pytest.approx(_ec50(BASE, "gaba", "insect_RDL") * 33.9)
    assert _ec50(lib, "fipronil", "insect_RDL") == pytest.approx(
        _ec50(BASE, "fipronil", "insect_RDL") * 1.1
    )
    # picrotoxin was not measured with this allele, so it must not move
    assert _ec50(lib, "picrotoxin", "insect_RDL") == _ec50(BASE, "picrotoxin", "insect_RDL")


def test_para_m918t_shifts_deltamethrin_but_not_ddt():
    """Vais 2000 / Usherwood 2005: 100x on deltamethrin, no meaningful DDT shift."""
    lib = apply_genotype(None, "para_M918T_superkdr")
    assert _ec50(lib, "deltamethrin", "insect_Nav") == pytest.approx(
        _ec50(BASE, "deltamethrin", "insect_Nav") * 100.0
    )
    assert _ec50(lib, "ddt", "insect_Nav") == pytest.approx(_ec50(BASE, "ddt", "insect_Nav"))


def test_super_kdr_states_the_efficacy_limitation():
    lib = apply_genotype(None, "para_M918T_superkdr")
    text = " ".join(lib["genotype"]["warnings"]).lower()
    assert "efficacy" in text and "binding sites" in text


def test_whole_animal_ratio_is_flagged_and_warned():
    lib = apply_genotype(None, "rdl_A302G_dsim")
    row = next(r for r in lib["genotype"]["shifted"] if r["compound"] == "fipronil")
    assert row["shift_kind"] == "whole_animal_RR"
    assert row["fold"] == pytest.approx(20000.0)
    spec = lib["compounds"]["fipronil"]["receptors"]["insect_RDL"]
    assert spec["genotype_shift"]["shift_kind"] == "whole_animal_RR"
    assert "metabolic" in spec["genotype_shift"]["warning"] or "metabolism" in spec["genotype_shift"]["warning"]
    assert any("OVERSTATES" in w for w in lib["genotype"]["warnings"])


def test_shifted_rows_keep_tier_and_source():
    lib = apply_genotype(None, "nachr_beta1_R81T_mper")
    spec = lib["compounds"]["imidacloprid"]["receptors"]["insect_nAChR"]
    assert spec["evidence_tier"] == "literature_order"
    assert "Bass" in spec["genotype_shift"]["source"]
    assert set(spec["genotype_shift"]) >= {"allele", "fold", "source", "assay_type"}
    assert "genotype nachr_beta1_R81T_mper" in spec["source"]


def test_one_shift_per_compound_even_with_two_published_rows():
    """R81T lists a 50x binding shift AND a 1679x whole-animal ratio for imidacloprid."""
    lib = apply_genotype(None, "nachr_beta1_R81T_mper")
    spec = lib["compounds"]["imidacloprid"]["receptors"]["insect_nAChR"]
    assert spec["genotype_shift"]["fold"] == pytest.approx(50.0)
    alternatives = spec["genotype_shift"]["alternatives_not_applied"]
    assert any(a["fold"] == pytest.approx(1679.0) for a in alternatives)


def test_genotype_occupancy_drops_with_a_resistance_allele():
    wt = genotype_occupancy("imidacloprid", 1e-7, WILD_TYPE)
    mut = genotype_occupancy("imidacloprid", 1e-7, "nachr_beta1_R81T_mper")
    wt_occ = {r["receptor"]: r["occupancy"] for r in wt["receptors"]}["insect_nAChR"]
    mut_row = {r["receptor"]: r for r in mut["receptors"]}["insect_nAChR"]
    assert mut_row["occupancy"] < wt_occ
    assert mut_row["delta_occupancy_vs_wt"] < 0
    assert mut["genotype"]["applies_to_this_compound"] is True


def test_genotype_occupancy_says_so_when_the_compound_is_not_shifted():
    result = genotype_occupancy("clothianidin", 1e-7, "nachr_beta1_R81T_mper")
    assert result["genotype"]["applies_to_this_compound"] is False
    assert any("NOT shifted" in w for w in result["warnings"])


def test_genotype_assay_carries_block_and_warnings():
    nb = genotype_assay("subgraph", "imidacloprid", 1e-7, "nachr_beta1_R81T_mper")
    assert nb["genotype"]["id"] == "nachr_beta1_R81T_mper"
    assert nb["readouts"]["mn9_hz"] is not None
    assert any("per compound" in w for w in nb["warnings"])


def test_genotype_assay_taste_uses_the_shifted_library():
    wt = genotype_assay("taste", "imidacloprid", 1e-7, WILD_TYPE)
    mut = genotype_assay("taste", "imidacloprid", 1e-7, "nachr_beta1_R81T_mper")
    assert wt["gains"]["g_ach"] != mut["gains"]["g_ach"]
    # and the swap is undone afterwards
    from flylab.assays.taste import run_taste_assay

    assert run_taste_assay("imidacloprid", 1e-7)["gains"]["g_ach"] == wt["gains"]["g_ach"]


def test_genotype_panel_for_fipronil_lists_every_rdl_allele():
    panel = genotype_panel("fipronil", 1e-7, assay=None)
    ids = [row["genotype"] for row in panel["rows"]]
    assert ids[0] == WILD_TYPE
    assert {"rdl_A302S_dmel", "rdl_A301S_nlug", "rdl_A302G_dsim"} <= set(ids)
    by_id = {row["genotype"]: row for row in panel["rows"]}
    # no sourced fipronil number for A302S: identical to wild type, and said so
    assert by_id["rdl_A302S_dmel"]["shifted"] is False
    assert by_id["rdl_A302S_dmel"]["occupancy"] == pytest.approx(by_id[WILD_TYPE]["occupancy"])
    assert "no sourced shift" in by_id["rdl_A302S_dmel"]["note"]
    # the 20000x whole-animal ratio essentially abolishes occupancy
    assert by_id["rdl_A302G_dsim"]["occupancy"] < 0.01


def test_genotype_panel_with_a_circuit_assay_reports_the_delta():
    panel = genotype_panel("imidacloprid", 1e-7, assay="subgraph")
    by_id = {row["genotype"]: row for row in panel["rows"]}
    mutant = by_id["nachr_beta1_R81T_mper"]
    assert mutant["readouts"]["mn9_hz"] is not None
    assert "mn9_hz" in mutant["circuit_effect_vs_wt"]


def test_unknown_genotype_raises():
    with pytest.raises(KeyError):
        apply_genotype(None, "not_an_allele")
