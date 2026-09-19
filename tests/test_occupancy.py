from flylab.pharm.occupancy import compare_compound, hill_occupancy


def test_zero_concentration():
    assert hill_occupancy(0.0, 1e-6) == 0.0


def test_ec50_is_half():
    assert abs(hill_occupancy(1e-6, 1e-6) - 0.5) < 1e-12


def test_imidacloprid_prefers_insect():
    result = compare_compound("imidacloprid", 1e-6)
    by = {row["receptor"]: row["occupancy"] for row in result["receptors"]}
    assert by["insect_nAChR"] > by["vertebrate_nAChR_a4b2"]
    assert by["insect_nAChR"] > 0.9


def test_diazepam_prefers_gaba():
    result = compare_compound("diazepam", 1e-7)
    by = {row["receptor"]: row["occupancy"] for row in result["receptors"]}
    assert by["vertebrate_GABA_A"] > by["insect_nAChR"]
