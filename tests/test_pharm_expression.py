"""Transcriptome-weighted patches: coverage, weights, and the uniform baseline."""
import pytest

from flylab.pharm.expression import (
    EXPRESSION_WARNING,
    LEVEL_MAP,
    coverage,
    expression_table,
    expression_weights,
    run_weighted_subgraph_assay,
    weighted_gains,
)
from flylab.pharm.mechanisms import gains_from_occupancy
from flylab.pharm.occupancy import compare_compound


def test_expression_table_rows_and_warning():
    table = expression_table()
    assert table["rows"]
    for row in table["rows"]:
        assert set(row) >= {
            "superclass", "receptor", "level", "weight", "known",
            "fraction_cells_expressing", "evidence_tier", "source_key", "mapping",
        }
        assert 0.0 <= row["weight"] <= 1.0
        assert row["level"] in set(LEVEL_MAP) | {"unknown"}
        assert row["mapping"] == "cross_atlas_inference"
    assert EXPRESSION_WARNING in table["warnings"]
    assert table["motor_neuron_gap"]["status"] == "not_sourced"
    public_lead = table["motor_neuron_gap"]["public_data_lead"]
    assert public_lead["accession"] == "GSE141807"
    assert public_lead["status"] == "reanalyzed_vnc_wide_cell_classes_unresolved"
    assert public_lead["summary_file"] == "vnc_receptor_detection_GSE141807.yaml"
    assert "not the brain subesophageal zone" in public_lead["scope"]


def test_coverage_is_reported_and_incomplete():
    cov = expression_table()["coverage"]
    pair_fraction = cov["fraction_cell_receptor_pairs_annotated"]
    assert 0.0 < pair_fraction < 1.0
    assert cov["coverage_unit"] == "graph_node_x_receptor_key_pair"
    assert cov["fraction_known_overall"] == pair_fraction  # legacy alias
    assert cov["n_annotated_cell_receptor_pairs"] / cov["n_cell_receptor_pairs"] == pytest.approx(pair_fraction)
    assert 0.0 < cov["fraction_graph_nodes_with_any_class_annotation"] < 1.0
    named = cov["graphs"]["named"]
    assert named["n_cells"] > 1000
    assert 0 <= named["n_cells_with_any_class_annotation"] <= named["n_cells"]
    for receptor in ("insect_nAChR", "insect_RDL", "insect_GluCl"):
        block = named[receptor]
        assert block["n_known"] + block["n_unknown"] == named["n_cells"]
        assert 0.0 <= block["fraction_known"] <= 1.0
    # the dataset has nothing at all for these two
    assert named["insect_OctR"]["fraction_known"] == 0.0
    assert named["insect_AChE"]["fraction_known"] == 0.0


def test_weights_are_in_range_and_unknown_defaults_to_one():
    weights = expression_weights("insect_nAChR")
    assert weights
    assert all(0.0 <= w <= 1.0 for w in weights.values())
    # most MaleCNS superclasses are 'unknown' in the dataset -> weight 1.0
    assert sum(1 for w in weights.values() if w == 1.0) > 0
    octr = expression_weights("insect_OctR")
    assert set(octr.values()) == {1.0}, "no OctR expression is sourced: everything defaults"


def test_default_unknown_is_honoured():
    weights = expression_weights("insect_OctR", default_unknown=0.0)
    assert set(weights.values()) == {0.0}
    half = expression_weights("insect_OctR", default_unknown=0.5)
    assert set(half.values()) == {0.5}
    with pytest.raises(ValueError):
        expression_weights("insect_nAChR", default_unknown=-1.0)


def test_weighted_gains_interpolate_between_one_and_the_full_gain():
    rows = compare_compound("imidacloprid", 1e-7)["receptors"]
    uniform = gains_from_occupancy(rows)
    node_gains = weighted_gains(rows, gains=uniform)
    weights = expression_weights("insect_nAChR")
    assert node_gains
    for body_id, gains in node_gains.items():
        assert set(gains) == {"g_ach", "g_gaba", "g_glu"}
        expected = 1.0 + (uniform["g_ach"] * uniform["ach_tone"] - 1.0) * weights[body_id]
        assert gains["g_ach"] == pytest.approx(expected)
        # a gain lies between "no receptor" (1.0) and the uniform value
        assert min(1.0, uniform["g_ach"]) - 1e-9 <= gains["g_ach"] <= max(1.0, uniform["g_ach"]) + 1e-9


def test_all_weights_one_reproduces_the_uniform_gains():
    rows = compare_compound("imidacloprid", 1e-7)["receptors"]
    uniform = gains_from_occupancy(rows)
    node_gains = weighted_gains(rows, gains=uniform, level_map={k: 1.0 for k in LEVEL_MAP})
    for gains in node_gains.values():
        assert gains["g_ach"] == pytest.approx(uniform["g_ach"] * uniform["ach_tone"])
        assert gains["g_gaba"] == pytest.approx(uniform["g_gaba"])


def test_weighted_assay_differs_from_uniform_by_a_finite_amount():
    nb = run_weighted_subgraph_assay("imidacloprid", 1e-7)
    readouts = nb["readouts"]
    delta = readouts["delta_vs_uniform"]
    assert set(delta) == {"mn9_hz", "dnp01_hz", "mean_hz", "max_hz"}
    assert all(v is not None for v in delta.values())
    assert 0.0 < readouts["max_abs_delta_vs_uniform"] < 1e3
    assert readouts["uniform"]["mn9_hz"] is not None
    assert readouts["mn9_hz"] != readouts["uniform"]["mn9_hz"]


def test_weighted_assay_is_identical_to_uniform_when_every_weight_is_one():
    """The per-node runtime mirrors RateNetwork.run exactly."""
    nb = run_weighted_subgraph_assay(
        "imidacloprid", 1e-7, level_map={k: 1.0 for k in LEVEL_MAP}, default_unknown=1.0
    )
    for value in nb["readouts"]["delta_vs_uniform"].values():
        assert value == pytest.approx(0.0, abs=1e-9)


def test_every_output_carries_the_sensitivity_warning():
    nb = run_weighted_subgraph_assay("imidacloprid", 1e-7)
    assert EXPRESSION_WARNING in nb["warnings"]
    assert any("motor" in w.lower() for w in nb["warnings"])
    assert any("uniform" in w for w in nb["warnings"])
    assert EXPRESSION_WARNING in expression_table()["warnings"]
    assert EXPRESSION_WARNING == coverage()["warning"]


def test_notebook_keeps_the_expression_provenance():
    nb = run_weighted_subgraph_assay("fipronil", 1e-7)
    block = nb["expression"]
    assert block["level_map"] == LEVEL_MAP
    assert block["default_unknown"] == 1.0
    assert block["aggregation"] == "mean"
    assert block["coverage"]["fraction_known_overall"] < 1.0
    assert nb["assay"] == "malecns_neighborhood_expression_weighted"
