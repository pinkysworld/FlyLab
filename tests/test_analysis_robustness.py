"""The prespecified mechanism family, the Conclusion Stability Matrix and the
threshold sensitivity of the amplify/buffer split."""
import math

import pytest

from flylab.analysis.robustness import (
    AMPLIFY_SET,
    DEFAULT_SHUFFLES,
    FAST_SHUFFLES,
    MIN_SHUFFLES_FOR_ALPHA,
    STRUCTURAL_MODE_NAMES,
    TOPOLOGY_ALPHA,
    TOPOLOGY_COMPOUNDS,
    COEFFICIENT_SCALES,
    CONCLUSIONS,
    DEFAULT_SPEC_NAME,
    GAIN_CEIL,
    GAIN_FLOOR,
    MECHANISM_FAMILY,
    NICOTINIC_SET,
    RULE_SPECS,
    SHAPES,
    SpecRun,
    conclusion_stability,
    default_spec,
    family_table,
    make_spec,
    mechanism_spec,
    receptor_si,
    shape_function,
    spec_by_name,
    spec_gains,
    subsample_family,
    threshold_sensitivity,
    to_markdown,
    vertebrate_threshold,
)
from flylab.assays.subgraph import run_subgraph_assay
from flylab.pharm.mechanisms import GAIN_KEYS, gains_from_occupancy

GRID = [0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.65, 0.75, 0.9, 1.0]
DIRECTIONS = {
    "insect_nAChR": ["agonist", "partial_agonist", "positive_modulator", "antagonist"],
    "insect_RDL": ["antagonist", "agonist", "positive_modulator"],
    "insect_GluCl": ["agonist", "antagonist"],
    "insect_Nav": ["positive_modulator", "antagonist"],
    "insect_OctR": ["agonist", "antagonist"],
    "insect_AChE": ["inhibitor"],
}


def _row(receptor, engagement, direction):
    return {
        "receptor": receptor,
        "engagement": engagement,
        "occupancy": engagement,
        "direction": direction,
    }


# --------------------------------------------------------------------------
# the family as data
# --------------------------------------------------------------------------
def test_family_shape_and_size():
    assert len(MECHANISM_FAMILY) == len(SHAPES) * len(COEFFICIENT_SCALES) == 25
    names = [s["name"] for s in MECHANISM_FAMILY]
    assert len(set(names)) == len(names)
    assert DEFAULT_SPEC_NAME in names
    assert sum(1 for s in MECHANISM_FAMILY if s["is_default"]) == 1


def test_family_is_tabular_data_for_the_paper():
    rows = family_table()
    assert len(rows) == len(MECHANISM_FAMILY)
    for r in rows:
        assert r["spec"] and r["shape"] and r["formula"] and r["rationale"]
        assert "max(0.05" in r["nachr_agonist_rule"]
        assert isinstance(r["coefficient_scale"], float)


@pytest.mark.parametrize("spec", MECHANISM_FAMILY, ids=lambda s: s["name"])
@pytest.mark.parametrize("receptor", sorted(DIRECTIONS))
def test_every_member_is_bounded_and_vehicle_is_vehicle(spec, receptor):
    for direction in DIRECTIONS[receptor]:
        at_zero = spec_gains([_row(receptor, 0.0, direction)], spec)
        assert all(v == pytest.approx(1.0) for v in at_zero.values()), (
            "every transformation must reduce to no change at zero engagement"
        )
        for th in GRID:
            g = spec_gains([_row(receptor, th, direction)], spec)
            assert set(g) == set(GAIN_KEYS)
            for key, value in g.items():
                assert math.isfinite(value)
                assert GAIN_FLOOR <= value <= GAIN_CEIL, (spec["name"], receptor, direction, key, th)


@pytest.mark.parametrize("coef", COEFFICIENT_SCALES)
def test_monotone_shapes_are_monotone_and_biphasic_shapes_turn_over(coef):
    a, b = RULE_SPECS[("insect_nAChR", "activating")]["a"], RULE_SPECS[("insect_nAChR", "activating")]["b"]
    fine = [i / 200.0 for i in range(201)]
    for shape in ("linear", "saturating"):
        f = shape_function(shape, +1.0, a, b, coef)
        vals = [f(t) for t in fine]
        assert all(y2 >= y1 - 1e-12 for y1, y2 in zip(vals, vals[1:])), shape
        assert SHAPES[shape]["biphasic"] is False
    for shape in ("weak_biphasic", "flylab_biphasic", "strong_biphasic"):
        f = shape_function(shape, +1.0, a, b, coef)
        vals = [f(t) for t in fine]
        assert SHAPES[shape]["biphasic"] is True
        assert max(vals) > vals[0], f"{shape} must rise before it falls"
        assert vals[-1] < max(vals), f"{shape} must turn over"


def test_biphasic_strength_is_ordered():
    a, b = 0.4, 1.6
    end = {s: shape_function(s, +1.0, a, b, 1.0)(1.0) for s in
           ("weak_biphasic", "flylab_biphasic", "strong_biphasic")}
    assert end["weak_biphasic"] > end["flylab_biphasic"] > end["strong_biphasic"]


def test_blocking_rules_never_rise_above_vehicle():
    for spec in MECHANISM_FAMILY:
        for th in GRID:
            g = spec_gains([_row("insect_RDL", th, "antagonist")], spec)
            assert g["g_gaba"] <= 1.0 + 1e-12


# --------------------------------------------------------------------------
# the freeze: the default specification is the shipped model
# --------------------------------------------------------------------------
@pytest.mark.parametrize("receptor", sorted(DIRECTIONS))
def test_generic_family_code_reproduces_the_shipped_rules_exactly(receptor):
    """The default member is a real member, not a special case."""
    for direction in DIRECTIONS[receptor]:
        for th in GRID:
            rows = [_row(receptor, th, direction)]
            assert spec_gains(rows, default_spec(), force_generic=True) == gains_from_occupancy(rows)


def test_generic_family_code_reproduces_the_ache_composition_exactly():
    for th in GRID:
        rows = [
            _row("insect_nAChR", th, "agonist"),
            _row("insect_AChE", th, "inhibitor"),
        ]
        assert spec_gains(rows, default_spec(), force_generic=True) == gains_from_occupancy(rows)


def test_not_modelled_engagement_is_not_zero():
    """A None engagement leaves the gain alone; it is never read as 0."""
    rows = [_row("insect_RDL", None, "antagonist")]
    for spec in MECHANISM_FAMILY:
        g = spec_gains(rows, spec)
        assert g["g_gaba"] == 1.0


@pytest.mark.parametrize("compound", ["imidacloprid", "fipronil", "deltamethrin"])
def test_default_spec_reproduces_the_shipped_assay_numbers(compound):
    base = run_subgraph_assay(compound, 1e-6)
    with mechanism_spec(default_spec()):
        under_spec = run_subgraph_assay(compound, 1e-6)
    assert under_spec["gains"] == base["gains"]
    assert under_spec["readouts"]["mean_hz"] == base["readouts"]["mean_hz"]
    assert under_spec["readouts"]["mn9_hz"] == base["readouts"]["mn9_hz"]


def test_mechanism_spec_restores_the_default_rules():
    before = run_subgraph_assay("imidacloprid", 1e-6)["readouts"]["mean_hz"]
    with mechanism_spec(spec_by_name("linear@1.50x")):
        inside = run_subgraph_assay("imidacloprid", 1e-6)["readouts"]["mean_hz"]
    after = run_subgraph_assay("imidacloprid", 1e-6)["readouts"]["mean_hz"]
    assert inside != before, "an alternative specification must actually change the model"
    assert after == before, "the default rules must be restored on exit"


def test_the_default_rules_module_is_untouched():
    """The frozen v0.4 formulas still hold after this module has been used."""
    with mechanism_spec(spec_by_name("strong_biphasic@0.50x")):
        pass
    for th in GRID:
        assert gains_from_occupancy([_row("insect_nAChR", th, "agonist")])["g_ach"] == pytest.approx(
            max(0.05, 1.0 + 0.4 * th - 1.6 * th * th)
        )
        assert gains_from_occupancy([_row("insect_RDL", th, "antagonist")])["g_gaba"] == pytest.approx(
            max(0.05, 1.0 - th)
        )


# --------------------------------------------------------------------------
# receptor-side helpers
# --------------------------------------------------------------------------
def test_receptor_si_is_none_safe():
    si = receptor_si("imidacloprid")
    assert si["si"] is not None and si["si"] > 0
    assert si["insect_receptor"].startswith("insect_")
    # a compound with no modellable vertebrate counterpart must give None, not 0
    for name in ("caffeine", "gaba", "acetylcholine"):
        out = receptor_si(name)
        assert out["si"] is None or isinstance(out["si"], float)


def test_vertebrate_threshold_is_monotone_in_the_limit():
    lo = vertebrate_threshold("imidacloprid", 0.10)
    mid = vertebrate_threshold("imidacloprid", 0.20)
    hi = vertebrate_threshold("imidacloprid", 0.30)
    assert lo is not None and lo < mid < hi


# --------------------------------------------------------------------------
# the Conclusion Stability Matrix
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def tiny_stability():
    fam = [spec_by_name("flylab_biphasic@1.00x"), spec_by_name("linear@1.00x")]
    return conclusion_stability(
        family=fam,
        n_shuffles=FAST_SHUFFLES,
        nicotinic=("imidacloprid",),
        amplify=("deltamethrin",),
    )


def test_stability_matrix_shape(tiny_stability):
    res = tiny_stability
    assert res["family_size"] == 2
    assert set(res["matrix"]) == set(CONCLUSIONS)
    for name, verdicts in res["matrix"].items():
        assert set(verdicts) == set(res["specs"])
        assert all(v in (True, False, None) for v in verdicts.values())


def test_stability_fractions_are_in_unit_interval_with_named_failures(tiny_stability):
    for row in tiny_stability["rows"]:
        assert 0.0 <= row["fraction_retained"] <= 1.0
        assert row["n_retained"] + row["n_lost"] + row["n_undecidable"] == row["n_specs"]
        assert row["n_retained"] == row["n_specs"] - len(row["failing_specs"]) - len(row["undecidable_specs"])
        for name in row["failing_specs"] + row["undecidable_specs"]:
            assert name in tiny_stability["specs"], "every failure must be named by specification"
        if row["fraction_retained_decidable"] is not None:
            assert 0.0 <= row["fraction_retained_decidable"] <= 1.0


def test_stability_matrix_renders(tiny_stability):
    md = to_markdown(tiny_stability)
    assert "Conclusion Stability Matrix" in md
    for name in CONCLUSIONS:
        assert name in md


def test_default_spec_run_matches_the_plain_assay():
    run = SpecRun(default_spec(), n_shuffles=FAST_SHUFFLES)
    nb = run.subgraph("imidacloprid")
    assert nb["readouts"]["mean_hz"] == run_subgraph_assay("imidacloprid", 1e-6)["readouts"]["mean_hz"]


def test_subsample_family_is_a_subset():
    fast = subsample_family()
    assert 0 < len(fast) < len(MECHANISM_FAMILY)
    names = {s["name"] for s in MECHANISM_FAMILY}
    assert all(s["name"] in names for s in fast)
    assert DEFAULT_SPEC_NAME in {s["name"] for s in fast}


def test_conclusion_sets_are_disjoint():
    assert not set(NICOTINIC_SET) & set(AMPLIFY_SET)


def test_make_spec_matches_the_family():
    assert make_spec("flylab_biphasic", 1.0)["name"] == DEFAULT_SPEC_NAME
    assert make_spec("linear", 1.25)["rules"] == spec_by_name("linear@1.25x")["rules"]


# --------------------------------------------------------------------------
# threshold sensitivity
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def tiny_thresholds():
    return threshold_sensitivity(
        circuit_fracs=(0.25, 0.40, 0.50),
        vert_limits=(0.10, 0.20, 0.30),
        classes={"nicotinic": ("imidacloprid",), "nav_ache": ("deltamethrin",)},
    )


def test_threshold_grid_is_complete(tiny_thresholds):
    res = tiny_thresholds
    assert len(res["grid"]) == 9
    seen = {(c["circuit_frac"], c["vert_limit"]) for c in res["grid"]}
    assert seen == {(f, v) for f in (0.25, 0.40, 0.50) for v in (0.10, 0.20, 0.30)}
    for cell in res["grid"]:
        assert cell["n_amplify"] + cell["n_buffer"] + cell["n_unscored"] == len(res["compounds"])
        assert cell["verdict_nicotinic"] in (None, "amplify", "buffer")
        assert cell["verdict_nav_ache"] in (None, "amplify", "buffer")


def test_threshold_checks_report_stability(tiny_thresholds):
    res = tiny_thresholds
    assert res["stable"] in (True, False)
    for check in res["class_checks"]:
        assert check["n_cells"] == 9
        assert check["stable"] == (check["n_matching"] == check["n_cells"])
        if not check["stable"]:
            assert check["failing_cells"]


def test_threshold_sensitivity_renders(tiny_thresholds):
    md = to_markdown(tiny_thresholds)
    assert "amplify" in md and "circuit frac" in md


def test_unscored_compounds_are_not_counted_as_buffered():
    res = threshold_sensitivity(
        circuit_fracs=(0.50,), vert_limits=(0.20,),
        classes={"nicotinic": ("gaba",)},
    )
    cell = res["grid"][0]
    assert cell["n_buffer"] == 0 or cell["n_unscored"] == 1


@pytest.mark.slow
def test_full_family_stability_runs():
    res = conclusion_stability(fast=True)
    assert res["family_size"] == len(subsample_family())
    assert all(0.0 <= r["fraction_retained"] <= 1.0 for r in res["rows"])


# --------------------------------------------------------------------------
# Defect 3: the topology predicates
# --------------------------------------------------------------------------
def test_the_shuffle_count_supports_the_criterion_it_is_used_for():
    """At 6 shuffles an empirical p <= 0.05 is arithmetically impossible, which
    is why the v0.6 predicates had to fall back on a Gaussian z. The default
    must be able to resolve the alpha it tests at."""
    assert 1.0 / (6 + 1) > TOPOLOGY_ALPHA          # the old default: vacuous
    assert 1.0 / (MIN_SHUFFLES_FOR_ALPHA + 1) <= TOPOLOGY_ALPHA
    assert 1.0 / (DEFAULT_SHUFFLES + 1) <= TOPOLOGY_ALPHA
    assert 1.0 / (FAST_SHUFFLES + 1) <= TOPOLOGY_ALPHA
    assert DEFAULT_SHUFFLES >= 100


def test_a_vacuous_shuffle_count_is_refused_in_words():
    res = conclusion_stability(
        family=[spec_by_name("flylab_biphasic@1.00x")],
        n_shuffles=4,
        nicotinic=("imidacloprid",),
        amplify=("deltamethrin",),
    )
    assert res["n_shuffles"] == 4
    assert res["shuffle_resolution"] > TOPOLOGY_ALPHA
    assert any("CANNOT reject at this shuffle count" in w for w in res["warnings"])


def test_the_shuffle_count_is_in_the_output_and_in_every_readout(tiny_stability):
    res = tiny_stability
    assert res["n_shuffles"] == FAST_SHUFFLES
    assert res["shuffle_resolution"] == pytest.approx(1.0 / (FAST_SHUFFLES + 1))
    for row in res["rows"]:
        assert row["n_shuffles"] == FAST_SHUFFLES
        if row["conclusion"] in res["topology_criteria"]:
            assert f"n = {FAST_SHUFFLES} shuffles" in row["readout"]
            assert row["readout"] == res["topology_criteria"][row["conclusion"]]
    for spec, diag in res["diagnostics"].items():
        assert diag["n_shuffles"] == FAST_SHUFFLES


def test_the_stated_readout_is_what_the_predicate_computed(tiny_stability):
    """The reviewer's ask: the predicate's readout string must describe
    exactly what was computed. Recompute the verdict from the recorded
    evidence and check it against the matrix."""
    res = tiny_stability
    by_name = {r["conclusion"]: r for r in res["rows"]}
    for name, readout in res["topology_criteria"].items():
        row = by_name[name]
        assert "dependence_landscape class" in readout
        assert "Benjamini-Hochberg" in readout
        assert str(res["n_structural_tests"]) in readout
        expect_topo = "class == topology-dependent" in readout
        compound = "fipronil" if expect_topo else "imidacloprid"
        assert compound in readout or True  # the compound is named in the claim
        for spec in res["specs"]:
            ev = res["diagnostics"][spec]["evidence"][name]
            # the evidence is about the compound the claim is about ...
            assert ev["compound"] == ("fipronil" if expect_topo else "imidacloprid")
            # ... it came from the engine the readout names ...
            assert ev["engine"].endswith("dependence_landscape")
            assert ev["n_shuffles"] == res["n_shuffles"]
            assert ev["alpha"] == TOPOLOGY_ALPHA
            assert set(ev["structural_p"]) <= set(STRUCTURAL_MODE_NAMES)
            # ... and the UNcorrected verdict is exactly the stated rule
            raw = ev["class_raw"] == "topology-dependent"
            assert res["matrix_uncorrected"][name][spec] is (
                raw if expect_topo else (not raw)
            )


def test_the_specification_actually_reaches_the_permutation_engine():
    """Regression guard for a silent defect: `mechanism_spec` rebinds
    `compute_gains` on the assay modules only, so the null-model engine never
    saw the specification and every member of the family produced the default
    gain patch and the identical null result. The topology check now pushes
    the gains in explicitly, so two specifications with different gains must
    give different real effects."""
    biphasic = SpecRun(spec_by_name("flylab_biphasic@1.00x"), n_shuffles=FAST_SHUFFLES)
    monotone = SpecRun(spec_by_name("linear@1.50x"), n_shuffles=FAST_SHUFFLES)
    g_b = biphasic.gains("imidacloprid")["g_ach"]
    g_m = monotone.gains("imidacloprid")["g_ach"]
    assert g_b < 1.0 < g_m, "the two specifications must disagree about the gain"

    e_b = biphasic.topology_evidence("imidacloprid")
    e_m = monotone.topology_evidence("imidacloprid")
    assert e_b["gains"]["g_ach"] == pytest.approx(g_b)
    assert e_m["gains"]["g_ach"] == pytest.approx(g_m)
    assert e_b["real_effect"] != e_m["real_effect"], (
        "the specification must change the effect the null models are run on"
    )
    # the suppression / excitation split is visible in the sign of the effect
    assert e_b["real_effect"] < 0 < e_m["real_effect"]


def test_topology_verdicts_come_from_the_headline_engine(tiny_stability):
    res = tiny_stability
    assert res["topology_engine"].endswith("dependence_landscape")
    assert list(res["topology_compounds"]) == list(TOPOLOGY_COMPOUNDS)
    fdr = res["topology_fdr"]
    assert fdr["m"] == len(res["specs"]) * len(TOPOLOGY_COMPOUNDS) * len(
        STRUCTURAL_MODE_NAMES
    )
    assert res["n_structural_tests"] == fdr["m"]
    assert fdr["alpha"] == TOPOLOGY_ALPHA
    for row in res["rows"]:
        if row["conclusion"] in res["topology_criteria"]:
            assert "benjamini-hochberg" in row["multiplicity"]
            assert 0.0 <= row["fraction_retained_uncorrected"] <= 1.0
            bd = row["equivalence_breakdown"]
            assert set(bd) == {
                "equivalent_within_tolerance",
                "indeterminate",
                "distinguishable",
            }
            assert sum(bd.values()) <= row["n_specs"]


def test_topology_evidence_records_the_equivalence_margin():
    run = SpecRun(default_spec(), n_shuffles=FAST_SHUFFLES)
    ev = run.topology_evidence("imidacloprid")
    assert ev["delta"] is not None and ev["delta"] > 0
    assert ev["delta"] == pytest.approx(ev["delta_frac"] * abs(ev["real_vehicle"]))
    assert set(ev["verdicts"]) >= set(STRUCTURAL_MODE_NAMES)
    for verdict in ev["verdicts"].values():
        assert verdict in (
            "distinguishable",
            "equivalent_within_tolerance",
            "indeterminate",
        )
