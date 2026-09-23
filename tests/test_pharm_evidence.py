"""Typed pharmacological evidence propagation (schema v3).

FlyLab must *refuse* transformations the evidence does not support. The two
peer-review objections this answers:

    "EC50 is a functional potency parameter, not a receptor binding constant.
     IC50 values are also used for some compounds. A Hill response derived from
     EC50 describes normalized functional response, not physical receptor
     occupancy."

    "Missing receptor values are represented by ec50_M = 0.01 even when
     direction: none, so the model reports things like diazepam/RDL = 10^-4.
     Missing evidence should not become a small quantitative response."
"""
import pytest

from flylab.pharm.evidence import (
    ALLOWED,
    BINDING_OCCUPANCY_RELATION,
    MODEL_STRENGTH,
    DISTANCE_LABELS,
    RELATION_DISTANCE,
    TRANSFORMATION_TABLE,
    EvidenceDistance,
    EngagementModel,
    EvidenceTypeError,
    ParameterType,
    SourceRelation,
    as_param_type,
    as_relation,
    check_transformation,
    describe,
    as_distance,
    is_modelled,
    model_for,
    provenance_warning,
    transformation_table_rows,
)
from flylab.pharm.occupancy import (
    compare_compound,
    engagement,
    evidence_distance_table,
    library_report,
    load_library,
    receptor_spec,
)

LIB = load_library()


# ---------------------------------------------------------------------------
# the type system itself
# ---------------------------------------------------------------------------
def test_every_parameter_type_has_exactly_one_allowed_model():
    assert set(ALLOWED) == set(ParameterType)
    assert ALLOWED[ParameterType.Kd] is EngagementModel.binding_occupancy
    assert ALLOWED[ParameterType.Ki] is EngagementModel.binding_occupancy
    for functional in (ParameterType.EC50, ParameterType.IC50, ParameterType.Kb):
        assert ALLOWED[functional] is EngagementModel.functional_engagement
    for none_of_it in (ParameterType.class_order, ParameterType.unknown,
                       ParameterType.relative_potency):
        assert ALLOWED[none_of_it] is EngagementModel.not_modelled


def test_binding_occupancy_is_only_legal_for_kd_and_ki():
    for pt in (ParameterType.Kd, ParameterType.Ki):
        check_transformation(pt, EngagementModel.binding_occupancy)  # no raise
    for pt in (ParameterType.EC50, ParameterType.IC50, ParameterType.Kb):
        with pytest.raises(EvidenceTypeError) as exc:
            check_transformation(pt, EngagementModel.binding_occupancy)
        assert "not a binding constant" in str(exc.value)
        assert "functional" in str(exc.value)


def test_class_order_can_never_produce_a_number():
    for pt in (ParameterType.class_order, ParameterType.unknown,
               ParameterType.relative_potency):
        for model in (EngagementModel.binding_occupancy, EngagementModel.functional_engagement):
            with pytest.raises(EvidenceTypeError):
                check_transformation(pt, model)
        # declining to model is always allowed
        check_transformation(pt, EngagementModel.not_modelled)


def test_model_for_uses_the_relation_only_to_weaken_the_claim():
    assert model_for("Kd", "exact_compound_exact_receptor_exact_species") is (
        EngagementModel.binding_occupancy
    )
    # a class-level statement cannot license any binding claim
    assert model_for("Kd", "class_extrapolation") is (
        EngagementModel.functional_engagement_proxy
    )
    # nothing survives an unsupported relation
    assert model_for("Kd", "unsupported") is EngagementModel.not_modelled
    assert model_for("EC50", "unsupported") is EngagementModel.not_modelled
    # a potency measured on the modelled target is not a proxy; a transferred
    # one is, and says so
    assert model_for("EC50", "exact_compound_exact_receptor_exact_species") is (
        EngagementModel.functional_engagement
    )
    assert model_for("EC50", "exact_compound_related_receptor") is (
        EngagementModel.functional_engagement_proxy
    )
    # distance can only ever weaken the claim, never strengthen it
    for pt in ParameterType:
        previous = MODEL_STRENGTH[ALLOWED[pt]]
        for dist in EvidenceDistance:
            strength = MODEL_STRENGTH[model_for(pt, dist)]
            assert strength <= previous
            previous = strength


# ---------------------------------------------------------------------------
# v0.6.1: the parameter type and the relation decide jointly
# ---------------------------------------------------------------------------
def test_binding_occupancy_needs_the_exact_species_relation():
    """A Kd measured in another species is not occupancy of the fly receptor."""
    assert BINDING_OCCUPANCY_RELATION is (
        SourceRelation.exact_compound_exact_receptor_exact_species
    )
    for pt in ("Kd", "Ki"):
        assert model_for(pt, BINDING_OCCUPANCY_RELATION) is EngagementModel.binding_occupancy
        for far in ("exact_compound_exact_receptor_other_species",
                    "exact_compound_related_receptor"):
            assert model_for(pt, far) is EngagementModel.binding_engagement_proxy
        # a class statement supports no binding claim at all
        assert model_for(pt, "class_extrapolation") is (
            EngagementModel.functional_engagement_proxy
        )
        assert model_for(pt, "unsupported") is EngagementModel.not_modelled


def test_a_cross_species_binding_constant_may_not_claim_occupancy():
    check_transformation("Kd", EngagementModel.binding_engagement_proxy,
                         "exact_compound_exact_receptor_other_species")  # no raise
    check_transformation("Kd", EngagementModel.binding_occupancy,
                         BINDING_OCCUPANCY_RELATION)  # no raise
    with pytest.raises(EvidenceTypeError) as exc:
        check_transformation("Kd", EngagementModel.binding_occupancy,
                             "exact_compound_exact_receptor_other_species")
    assert "exact_compound_exact_receptor_exact_species" in str(exc.value)
    assert "binding_engagement_proxy" in str(exc.value)
    with pytest.raises(EvidenceTypeError):
        check_transformation("Ki", EngagementModel.binding_occupancy,
                             "exact_compound_related_receptor")
    # a functional potency may not claim EITHER binding model
    for model in (EngagementModel.binding_occupancy,
                  EngagementModel.binding_engagement_proxy):
        with pytest.raises(EvidenceTypeError):
            check_transformation("EC50", model, BINDING_OCCUPANCY_RELATION)
    # nor may a transferred potency drop the word "proxy"
    with pytest.raises(EvidenceTypeError) as exc:
        check_transformation("EC50", EngagementModel.functional_engagement, "E1")
    assert "functional_engagement_proxy" in str(exc.value)
    check_transformation("EC50", EngagementModel.functional_engagement_proxy, "E1")
    # and nothing numeric survives an unsupported relation
    with pytest.raises(EvidenceTypeError):
        check_transformation("EC50", EngagementModel.functional_engagement, "unsupported")


def test_a_binding_derived_row_carries_a_warning_naming_the_gap():
    warning = provenance_warning(
        EngagementModel.binding_engagement_proxy,
        "exact_compound_exact_receptor_other_species",
        "Myzus persicae",
    )
    assert "Myzus persicae" in warning and "not in Drosophila melanogaster" in warning
    assert provenance_warning(EngagementModel.binding_occupancy,
                              BINDING_OCCUPANCY_RELATION, "Drosophila") is None
    assert provenance_warning(EngagementModel.functional_engagement,
                              "exact_compound_related_receptor", "x") is None


def test_unknown_names_are_rejected_loudly():
    with pytest.raises(EvidenceTypeError):
        as_param_type("pIC50ish")
    with pytest.raises(EvidenceTypeError):
        as_relation("close enough")
    assert as_param_type("ec50") is ParameterType.EC50
    assert as_relation("class_extrapolation") is SourceRelation.class_extrapolation
    assert as_param_type(None) is ParameterType.unknown


def test_is_modelled_needs_both_a_type_and_a_value():
    assert is_modelled("EC50", "exact_compound_related_receptor", 1e-8)
    assert not is_modelled("EC50", "exact_compound_related_receptor", None)
    assert not is_modelled("unknown", "unsupported", None)


def test_engagement_refuses_to_transform_a_placeholder():
    with pytest.raises(EvidenceTypeError) as exc:
        engagement(1e-6, 0.01, param_type="unknown", relation="unsupported")
    assert "not modelled" in str(exc.value)
    with pytest.raises(EvidenceTypeError):
        engagement(1e-6, 1e-3, param_type="class_order")


def test_describe_is_a_provenance_record():
    spec = LIB["compounds"]["fipronil"]["receptors"]["vertebrate_GABA_A"]
    record = describe({"receptor": "vertebrate_GABA_A", **spec})
    assert record["param_type"] == "IC50"
    assert record["units"] == "M"
    assert record["engagement_model"] == "functional_engagement_proxy"
    assert record["evidence_distance"] == "E1"
    assert "NOT fractional receptor occupancy" in record["engagement_model_note"] or (
        "not fractional receptor occupancy" in record["engagement_model_note"].lower()
    )
    assert record["relation"].startswith("exact_compound_exact_receptor")
    assert record["species"]
    assert record["evidence_tier"] == "literature_order"
    assert "Hainzl" in record["source"]
    assert "multi-species binding-site proxy" in record["source"]
    assert record["modelled"] is True

    # the aphid Kd: a genuine binding constant, but not measured in the fly
    kd = describe({"receptor": "insect_nAChR_beta1",
                   **LIB["compounds"]["imidacloprid"]["receptors"]["insect_nAChR_beta1"]})
    assert kd["param_type"] == "Kd"
    assert kd["engagement_model"] == "binding_engagement_proxy"
    assert kd["evidence_distance"] == "E1"
    assert kd["evidence_distance_label"] == "E1 (cross-species)"
    assert "Myzus persicae" in kd["provenance_warning"]
    assert kd["doi"] == "10.1186/1471-2202-12-51"

    # the Drosophila Kd: the one row that earns the word "occupancy"
    dmel = describe({"receptor": "insect_nAChR_native_dmel",
                     **LIB["compounds"]["imidacloprid"]["receptors"]["insect_nAChR_native_dmel"]})
    assert dmel["engagement_model"] == "binding_occupancy"
    assert dmel["provenance_warning"] is None
    assert dmel["evidence_distance"] == "E0"
    assert dmel["relation"] == "exact_compound_exact_receptor_exact_species"
    assert dmel["doi"] == "10.1046/j.1471-4159.1996.67041669.x"

    placeholder = describe({"receptor": "insect_RDL",
                            **receptor_spec(LIB["compounds"]["diazepam"], "insect_RDL")})
    assert placeholder["modelled"] is False
    assert placeholder["param_value_M"] is None
    assert placeholder["engagement_model"] == "not_modelled"


# ---------------------------------------------------------------------------
# the objection, end to end
# ---------------------------------------------------------------------------
def test_diazepam_at_insect_rdl_is_not_modelled_not_1e_minus_4():
    result = compare_compound("diazepam", 1e-6)
    rows = {r["receptor"]: r for r in result["receptors"]}
    rdl = rows["insect_RDL"]
    assert rdl["engagement"] is None
    assert rdl["occupancy"] is None           # deprecated alias agrees
    assert rdl["param_value_M"] is None       # no 0.01 M hiding anywhere
    assert rdl["ec50_M"] is None
    assert rdl["param_type"] == "unknown"
    assert rdl["engagement_model"] == "not_modelled"
    assert "placeholder" in rdl["not_modelled_reason"]


def test_diazepam_is_absent_from_the_selectivity_numbers():
    pairs = compare_compound("diazepam", 1e-6)["selectivity"]
    gaba = pairs["GABA_A"]
    assert gaba["placeholder"] is True
    assert gaba["ratio"] is None and gaba["reason"] == "placeholder"
    assert gaba["log10_ec50_ratio_vert_over_insect"] is None
    assert gaba["occupancy_difference"] is None
    numeric = [p for p in pairs.values() if p["ratio"] is not None]
    assert numeric == []  # nothing about diazepam is comparable across organisms


def test_every_row_of_every_compound_is_typed_and_none_safe():
    for key in LIB["compounds"]:
        result = compare_compound(key, 1e-6)
        assert result["engagement_is_not_occupancy"]
        for row in result["receptors"]:
            assert row["param_type"] in {p.value for p in ParameterType}
            assert row["relation"] in {r.value for r in SourceRelation}
            if row["engagement"] is None:
                assert row["param_value_M"] is None
                assert row["engagement_model"] == "not_modelled"
            else:
                assert 0.0 <= row["engagement"] <= 1.0
                assert row["engagement_model"] in {
                    "binding_occupancy", "binding_engagement_proxy",
                    "functional_engagement", "functional_engagement_proxy",
                }
                if row["engagement_model"].endswith("_proxy"):
                    assert row["provenance_warning"]
                    assert row["evidence_distance"] in {"E1", "E2", "E3"}
                else:
                    assert "provenance_warning" not in row
                    assert row["evidence_distance"] == "E0"


def _rows_with_model(model):
    return [
        (key, row["receptor"])
        for key in LIB["compounds"]
        for row in compare_compound(key, 1e-6)["receptors"]
        if row["engagement_model"] == model
    ]


def test_only_a_drosophila_kd_row_reports_binding_occupancy():
    """v0.6.1: a Kd is necessary but not sufficient - the species matters too."""
    assert _rows_with_model("binding_occupancy") == [
        ("imidacloprid", "insect_nAChR_native_dmel")
    ]
    # the aphid Kd is still modelled, but it is no longer called an occupancy
    assert _rows_with_model("binding_engagement_proxy") == [
        ("imidacloprid", "insect_nAChR_beta1")
    ]
    beta1 = next(
        r for r in compare_compound("imidacloprid", 1e-6)["receptors"]
        if r["receptor"] == "insect_nAChR_beta1"
    )
    assert beta1["param_type"] == "Kd"
    assert "Myzus persicae" in beta1["provenance_warning"]
    assert "not in Drosophila melanogaster" in beta1["provenance_warning"]


def test_library_report_separates_evidence_from_missing_evidence():
    report = library_report()
    assert report["by_engagement_model"]["not_modelled"] == report["n_rows_not_modelled"]
    by_model = report["by_engagement_model"]
    assert by_model["binding_occupancy"] == 1          # the Drosophila Kd only
    assert by_model["binding_engagement_proxy"] == 1   # the aphid Kd
    assert by_model["binding_occupancy"] + by_model["binding_engagement_proxy"] == (
        report["by_param_type"].get("Kd", 0) + report["by_param_type"].get("Ki", 0)
    )
    assert sum(by_model.values()) == report["n_rows"]
    assert by_model["not_modelled"] == report["by_evidence_distance"]["E4"]
    assert report["note"]


# ---------------------------------------------------------------------------
# the evidence-distance hierarchy as a reportable artefact
# ---------------------------------------------------------------------------
def test_the_transformation_table_is_the_whole_rule():
    assert len(TRANSFORMATION_TABLE) == len(ParameterType) * len(EvidenceDistance)
    for (pt, dist), model in TRANSFORMATION_TABLE.items():
        assert model_for(pt, dist) is model
        assert MODEL_STRENGTH[model] <= MODEL_STRENGTH[ALLOWED[pt]]
        if dist is EvidenceDistance.E4:
            assert model is EngagementModel.not_modelled
        if model is EngagementModel.binding_occupancy:
            assert pt in (ParameterType.Kd, ParameterType.Ki)
            assert dist is EvidenceDistance.E0
    rows = transformation_table_rows()
    assert len(rows) == len(TRANSFORMATION_TABLE)
    assert {r["evidence_distance"] for r in rows} == {d.value for d in EvidenceDistance}
    assert all(r["engagement_model_note"] for r in rows)


def test_distance_is_a_first_class_field_mapped_from_the_relation():
    assert as_distance("exact_compound_exact_receptor_other_species") is EvidenceDistance.E1
    assert as_distance("E2") is EvidenceDistance.E2
    assert as_distance(None) is EvidenceDistance.E4
    assert [d.rank for d in EvidenceDistance] == [0, 1, 2, 3, 4]
    assert set(RELATION_DISTANCE) == set(SourceRelation)
    assert len(set(RELATION_DISTANCE.values())) == len(EvidenceDistance)
    assert DISTANCE_LABELS[EvidenceDistance.E1] == "E1 (cross-species)"


def test_evidence_distance_table_is_the_paper_census():
    table = evidence_distance_table()
    report = library_report()
    assert table["n_rows"] == report["n_rows"]
    for grade, block in table["grades"].items():
        assert block["n_rows"] == report["by_evidence_distance"].get(grade, 0)
        assert sum(block["by_param_type"].values()) == block["n_rows"]
        assert len(block["rows"]) == block["n_rows"]
    # only E0 may carry a physical occupancy, and the census says how many
    for grade, block in table["grades"].items():
        if grade != "E0":
            assert "binding_occupancy" not in block["by_engagement_model"]
    assert table["n_binding_occupancy"] == report["by_engagement_model"].get(
        "binding_occupancy", 0
    )


# ---------------------------------------------------------------------------
# the soundness invariant
# ---------------------------------------------------------------------------
#
# The type system is only a contribution if it holds across the pipeline, not
# just at the call site that happens to check.  The invariant, stated once:
#
#   NO numeric engagement that reaches a readout may originate from a library
#   row whose (param_type, evidence_distance) pair does not admit one, and the
#   model label it travels under must be exactly the one that pair admits.
#
# This test enumerates the whole library and every entry point that turns a
# library row into a number a readout can see.
def _admitted() -> dict[tuple[str, str], EngagementModel]:
    """(compound, receptor) -> the one model the evidence admits."""
    from flylab.pharm.occupancy import spec_value_M

    admitted: dict[tuple[str, str], EngagementModel] = {}
    for key, entry in LIB["compounds"].items():
        for receptor, spec in entry["receptors"].items():
            value = spec_value_M(spec)
            admitted[(key, receptor)] = (
                model_for(spec.get("param_type"), spec.get("relation"))
                if value is not None else EngagementModel.not_modelled
            )
    return admitted


ADMITTED = _admitted()


def _assert_sound(compound: str, receptor: str, value, model, where: str) -> None:
    """One number, one label, checked against what the evidence admits."""
    admits = ADMITTED[(compound, receptor)]
    if value is None:
        return
    assert admits is not EngagementModel.not_modelled, (
        f"{where}: {compound}/{receptor} produced the number {value!r} although "
        f"its evidence admits no numeric engagement"
    )
    if model is not None:
        assert str(model) == admits.value, (
            f"{where}: {compound}/{receptor} reported engagement_model={model!r}, "
            f"but its evidence admits only {admits.value}"
        )
        spec = LIB["compounds"][compound]["receptors"][receptor]
        # the refusal itself must agree, given both facts
        check_transformation(spec.get("param_type"), model, spec.get("relation"))


def test_no_number_reaching_a_readout_comes_from_inadmissible_evidence():
    """The soundness invariant, over the whole library and every entry point."""
    from flylab.assays.subgraph import run_subgraph_assay
    from flylab.assays.taste import run_taste_assay
    from flylab.browser.bridge import build_dashboard
    from flylab.pharm.mechanisms import gains_from_occupancy
    from flylab.pharm.occupancy import engagement_curve, occupancy_curve
    from flylab.pharm.uncertainty import occupancy_ci

    checked = 0
    for compound in sorted(LIB["compounds"]):
        conc = 1e-6

        # 1. the typed comparison itself
        result = compare_compound(compound, conc)
        for row in result["receptors"]:
            _assert_sound(compound, row["receptor"], row["engagement"],
                          row["engagement_model"], "compare_compound")
            checked += 1

        # 2. every point of every dose-response curve
        for point in occupancy_curve(compound, [1e-9, 1e-6, 1e-3]):
            for receptor, value in point["receptors"].items():
                _assert_sound(compound, receptor, value, None, "occupancy_curve")
        for point in engagement_curve(compound, [1e-6]):
            for receptor, value in point["receptors"].items():
                _assert_sound(compound, receptor, value, None, "engagement_curve")

        # 3. the Monte-Carlo interval (point AND both tails)
        ci = occupancy_ci(compound, conc, n=8, seed=0)
        for receptor, block in ci["receptors"].items():
            for field in ("point", "mean", "p2_5", "p97_5"):
                _assert_sound(compound, receptor, block.get(field),
                              block["engagement_model"], f"occupancy_ci.{field}")

        # 4. the circuit assays -- the numbers that actually reach a readout
        for name, assay in (
            ("run_subgraph_assay", run_subgraph_assay(compound, conc)),
            ("run_taste_assay", run_taste_assay(compound, conc)),
        ):
            rows = assay["occupancy"]
            for row in rows:
                _assert_sound(compound, row["receptor"], row.get("engagement"),
                              row.get("engagement_model"), name)
            # and the gain patch, which is the bridge from evidence to circuit:
            # dropping every inadmissible row must not change a single gain
            admissible_rows = [
                r for r in rows
                if ADMITTED[(compound, r["receptor"])] is not EngagementModel.not_modelled
            ]
            assert gains_from_occupancy(rows) == gains_from_occupancy(admissible_rows), (
                f"{name}: {compound} gain patch depends on a row the evidence "
                "does not admit"
            )
            not_modelled_rows = [r for r in rows if r.get("engagement") is None]
            assert gains_from_occupancy(not_modelled_rows) == gains_from_occupancy([]), (
                f"{name}: {compound} not-modelled rows moved a gain off its default"
            )

        # 5. the dashboard, including its headline tiles
        dash = build_dashboard(compound, conc, include_dependence=False,
                               include_ladder=False)
        for row in dash["evidence"]:
            _assert_sound(compound, row["receptor"], row.get("engagement"),
                          row.get("engagement_model"), "build_dashboard")
        for tile in ("insect_engagement", "vertebrate_engagement"):
            block = dash["headline"][tile]
            if block.get("receptor") is not None:
                _assert_sound(compound, block["receptor"], block["value"],
                              block.get("engagement_model"), f"headline.{tile}")

    assert checked == sum(len(e["receptors"]) for e in LIB["compounds"].values())


def test_the_invariant_would_catch_a_violation():
    """The test above is only worth having if it can fail: prove it can.

    A hand-forged row that claims occupancy from a cross-species Kd must be
    refused by the same check the pipeline runs.
    """
    spec = LIB["compounds"]["imidacloprid"]["receptors"]["insect_nAChR_beta1"]
    with pytest.raises(EvidenceTypeError):
        check_transformation(spec["param_type"], EngagementModel.binding_occupancy,
                             spec["relation"])
    placeholder = LIB["compounds"]["diazepam"]["receptors"]["insect_RDL"]
    with pytest.raises(EvidenceTypeError):
        check_transformation(placeholder["param_type"],
                             EngagementModel.functional_engagement_proxy,
                             placeholder["relation"])


def test_the_census_states_whether_the_binding_branch_is_exercised():
    """The paper must not claim an example the library no longer contains."""
    report = library_report()
    branch = report["binding_branch"]
    assert branch["n_rows"] == report["n_rows_binding_occupancy"]
    assert branch["exercised"] is (branch["n_rows"] > 0)
    assert branch["n_rows"] == report["by_engagement_model"].get("binding_occupancy", 0)
    assert report["n_rows_binding_engagement_proxy"] == report[
        "by_engagement_model"
    ].get("binding_engagement_proxy", 0)
    for row in branch["rows"]:
        assert row["evidence_distance"] == "E0"
        assert row["param_type"] in ("Kd", "Ki")
        assert "Drosophila melanogaster" in str(row["species"])
    for row in branch["proxy_rows"]:
        assert row["evidence_distance"] in ("E1", "E2")
    assert branch["note"]
    if not branch["exercised"]:
        assert "EMPTY" in branch["note"]
