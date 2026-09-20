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
    EngagementModel,
    EvidenceTypeError,
    ParameterType,
    SourceRelation,
    as_param_type,
    as_relation,
    check_transformation,
    describe,
    is_modelled,
    model_for,
    provenance_warning,
)
from flylab.pharm.occupancy import (
    compare_compound,
    engagement,
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
    # a class-level statement cannot license a physical occupancy claim
    assert model_for("Kd", "class_extrapolation") is EngagementModel.functional_engagement
    # nothing survives an unsupported relation
    assert model_for("Kd", "unsupported") is EngagementModel.not_modelled
    assert model_for("EC50", "unsupported") is EngagementModel.not_modelled
    assert model_for("EC50", "exact_compound_related_receptor") is (
        EngagementModel.functional_engagement
    )


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
            assert model_for(pt, far) is EngagementModel.binding_derived_engagement
        # unchanged by v0.6.1: a class statement supports no binding claim at all
        assert model_for(pt, "class_extrapolation") is EngagementModel.functional_engagement
        assert model_for(pt, "unsupported") is EngagementModel.not_modelled


def test_a_cross_species_binding_constant_may_not_claim_occupancy():
    check_transformation("Kd", EngagementModel.binding_derived_engagement,
                         "exact_compound_exact_receptor_other_species")  # no raise
    check_transformation("Kd", EngagementModel.binding_occupancy,
                         BINDING_OCCUPANCY_RELATION)  # no raise
    with pytest.raises(EvidenceTypeError) as exc:
        check_transformation("Kd", EngagementModel.binding_occupancy,
                             "exact_compound_exact_receptor_other_species")
    assert "exact_compound_exact_receptor_exact_species" in str(exc.value)
    assert "binding_derived_engagement" in str(exc.value)
    with pytest.raises(EvidenceTypeError):
        check_transformation("Ki", EngagementModel.binding_occupancy,
                             "exact_compound_related_receptor")
    # a functional potency may not claim EITHER binding model
    for model in (EngagementModel.binding_occupancy,
                  EngagementModel.binding_derived_engagement):
        with pytest.raises(EvidenceTypeError):
            check_transformation("EC50", model, BINDING_OCCUPANCY_RELATION)
    # and nothing numeric survives an unsupported relation
    with pytest.raises(EvidenceTypeError):
        check_transformation("EC50", EngagementModel.functional_engagement, "unsupported")


def test_a_binding_derived_row_carries_a_warning_naming_the_gap():
    warning = provenance_warning(
        EngagementModel.binding_derived_engagement,
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
    assert record["engagement_model"] == "functional_engagement"
    assert "NOT fractional receptor occupancy" in record["engagement_model_note"] or (
        "not fractional receptor occupancy" in record["engagement_model_note"].lower()
    )
    assert record["relation"].startswith("exact_compound_exact_receptor")
    assert record["species"]
    assert record["evidence_tier"] == "literature_order"
    assert "Ratra" in record["source"]
    assert record["modelled"] is True

    # the aphid Kd: a genuine binding constant, but not measured in the fly
    kd = describe({"receptor": "insect_nAChR_beta1",
                   **LIB["compounds"]["imidacloprid"]["receptors"]["insect_nAChR_beta1"]})
    assert kd["param_type"] == "Kd"
    assert kd["engagement_model"] == "binding_derived_engagement"
    assert "Myzus persicae" in kd["provenance_warning"]
    assert kd["doi"] == "10.1186/1471-2202-12-51"

    # the Drosophila Kd: the one row that earns the word "occupancy"
    dmel = describe({"receptor": "insect_nAChR_native_dmel",
                     **LIB["compounds"]["imidacloprid"]["receptors"]["insect_nAChR_native_dmel"]})
    assert dmel["engagement_model"] == "binding_occupancy"
    assert dmel["provenance_warning"] is None
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
                    "binding_occupancy", "binding_derived_engagement",
                    "functional_engagement",
                }
                if row["engagement_model"] == "binding_derived_engagement":
                    assert row["provenance_warning"]
                else:
                    assert "provenance_warning" not in row


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
    assert _rows_with_model("binding_derived_engagement") == [
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
    assert by_model["binding_derived_engagement"] == 1  # the aphid Kd
    assert by_model["functional_engagement"] == report["n_rows_modelled"] - 2
    assert by_model["binding_occupancy"] + by_model["binding_derived_engagement"] == (
        report["by_param_type"].get("Kd", 0) + report["by_param_type"].get("Ki", 0)
    )
    assert report["note"]
