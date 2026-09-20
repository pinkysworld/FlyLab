"""Typed pharmacological evidence propagation (schema v3).

Why this module exists (peer review, v0.5 -> v0.6), quoted so the intent
survives every later refactor:

    "The manuscript defines theta = C^n/(EC50^n + C^n) as fractional receptor
    occupancy. EC50 is a functional potency parameter, not a receptor binding
    constant. IC50 values are also used for some compounds. A Hill response
    derived from EC50 describes normalized functional response, not physical
    receptor occupancy."

    "Missing receptor values are represented by ec50_M = 0.01 even when
    direction: none, so the model reports things like diazepam/RDL = 10^-4.
    Missing evidence should not become a small quantitative response."

The answer to the first objection is *typing*: every library row records the
kind of parameter its source actually measured (:class:`ParameterType`) and how
far that source is from "this compound, this receptor, this species"
(:class:`SourceRelation`).  From those two facts, and only from them, FlyLab
derives which transformation it is allowed to apply (:class:`EngagementModel`):

* ``binding_occupancy`` -- C^n/(Kd^n + C^n) is fractional occupancy of the
  binding site *of the receptor FlyLab is modelling*.
* ``binding_engagement_proxy`` -- the same Hill expression on a binding
  constant that was not measured on the modelled target.
* ``functional_engagement`` -- the same algebra on a functional potency
  (``EC50``, ``IC50``, ``Kb``) measured on the modelled target: a *normalised
  functional response*, not a physical occupancy.
* ``functional_engagement_proxy`` -- a functional potency transferred from
  another species, a related preparation or a chemical class.
* ``not_modelled`` -- no number may be produced at all.  This is the answer to
  the second objection: a ``class_order``/``unknown`` row, or a row whose
  relation is ``unsupported``, yields ``None``, never 1e-4.

The evidence-distance hierarchy (v0.6.1)
----------------------------------------

The two facts decide the model *jointly*, and the second of them is now a
first-class, ordered quantity rather than a string: :class:`EvidenceDistance`
grades how far the source sits from "this compound, at this receptor, in
*Drosophila melanogaster*".  :data:`RELATION_DISTANCE` maps the existing
:class:`SourceRelation` vocabulary onto it, so nothing downstream breaks and
every row can report ``evidence_distance: "E1"`` beside its relation string:

===  ==========================================  ===================================
E    ``SourceRelation``                          meaning
===  ==========================================  ===================================
E0   exact_compound_exact_receptor_exact_species this compound, this receptor, this
                                                 species -- nothing is transferred
E1   exact_compound_exact_receptor_other_species the right receptor, another organism
E2   exact_compound_related_receptor             a related/hybrid/native-mixed
                                                 preparation
E3   class_extrapolation                         a statement about the chemical class
E4   unsupported                                 nothing supports a value here
===  ==========================================  ===================================

:data:`TRANSFORMATION_TABLE` is then the whole of the rule -- the permitted
transformation as a function of (parameter type, evidence distance):

=================  ====  ==============================
``param_type``     dist  ``model_for`` result
=================  ====  ==============================
Kd / Ki            E0    ``binding_occupancy``
Kd / Ki            E1    ``binding_engagement_proxy``
Kd / Ki            E2    ``binding_engagement_proxy``
Kd / Ki            E3    ``functional_engagement_proxy``
EC50 / IC50 / Kb   E0    ``functional_engagement``
EC50 / IC50 / Kb   E1-E3 ``functional_engagement_proxy``
anything           E4    ``not_modelled``
relative_potency / class_order / unknown, any distance ``not_modelled``
=================  ====  ==============================

Distance alone never *raises* a claim: it can only move a row down the ladder
of :data:`MODEL_STRENGTH`.  E1 and E2 differ in what the proxy's provenance
warning has to name -- the species gap, or the preparation gap -- not in the
arithmetic (:func:`provenance_warning`).

The novel part is that FlyLab **refuses** invalid transformations rather than
silently performing them: :func:`check_transformation` raises
:class:`EvidenceTypeError` when asked for binding occupancy from an EC50, from
a cross-species Kd (the v0.6.1 fix: "a Kd measured in an aphid still yielded a
model the paper described as physical receptor occupancy"), or for any numeric
engagement from a placeholder row.
"""

from __future__ import annotations

from enum import Enum
from typing import Any, Mapping

__all__ = [
    "ParameterType",
    "SourceRelation",
    "EvidenceDistance",
    "EngagementModel",
    "EvidenceTypeError",
    "ALLOWED",
    "BINDING_MODELS",
    "PROXY_MODELS",
    "BINDING_OCCUPANCY_RELATION",
    "BINDING_OCCUPANCY_DISTANCE",
    "MODEL_STRENGTH",
    "RELATION_DISTANCE",
    "DISTANCE_RELATION",
    "DISTANCE_NOTES",
    "DISTANCE_LABELS",
    "TRANSFORMATION_TABLE",
    "BINDING_PARAM_TYPES",
    "FUNCTIONAL_PARAM_TYPES",
    "NOT_MODELLED_PARAM_TYPES",
    "PARAM_TYPE_NOTES",
    "RELATION_NOTES",
    "as_param_type",
    "as_relation",
    "as_distance",
    "distance_for",
    "model_for",
    "check_transformation",
    "is_modelled",
    "describe",
    "engagement_note",
    "provenance_warning",
    "transformation_table_rows",
]


class ParameterType(str, Enum):
    """What the cited source actually measured."""

    Kd = "Kd"                                # equilibrium dissociation constant (binding)
    Ki = "Ki"                                # inhibition constant (binding, competition)
    EC50 = "EC50"                            # half-maximal *effective* concentration (function)
    IC50 = "IC50"                            # half-maximal *inhibitory* concentration (function)
    Kb = "Kb"                                # antagonist dissociation constant from a Schild/Gaddum fit
    relative_potency = "relative_potency"    # "n-fold weaker than X"; no concentration scale of its own
    class_order = "class_order"              # "compounds of this class act in the micromolar range"
    unknown = "unknown"                      # no usable number exists (placeholder row)

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return self.value


class SourceRelation(str, Enum):
    """How far the source is from *this compound, this receptor, this species*."""

    exact_compound_exact_receptor_exact_species = "exact_compound_exact_receptor_exact_species"
    exact_compound_exact_receptor_other_species = "exact_compound_exact_receptor_other_species"
    exact_compound_related_receptor = "exact_compound_related_receptor"
    class_extrapolation = "class_extrapolation"
    unsupported = "unsupported"

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return self.value


class EvidenceDistance(str, Enum):
    """How far the source sits from the modelled target, as an ordered grade.

    The same fact as :class:`SourceRelation`, but ordered and short enough to
    print: a row reports ``relation`` (the vocabulary) *and*
    ``evidence_distance`` (the grade), and the transformation rule is written
    against the grade.  ``E0`` transfers nothing; every step away transfers one
    more assumption.
    """

    E0 = "E0"   # this compound, this receptor, this species
    E1 = "E1"   # this compound and receptor, another species
    E2 = "E2"   # this compound, a related / hybrid / native-mixed preparation
    E3 = "E3"   # a chemical-class statement, not a measurement here
    E4 = "E4"   # nothing supports a value at this receptor

    @property
    def rank(self) -> int:
        """0 for ``E0`` .. 4 for ``E4``: bigger means further away."""
        return int(self.value[1:])

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return self.value


class EngagementModel(str, Enum):
    """What FlyLab is allowed to compute from a row."""

    binding_occupancy = "binding_occupancy"
    binding_engagement_proxy = "binding_engagement_proxy"
    functional_engagement = "functional_engagement"
    functional_engagement_proxy = "functional_engagement_proxy"
    not_modelled = "not_modelled"

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return self.value


class EvidenceTypeError(ValueError):
    """Raised when a transformation is not supported by the evidence type."""


#: The **strongest** model each parameter type may drive, before the relation
#: is taken into account.  Binding occupancy is legal for Kd/Ki alone -- and
#: even for those it is only reached when the relation is
#: :data:`BINDING_OCCUPANCY_RELATION`; see :func:`model_for`.  Read this table
#: as a ceiling on the claim, never as the claim itself.
ALLOWED: dict[ParameterType, EngagementModel] = {
    ParameterType.Kd: EngagementModel.binding_occupancy,
    ParameterType.Ki: EngagementModel.binding_occupancy,
    ParameterType.EC50: EngagementModel.functional_engagement,
    ParameterType.IC50: EngagementModel.functional_engagement,
    ParameterType.Kb: EngagementModel.functional_engagement,
    ParameterType.relative_potency: EngagementModel.not_modelled,
    ParameterType.class_order: EngagementModel.not_modelled,
    ParameterType.unknown: EngagementModel.not_modelled,
}

#: The evidence-distance hierarchy, as a map from the existing relation
#: vocabulary.  The relation names are unchanged (nothing downstream breaks);
#: the distance is the ordered grade the rule is written against.
RELATION_DISTANCE: dict[SourceRelation, EvidenceDistance] = {
    SourceRelation.exact_compound_exact_receptor_exact_species: EvidenceDistance.E0,
    SourceRelation.exact_compound_exact_receptor_other_species: EvidenceDistance.E1,
    SourceRelation.exact_compound_related_receptor: EvidenceDistance.E2,
    SourceRelation.class_extrapolation: EvidenceDistance.E3,
    SourceRelation.unsupported: EvidenceDistance.E4,
}

#: Inverse of :data:`RELATION_DISTANCE` (the map is one-to-one).
DISTANCE_RELATION: dict[EvidenceDistance, SourceRelation] = {
    d: r for r, d in RELATION_DISTANCE.items()
}

#: What each grade means, in one clause, for a dashboard chip or a table cell.
DISTANCE_NOTES: dict[EvidenceDistance, str] = {
    EvidenceDistance.E0: "on-target: this compound, this receptor, this species",
    EvidenceDistance.E1: "cross-species: the right receptor, another organism",
    EvidenceDistance.E2: "related receptor: a hybrid, chimeric, native-mixed or "
                         "orthologous preparation",
    EvidenceDistance.E3: "class extrapolation: a statement about the chemical class, "
                         "not a measurement at this receptor",
    EvidenceDistance.E4: "unsupported: no source supports any value here",
}

#: Two-word labels, e.g. ``"E1 (cross-species)"``.
DISTANCE_LABELS: dict[EvidenceDistance, str] = {
    EvidenceDistance.E0: "E0 (on-target)",
    EvidenceDistance.E1: "E1 (cross-species)",
    EvidenceDistance.E2: "E2 (related receptor)",
    EvidenceDistance.E3: "E3 (class extrapolation)",
    EvidenceDistance.E4: "E4 (unsupported)",
}

#: The one relation, and the one distance, that license a *physical occupancy*
#: claim: the source must have measured this compound, at this receptor, in
#: this species.  Anything further away gives a proxy model.
BINDING_OCCUPANCY_RELATION = SourceRelation.exact_compound_exact_receptor_exact_species
BINDING_OCCUPANCY_DISTANCE = EvidenceDistance.E0

#: The two models that may only ever be derived from a Kd/Ki.
BINDING_MODELS = frozenset(
    {EngagementModel.binding_occupancy, EngagementModel.binding_engagement_proxy}
)

#: The two models that mean "transferred from somewhere else": every row using
#: one carries a :func:`provenance_warning`.
PROXY_MODELS = frozenset(
    {
        EngagementModel.binding_engagement_proxy,
        EngagementModel.functional_engagement_proxy,
    }
)

#: How strong a claim each model makes.  Evidence distance may only ever move a
#: row *down* this ladder (see :func:`model_for`), and
#: :func:`check_transformation` refuses any request above the row's ceiling.
MODEL_STRENGTH: dict[EngagementModel, int] = {
    EngagementModel.not_modelled: 0,
    EngagementModel.functional_engagement_proxy: 1,
    EngagementModel.functional_engagement: 2,
    EngagementModel.binding_engagement_proxy: 3,
    EngagementModel.binding_occupancy: 4,
}


def _transformation(param_type: ParameterType, distance: EvidenceDistance) -> EngagementModel:
    """The rule itself, written once, as (parameter type, distance) -> model."""
    ceiling = ALLOWED[param_type]
    if distance is EvidenceDistance.E4 or ceiling is EngagementModel.not_modelled:
        # No evidence to transform, or a parameter that carries no scale.
        return EngagementModel.not_modelled
    if ceiling is EngagementModel.binding_occupancy:        # Kd / Ki
        if distance is EvidenceDistance.E0:
            return EngagementModel.binding_occupancy
        if distance in (EvidenceDistance.E1, EvidenceDistance.E2):
            return EngagementModel.binding_engagement_proxy
        # E3: a class-level statement supports no binding claim at all
        return EngagementModel.functional_engagement_proxy
    if distance is EvidenceDistance.E0:                     # EC50 / IC50 / Kb
        return EngagementModel.functional_engagement
    return EngagementModel.functional_engagement_proxy


#: The whole rule, materialised: ``TRANSFORMATION_TABLE[(param_type, distance)]``.
#: The paper's methods section prints this table; :func:`model_for` reads it.
TRANSFORMATION_TABLE: dict[tuple[ParameterType, EvidenceDistance], EngagementModel] = {
    (pt, d): _transformation(pt, d) for pt in ParameterType for d in EvidenceDistance
}

BINDING_PARAM_TYPES = frozenset(
    p for p, m in ALLOWED.items() if m is EngagementModel.binding_occupancy
)
FUNCTIONAL_PARAM_TYPES = frozenset(
    p for p, m in ALLOWED.items() if m is EngagementModel.functional_engagement
)
NOT_MODELLED_PARAM_TYPES = frozenset(
    p for p, m in ALLOWED.items() if m is EngagementModel.not_modelled
)

PARAM_TYPE_NOTES: dict[ParameterType, str] = {
    ParameterType.Kd: "equilibrium dissociation constant from a binding assay; "
                      "C^n/(Kd^n+C^n) is physical fractional occupancy OF THE "
                      "RECEPTOR IT WAS MEASURED ON, so it is an occupancy of the "
                      "modelled target only at evidence distance E0",
    ParameterType.Ki: "inhibition constant from a competition binding assay; "
                      "treated as a binding constant, with the same E0 condition "
                      "as Kd before it may be called an occupancy here",
    ParameterType.EC50: "half-maximal effective concentration of a functional response; "
                        "the Hill curve is a normalised response, NOT occupancy",
    ParameterType.IC50: "half-maximal inhibitory concentration of a functional or "
                        "displacement assay; assay-condition dependent, NOT occupancy",
    ParameterType.Kb: "antagonist dissociation constant from a Schild/Gaddum analysis; "
                      "a functional estimate, reported as an engagement",
    ParameterType.relative_potency: "a fold-difference with no concentration scale; "
                                    "cannot be put through a Hill equation",
    ParameterType.class_order: "a statement about a chemical class, not a measurement "
                               "at this receptor; no number is produced",
    ParameterType.unknown: "no usable number exists for this compound at this receptor",
}

RELATION_NOTES: dict[SourceRelation, str] = {
    SourceRelation.exact_compound_exact_receptor_exact_species:
        "this compound, this receptor, Drosophila melanogaster",
    SourceRelation.exact_compound_exact_receptor_other_species:
        "this compound and this receptor, but in another species",
    SourceRelation.exact_compound_related_receptor:
        "this compound on a hybrid, chimeric, native-mixed or orthologous receptor "
        "preparation that is not exactly this key",
    SourceRelation.class_extrapolation:
        "an order-of-magnitude statement about the chemical class or a contrast "
        "derived from a different receptor; never a measured constant here",
    SourceRelation.unsupported:
        "no source supports any value at this receptor",
}


def as_param_type(value: Any) -> ParameterType:
    """Coerce a string / enum to :class:`ParameterType` (case-insensitive)."""
    if isinstance(value, ParameterType):
        return value
    if value is None:
        return ParameterType.unknown
    text = str(value).strip()
    for member in ParameterType:
        if member.value.lower() == text.lower():
            return member
    raise EvidenceTypeError(
        f"unknown parameter type {value!r}; allowed: "
        + ", ".join(m.value for m in ParameterType)
    )


def as_relation(value: Any) -> SourceRelation:
    """Coerce a string / enum to :class:`SourceRelation` (case-insensitive)."""
    if isinstance(value, SourceRelation):
        return value
    if value is None:
        return SourceRelation.unsupported
    text = str(value).strip()
    for member in SourceRelation:
        if member.value.lower() == text.lower():
            return member
    raise EvidenceTypeError(
        f"unknown source relation {value!r}; allowed: "
        + ", ".join(m.value for m in SourceRelation)
    )


def as_distance(value: Any) -> EvidenceDistance:
    """Coerce a relation name, an ``"E1"``-style grade or an enum to a distance."""
    if isinstance(value, EvidenceDistance):
        return value
    if value is None:
        return EvidenceDistance.E4
    text = str(value).strip()
    for member in EvidenceDistance:
        if member.value.lower() == text.lower():
            return member
    return RELATION_DISTANCE[as_relation(text)]


def distance_for(relation: Any) -> EvidenceDistance:
    """The evidence distance of a ``relation`` (alias of :func:`as_distance`)."""
    return as_distance(relation)


def model_for(param_type: Any, relation: Any = None) -> EngagementModel:
    """The engagement model a row is allowed to drive.

    This is :data:`TRANSFORMATION_TABLE` looked up: the parameter type sets a
    *ceiling* (``ALLOWED[pt]``) and the evidence distance can only move the row
    down :data:`MODEL_STRENGTH` from there.  ``relation`` may be a
    :class:`SourceRelation`, its string, or an :class:`EvidenceDistance`
    (``"E1"``).

    * ``E4`` / ``unsupported`` -> ``not_modelled`` whatever the parameter type
      (there is no evidence to transform).
    * a ``Kd``/``Ki`` reaches ``binding_occupancy`` **only** at ``E0``.  This
      compound, this receptor, this species is what makes "fractional occupancy
      of the modelled receptor" a statement about the modelled receptor.
    * a ``Kd``/``Ki`` at ``E1``/``E2`` -> ``binding_engagement_proxy``: the same
      Hill number, labelled as a proxy derived from a binding constant measured
      elsewhere, and carrying a provenance warning naming the gap.
    * a ``Kd``/``Ki`` at ``E3`` -> ``functional_engagement_proxy``: a
      class-level statement licenses no binding claim at all.
    * a functional potency at ``E0`` -> ``functional_engagement``; at
      ``E1``-``E3`` -> ``functional_engagement_proxy``.

    Called without a relation it returns the ceiling, i.e. ``ALLOWED[pt]``.
    """
    pt = as_param_type(param_type)
    if relation is None:
        return ALLOWED[pt]
    return TRANSFORMATION_TABLE[(pt, as_distance(relation))]


def check_transformation(param_type: Any, model: Any, relation: Any = None) -> None:
    """Raise :class:`EvidenceTypeError` if ``model`` may not be derived from
    ``param_type`` (and, when given, from ``relation`` / its evidence distance).

    This is the refusal the peer review asked for: FlyLab does not quietly
    relabel a functional potency as an occupancy, it does not call a Kd
    measured in another species an occupancy of the *Drosophila* receptor, and
    it does not turn a placeholder into a small number.

    Without ``relation`` the check is the type-level ceiling ``ALLOWED[pt]``;
    with it, the ceiling is the row's actual model, :func:`model_for`.  A
    *weaker* model than the ceiling is always allowed -- declining to claim as
    much as the evidence permits is never an error.
    """
    pt = as_param_type(param_type)
    want = model if isinstance(model, EngagementModel) else EngagementModel(str(model))
    dist = as_distance(relation) if relation is not None else None
    ceiling = model_for(pt, dist) if dist is not None else ALLOWED[pt]
    if MODEL_STRENGTH[want] <= MODEL_STRENGTH[ceiling]:
        return
    where = f" at evidence distance {DISTANCE_LABELS[dist]}" if dist is not None else ""
    reason = DISTANCE_NOTES[dist] if dist is not None else PARAM_TYPE_NOTES[pt]
    if ceiling is EngagementModel.not_modelled:
        raise EvidenceTypeError(
            f"{pt.value}{where} cannot drive engagement_model={want.value}: "
            f"{reason}. A row of this kind reports N/A (not modelled), never a "
            "number."
        )
    if want in BINDING_MODELS and pt not in BINDING_PARAM_TYPES:
        raise EvidenceTypeError(
            f"{pt.value} is not a binding constant: "
            f"{PARAM_TYPE_NOTES[pt]}. A Hill curve built from a {pt.value} is a "
            "normalised functional response (engagement_model="
            f"{ceiling.value}), not physical receptor occupancy and not a "
            "binding-derived proxy. "
            "Use a Kd or Ki row if you need a binding-based model."
        )
    if want is EngagementModel.binding_occupancy:
        # A genuine binding constant, but not measured on the modelled target.
        raise EvidenceTypeError(
            f"{pt.value}{where} cannot drive engagement_model=binding_occupancy: "
            f"{reason}. Physical fractional occupancy of the modelled receptor "
            f"requires evidence distance {DISTANCE_LABELS[BINDING_OCCUPANCY_DISTANCE]}, "
            f"i.e. relation={BINDING_OCCUPANCY_RELATION.value}; this row supports "
            f"at most engagement_model={ceiling.value}, which is a proxy, not an "
            "occupancy of this receptor."
        )
    raise EvidenceTypeError(
        f"{pt.value}{where} cannot drive engagement_model={want.value}: "
        f"{reason}. At most {ceiling.value} is allowed for this row, because a "
        "transferred number is a proxy and must say so."
    )


def is_modelled(param_type: Any, relation: Any = None, value: Any = None) -> bool:
    """True when a numeric engagement may be reported for this row."""
    if value is None:
        return False
    return model_for(param_type, relation) is not EngagementModel.not_modelled


def describe(row: Mapping[str, Any]) -> dict[str, Any]:
    """UI / paper-ready provenance record for one library row or output row.

    Accepts either a raw library spec (``value_M``/``param_type``/...) or a row
    produced by :func:`flylab.pharm.occupancy.compare_compound`.

    Returns parameter type, value, units, the engagement model actually used,
    the source relation, species, evidence tier, the source string and any
    DOI / PMID found in it.
    """
    param_type = as_param_type(row.get("param_type"))
    relation = as_relation(row.get("relation"))
    distance = RELATION_DISTANCE[relation]
    value = row.get("param_value_M", row.get("value_M", row.get("ec50_M")))
    value = None if value is None else float(value)
    model = model_for(param_type, relation) if value is not None else EngagementModel.not_modelled
    source = str(row.get("source") or "")
    warning = provenance_warning(model, relation, row.get("species"))
    return {
        "receptor": row.get("receptor"),
        "param_type": param_type.value,
        "param_value_M": value,
        "units": "M" if value is not None else None,
        "n": float(row.get("n", 1.0) or 1.0),
        "direction": row.get("direction") or "none",
        "engagement_model": model.value,
        "engagement_model_note": engagement_note(model),
        "provenance_warning": warning,
        "relation": relation.value,
        "relation_note": RELATION_NOTES[relation],
        "evidence_distance": distance.value,
        "evidence_distance_label": DISTANCE_LABELS[distance],
        "evidence_distance_note": DISTANCE_NOTES[distance],
        "param_type_note": PARAM_TYPE_NOTES[param_type],
        "species": row.get("species"),
        "evidence_tier": row.get("evidence_tier", "class_placeholder"),
        "source": source,
        "doi": _find(source, "doi:", "10."),
        "pmid": _find(source, "pmid:", "PMID"),
        "modelled": model is not EngagementModel.not_modelled,
    }


def engagement_note(model: EngagementModel | str) -> str:
    """One sentence explaining what a given engagement model reports."""
    m = model if isinstance(model, EngagementModel) else EngagementModel(str(model))
    if m is EngagementModel.binding_occupancy:
        return (
            "fractional occupancy of the binding site, from a Kd/Ki measured on "
            "this compound at this receptor in this species (evidence distance "
            "E0): a physical occupancy claim about the modelled receptor"
        )
    if m is EngagementModel.binding_engagement_proxy:
        return (
            "engagement proxy derived from a binding constant (Kd/Ki) measured "
            "at evidence distance E1/E2 -- another species, or a related "
            "preparation: the same Hill expression, but it is NOT fractional "
            "receptor occupancy of the modelled target"
        )
    if m is EngagementModel.functional_engagement:
        return (
            "normalised functional engagement from a potency parameter "
            "(EC50/IC50/Kb) measured on this target (evidence distance E0); it "
            "is NOT fractional receptor occupancy"
        )
    if m is EngagementModel.functional_engagement_proxy:
        return (
            "normalised functional engagement from a potency parameter "
            "transferred from another species, a related preparation or a "
            "chemical class (evidence distance E1-E3); it is NOT fractional "
            "receptor occupancy and it is not a measurement on this target"
        )
    return "not modelled: the evidence does not support any number at this receptor"


def provenance_warning(
    model: EngagementModel | str,
    relation: Any = None,
    species: Any = None,
) -> str | None:
    """The warning a row must carry, or ``None`` when it needs none.

    Every proxy model carries one, and it names the gap explicitly -- the
    species the number was measured in (E1), the preparation it was measured on
    (E2), or the fact that it is a class statement (E3).  A reader must never
    see a transferred number without being told what was transferred.
    """
    m = model if isinstance(model, EngagementModel) else EngagementModel(str(model))
    if m not in PROXY_MODELS:
        return None
    dist = as_distance(relation) if relation is not None else None
    where = str(species).strip() if species else ""
    kind = ("binding constant" if m is EngagementModel.binding_engagement_proxy
            else "potency value")
    if dist is EvidenceDistance.E2:
        gap = (
            f"the {kind} was measured on a related, hybrid or native-mixed "
            "preparation" + (f" ({where})" if where else "")
            + ", not on the receptor this key names"
        )
    elif dist is EvidenceDistance.E3:
        gap = (
            f"the {kind} is a chemical-class statement"
            + (f" ({where})" if where else "")
            + ", not a measurement of this compound at this receptor"
        )
    else:  # E1, or an unstated distance: assume the species gap
        gap = (
            f"the {kind} was measured in "
            + (where or "another species")
            + ", not in Drosophila melanogaster"
        )
    label = DISTANCE_LABELS[dist] if dist is not None else "E1-E3"
    if m is EngagementModel.binding_engagement_proxy:
        reported = (
            "a binding-derived engagement proxy, not as physical receptor "
            "occupancy of the modelled target"
        )
    else:
        reported = (
            "a functional engagement proxy, not as a measurement on the "
            "modelled target"
        )
    return (
        f"EVIDENCE DISTANCE {label}: {gap}. The number is therefore reported as "
        f"{reported}; transferring it here is an assumption, not a measurement."
    )


def transformation_table_rows() -> list[dict[str, Any]]:
    """:data:`TRANSFORMATION_TABLE` as printable rows (the methods table).

    One row per (parameter type, evidence distance) pair, with the model the
    rule permits and one sentence about what that model reports.
    """
    rows: list[dict[str, Any]] = []
    for pt in ParameterType:
        for dist in EvidenceDistance:
            model = TRANSFORMATION_TABLE[(pt, dist)]
            rows.append(
                {
                    "param_type": pt.value,
                    "param_type_note": PARAM_TYPE_NOTES[pt],
                    "evidence_distance": dist.value,
                    "evidence_distance_label": DISTANCE_LABELS[dist],
                    "relation": DISTANCE_RELATION[dist].value,
                    "engagement_model": model.value,
                    "engagement_model_note": engagement_note(model),
                    "model_strength": MODEL_STRENGTH[model],
                }
            )
    return rows


def _find(source: str, *markers: str) -> str | None:
    """Pull a DOI / PMID out of a free-text source string, if one is there.

    The marker itself ("doi:", "PMID") is stripped, so the result is the
    identifier alone.
    """
    low = source.lower()
    for marker in markers:
        idx = low.find(marker.lower())
        if idx < 0:
            continue
        token = source[idx:].split()[0].strip().strip(".,;)")
        if marker.lower() in ("doi:", "pmid:", "pmid"):
            token = token[len(marker):].lstrip(":").strip()
        if not token:
            continue
        if marker == "10." and not token.startswith("10."):
            continue
        return token
    return None
