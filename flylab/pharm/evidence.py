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
  binding site *of the receptor FlyLab is modelling*.  That claim needs two
  things at once, and v0.6.1 enforces both: the parameter must be a
  dissociation constant (``Kd``) or an inhibition constant (``Ki``) from a
  binding assay, **and** the source relation must be
  ``exact_compound_exact_receptor_exact_species``.  A binding constant measured
  on another species' receptor, or on a related/hybrid preparation, does not
  describe occupancy of the modelled *Drosophila* target, however physical the
  measurement was in its own organism.
* ``binding_derived_engagement`` -- the same Hill expression applied to a
  ``Kd``/``Ki`` whose source is one step away (another species, or a related
  receptor).  It is an *engagement proxy derived from a binding constant
  measured elsewhere*, not an occupancy of this receptor, and every row that
  uses it carries a provenance warning naming the gap (see
  :func:`provenance_warning`).  This is the model the second peer review asked
  for: "a Kd measured in an aphid still yielded a model the paper described as
  physical receptor occupancy".
* ``functional_engagement`` -- the same algebra applied to a functional
  potency (``EC50``, ``IC50``, ``Kb``) gives a *normalised functional
  response*, not a physical occupancy.  FlyLab uses that number as the input
  to its gain rules and says so in every output.
* ``not_modelled`` -- no number may be produced at all.  This is the answer to
  the second objection: a ``class_order``/``unknown`` row, or a row whose
  relation is ``unsupported``, yields ``None``, never 1e-4.

The novel part is that FlyLab **refuses** invalid transformations rather than
silently performing them: :func:`check_transformation` raises
:class:`EvidenceTypeError` when asked for binding occupancy from an EC50, from
a cross-species Kd, or for any numeric engagement from a placeholder row.

The two facts therefore decide the model *jointly*; neither alone is enough:

===============  ==========================================  =====================================
``param_type``   ``relation``                                ``model_for`` result
===============  ==========================================  =====================================
Kd / Ki          exact_compound_exact_receptor_exact_species ``binding_occupancy``
Kd / Ki          exact_compound_exact_receptor_other_species ``binding_derived_engagement``
Kd / Ki          exact_compound_related_receptor             ``binding_derived_engagement``
Kd / Ki          class_extrapolation                         ``functional_engagement``
EC50 / IC50 / Kb any relation except ``unsupported``          ``functional_engagement``
any              unsupported                                 ``not_modelled``
relative_potency / class_order / unknown, any relation        ``not_modelled``
===============  ==========================================  =====================================
"""

from __future__ import annotations

from enum import Enum
from typing import Any, Mapping

__all__ = [
    "ParameterType",
    "SourceRelation",
    "EngagementModel",
    "EvidenceTypeError",
    "ALLOWED",
    "BINDING_MODELS",
    "BINDING_OCCUPANCY_RELATION",
    "MODEL_STRENGTH",
    "BINDING_PARAM_TYPES",
    "FUNCTIONAL_PARAM_TYPES",
    "NOT_MODELLED_PARAM_TYPES",
    "PARAM_TYPE_NOTES",
    "RELATION_NOTES",
    "as_param_type",
    "as_relation",
    "model_for",
    "check_transformation",
    "is_modelled",
    "describe",
    "engagement_note",
    "provenance_warning",
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


class EngagementModel(str, Enum):
    """What FlyLab is allowed to compute from a row."""

    binding_occupancy = "binding_occupancy"
    binding_derived_engagement = "binding_derived_engagement"
    functional_engagement = "functional_engagement"
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

#: The one relation that licenses a *physical occupancy* claim: the source must
#: have measured this compound, at this receptor, in this species.  Anything
#: further away gives :attr:`EngagementModel.binding_derived_engagement`.
BINDING_OCCUPANCY_RELATION = SourceRelation.exact_compound_exact_receptor_exact_species

#: The two models that may only ever be derived from a Kd/Ki.
BINDING_MODELS = frozenset(
    {EngagementModel.binding_occupancy, EngagementModel.binding_derived_engagement}
)

#: How strong a claim each model makes.  A relation may only ever move a row
#: *down* this ladder (see :func:`model_for`), and
#: :func:`check_transformation` refuses any request above the row's ceiling.
MODEL_STRENGTH: dict[EngagementModel, int] = {
    EngagementModel.not_modelled: 0,
    EngagementModel.functional_engagement: 1,
    EngagementModel.binding_derived_engagement: 2,
    EngagementModel.binding_occupancy: 3,
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
                      "C^n/(Kd^n+C^n) is physical fractional occupancy",
    ParameterType.Ki: "inhibition constant from a competition binding assay; "
                      "treated as a binding constant",
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


def model_for(param_type: Any, relation: Any = None) -> EngagementModel:
    """The engagement model a row is allowed to drive.

    ``ALLOWED`` sets a *ceiling* from the parameter type alone; the relation can
    only ever weaken the claim, and for a binding constant it must actively
    license the strongest one:

    * ``relation=unsupported`` -> ``not_modelled`` whatever the parameter type
      (there is no evidence to transform).
    * a ``Kd``/``Ki`` reaches ``binding_occupancy`` **only** with
      ``relation=exact_compound_exact_receptor_exact_species``.  This compound,
      this receptor, this species is what makes "fractional occupancy of the
      modelled receptor" a statement about the modelled receptor.
    * a ``Kd``/``Ki`` from another species or a related/hybrid receptor ->
      ``binding_derived_engagement``: the same Hill number, labelled as an
      engagement proxy derived from a binding constant measured elsewhere.
    * ``relation=class_extrapolation`` -> a binding constant is demoted all the
      way to ``functional_engagement``: a class-level statement cannot license
      any binding claim about this particular receptor.

    Called without a relation it returns the ceiling, i.e. ``ALLOWED[pt]``.
    """
    pt = as_param_type(param_type)
    rel = as_relation(relation) if relation is not None else None
    model = ALLOWED[pt]
    if rel is None:
        return model
    if rel is SourceRelation.unsupported:
        return EngagementModel.not_modelled
    if model is EngagementModel.binding_occupancy:
        if rel is BINDING_OCCUPANCY_RELATION:
            return EngagementModel.binding_occupancy
        if rel is SourceRelation.class_extrapolation:
            return EngagementModel.functional_engagement
        # exact receptor in another species, or a related/hybrid preparation
        return EngagementModel.binding_derived_engagement
    return model


def check_transformation(param_type: Any, model: Any, relation: Any = None) -> None:
    """Raise :class:`EvidenceTypeError` if ``model`` may not be derived from
    ``param_type`` (and, when given, from ``relation``).

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
    rel = as_relation(relation) if relation is not None else None
    ceiling = model_for(pt, rel) if rel is not None else ALLOWED[pt]
    if MODEL_STRENGTH[want] <= MODEL_STRENGTH[ceiling]:
        return
    where = f" with relation={rel.value}" if rel is not None else ""
    if ceiling is EngagementModel.not_modelled:
        reason = RELATION_NOTES[rel] if rel is not None else PARAM_TYPE_NOTES[pt]
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
            "binding-derived engagement. "
            "Use a Kd or Ki row if you need a binding-based model."
        )
    if want is EngagementModel.binding_occupancy:
        # A genuine binding constant, but not measured on the modelled target.
        raise EvidenceTypeError(
            f"{pt.value}{where} cannot drive engagement_model=binding_occupancy: "
            f"{RELATION_NOTES[rel] if rel is not None else PARAM_TYPE_NOTES[pt]}. "
            "Physical fractional occupancy of the modelled receptor requires "
            f"relation={BINDING_OCCUPANCY_RELATION.value}; this row supports at "
            f"most engagement_model={ceiling.value}, which is a binding-derived "
            "engagement proxy, not an occupancy of this receptor."
        )
    raise EvidenceTypeError(
        f"{pt.value}{where} cannot drive engagement_model={want.value}: "
        f"{PARAM_TYPE_NOTES[pt]}. At most {ceiling.value} is allowed for this "
        "row."
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
            "this compound at this receptor in this species (a physical "
            "occupancy claim about the modelled receptor)"
        )
    if m is EngagementModel.binding_derived_engagement:
        return (
            "engagement proxy derived from a binding constant (Kd/Ki) that was "
            "NOT measured on the modelled receptor in Drosophila: the same Hill "
            "expression, but it is NOT fractional occupancy of this receptor, "
            "because the constant comes from another species or a related "
            "preparation"
        )
    if m is EngagementModel.functional_engagement:
        return (
            "normalised functional engagement from a potency parameter "
            "(EC50/IC50/Kb); it is NOT fractional receptor occupancy"
        )
    return "not modelled: the evidence does not support any number at this receptor"


def provenance_warning(
    model: EngagementModel | str,
    relation: Any = None,
    species: Any = None,
) -> str | None:
    """The warning a row must carry, or ``None`` when it needs none.

    A ``binding_derived_engagement`` row always carries one, and it names the
    gap explicitly -- the species the constant was measured in, or the fact
    that the preparation was not this receptor.  A reader must never see such a
    number without being told that it is not an occupancy of the *Drosophila*
    target FlyLab is modelling.
    """
    m = model if isinstance(model, EngagementModel) else EngagementModel(str(model))
    if m is not EngagementModel.binding_derived_engagement:
        return None
    rel = as_relation(relation) if relation is not None else None
    where = str(species).strip() if species else ""
    if rel is SourceRelation.exact_compound_related_receptor:
        gap = (
            "the binding constant was measured on a related, hybrid or "
            "native-mixed preparation" + (f" ({where})" if where else "")
            + ", not on the receptor this key names"
        )
    else:
        gap = (
            "the binding constant was measured in "
            + (where or "another species")
            + ", not in Drosophila melanogaster"
        )
    return (
        "SPECIES/PREPARATION GAP: " + gap + ". The number is therefore reported "
        "as a binding-derived engagement proxy, not as physical receptor "
        "occupancy of the modelled target; cross-species transfer of an "
        "affinity is an assumption, not a measurement."
    )


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
