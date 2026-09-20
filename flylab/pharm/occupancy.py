"""Typed receptor engagement over the sourced compound library (schema v3).

Peer review of v0.5 raised two blocking objections, quoted here so the intent
survives:

    "The manuscript defines theta = C^n/(EC50^n + C^n) as fractional receptor
    occupancy. EC50 is a functional potency parameter, not a receptor binding
    constant. IC50 values are also used for some compounds. A Hill response
    derived from EC50 describes normalized functional response, not physical
    receptor occupancy."

    "Missing receptor values are represented by ec50_M = 0.01 even when
    direction: none, so the model reports things like diazepam/RDL = 10^-4.
    Missing evidence should not become a small quantitative response."

Consequences, and the whole point of this module:

* The quantity is called **engagement**, not occupancy, unless the row's
  ``param_type`` is a genuine binding constant (``Kd``/``Ki``), in which case
  ``engagement_model`` is ``binding_occupancy`` and the word is earned.  The
  legal transformations live in :mod:`flylab.pharm.evidence`, which *refuses*
  the illegal ones instead of performing them quietly.
* A row with no sourced value reports ``engagement: None`` -- N/A, not
  modelled -- and is excluded from every numeric aggregation (selectivity
  ratios, curves, Monte-Carlo intervals, gain patches).  ``1e-4`` is never
  produced from missing evidence again.

Backward compatibility: :func:`hill_occupancy`, :func:`load_library`,
:func:`compare_compound`, :func:`occupancy_curve` keep their signatures; rows
keep an ``occupancy`` alias key (same value as ``engagement``, possibly
``None``) and a deprecated ``ec50_M`` alias (which may hold an IC50 or a Kd and
must not be read as an EC50).
"""

from __future__ import annotations

import copy
import hashlib
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

from flylab.pharm.evidence import (
    ALLOWED,
    EngagementModel,
    EvidenceTypeError,
    ParameterType,
    SourceRelation,
    as_param_type,
    as_relation,
    check_transformation,
    describe,
    engagement_note,
    model_for,
)

LIBRARY_PATH = Path(__file__).with_name("library.yaml")

SCHEMA_VERSION = 3

#: One sentence, carried in every :func:`compare_compound` result.
ENGAGEMENT_IS_NOT_OCCUPANCY = (
    "engagement is a normalised Hill response computed from the sourced potency "
    "parameter named in param_type; it is physical receptor occupancy only when "
    "engagement_model is 'binding_occupancy' (a measured Kd/Ki) and is otherwise "
    "a functional engagement, while rows with no sourced value report None "
    "(not modelled) rather than a small number."
)

#: Row used when a compound does not list a receptor at all. No number: the
#: peer review's second objection is that missing evidence must not become a
#: small quantitative response.
PLACEHOLDER_ROW: dict[str, Any] = {
    "param_type": ParameterType.unknown.value,
    "value_M": None,
    "n": 1.0,
    "direction": "none",
    "efficacy": 0.0,
    "relation": SourceRelation.unsupported.value,
    "source": "class-order placeholder: receptor not listed for this compound",
    "evidence_tier": "class_placeholder",
}

ALLOWED_DIRECTIONS = (
    "agonist",
    "partial_agonist",
    "antagonist",
    "positive_modulator",
    "negative_modulator",
    "inhibitor",
    "none",
)
ALLOWED_TIERS = ("literature_order", "class_placeholder", "measured_fit")
ALLOWED_PARAM_TYPES = tuple(p.value for p in ParameterType)
ALLOWED_RELATIONS = tuple(r.value for r in SourceRelation)

#: Fallback pairs if a library predates the ``selectivity_pairs`` section.
DEFAULT_SELECTIVITY_PAIRS: dict[str, dict[str, str]] = {
    "nAChR": {"insect": "insect_nAChR", "vertebrate": "vertebrate_nAChR_a4b2"},
    "GABA_A": {"insect": "insect_RDL", "vertebrate": "vertebrate_GABA_A"},
    "GluCl": {"insect": "insect_GluCl", "vertebrate": "vertebrate_GlyR"},
    "AChE": {"insect": "insect_AChE", "vertebrate": "vertebrate_AChE"},
    "Nav": {"insect": "insect_Nav", "vertebrate": "vertebrate_Nav1_x"},
}

#: Half-decade grid used by :func:`occupancy_curve` (1e-11 .. 1e-3 M).
DEFAULT_CURVE_CONCS: tuple[float, ...] = tuple(
    round(10 ** (-11 + 0.5 * i), 15) for i in range(17)
)


# ---------------------------------------------------------------------------
# the Hill transformation, typed
# ---------------------------------------------------------------------------
def engagement(
    conc_M: float,
    value_M: float | None,
    n: float = 1.0,
    param_type: Any = ParameterType.EC50,
    relation: Any = None,
) -> float | None:
    """Hill engagement ``C^n / (value^n + C^n)`` for a **typed** parameter.

    ``param_type`` says what the sourced number is.  For ``Kd``/``Ki`` the
    result is fractional receptor occupancy; for ``EC50``/``IC50``/``Kb`` it is
    a normalised functional response and must not be called occupancy (the
    peer review: "A Hill response derived from EC50 describes normalized
    functional response, not physical receptor occupancy").

    Returns ``None`` when the row has no value.  Raises
    :class:`~flylab.pharm.evidence.EvidenceTypeError` when the parameter type
    cannot drive a number at all (``class_order``, ``unknown``,
    ``relative_potency``, or ``relation='unsupported'``) -- FlyLab refuses the
    transformation rather than inventing a response.
    """
    pt = as_param_type(param_type)
    model = model_for(pt, relation)
    if model is EngagementModel.not_modelled:
        raise EvidenceTypeError(
            f"cannot compute an engagement from param_type={pt.value}"
            + (f", relation={as_relation(relation).value}" if relation is not None else "")
            + ": this row carries no modellable evidence, so it reports None "
            "(not modelled), never a number."
        )
    if value_M is None:
        return None
    value = float(value_M)
    if conc_M < 0 or value <= 0:
        raise ValueError("concentration must be >= 0 and the potency/affinity value > 0")
    if conc_M == 0:
        return 0.0
    cn = conc_M ** float(n)
    return cn / (value ** float(n) + cn)


def hill_occupancy(conc_M: float, ec50_M: float, n: float = 1.0) -> float:
    """Deprecated alias of :func:`engagement` (kept: callers and tests use it).

    The name is wrong for every row whose ``param_type`` is not ``Kd``/``Ki``,
    which is why the library and :func:`compare_compound` use ``engagement``.
    Assumes the value is an EC50, i.e. returns a *functional engagement*.
    """
    if ec50_M is None:
        raise ValueError(
            "hill_occupancy needs a value; a row with no sourced value is not "
            "modelled (use compare_compound, which reports None for it)"
        )
    value = engagement(conc_M, ec50_M, n, param_type=ParameterType.EC50)
    assert value is not None  # ec50_M is not None on this path
    return value


def binding_occupancy(conc_M: float, kd_M: float, n: float = 1.0) -> float:
    """Fractional occupancy from a genuine binding constant (``Kd``/``Ki``)."""
    if kd_M is None:
        raise ValueError("binding_occupancy needs a Kd/Ki value")
    value = engagement(conc_M, kd_M, n, param_type=ParameterType.Kd)
    assert value is not None
    return value


# ---------------------------------------------------------------------------
# library access
# ---------------------------------------------------------------------------
def _normalise_spec(spec: dict[str, Any]) -> dict[str, Any]:
    """v3 spec, with the deprecated ``ec50_M`` alias mirrored for old callers.

    Accepts a v2 row (``ec50_M`` only): it is read as an untyped EC50 so old
    libraries keep loading, but the type is then explicit in the result.
    """
    out = dict(spec)
    # ``ec50_M`` wins when both are present: callers that perturb a loaded
    # library (genotype shifts, ensembles, sensitivity sweeps) write to the
    # deprecated alias, and their edit must not be overwritten.
    value = out.get("ec50_M") if out.get("ec50_M") is not None else out.get("value_M")
    out["value_M"] = None if value is None else float(value)
    if "param_type" not in out:
        out["param_type"] = (
            ParameterType.EC50.value if out["value_M"] is not None else ParameterType.unknown.value
        )
    if "relation" not in out:
        out["relation"] = (
            SourceRelation.class_extrapolation.value
            if out["value_M"] is not None
            else SourceRelation.unsupported.value
        )
    out["n"] = float(out.get("n", 1.0))
    out["direction"] = out.get("direction", "none")
    out["evidence_tier"] = out.get("evidence_tier", "class_placeholder")
    if out["value_M"] is None:
        out.pop("ec50_M", None)
    else:
        # Deprecated alias: may hold an IC50/Ki/Kd. Do not read it as an EC50.
        out["ec50_M"] = out["value_M"]
    return out


@lru_cache(maxsize=4)
def _read_library(path_str: str) -> dict[str, Any]:
    with Path(path_str).open() as fh:
        lib = yaml.safe_load(fh)
    for entry in (lib.get("compounds") or {}).values():
        rows = entry.get("receptors") or {}
        for receptor, spec in list(rows.items()):
            rows[receptor] = _normalise_spec(spec)
    return lib


def load_library(path: Path | None = None) -> dict[str, Any]:
    """Load (and cache) the YAML library. Returns a private deep copy."""
    return copy.deepcopy(_read_library(str(path or LIBRARY_PATH)))


def library_sha256(path: Path | None = None) -> str:
    """SHA-256 of the library file on disk, for notebook provenance."""
    return hashlib.sha256((path or LIBRARY_PATH).read_bytes()).hexdigest()


def receptor_table(library: dict[str, Any] | None = None) -> dict[str, Any]:
    """The ``receptors`` section (organism / transmitter / family per target)."""
    lib = library or load_library()
    return lib.get("receptors", {})


def selectivity_pairs(library: dict[str, Any] | None = None) -> dict[str, dict[str, str]]:
    """Insect/vertebrate counterpart pairs used by the scorecard."""
    lib = library or load_library()
    return lib.get("selectivity_pairs") or copy.deepcopy(DEFAULT_SELECTIVITY_PAIRS)


def list_compounds(library: dict[str, Any] | None = None) -> list[str]:
    """Sorted library keys."""
    lib = library or load_library()
    return sorted(lib["compounds"])


def _compound(name: str, lib: dict[str, Any]) -> dict[str, Any]:
    key = name.lower().strip()
    if key not in lib["compounds"]:
        raise KeyError(f"unknown compound {name!r}. known: {', '.join(sorted(lib['compounds']))}")
    return lib["compounds"][key]


def receptor_spec(compound: dict[str, Any], receptor: str) -> dict[str, Any]:
    """Declared spec for ``receptor``, or the inert placeholder row (no value)."""
    spec = compound.get("receptors", {}).get(receptor)
    return _normalise_spec(spec) if spec else dict(PLACEHOLDER_ROW)


def spec_value_M(spec: dict[str, Any]) -> float | None:
    """Sourced parameter value in M, or ``None`` for a placeholder row.

    Prefers the deprecated ``ec50_M`` alias when present, because callers that
    perturb a library (genotype shifts, Monte-Carlo ensembles) write to it.
    """
    for key in ("ec50_M", "value_M", "value"):
        value = spec.get(key)
        if value is not None:
            return float(value)
    return None


def is_placeholder(spec: dict[str, Any]) -> bool:
    """True when a row carries no modellable evidence."""
    if spec_value_M(spec) is None:
        return True
    return model_for(spec.get("param_type"), spec.get("relation")) is EngagementModel.not_modelled


# ---------------------------------------------------------------------------
# rows
# ---------------------------------------------------------------------------
def _row(receptor: str, spec: dict[str, Any], conc_M: float) -> dict[str, Any]:
    spec = _normalise_spec(spec)
    param_type = as_param_type(spec.get("param_type"))
    relation = as_relation(spec.get("relation"))
    value = spec_value_M(spec)
    model = model_for(param_type, relation) if value is not None else EngagementModel.not_modelled
    n = float(spec.get("n", 1.0))
    if model is EngagementModel.not_modelled:
        value_engagement: float | None = None
    else:
        check_transformation(param_type, model)
        value_engagement = engagement(conc_M, value, n, param_type=param_type, relation=relation)
    row = {
        "receptor": receptor,
        # v0.6 canonical keys
        "engagement": value_engagement,
        "engagement_model": model.value,
        "engagement_note": engagement_note(model),
        "param_type": param_type.value,
        "param_value_M": value,
        "relation": relation.value,
        "species": spec.get("species"),
        # unchanged keys
        "direction": spec.get("direction", "unknown"),
        "source": spec.get("source", ""),
        "evidence_tier": spec.get("evidence_tier", "class_placeholder"),
        "n": n,
        # deprecated aliases, kept so nothing downstream breaks immediately
        "occupancy": value_engagement,
        "ec50_M": value,
    }
    if model is EngagementModel.not_modelled:
        row["not_modelled_reason"] = (
            "placeholder: no sourced value for this compound at this receptor "
            f"(param_type={param_type.value}, relation={relation.value})"
        )
    if spec.get("efficacy") is not None:
        row["efficacy"] = float(spec["efficacy"])
    return row


def row_provenance(row: dict[str, Any]) -> dict[str, Any]:
    """Provenance record for one row (see :func:`flylab.pharm.evidence.describe`)."""
    return describe(row)


# ---------------------------------------------------------------------------
# selectivity
# ---------------------------------------------------------------------------
def _selectivity(compound: dict[str, Any], conc_M: float, pairs: dict[str, dict[str, str]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for pair_name, pair in pairs.items():
        ins_name, vert_name = pair["insect"], pair["vertebrate"]
        ins_row = _row(ins_name, receptor_spec(compound, ins_name), conc_M)
        vert_row = _row(vert_name, receptor_spec(compound, vert_name), conc_M)
        ins_value, vert_value = ins_row["param_value_M"], vert_row["param_value_M"]
        placeholder = ins_value is None or vert_value is None
        tiers = {ins_row["evidence_tier"], vert_row["evidence_tier"]}
        block: dict[str, Any] = {
            "insect_receptor": ins_name,
            "vertebrate_receptor": vert_name,
            "insect_param_type": ins_row["param_type"],
            "vertebrate_param_type": vert_row["param_type"],
            "insect_param_value_M": ins_value,
            "vertebrate_param_value_M": vert_value,
            # deprecated aliases (may hold an IC50/Ki/Kd, never read as EC50)
            "insect_ec50_M": ins_value,
            "vertebrate_ec50_M": vert_value,
            "insect_engagement": ins_row["engagement"],
            "vertebrate_engagement": vert_row["engagement"],
            "insect_occupancy": ins_row["engagement"],
            "vertebrate_occupancy": vert_row["engagement"],
            "evidence_tier": "class_placeholder" if "class_placeholder" in tiers else sorted(tiers)[0],
            "placeholder": placeholder,
        }
        if placeholder:
            # The peer review: missing evidence must not become a number. A
            # pair with a placeholder on either side has no ratio at all - not
            # a log-ratio of 0.00.
            block.update(
                {
                    "ratio": None,
                    "reason": "placeholder",
                    "ec50_ratio_vert_over_insect": None,
                    "log10_ratio_vert_over_insect": None,
                    "log10_ec50_ratio_vert_over_insect": None,
                    "occupancy_difference": None,
                    "engagement_difference": None,
                    "comparable": False,
                }
            )
        else:
            ratio = vert_value / ins_value
            diff = ins_row["engagement"] - vert_row["engagement"]
            comparable = ins_row["param_type"] == vert_row["param_type"]
            block.update(
                {
                    # >1 means the insect target is the more potent one.
                    "ratio": ratio,
                    "reason": None if comparable else "mixed_param_types",
                    "ec50_ratio_vert_over_insect": ratio,
                    "log10_ratio_vert_over_insect": _log10(ratio),
                    "log10_ec50_ratio_vert_over_insect": _log10(ratio),
                    "occupancy_difference": diff,
                    "engagement_difference": diff,
                    # False when the two sides are different KINDS of parameter
                    # (e.g. an insect IC50 against a vertebrate Kd).
                    "comparable": comparable,
                }
            )
        out[pair_name] = block
    return out


def _log10(x: float) -> float:
    from math import log10

    return log10(x) if x > 0 else float("-inf")


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------
def compare_compound(name: str, conc_M: float, library: dict[str, Any] | None = None) -> dict[str, Any]:
    """Typed engagement of every listed receptor, plus the selectivity scorecard.

    Output keys from v0.4/v0.5 are unchanged (``compound``, ``class``,
    ``concentration_M``, ``receptors``, ``selectivity``, ``disclaimer``, ...).
    Each row gains ``engagement`` (float **or None** for a placeholder),
    ``engagement_model``, ``param_type``, ``param_value_M``, ``relation`` and
    ``species``; ``occupancy`` and ``ec50_M`` remain as deprecated aliases
    holding the same values.
    """
    lib = library or load_library()
    key = name.lower().strip()
    compound = _compound(key, lib)
    rows = [_row(receptor, spec, conc_M) for receptor, spec in compound["receptors"].items()]
    n_not_modelled = sum(1 for r in rows if r["engagement"] is None)
    return {
        "compound": compound["name"],
        "key": key,
        "class": compound.get("class"),
        "cas": compound.get("cas"),
        "concentration_M": conc_M,
        "receptors": rows,
        "selectivity": _selectivity(compound, conc_M, selectivity_pairs(lib)),
        "library_version": lib.get("library_version"),
        "schema_version": lib.get("schema_version"),
        "notes": compound.get("notes"),
        "engagement_is_not_occupancy": ENGAGEMENT_IS_NOT_OCCUPANCY,
        "n_receptors_not_modelled": n_not_modelled,
        "disclaimer": (
            "Library values are literature-order teaching defaults of the kind named in each "
            "row's param_type (EC50/IC50/Kd/Ki), not measured constants for this receptor and "
            "species unless relation says so. Rows with no sourced value report engagement None "
            "(not modelled) and are excluded from every numeric aggregation. Circuit injection "
            "currently uses reduced_taste_v0, not a FlyWire/MaleCNS extract."
        ),
    }


def occupancy_curve(
    compound: str,
    concs: list[float] | None = None,
    library: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Per-receptor engagement across concentrations (default 1e-11..1e-3 M).

    Returns ``[{"conc_M": float, "receptors": {name: engagement | None}}, ...]``.
    Placeholder receptors are ``None`` at every concentration: they are not
    modelled, so they have no curve.
    """
    lib = library or load_library()
    entry = _compound(compound, lib)
    grid = list(concs) if concs else list(DEFAULT_CURVE_CONCS)
    curve = []
    for c in grid:
        values: dict[str, float | None] = {}
        for receptor, spec in entry["receptors"].items():
            row = _row(receptor, spec, float(c))
            values[receptor] = row["engagement"]
        curve.append({"conc_M": float(c), "receptors": values})
    return curve


def engagement_curve(
    compound: str,
    concs: list[float] | None = None,
    library: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Preferred name for :func:`occupancy_curve` (same output)."""
    return occupancy_curve(compound, concs, library)


# ---------------------------------------------------------------------------
# reporting / validation
# ---------------------------------------------------------------------------
def validate_library(library: dict[str, Any] | None = None) -> list[str]:
    """Problems with the library under the v3 schema; empty list means valid."""
    lib = library or load_library()
    problems: list[str] = []
    receptors = lib.get("receptors") or {}
    if int(lib.get("schema_version", 0)) != SCHEMA_VERSION:
        problems.append(f"schema_version is {lib.get('schema_version')!r}, expected {SCHEMA_VERSION}")
    for key, entry in (lib.get("compounds") or {}).items():
        for receptor, spec in (entry.get("receptors") or {}).items():
            where = f"{key}:{receptor}"
            if receptor not in receptors:
                problems.append(f"{where}: receptor not declared in the receptors section")
            try:
                param_type = as_param_type(spec.get("param_type"))
                relation = as_relation(spec.get("relation"))
            except EvidenceTypeError as exc:
                problems.append(f"{where}: {exc}")
                continue
            value = spec_value_M(spec)
            model = model_for(param_type, relation)
            if spec.get("direction") not in ALLOWED_DIRECTIONS:
                problems.append(f"{where}: bad direction {spec.get('direction')!r}")
            if spec.get("evidence_tier") not in ALLOWED_TIERS:
                problems.append(f"{where}: bad evidence_tier {spec.get('evidence_tier')!r}")
            if not str(spec.get("source", "")).strip():
                problems.append(f"{where}: no source string")
            if model is EngagementModel.not_modelled and value is not None:
                problems.append(
                    f"{where}: param_type={param_type.value}/relation={relation.value} is not "
                    f"modellable but carries value_M={value!r}"
                )
            if param_type is ParameterType.unknown and value is not None:
                problems.append(f"{where}: param_type 'unknown' must not carry a value")
            if value is not None and value <= 0:
                problems.append(f"{where}: value_M must be > 0, got {value!r}")
            if value is None and spec.get("evidence_tier") != "class_placeholder":
                problems.append(f"{where}: no value but tier is {spec.get('evidence_tier')!r}")
            if spec.get("evidence_tier") == "class_placeholder" and value is not None:
                problems.append(f"{where}: class_placeholder rows must carry no value")
    return problems


def library_report(library: dict[str, Any] | None = None) -> dict[str, Any]:
    """Census of the library by param_type / evidence tier / relation.

    This is the methods table the paper needs: how many rows rest on a genuine
    binding constant, how many on a functional potency, and how many on nothing
    at all.
    """
    lib = library or load_library()
    by_param: Counter[str] = Counter()
    by_tier: Counter[str] = Counter()
    by_relation: Counter[str] = Counter()
    by_model: Counter[str] = Counter()
    by_receptor: Counter[str] = Counter()
    by_species: Counter[str] = Counter()
    n_rows = 0
    modelled = 0
    for entry in (lib.get("compounds") or {}).values():
        for receptor, spec in (entry.get("receptors") or {}).items():
            n_rows += 1
            param_type = as_param_type(spec.get("param_type"))
            relation = as_relation(spec.get("relation"))
            value = spec_value_M(spec)
            model = model_for(param_type, relation) if value is not None else EngagementModel.not_modelled
            by_param[param_type.value] += 1
            by_tier[str(spec.get("evidence_tier"))] += 1
            by_relation[relation.value] += 1
            by_model[model.value] += 1
            by_receptor[receptor] += 1
            by_species[str(spec.get("species") or "Drosophila melanogaster (or not applicable)")] += 1
            if model is not EngagementModel.not_modelled:
                modelled += 1
    return {
        "schema_version": lib.get("schema_version"),
        "library_version": lib.get("library_version"),
        "n_compounds": len(lib.get("compounds") or {}),
        "n_receptor_keys": len(lib.get("receptors") or {}),
        "n_rows": n_rows,
        "n_rows_modelled": modelled,
        "n_rows_not_modelled": n_rows - modelled,
        "by_param_type": dict(sorted(by_param.items())),
        "by_evidence_tier": dict(sorted(by_tier.items())),
        "by_relation": dict(sorted(by_relation.items())),
        "by_engagement_model": dict(sorted(by_model.items())),
        "by_receptor": dict(sorted(by_receptor.items())),
        "by_species": dict(sorted(by_species.items())),
        "note": ENGAGEMENT_IS_NOT_OCCUPANCY,
    }


__all__ = [
    "ALLOWED",
    "ALLOWED_DIRECTIONS",
    "ALLOWED_PARAM_TYPES",
    "ALLOWED_RELATIONS",
    "ALLOWED_TIERS",
    "ENGAGEMENT_IS_NOT_OCCUPANCY",
    "PLACEHOLDER_ROW",
    "SCHEMA_VERSION",
    "binding_occupancy",
    "compare_compound",
    "engagement",
    "engagement_curve",
    "hill_occupancy",
    "is_placeholder",
    "library_report",
    "library_sha256",
    "list_compounds",
    "load_library",
    "occupancy_curve",
    "receptor_spec",
    "receptor_table",
    "row_provenance",
    "selectivity_pairs",
    "spec_value_M",
    "validate_library",
]
