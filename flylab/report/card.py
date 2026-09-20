"""Claim cards: one result, nine sections, nothing left implicit.

``flylab.analysis.claims`` walks the nine-link dependency chain behind a
readout.  A **claim card** is the reader-facing form of that chain: a fixed set
of sections, always in the same order, always all present, so that a claim can
never quietly omit the section that would have undermined it.

    1. claim                              what is being claimed, in one sentence
    2. status                             what kind of statement it is
    3. evidence_inputs                    the typed parameters, with their evidence distance
    4. measured_not_assumed               what was measured rather than asserted
    5. assumptions_introduced             transmitter sign, gain mapping, free concentration
    6. connectome_dependence              which nulls it is distinguishable from, and the p
    7. specification_stability            how much of the gain-rule family retains it
    8. independent_biological_validation  none
    9. named_unknowns                     what nothing in the chain can supply

Sections 6 and 7 depend on analyses a run may not have requested.  They are
then filled in with ``assessed: false`` and the sentence *"not assessed in this
run"* -- they are never dropped, because a missing section reads as an absent
concern rather than an unanswered question.

Public API: :func:`claim_card`, :func:`to_markdown`, :func:`cards_for_run`.
Every import is function-local, so this module costs nothing to import and does
not drag the analysis layer in.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping, Sequence

__all__ = [
    "CARD_VERSION",
    "SECTIONS",
    "SECTION_TITLES",
    "claim_card",
    "to_markdown",
    "cards_for_run",
]

CARD_VERSION = "1.0"

#: the sections of a card, in order.  A card has exactly these keys.
SECTIONS: tuple[str, ...] = (
    "claim",
    "status",
    "evidence_inputs",
    "measured_not_assumed",
    "assumptions_introduced",
    "connectome_dependence",
    "specification_stability",
    "independent_biological_validation",
    "named_unknowns",
)

SECTION_TITLES: dict[str, str] = {
    "claim": "Claim",
    "status": "Status",
    "evidence_inputs": "Evidence inputs",
    "measured_not_assumed": "Measured, not assumed",
    "assumptions_introduced": "Assumptions introduced",
    "connectome_dependence": "Connectome dependence",
    "specification_stability": "Specification stability",
    "independent_biological_validation": "Independent biological validation",
    "named_unknowns": "Named unknowns",
}

NOT_ASSESSED = "not assessed in this run"

#: the three assumptions every FlyLab circuit claim carries, whatever the
#: compound.  They are stated explicitly rather than inferred, because they are
#: the ones a reader is most likely to mistake for measurements.
CORE_ASSUMPTIONS: tuple[dict[str, str], ...] = (
    {
        "id": "transmitter_sign",
        "statement": (
            "Every edge sign comes from the MaleCNS *predicted* consensus transmitter, not "
            "from immunostaining: acetylcholine excites, GABA and glutamate inhibit."
        ),
        "consequence": (
            "A wrong prediction flips a synapse, and the drug patch follows the label rather "
            "than the biology."
        ),
    },
    {
        "id": "gain_mapping",
        "statement": (
            "Receptor engagement becomes a synaptic gain through a rule FlyLab asserts; the "
            "shape of that rule has never been measured for these receptors."
        ),
        "consequence": (
            "The variance budget names this as the single assumption worth an experiment, and "
            "the specification-stability section says how much of this claim survives changing it."
        ),
    },
    {
        "id": "free_concentration",
        "statement": (
            "The concentration is a free concentration at the receptor, applied directly; no "
            "model maps an administered dose to it."
        ),
        "consequence": (
            "Cuticular penetration, haemolymph binding, the blood-brain barrier and metabolism "
            "are all outside the model."
        ),
    },
    {
        "id": "uniform_expression",
        "statement": (
            "The gain is applied uniformly to every cell of a transmitter class, because "
            "per-cell receptor expression is largely unmapped."
        ),
        "consequence": "A cell that does not express the receptor is patched exactly like one that does.",
    },
)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _num(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if out == out else None


def _mapping(value: Any) -> dict[str, Any]:
    return dict(value) if isinstance(value, Mapping) else {}


def _close(a: Any, b: Any, rel: float = 1e-6) -> bool:
    x, y = _num(a), _num(b)
    if x is None or y is None:
        return False
    if x == y:
        return True
    scale = max(abs(x), abs(y))
    return scale > 0 and abs(x - y) / scale <= rel


def _readout_value(readouts: Mapping[str, Any], readout: str) -> float | None:
    value = _num(readouts.get(readout))
    if value is not None:
        return value
    gains = _mapping(readouts.get("gains"))
    return _num(gains.get(readout))


def _pick_readout(readouts: Mapping[str, Any], preferred: str | None) -> str:
    if preferred and _readout_value(readouts, preferred) is not None:
        return preferred
    for key in ("mean_hz", "mn9_hz", "dnp01_hz"):
        if _readout_value(readouts, key) is not None:
            return key
    for key, value in readouts.items():
        if isinstance(key, str) and key.endswith("_hz") and _num(value) is not None:
            return key
    return preferred or "mean_hz"


# --------------------------------------------------------------------------
# sections
# --------------------------------------------------------------------------
def _claim_section(
    subject: Mapping[str, Any],
    readout: str,
    value: float | None,
    vehicle: float | None,
    engine: str,
    graph: str,
    graph_facts: Mapping[str, Any],
) -> dict[str, Any]:
    compound = subject.get("compound") or "vehicle"
    conc = _num(subject.get("concentration_M"))
    change = None
    if value is not None and vehicle is not None and abs(vehicle) > 1e-12:
        change = (value - vehicle) / vehicle
    floor = graph_facts.get("min_weight")
    floor_text = "" if floor is None else f", >= {floor} synapses"
    where = (
        f"the MaleCNS v1.0 `{graph}` cut "
        f"({graph_facts.get('n_nodes')} cells, {graph_facts.get('n_edges')} edges{floor_text})"
    )
    dose = f" at {conc:.2e} M" if conc else ""
    if value is None:
        sentence = (
            f"FlyLab's {engine} engine returns no value for {readout} with {compound}{dose} on {where}."
        )
    elif vehicle is None:
        sentence = (
            f"On {where}, FlyLab's {engine} engine predicts {readout} = {value:.3f} "
            f"for {compound}{dose}."
        )
    else:
        direction = "raises" if value > vehicle else "lowers" if value < vehicle else "leaves unchanged"
        sentence = (
            f"On {where}, FlyLab's {engine} engine predicts that {compound}{dose} {direction} "
            f"{readout} from {vehicle:.3f} under vehicle to {value:.3f}"
            + (f" ({change * 100:+.1f}%)" if change is not None else "")
            + "."
        )
    return {
        "sentence": sentence,
        "compound": subject.get("compound"),
        "concentration_M": conc,
        "engine": engine,
        "graph": graph,
        "readout": readout,
        "value": value,
        "vehicle": vehicle,
        "relative_change": change,
    }


def _status_section() -> dict[str, Any]:
    return {
        "status": "model-derived",
        "label": "model_derived",
        "statement": (
            "This is a prediction of a simulation, not an observation. It is model-derived: "
            "a typed literature parameter, an asserted gain rule and a measured connectome cut "
            "were combined by a deterministic engine. No part of it was measured in a living fly, "
            "and the notebook's live_lab field is null."
        ),
        "assessed": True,
    }


def _evidence_section(audit: Mapping[str, Any]) -> dict[str, Any]:
    by_step = {link.get("step"): link for link in audit.get("chain") or []}
    sources = list((by_step.get("parameter_source") or {}).get("sources") or [])
    per_receptor = {
        str(row.get("receptor")): row
        for row in ((by_step.get("engagement_transformation") or {}).get("detail") or {}).get(
            "per_receptor"
        )
        or []
    }
    rows: list[dict[str, Any]] = []
    for row in sources:
        receptor = str(row.get("receptor"))
        extra = _mapping(per_receptor.get(receptor))
        rows.append(
            {
                "receptor": receptor,
                "parameter_type": row.get("param_type"),
                "parameter_value_M": row.get("param_value_M"),
                "evidence_distance": row.get("relation"),
                "evidence_distance_note": row.get("relation_note"),
                "evidence_tier": row.get("evidence_tier"),
                "species": row.get("species"),
                "source": row.get("source"),
                "doi": row.get("doi"),
                "pmid": row.get("pmid"),
                "engagement": extra.get("engagement"),
                "engagement_model": extra.get("engagement_model"),
                "modelled": row.get("classification") != "NOT MODELLED",
            }
        )
    modelled = [r for r in rows if r["modelled"]]
    by_type: dict[str, int] = {}
    for r in modelled:
        key = str(r["parameter_type"] or "unknown")
        by_type[key] = by_type.get(key, 0) + 1
    return {
        "assessed": True,
        "n_rows": len(rows),
        "n_modelled": len(modelled),
        "n_not_modelled": len(rows) - len(modelled),
        "by_parameter_type": dict(sorted(by_type.items())),
        "rows": rows,
        "statement": (
            f"{len(modelled)} of {len(rows)} receptor rows carry a sourced, typed parameter. "
            "The parameter type decides the transformation: only a measured Kd or Ki yields a "
            "binding occupancy; an EC50/IC50 yields a functional engagement, which is a "
            "different quantity. A row with no sourced value returns N/A and is excluded, "
            "never replaced by a small number."
        ),
    }


def _measured_section(audit: Mapping[str, Any], graph_facts: Mapping[str, Any]) -> dict[str, Any]:
    by_step = {link.get("step"): link for link in audit.get("chain") or []}
    link = by_step.get("malecns_edges") or {}
    detail = _mapping(link.get("detail")) or dict(graph_facts)
    return {
        "assessed": True,
        "what": "the MaleCNS v1.0 edges",
        "n_nodes": detail.get("n_nodes"),
        "n_edges": detail.get("n_edges"),
        "min_weight": detail.get("min_weight"),
        "hops": detail.get("hops"),
        "map_id": detail.get("map"),
        "citation": detail.get("citation"),
        "statement": (
            "Exactly one link of the nine-link chain behind this number is a measurement of the "
            "system being simulated: the synaptic wiring. It is a hops-limited neighbourhood of "
            "the public MaleCNS v1.0 reconstruction with a synapse-count floor -- "
            f"{detail.get('n_nodes')} cells and {detail.get('n_edges')} edges. Everything "
            "downstream of it is assertion or computation."
        ),
    }


def _assumptions_section(audit: Mapping[str, Any]) -> dict[str, Any]:
    items = [dict(item) for item in CORE_ASSUMPTIONS]
    known = {item["statement"] for item in items}
    counts: dict[str, int] = {}
    for link in audit.get("chain") or []:
        step = str(link.get("step"))
        for text in link.get("assumptions") or []:
            if text in known:
                continue
            known.add(text)
            counts[step] = counts.get(step, 0) + 1
            items.append(
                {
                    "id": f"{step}_{counts[step]}",
                    "statement": text,
                    "consequence": None,
                    "from_step": step,
                }
            )
    return {
        "assessed": True,
        "n_assumptions": len(items),
        "items": items,
        "statement": (
            "These are asserted by FlyLab, not measured. They are listed here because each of "
            "them can change the sign or the size of the claim above."
        ),
    }


def _dependence_cell(payload: Any, compound: str | None, conc_M: float | None) -> dict[str, Any] | None:
    """The dependence profile matching this card, from a run payload."""
    cells: list[Mapping[str, Any]] = []
    if isinstance(payload, Mapping):
        if payload.get("cells"):
            cells = [c for c in payload["cells"] if isinstance(c, Mapping)]
        elif payload.get("modes"):
            cells = [payload]
    elif isinstance(payload, Sequence) and not isinstance(payload, (str, bytes)):
        cells = [c for c in payload if isinstance(c, Mapping)]
    if not cells:
        return None
    exact = [
        c
        for c in cells
        if (compound is None or str(c.get("compound")) == str(compound))
        and _close(c.get("conc_M"), conc_M)
    ]
    if exact:
        return dict(exact[0])
    same_compound = [c for c in cells if str(c.get("compound")) == str(compound)]
    if len(same_compound) == 1:
        out = dict(same_compound[0])
        out["_concentration_mismatch"] = True
        return out
    return None


def _dependence_section(payload: Any, compound: str | None, conc_M: float | None, requested: bool) -> dict[str, Any]:
    if not requested and payload is None:
        return {
            "assessed": False,
            "statement": (
                f"Connectome dependence was {NOT_ASSESSED}. Nothing here says whether this "
                "prediction needs the measured wiring or would survive a degraded graph."
            ),
            "how_to_assess": "add `dependence: {permutations: 1000}` to the spec's `analyses:` block",
            "distinguishable_from": [],
            "not_distinguishable_from": [],
        }
    cell = _dependence_cell(payload, compound, conc_M)
    if cell is None:
        return {
            "assessed": False,
            "statement": (
                f"Connectome dependence was requested but {NOT_ASSESSED} for this compound and "
                "concentration; the analysis covered other cells."
            ),
            "how_to_assess": "set `dependence: {concentrations: all}` to cover every dose in the spec",
            "distinguishable_from": [],
            "not_distinguishable_from": [],
        }
    modes = [dict(m) for m in cell.get("modes") or [] if isinstance(m, Mapping)]
    def _p(mode: Mapping[str, Any]) -> float | None:
        return _num(mode.get("p_two_sided", mode.get("p")))

    beat = [m for m in modes if m.get("beats_null") is True]
    not_beat = [m for m in modes if m.get("beats_null") is not True]
    slim = lambda m: {  # noqa: E731 - a local projection, not a policy
        "null_model": m.get("mode"),
        "information_kept": m.get("information_kept"),
        "p_two_sided": _p(m),
        "p_resolution": _num(m.get("p_resolution")),
        "z": _num(m.get("z")),
        "permutations": m.get("n"),
    }
    if beat:
        names = ", ".join(f"{m.get('mode')} (p = {_fmt_p(_p(m))})" for m in beat)
        statement = (
            f"The effect is distinguishable from {len(beat)} of {len(modes)} degraded graph models: "
            f"{names}. It is not distinguishable from "
            + (", ".join(str(m.get("mode")) for m in not_beat) or "none of the others")
            + "."
        )
    else:
        statement = (
            f"The effect is not distinguishable from any of the {len(modes)} degraded graph models "
            "tested: a graph that keeps only composition-level information reproduces it. "
            "This prediction does not need the measured wiring."
        )
    return {
        "assessed": True,
        "statement": statement,
        "class": cell.get("class"),
        "necessary_information_level": _level_name(cell.get("necessary_information_level")),
        "necessary_information_level_detail": cell.get("necessary_information_level"),
        "permutations": cell.get("n"),
        "real_effect": _num(cell.get("real_effect")),
        "readout": cell.get("readout"),
        "concentration_M": _num(cell.get("conc_M")),
        "concentration_mismatch": bool(cell.get("_concentration_mismatch")),
        "distinguishable_from": [slim(m) for m in beat],
        "not_distinguishable_from": [slim(m) for m in not_beat],
        "note": (
            "p is the empirical two-sided permutation probability with resolution 1/(n+1); "
            "z is secondary and a large |z| with a coarse p means only that n was small."
        ),
    }


def _level_name(level: Any) -> str | None:
    """The weakest graph model that still reproduces the effect, as a name.

    The analysis layer returns either a bare string or a block describing the
    level; the card shows the name and keeps the block beside it.
    """
    if isinstance(level, Mapping):
        return level.get("level") or level.get("mode")
    return str(level) if level else None


def _fmt_p(p: float | None) -> str:
    if p is None:
        return "n/a"
    return f"{p:.4f}" if p >= 1e-4 else f"{p:.1e}"


def _stability_section(payload: Any, compound: str | None, requested: bool) -> dict[str, Any]:
    result = None
    if isinstance(payload, Mapping):
        result = payload.get("result") if isinstance(payload.get("result"), Mapping) else payload
    if not isinstance(result, Mapping) or not result.get("rows"):
        return {
            "assessed": False,
            "statement": (
                f"Specification stability was {NOT_ASSESSED}. Nothing here says whether this "
                "claim survives the alternative engagement-to-gain rules."
            )
            if not requested
            else (
                f"Specification stability was requested but {NOT_ASSESSED}: the analysis "
                "produced no conclusion rows."
            ),
            "how_to_assess": "add `robustness: {specifications: default_family}` to the spec's `analyses:` block",
            "conclusions": [],
        }
    family_size = result.get("family_size")
    conclusions = []
    for row in result.get("rows") or []:
        if not isinstance(row, Mapping):
            continue
        claim_text = str(row.get("claim") or "")
        conclusions.append(
            {
                "conclusion": row.get("conclusion"),
                "claim": claim_text,
                "retained": row.get("n_retained"),
                "of": row.get("n_specs") or family_size,
                "fraction_retained": _num(row.get("fraction_retained")),
                "fragile": row.get("fragile"),
                "failing_specifications": list(row.get("failing_specs") or []),
                "mentions_this_compound": bool(compound and compound.lower() in claim_text.lower()),
            }
        )
    mine = [c for c in conclusions if c["mentions_this_compound"]] or conclusions
    parts = [
        f"{c['conclusion']}: retained in {c['retained']}/{c['of']} specifications"
        for c in mine[:4]
        if c.get("retained") is not None
    ]
    return {
        "assessed": True,
        "family_size": family_size,
        "default_specification": result.get("default_spec"),
        "statement": (
            "Across the prespecified family of engagement-to-gain specifications, "
            + ("; ".join(parts) if parts else "no conclusion could be scored")
            + ". A conclusion that is not retained by the whole family is a property of the "
            "gain rule, not of the pharmacology-connectome integration."
        ),
        "conclusions": conclusions,
        "conclusions_mentioning_compound": [c["conclusion"] for c in conclusions if c["mentions_this_compound"]],
    }


def _validation_section() -> dict[str, Any]:
    return {
        "assessed": True,
        "status": "none",
        "statement": (
            "None. There is no independent, out-of-sample biological validation of this claim. "
            "Rank comparisons against published orderings are literature *concordance*, and some "
            "of them share a source with the compound library, so they are not out-of-sample. "
            "Agreement between the rate and LIF engines is implementation consistency, not "
            "cross-model validation: the two agree on direction only."
        ),
        "what_would_change_it": (
            "An electrophysiological or behavioural measurement in adult Drosophila, obtained "
            "after this prediction was recorded, at a dose and readout the model actually names."
        ),
    }


def _unknowns_section(audit: Mapping[str, Any]) -> dict[str, Any]:
    fiu = _mapping(audit.get("fact_inference_unknown"))
    items = [
        {
            "id": item.get("id"),
            "statement": item.get("statement"),
            "why": item.get("why"),
            "would_resolve": item.get("would_resolve"),
        }
        for item in fiu.get("unknown") or []
        if isinstance(item, Mapping)
    ]
    return {
        "assessed": True,
        "n_unknowns": len(items),
        "items": items,
        "statement": (
            "These are gaps no part of the chain can fill. They are named rather than "
            "parameterised: FlyLab does not supply a plausible number in their place."
        ),
    }


# --------------------------------------------------------------------------
# the card
# --------------------------------------------------------------------------
def claim_card(result_or_notebook: Mapping[str, Any] | None = None, **context: Any) -> dict[str, Any]:
    """One claim card: nine sections, as a plain dict.

    ``result_or_notebook`` is a FlyLab notebook (schema 0.3) or any mapping
    carrying ``compound`` / ``concentration_M``.  Recognised context:

    ``compound``, ``conc_M``, ``engine``, ``graph``, ``readout``
        identify the claim; taken from the notebook when omitted.
    ``vehicle``
        the matching vehicle notebook or readouts block, so the claim sentence
        can quote a change rather than a bare value.
    ``dependence`` / ``robustness``
        the payloads written by :func:`flylab.spec.run_spec`.  When absent, the
        corresponding section says *"not assessed in this run"* rather than
        disappearing.
    ``requested_analyses``
        which analyses the run asked for, so a requested-but-failed analysis
        reads differently from one that was never requested.
    """
    from flylab.analysis.claims import claim_audit

    payload = dict(result_or_notebook or {})
    compound = context.get("compound") or payload.get("compound")
    conc = _num(context.get("conc_M", payload.get("concentration_M")))
    readouts = payload.get("readouts") if isinstance(payload.get("readouts"), Mapping) else {}
    engine = str(
        context.get("engine")
        or readouts.get("engine")
        or ("lif" if "spik" in str(payload.get("assay") or "") else "rate")
    )
    graph = str(context.get("graph") or "named")

    audit = claim_audit(payload, compound=compound, conc_M=conc, graph=graph)
    graph_facts = _mapping(
        next(
            (
                link.get("detail")
                for link in audit.get("chain") or []
                if link.get("step") == "malecns_edges"
            ),
            {},
        )
    )

    readout = _pick_readout(readouts, context.get("readout"))
    value = _readout_value(readouts, readout)
    vehicle_src = context.get("vehicle")
    vehicle_readouts = {}
    if isinstance(vehicle_src, Mapping):
        vehicle_readouts = (
            vehicle_src["readouts"]
            if isinstance(vehicle_src.get("readouts"), Mapping)
            else vehicle_src
        )
    vehicle = _readout_value(vehicle_readouts, readout) if vehicle_readouts else None

    requested = {str(a) for a in context.get("requested_analyses") or []}
    card: dict[str, Any] = {
        "flylab_card_version": CARD_VERSION,
        "kind": "flylab-claim-card",
        "run": context.get("run_name"),
        "notebook": context.get("notebook_path"),
        "subject": {
            "compound": compound,
            "concentration_M": conc,
            "engine": engine,
            "graph": graph,
            "readout": readout,
            "assay": payload.get("assay"),
        },
        "sections": list(SECTIONS),
        "claim": _claim_section(
            audit.get("subject") or {"compound": compound, "concentration_M": conc},
            readout,
            value,
            vehicle,
            engine,
            graph,
            graph_facts,
        ),
        "status": _status_section(),
        "evidence_inputs": _evidence_section(audit),
        "measured_not_assumed": _measured_section(audit, graph_facts),
        "assumptions_introduced": _assumptions_section(audit),
        "connectome_dependence": _dependence_section(
            context.get("dependence"), compound, conc, "dependence" in requested
        ),
        "specification_stability": _stability_section(
            context.get("robustness"), compound, "robustness" in requested
        ),
        "independent_biological_validation": _validation_section(),
        "named_unknowns": _unknowns_section(audit),
        "provenance": dict(payload.get("provenance") or {}),
        "notebook_warnings": [w for w in (payload.get("warnings") or []) if isinstance(w, str)],
        "disclaimer": (
            "FlyLab is a simulation on a public connectome. Nothing on this card is a measured "
            "drug effect in a living fly."
        ),
    }
    return card


def to_markdown(card: Mapping[str, Any]) -> str:
    """The card as readable Markdown, sections in :data:`SECTIONS` order."""
    subject = _mapping(card.get("subject"))
    conc = _num(subject.get("concentration_M"))
    title = f"{subject.get('compound') or 'vehicle'}" + (f" at {conc:.2e} M" if conc else "")
    lines = [
        f"# Claim card -- {title}",
        "",
        f"*{subject.get('engine')} engine on the `{subject.get('graph')}` MaleCNS cut; "
        f"readout {subject.get('readout')}.*",
        "",
    ]
    for key in SECTIONS:
        section = _mapping(card.get(key))
        lines.append(f"## {SECTION_TITLES.get(key, key)}")
        lines.append("")
        if key == "claim":
            lines += [f"> {section.get('sentence')}", ""]
            continue
        if not section.get("assessed", True):
            lines += [section.get("statement") or NOT_ASSESSED, ""]
            if section.get("how_to_assess"):
                lines += [f"*How to assess it: {section['how_to_assess']}.*", ""]
            continue
        if section.get("statement"):
            lines += [str(section["statement"]), ""]
        lines += _section_body(key, section)
    lines += [
        "---",
        "",
        str(card.get("disclaimer") or ""),
        "",
    ]
    return "\n".join(lines).rstrip() + "\n"


def _short(text: Any, limit: int = 110) -> str:
    """A one-cell rendering of a long source string; the JSON keeps it whole."""
    out = str(text or "-").replace("|", "/").replace("\n", " ")
    return out if len(out) <= limit else out[: limit - 1].rstrip() + "\u2026"


def _section_body(key: str, section: Mapping[str, Any]) -> list[str]:
    out: list[str] = []
    if key == "evidence_inputs":
        rows = [r for r in section.get("rows") or [] if isinstance(r, Mapping)]
        if rows:
            out += [
                "| receptor | parameter type | value (M) | evidence distance | engagement | source |",
                "|---|---|---|---|---|---|",
            ]
            for r in rows:
                value = r.get("parameter_value_M")
                engagement = r.get("engagement")
                out.append(
                    "| {receptor} | {ptype} | {value} | {dist} | {eng} | {src} |".format(
                        receptor=r.get("receptor"),
                        ptype=r.get("parameter_type") or "-",
                        value="not modelled" if value is None else f"{float(value):.2e}",
                        dist=r.get("evidence_distance") or "-",
                        eng="N/A" if engagement is None else f"{float(engagement):.3f}",
                        src=_short(r.get("source")),
                    )
                )
            out.append("")
    elif key == "assumptions_introduced":
        for item in section.get("items") or []:
            if not isinstance(item, Mapping):
                continue
            line = f"- **{item.get('id')}** -- {item.get('statement')}"
            if item.get("consequence"):
                line += f" *{item['consequence']}*"
            out.append(line)
        out.append("")
    elif key == "connectome_dependence":
        for label, field in (
            ("Distinguishable from", "distinguishable_from"),
            ("Not distinguishable from", "not_distinguishable_from"),
        ):
            items = [m for m in section.get(field) or [] if isinstance(m, Mapping)]
            if not items:
                continue
            out.append(f"- {label}:")
            for m in items:
                out.append(
                    f"    - `{m.get('null_model')}` -- p = {_fmt_p(_num(m.get('p_two_sided')))}"
                    + (f", n = {m.get('permutations')}" if m.get("permutations") else "")
                    + (f" ({m.get('information_kept')})" if m.get("information_kept") else "")
                )
        if section.get("class"):
            out.append(f"- Dependence class: **{section['class']}**")
        if section.get("necessary_information_level"):
            out.append(f"- Necessary information level: `{section['necessary_information_level']}`")
        out.append("")
    elif key == "specification_stability":
        rows = [c for c in section.get("conclusions") or [] if isinstance(c, Mapping)]
        if rows:
            out += ["| conclusion | retained | fragile |", "|---|---|---|"]
            for c in rows:
                out.append(
                    f"| {c.get('conclusion')} | {c.get('retained')}/{c.get('of')} | "
                    f"{'yes' if c.get('fragile') else 'no'} |"
                )
            out.append("")
    elif key == "independent_biological_validation":
        if section.get("what_would_change_it"):
            out += [f"*What would change this: {section['what_would_change_it']}*", ""]
    elif key == "named_unknowns":
        for item in section.get("items") or []:
            if not isinstance(item, Mapping):
                continue
            out.append(f"- **{item.get('id')}** -- {item.get('statement')}")
            if item.get("would_resolve"):
                out.append(f"    - *Would resolve it:* {item['would_resolve']}")
        out.append("")
    elif key == "measured_not_assumed":
        if section.get("citation"):
            out += [f"*Source: {section['citation']}*", ""]
    return out


def cards_for_run(run_dir: str | Path) -> list[dict[str, Any]]:
    """Every claim card of a run directory.

    Reads ``cards/*.json`` when :func:`flylab.spec.run_spec` already wrote them;
    otherwise rebuilds them from the run's notebooks and analysis artifacts, so
    a run made before cards existed still yields cards.
    """
    root = Path(run_dir)
    if not root.exists():
        raise FileNotFoundError(f"run directory not found: {root}")
    card_dir = root / "cards"
    if card_dir.is_dir():
        found = sorted(card_dir.glob("*.json"))
        if found:
            return [json.loads(p.read_text()) for p in found]

    manifest = {}
    manifest_path = root / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    spec = _mapping(manifest.get("spec"))
    requested = list(manifest.get("analyses_requested") or spec.get("analyses") or [])
    dependence = _load_json(root / "dependence.json")
    robustness = _load_json(root / "robustness.json")

    nb_dir = root / "notebooks"
    if not nb_dir.is_dir():
        return []
    notebooks = [json.loads(p.read_text()) for p in sorted(nb_dir.glob("*.json"))]
    vehicles = {
        str(nb.get("assay")): nb for nb in notebooks if nb.get("compound") in (None, "", "vehicle")
    }
    concs = [
        _num(nb.get("concentration_M"))
        for nb in notebooks
        if nb.get("compound") and _num(nb.get("concentration_M")) is not None
    ]
    headline = max(concs) if concs else None

    cards: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for nb in notebooks:
        compound = nb.get("compound")
        conc = _num(nb.get("concentration_M"))
        if not compound or conc is None or (headline is not None and not _close(conc, headline)):
            continue
        assay = str(nb.get("assay") or "")
        engine = "lif" if "spik" in assay else "rate"
        key = (str(compound), engine)
        if key in seen:
            continue
        seen.add(key)
        cards.append(
            claim_card(
                nb,
                compound=compound,
                conc_M=conc,
                engine=engine,
                graph=str(spec.get("graph") or "named"),
                readout=(spec.get("readouts") or [None])[0],
                vehicle=vehicles.get(assay),
                dependence=dependence,
                robustness=robustness,
                requested_analyses=requested,
                run_name=manifest.get("name"),
            )
        )
    return cards


def _load_json(path: Path) -> Any:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:  # pragma: no cover - a truncated artifact
        return None
