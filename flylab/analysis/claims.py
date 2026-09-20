"""Automated claim provenance: what every FlyLab number actually rests on.

A FlyLab readout is the end of a chain.  A cited potency parameter becomes a
typed *engagement*; engagement becomes a *gain* through a mechanism rule the
project asserts rather than measures; the gain is applied uniformly because
per-cell receptor expression is not known; the patched network is the real
MaleCNS wiring but with *predicted* transmitters; and a deterministic engine
turns that into a firing rate.  Only the first and the fifth of those links are
measurements, and they are measurements of different things.

This module writes that chain down, per result, in a machine-readable form, so
a reader never has to reconstruct it from prose:

* :func:`claim_audit` -- the full dependency chain, each link labelled
  ``OBSERVED`` / ``LITERATURE-DERIVED`` / ``MODEL-ASSUMPTION`` / ``COMPUTED`` /
  ``UNKNOWN`` and carrying the sources, assumptions and unknowns it introduces.
* :func:`fact_inference_unknown` -- the short "Fact / Inference / Unknown"
  summary that should accompany every experiment.

Nothing here computes pharmacology or circuit maths: it reads a notebook (or a
compound + concentration) and reports provenance.  Imports are function-local
so ``import flylab.analysis.claims`` stays cheap in the browser build, and no
web framework is touched.
"""

from __future__ import annotations

from typing import Any, Mapping

__all__ = [
    "LABELS",
    "CLASSIFICATIONS",
    "CHAIN_STEPS",
    "CLAIMS_VERSION",
    "claim_audit",
    "fact_inference_unknown",
    "label_counts",
    "to_markdown",
]

CLAIMS_VERSION = "1.0"

#: provenance labels, strongest evidence first.  ``OBSERVED`` is a direct
#: measurement of *this* system (the MaleCNS reconstruction);
#: ``LITERATURE-DERIVED`` is a measurement of something else that FlyLab reuses.
LABELS: tuple[str, ...] = (
    "OBSERVED",
    "LITERATURE-DERIVED",
    "MODEL-ASSUMPTION",
    "COMPUTED",
    "UNKNOWN",
)

#: the UI-facing chip vocabulary.  It is deliberately *not* the same list: a
#: reader needs to know whether a number was measured (LITERATURE), derived
#: from one by a transformation FlyLab chose (MODEL-DERIVED), asserted
#: (MODEL-ASSUMPTION), simulated (PREDICTION) or simply absent (NOT MODELLED).
CLASSIFICATIONS: tuple[str, ...] = (
    "LITERATURE",
    "MODEL-DERIVED",
    "MODEL-ASSUMPTION",
    "PREDICTION",
    "NOT MODELLED",
)

#: the dependency chain, in the order a value travels it
CHAIN_STEPS: tuple[str, ...] = (
    "result",
    "parameter_source",
    "engagement_transformation",
    "mechanism_rule",
    "expression_assumption",
    "malecns_edges",
    "transmitter_predictions",
    "engine",
    "readout",
)

#: unknowns that apply to every FlyLab run, whatever the compound
UNIVERSAL_UNKNOWNS: tuple[dict[str, str], ...] = (
    {
        "id": "cns_free_concentration",
        "statement": (
            "The free concentration of the compound in the fly CNS after a dose is unknown."
        ),
        "why": (
            "FlyLab is dosed with a free concentration at the receptor. Nothing maps an "
            "applied dose to that number: cuticular penetration, haemolymph binding, the "
            "blood-brain barrier and metabolism are all unmeasured here, and the "
            "one-compartment exposure model uses placeholder constants."
        ),
        "would_resolve": "Measured CNS/haemolymph concentration-time curves for this compound in adult Drosophila.",
    },
    {
        "id": "mn9_receptor_expression",
        "statement": "The receptor complement of MN9 -- the main motor readout -- is not measured.",
        "why": (
            "No cited atlas reports Rdl or nAChR expression for identified adult leg or "
            "labellar motor neurons, so MN9 carries the 'unknown' default weight and the "
            "model patches it with the same uniform gain as every other cell."
        ),
        "would_resolve": (
            "Adult VNC snRNA-seq restricted to annotated motor-neuron types, reporting "
            "nAChR / Rdl / GluCl / Ace expression."
        ),
    },
    {
        "id": "biological_variance",
        "statement": "Biological variance -- between flies, between cells, across time -- is not represented.",
        "why": (
            "Every run is deterministic given its seed. The ensemble and Monte-Carlo layers "
            "perturb FlyLab's own parameters, which is model sensitivity, not the spread you "
            "would see across animals."
        ),
        "would_resolve": "Replicated recordings from multiple animals with the same stimulus and dose.",
    },
)


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def _num(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if out == out else None  # drop NaN


def _link(
    step: str,
    label: str,
    classification: str,
    statement: str,
    *,
    detail: Any = None,
    sources: list[dict[str, Any]] | None = None,
    assumptions: list[str] | None = None,
    unknowns: list[str] | None = None,
) -> dict[str, Any]:
    if label not in LABELS:  # pragma: no cover - guards a typo in this module
        raise ValueError(f"unknown provenance label {label!r}")
    if classification not in CLASSIFICATIONS:  # pragma: no cover
        raise ValueError(f"unknown classification {classification!r}")
    return {
        "step": step,
        "order": CHAIN_STEPS.index(step) if step in CHAIN_STEPS else len(CHAIN_STEPS),
        "label": label,
        "classification": classification,
        "statement": statement,
        "detail": detail if detail is not None else {},
        "sources": sources or [],
        "assumptions": assumptions or [],
        "unknowns": unknowns or [],
    }


#: readout keys that are inputs or bookkeeping, not predictions
_NON_READOUT = frozenset(
    {"drive_hz", "steps", "t_ms", "seed", "n_nodes", "n_edges", "runtime_ms", "runtime_s"}
)


def _rate_values(readouts: Mapping[str, Any] | None) -> dict[str, float]:
    """The predicted quantities in a readouts block, without the design fields."""
    out: dict[str, float] = {}
    for key, value in (readouts or {}).items():
        if key in _NON_READOUT or not isinstance(key, str):
            continue
        number = _num(value)
        if number is None:
            continue
        if key.endswith("_hz") or key.startswith("g_") or key.endswith("_index") or key.endswith("_ratio"):
            out[key] = number
    return out


def _subject(payload: Mapping[str, Any] | None) -> dict[str, Any]:
    """Compound / concentration / assay, from a notebook or a bare result."""
    p = dict(payload or {})
    readouts = p.get("readouts") if isinstance(p.get("readouts"), Mapping) else {}
    conc = p.get("concentration_M", p.get("conc_M"))
    graph = p.get("graph")
    if graph is None and isinstance(p.get("map"), Mapping):
        graph = p["map"].get("version")
    return {
        "compound": p.get("compound") or p.get("key"),
        "concentration_M": _num(conc),
        "assay": p.get("assay") or "occupancy",
        "graph": graph,
        "engine": p.get("engine") or (readouts.get("engine") if readouts else None),
        "readouts": sorted(_rate_values(readouts)),
    }


def _graph_name(subject: Mapping[str, Any], fallback: str = "named") -> str:
    name = subject.get("graph")
    if isinstance(name, str):
        if "taste_motor" in name:
            return "taste_motor"
        if name in ("named", "taste_motor"):
            return name
    return fallback


def _occupancy_rows(payload: Mapping[str, Any], subject: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = payload.get("occupancy")
    if isinstance(rows, list) and rows:
        return [dict(r) for r in rows if isinstance(r, Mapping)]
    compound = subject.get("compound")
    if not compound:
        return []
    from flylab.pharm.occupancy import compare_compound

    conc = subject.get("concentration_M")
    try:
        return list(compare_compound(str(compound), float(conc or 0.0))["receptors"])
    except (KeyError, ValueError):
        return []


def _described(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    from flylab.pharm.evidence import describe

    out = []
    for row in rows:
        try:
            rec = describe(row)
        except Exception:  # pragma: no cover - a malformed row must not break an audit
            rec = {"receptor": row.get("receptor"), "modelled": False}
        rec = dict(rec)
        eng = row.get("engagement", row.get("occupancy"))
        rec["engagement"] = _num(eng)
        rec["classification"] = (
            "NOT MODELLED"
            if not rec.get("modelled")
            else "LITERATURE"
            if rec.get("evidence_tier") == "literature_order"
            else "MODEL-ASSUMPTION"
        )
        out.append(rec)
    return out


# --------------------------------------------------------------------------
# claim audit
# --------------------------------------------------------------------------
def claim_audit(notebook_or_result: Mapping[str, Any] | None = None, **kw: Any) -> dict[str, Any]:
    """The full dependency chain behind one result.

    ``notebook_or_result`` is a FlyLab notebook (schema 0.3) or any dict that
    carries ``compound`` / ``concentration_M``; ``compound=`` and ``conc_M=``
    may be passed instead.  Every link in the returned ``chain`` says what it
    contributes, how strong that contribution is (``label``), how the UI should
    chip it (``classification``), which sources back it, which assumptions it
    introduces and which unknowns it leaves open.
    """
    payload: dict[str, Any] = dict(notebook_or_result or {})
    if kw.get("compound"):
        payload.setdefault("compound", kw["compound"])
    if kw.get("conc_M") is not None:
        payload.setdefault("concentration_M", kw["conc_M"])
    if kw.get("assay"):
        payload.setdefault("assay", kw["assay"])
    if kw.get("graph"):
        payload.setdefault("graph", kw["graph"])

    subject = _subject(payload)
    rows = _occupancy_rows(payload, subject)
    evidence = _described(rows)
    modelled = [e for e in evidence if e.get("modelled")]
    not_modelled = [e for e in evidence if not e.get("modelled")]
    readouts = payload.get("readouts") if isinstance(payload.get("readouts"), Mapping) else {}
    gains = payload.get("gains") if isinstance(payload.get("gains"), Mapping) else {}

    chain: list[dict[str, Any]] = []

    # 1. the result itself ---------------------------------------------------
    named = _rate_values(readouts)
    chain.append(
        _link(
            "result",
            "COMPUTED",
            "PREDICTION",
            (
                f"{subject['assay']} readout for "
                f"{subject['compound'] or 'vehicle'}"
                + (
                    f" at {subject['concentration_M']:.2e} M"
                    if subject["concentration_M"]
                    else ""
                )
                + ". It is a simulation output, not a measurement from a fly."
            ),
            detail={
                "assay": subject["assay"],
                "compound": subject["compound"],
                "concentration_M": subject["concentration_M"],
                "values": dict(sorted(named.items())),
            },
            assumptions=["Every number downstream of this point is model output."],
        )
    )

    # 2. parameter sources ---------------------------------------------------
    sources = [
        {
            "receptor": e.get("receptor"),
            "param_type": e.get("param_type"),
            "param_value_M": e.get("param_value_M"),
            "relation": e.get("relation"),
            "species": e.get("species"),
            "evidence_tier": e.get("evidence_tier"),
            "source": e.get("source"),
            "doi": e.get("doi"),
            "pmid": e.get("pmid"),
            "classification": e.get("classification"),
        }
        for e in evidence
    ]
    by_type: dict[str, int] = {}
    for e in modelled:
        key = str(e.get("param_type") or "unknown")
        by_type[key] = by_type.get(key, 0) + 1
    chain.append(
        _link(
            "parameter_source",
            "LITERATURE-DERIVED" if modelled else "UNKNOWN",
            "LITERATURE" if modelled else "NOT MODELLED",
            (
                f"{len(modelled)} of {len(evidence)} receptor rows carry a sourced, typed "
                "potency or affinity parameter; the rest are class placeholders and are "
                "reported as not modelled, never as zero."
            ),
            detail={
                "n_rows": len(evidence),
                "n_sourced": len(modelled),
                "n_not_modelled": len(not_modelled),
                "by_param_type": by_type,
            },
            sources=sources,
            unknowns=(
                [f"no sourced value at {e.get('receptor')}" for e in not_modelled] or []
            ),
        )
    )

    # 3. engagement transformation ------------------------------------------
    models = sorted({str(e.get("engagement_model")) for e in evidence if e.get("engagement_model")})
    chain.append(
        _link(
            "engagement_transformation",
            "COMPUTED",
            "MODEL-DERIVED",
            (
                "Engagement is the Hill expression C^n / (value^n + C^n) applied to that "
                "typed parameter. Only a Kd/Ki row yields fractional receptor occupancy; "
                "an EC50/IC50 row yields normalised functional engagement, which is not "
                "the same quantity."
            ),
            detail={
                "formula": "theta = C^n / (value^n + C^n)",
                "engagement_models": models,
                "per_receptor": [
                    {
                        "receptor": e.get("receptor"),
                        "engagement": e.get("engagement"),
                        "engagement_model": e.get("engagement_model"),
                        "param_type": e.get("param_type"),
                        "n": e.get("n"),
                    }
                    for e in evidence
                ],
            },
            assumptions=[
                "A single-site Hill relation with the library's n is assumed at every receptor.",
                "Assay conditions behind an EC50/IC50 are not carried through the transformation.",
            ],
        )
    )

    # 4. mechanism rule ------------------------------------------------------
    from flylab.pharm.mechanisms import GAIN_KEYS, mechanism_table_rows

    active_keys = {
        f"{e.get('receptor')}:{r.get('direction')}"
        for e, r in zip(evidence, rows)
        if e.get("modelled") and r.get("direction") not in (None, "none")
    }
    rules = [
        r
        for r in mechanism_table_rows()
        if f"{r['receptor']}:{r['direction']}" in active_keys
    ]
    chain.append(
        _link(
            "mechanism_rule",
            "MODEL-ASSUMPTION",
            "MODEL-ASSUMPTION",
            (
                "Engagement becomes a synaptic gain through a rule FlyLab asserts. The "
                "shape of that rule -- monotone, saturating or rise-then-fall -- has never "
                "been measured for these receptors; it is the project's central modelling "
                "choice, not a result."
            ),
            detail={
                "rules_applied": rules,
                "gains": {k: _num(gains.get(k)) for k in GAIN_KEYS} if gains else {},
                "gain_floor": 0.05,
            },
            assumptions=[
                "The occupancy-to-gain transformation is asserted, not fitted to any animal data.",
                "The same rule is applied to every synapse of a transmitter class.",
            ],
        )
    )

    # 5. expression assumption ----------------------------------------------
    coverage, mn9_gap = _expression_facts()
    chain.append(
        _link(
            "expression_assumption",
            "MODEL-ASSUMPTION",
            "MODEL-ASSUMPTION",
            (
                "The gain is applied uniformly to every cell of a transmitter class, because "
                "per-cell receptor expression is largely unmapped"
                + (f" ({coverage * 100:.1f}% of cells mapped)" if coverage is not None else "")
                + ". A cell that does not express the receptor is patched exactly like one that does."
            ),
            detail={
                "fraction_cells_mapped": coverage,
                "motor_neuron_gap": mn9_gap,
                "primary_model": "uniform gains",
            },
            assumptions=["Uniform receptor expression across the patched transmitter class."],
            unknowns=["receptor expression in MN9 and in every other adult motor neuron"],
        )
    )

    # 6. MaleCNS edges -------------------------------------------------------
    graph_name = _graph_name(subject)
    graph_detail = _graph_facts(graph_name)
    chain.append(
        _link(
            "malecns_edges",
            "OBSERVED",
            "LITERATURE",
            (
                "The wiring is measured: a hops-limited neighbourhood of the public MaleCNS "
                "v1.0 reconstruction, with a synapse-count floor. This is the one layer of "
                "the chain that is a direct observation of the system being modelled."
            ),
            detail=graph_detail,
            sources=(
                [{"source": graph_detail.get("citation"), "relation": "exact_system", "classification": "LITERATURE"}]
                if graph_detail.get("citation")
                else []
            ),
            assumptions=[
                "The cut is a neighbourhood with a weight floor, not the whole CNS.",
            ],
        )
    )

    # 7. transmitter predictions --------------------------------------------
    chain.append(
        _link(
            "transmitter_predictions",
            "MODEL-ASSUMPTION",
            "PREDICTION",
            (
                "Edge signs come from MaleCNS *predicted* consensus transmitters, not from "
                "immunostaining. A wrong prediction flips a synapse from excitatory to "
                "inhibitory, and the drug patch follows the label, not the biology."
            ),
            detail={
                "source": "MaleCNS v1.0 predicted consensus neurotransmitter per body",
                "transmitter_census": graph_detail.get("transmitters"),
                "sign_rule": "acetylcholine +1, gaba/glutamate -1, others 0 unless typed",
            },
            assumptions=["Predicted transmitter labels are taken as ground truth for the sign of every edge."],
            unknowns=["ground-truth transmitter identity for the cells in this cut"],
        )
    )

    # 8. engine --------------------------------------------------------------
    engine = subject.get("engine") or ("lif" if "spik" in str(subject.get("assay")) else "rate")
    chain.append(
        _link(
            "engine",
            "MODEL-ASSUMPTION",
            "MODEL-ASSUMPTION",
            (
                f"The {engine} engine turns patched weights into activity. Its time constants, "
                "drive and settling are teaching choices; they are not fitted to any recording."
            ),
            detail={
                "engine": engine,
                "drive_hz": _num((readouts or {}).get("drive_hz")),
                "steps": (readouts or {}).get("steps"),
                "t_ms": (readouts or {}).get("t_ms"),
                "seed": ((payload.get("provenance") or {}) or {}).get("rng_seed"),
            },
            assumptions=[
                "Sensory drive amplitude is set by hand, not measured.",
                "No fitted parameter separates the engine from the data.",
            ],
        )
    )

    # 9. readout -------------------------------------------------------------
    chain.append(
        _link(
            "readout",
            "COMPUTED",
            "PREDICTION",
            (
                "The reported rates are the engine's output on the patched graph. They are "
                "predictions about a simulation, and carry no animal units."
            ),
            detail={
                "values": dict(sorted(named.items())),
                "vehicle": (readouts or {}).get("vehicle"),
            },
            unknowns=["how any of these numbers maps to a behavioural or electrophysiological measurement"],
        )
    )

    warnings = list(payload.get("warnings") or [])
    warnings = [w for w in warnings if isinstance(w, str)]
    audit = {
        "flylab_claims_version": CLAIMS_VERSION,
        "subject": subject,
        "labels": list(LABELS),
        "classifications": list(CLASSIFICATIONS),
        "chain": chain,
        "label_counts": label_counts(chain),
        "provenance": dict(payload.get("provenance") or {}),
        "notebook_warnings": warnings,
        "warnings": [
            "This audit describes provenance, not correctness: a well-sourced parameter can "
            "still be used by a rule that is wrong.",
            "Exactly one link in the chain (the MaleCNS edges) is an observation of the "
            "system being modelled. Everything after it is assertion or computation.",
            "No link is a live-animal measurement. live_lab stays null unless a human imports a table.",
        ],
        "disclaimer": (
            "FlyLab is a simulation on a public connectome. Nothing in this chain is a "
            "measured drug effect in a living fly."
        ),
    }
    audit["fact_inference_unknown"] = fact_inference_unknown(audit=audit)
    return audit


def label_counts(chain: list[dict[str, Any]] | Mapping[str, Any]) -> dict[str, int]:
    """How many links carry each provenance label."""
    links = chain["chain"] if isinstance(chain, Mapping) else chain
    out = {label: 0 for label in LABELS}
    for link in links or []:
        label = str(link.get("label"))
        if label in out:
            out[label] += 1
    return out


# --------------------------------------------------------------------------
# fact / inference / unknown
# --------------------------------------------------------------------------
def fact_inference_unknown(
    notebook_or_result: Mapping[str, Any] | None = None,
    *,
    audit: Mapping[str, Any] | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """The three-column summary that should accompany every experiment.

    ``facts`` are things somebody measured (the MaleCNS topology and the
    sourced library parameters); ``model_inference`` is what FlyLab computed
    from them (the predicted rate changes and the circuit indices); ``unknown``
    is what neither the model nor its sources can supply -- the CNS free
    concentration after dosing, receptor expression in MN9, and biological
    variance.
    """
    if audit is None:
        audit = claim_audit(notebook_or_result, **kw)
    by_step = {link["step"]: link for link in audit.get("chain") or []}
    subject = dict(audit.get("subject") or {})

    facts: list[dict[str, Any]] = []
    edges = by_step.get("malecns_edges") or {}
    edge_detail = dict(edges.get("detail") or {})
    facts.append(
        {
            "kind": "connectome",
            "label": "OBSERVED",
            "classification": "LITERATURE",
            "statement": (
                f"The MaleCNS v1.0 cut used here has {edge_detail.get('n_nodes')} cells and "
                f"{edge_detail.get('n_edges')} edges above a {edge_detail.get('min_weight')}-synapse floor."
            ),
            "detail": edge_detail,
            "source": edge_detail.get("citation"),
        }
    )
    for row in (by_step.get("parameter_source") or {}).get("sources") or []:
        if row.get("classification") == "NOT MODELLED":
            continue
        facts.append(
            {
                "kind": "parameter",
                "label": "LITERATURE-DERIVED",
                "classification": "LITERATURE",
                "statement": (
                    f"{row.get('receptor')}: {row.get('param_type')} = "
                    f"{row.get('param_value_M')} M ({row.get('species')})."
                ),
                "detail": row,
                "source": row.get("source"),
            }
        )

    inference: list[dict[str, Any]] = []
    readouts = dict((by_step.get("readout") or {}).get("detail") or {})
    values = dict(readouts.get("values") or {})
    vehicle = readouts.get("vehicle") if isinstance(readouts.get("vehicle"), Mapping) else {}
    for key in sorted(values):
        if not key.endswith("_hz"):
            continue  # gains get their own entry below
        base = _num((vehicle or {}).get(key))
        treated = _num(values.get(key))
        if treated is None:
            continue
        change = None
        if base is not None and abs(base) > 1e-12:
            change = (treated - base) / base
        inference.append(
            {
                "kind": "rate_change",
                "label": "COMPUTED",
                "classification": "PREDICTION",
                "readout": key,
                "vehicle": base,
                "treated": treated,
                "relative_change": change,
                "statement": (
                    f"Predicted {key}: {treated:.3f}"
                    + (f" against a vehicle {base:.3f}" if base is not None else "")
                    + (f" ({change * 100:+.1f}%)" if change is not None else "")
                    + "."
                ),
            }
        )
    for name, block in (("mechanism_rule", "gain"), ("engagement_transformation", "engagement")):
        detail = dict((by_step.get(name) or {}).get("detail") or {})
        if block == "gain" and detail.get("gains"):
            inference.append(
                {
                    "kind": "gain",
                    "label": "MODEL-ASSUMPTION",
                    "classification": "MODEL-ASSUMPTION",
                    "statement": (
                        "Synaptic gains under the asserted mechanism rule: "
                        + ", ".join(
                            f"{k}={v:.3f}" for k, v in sorted(detail["gains"].items()) if v is not None
                        )
                        + "."
                    ),
                    "detail": detail.get("gains"),
                }
            )
        if block == "engagement" and detail.get("per_receptor"):
            inference.append(
                {
                    "kind": "engagement",
                    "label": "COMPUTED",
                    "classification": "MODEL-DERIVED",
                    "statement": (
                        "Receptor engagement is a Hill transform of the cited parameters, "
                        "not a measured occupancy."
                    ),
                    "detail": detail.get("per_receptor"),
                }
            )
    indices = kw.get("indices") or (
        notebook_or_result.get("indices") if isinstance(notebook_or_result, Mapping) else None
    )
    for key, value in sorted((indices or {}).items()):
        if _num(value) is None:
            continue
        inference.append(
            {
                "kind": "index",
                "label": "COMPUTED",
                "classification": "MODEL-DERIVED",
                "statement": f"Circuit index {key} = {_num(value):.3f}.",
                "detail": {key: _num(value)},
            }
        )

    unknown = [
        {
            "kind": "unknown",
            "label": "UNKNOWN",
            "classification": "NOT MODELLED",
            "id": item["id"],
            "statement": item["statement"],
            "why": item["why"],
            "would_resolve": item["would_resolve"],
        }
        for item in UNIVERSAL_UNKNOWNS
    ]
    for row in (by_step.get("parameter_source") or {}).get("sources") or []:
        if row.get("classification") != "NOT MODELLED":
            continue
        unknown.append(
            {
                "kind": "unknown",
                "label": "UNKNOWN",
                "classification": "NOT MODELLED",
                "id": f"no_value_{row.get('receptor')}",
                "statement": f"No sourced value for this compound at {row.get('receptor')}.",
                "why": "The row is a class placeholder; FlyLab reports it as not modelled rather than as zero.",
                "would_resolve": f"A published potency or affinity measurement at {row.get('receptor')}.",
            }
        )

    return {
        "compound": subject.get("compound"),
        "concentration_M": subject.get("concentration_M"),
        "assay": subject.get("assay"),
        "facts": facts,
        "model_inference": inference,
        "unknown": unknown,
        "counts": {
            "facts": len(facts),
            "model_inference": len(inference),
            "unknown": len(unknown),
        },
        "note": (
            "Facts are measurements somebody made; inference is what this model computed "
            "from them; unknowns are gaps no part of the chain can fill."
        ),
    }


# --------------------------------------------------------------------------
# background facts (cheap, cached by the callers that need them)
# --------------------------------------------------------------------------
def _graph_facts(name: str) -> dict[str, Any]:
    try:
        from flylab.circuit.rate import load_graph

        g = load_graph(name)
    except Exception:  # pragma: no cover - a missing graph must not break an audit
        return {"graph": name, "available": False}
    nodes = g.get("nodes") or []
    census: dict[str, int] = {}
    for node in nodes:
        nt = str((node or {}).get("nt") or "unclear")
        census[nt] = census.get(nt, 0) + 1
    return {
        "graph": name,
        "available": True,
        "map": g.get("map"),
        "citation": g.get("citation"),
        "n_nodes": g.get("n_nodes") or len(nodes),
        "n_edges": g.get("n_edges"),
        "min_weight": g.get("min_weight"),
        "hops": g.get("hops"),
        "transmitters": dict(sorted(census.items(), key=lambda kv: -kv[1])),
    }


def _expression_facts() -> tuple[float | None, dict[str, Any] | None]:
    try:
        from flylab.pharm.expression import expression_table

        table = expression_table()
    except Exception:  # pragma: no cover - optional module
        return None, None
    coverage = _num((table.get("coverage") or {}).get("fraction_known_overall"))
    gap = table.get("motor_neuron_gap")
    slim = None
    if isinstance(gap, Mapping):
        slim = {
            "status": gap.get("status"),
            "description": gap.get("description"),
            "consequence_for_flylab": gap.get("consequence_for_flylab"),
        }
    return coverage, slim


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------
def to_markdown(audit: Mapping[str, Any]) -> str:
    """The audit as a short markdown table, for a supplement or a terminal."""
    subject = dict(audit.get("subject") or {})
    lines = [
        f"# Claim provenance -- {subject.get('compound') or 'vehicle'}"
        + (f" at {subject['concentration_M']:.2e} M" if subject.get("concentration_M") else ""),
        "",
        "| step | label | chip | what it contributes |",
        "|---|---|---|---|",
    ]
    for link in audit.get("chain") or []:
        lines.append(
            f"| {link['step']} | {link['label']} | {link['classification']} | "
            + str(link["statement"]).replace("|", "/")
            + " |"
        )
    fiu = audit.get("fact_inference_unknown") or {}
    for key, title in (("facts", "Fact"), ("model_inference", "Inference"), ("unknown", "Unknown")):
        lines += ["", f"## {title}", ""]
        for item in fiu.get(key) or []:
            lines.append(f"- {item.get('statement')}")
    return "\n".join(lines) + "\n"
