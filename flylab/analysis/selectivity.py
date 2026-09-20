"""Circuit-level selectivity: does the *circuit* effect arrive before the
vertebrate receptor does?

The receptor scorecard FlyLab already ships answers one question: at a given
concentration, how much insect target is occupied versus vertebrate target.
That is a **receptor** selectivity index -- essentially ``log10(EC50_vert /
EC50_insect)`` -- and it knows nothing about wiring.

This module adds the circuit counterpart:

``C_circuit``
    the lowest concentration at which the simulated network readout moves by
    ``threshold_frac`` of its vehicle value (monotone-interpolated on log
    concentration);
``C_vert``
    the concentration at which the most potent *vertebrate* receptor in the
    teaching library reaches ``vert_occ_limit`` occupancy (closed form from the
    Hill equation);
``circuit SI``
    ``log10(C_vert) - log10(C_circuit)``: the log10 width of the window in
    which the fly circuit is already perturbed and the vertebrate receptor is
    not.  Positive = circuit effect precedes vertebrate occupancy.

Both indices are reported side by side so the paper can ask whether the
neighborhood amplifies or buffers what the receptor numbers promise.  Every
number here is ``model_derived``: it is a property of the teaching EC50
library, the gain patch rules and one hops-limited MaleCNS neighborhood.
"""
from __future__ import annotations

import math
import time
from typing import Any, Iterable, Sequence

from flylab.analysis.nullmodels import DEFAULT_GRAPH, DEFAULT_READOUT, drug_effect

__all__ = [
    "DEFAULT_CONCS",
    "receptor_selectivity_table",
    "circuit_threshold_conc",
    "vertebrate_threshold_conc",
    "circuit_selectivity_index",
    "graph_composition",
    "selectivity_landscape",
]

#: 13-point half-decade ladder, 1e-10 .. 1e-4 M
DEFAULT_CONCS: tuple[float, ...] = tuple(10.0 ** (-10 + 0.5 * i) for i in range(13))

BASE_WARNINGS = [
    "model_derived: every index here comes from the teaching EC50 library, the "
    "gain patch rules and one hops-limited MaleCNS neighborhood. None of it is "
    "a measured selectivity margin and none of it is a safety statement.",
    "Receptors with evidence_tier 'class_placeholder' carry no number at all "
    "(schema v3): they report engagement None, are excluded from every ratio, "
    "and any row that would have depended on one is flagged or skipped.",
]


def _is_placeholder(row: dict[str, Any]) -> bool:
    """True when a row carries no sourced value (schema v3: engagement None)."""
    return (
        row.get("engagement", row.get("occupancy")) is None
        or row.get("param_value_M", row.get("ec50_M")) is None
        or row.get("evidence_tier") == "class_placeholder"
        or str(row.get("direction") or "none") == "none"
    )


def _engagement(row: dict[str, Any]) -> float | None:
    value = row.get("engagement", row.get("occupancy"))
    return None if value is None else float(value)


# --------------------------------------------------------------------------
# receptor level
# --------------------------------------------------------------------------
def _compound_rows(compound: str, conc_M: float, library: dict[str, Any] | None = None):
    from flylab.pharm.occupancy import compare_compound

    occ = compare_compound(compound, conc_M, library=library) if library else compare_compound(compound, conc_M)
    return occ


def _best_pair(occ: dict[str, Any]) -> tuple[str | None, dict[str, Any] | None]:
    """Insect/vertebrate pair with the largest log10 EC50 ratio.

    Pairs that rest on a class placeholder lose to any real pair, whatever
    their ratio.
    """
    pairs = occ.get("selectivity") or {}
    if not pairs:
        return None, None
    # Schema v3: a pair with a placeholder on either side has ratio None and is
    # excluded from the numbers entirely (it is not a log-ratio of 0.00).
    real = {
        k: v for k, v in pairs.items()
        if not v.get("placeholder") and v.get("log10_ec50_ratio_vert_over_insect") is not None
    }
    if not real:
        return None, None
    name = max(real, key=lambda k: real[k]["log10_ec50_ratio_vert_over_insect"])
    return name, real[name]


def receptor_selectivity_table(
    conc_M: float = 1e-6,
    compounds: Sequence[str] | None = None,
    library: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Receptor scorecard for every compound in the library.

    One row per compound: the highest insect occupancy and the highest
    vertebrate occupancy at ``conc_M``, the best insect/vertebrate pair's
    ``log10(EC50_vert / EC50_insect)`` (the receptor SI), and whether that row
    rests on a class placeholder.
    """
    from flylab.pharm.occupancy import list_compounds

    names = list(compounds) if compounds is not None else list_compounds(library)
    rows: list[dict[str, Any]] = []
    for name in names:
        occ = _compound_rows(name, conc_M, library)
        recs = occ["receptors"]
        ins = [r for r in recs if str(r["receptor"]).startswith("insect_")]
        vert = [r for r in recs if str(r["receptor"]).startswith("vertebrate_")]
        ins_real = [r for r in ins if not _is_placeholder(r)]
        vert_real = [r for r in vert if not _is_placeholder(r)]
        # not-modelled rows never win a max and never become a 0
        ins_num = [r for r in ins if _engagement(r) is not None]
        vert_num = [r for r in vert if _engagement(r) is not None]
        top_ins = max(ins_num, key=_engagement) if ins_num else None
        top_vert = max(vert_num, key=_engagement) if vert_num else None
        pair_name, pair = _best_pair(occ)
        rows.append(
            {
                "compound": occ["key"],
                "name": occ["compound"],
                "class": occ.get("class"),
                "conc_M": float(conc_M),
                "max_insect_engagement": _engagement(top_ins) if top_ins else None,
                "max_vertebrate_engagement": _engagement(top_vert) if top_vert else None,
                # deprecated aliases (same values; None means "not modelled")
                "max_insect_occupancy": _engagement(top_ins) if top_ins else None,
                "insect_receptor": top_ins["receptor"] if top_ins else None,
                "max_vertebrate_occupancy": _engagement(top_vert) if top_vert else None,
                "vertebrate_receptor": top_vert["receptor"] if top_vert else None,
                "best_pair": pair_name,
                "receptor_si_log10": (
                    float(pair["log10_ec50_ratio_vert_over_insect"]) if pair else None
                ),
                "insect_param_type": pair.get("insect_param_type") if pair else None,
                "vertebrate_param_type": pair.get("vertebrate_param_type") if pair else None,
                "insect_ec50_M": float(pair["insect_ec50_M"]) if pair else None,
                "vertebrate_ec50_M": float(pair["vertebrate_ec50_M"]) if pair else None,
                "placeholder": bool(pair.get("placeholder")) if pair else True,
                "skipped_reason": None if pair else (
                    "every insect/vertebrate pair for this compound rests on a placeholder row, "
                    "so no receptor selectivity index is defined"
                ),
                "insect_all_placeholder": not ins_real,
                "vertebrate_all_placeholder": not vert_real,
                "evidence_tier": pair.get("evidence_tier") if pair else None,
            }
        )
    warnings = list(BASE_WARNINGS)
    flagged = [r["compound"] for r in rows if r["placeholder"]]
    if flagged:
        warnings.append(
            "best pair rests on a class placeholder for: " + ", ".join(sorted(flagged))
        )
    return {
        "conc_M": float(conc_M),
        "n_compounds": len(rows),
        "rows": rows,
        "label": "model_derived",
        "warnings": warnings,
    }


def vertebrate_threshold_conc(
    compound: str,
    occ_limit: float = 0.2,
    library: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Concentration at which the most potent vertebrate target hits ``occ_limit``.

    Closed form from the Hill equation: ``C = EC50 * (th / (1 - th)) ** (1/n)``.
    Real (non-placeholder) vertebrate rows win over placeholders.
    """
    occ = _compound_rows(compound, 1e-6, library)
    vert = [r for r in occ["receptors"] if str(r["receptor"]).startswith("vertebrate_")]
    if not vert:
        return {"conc_M": None, "receptor": None, "placeholder": True, "occ_limit": float(occ_limit)}
    th = min(max(float(occ_limit), 1e-9), 1 - 1e-9)
    cands = []
    for r in vert:
        value = r.get("param_value_M", r.get("ec50_M"))
        if value is None:
            continue  # not modelled: no threshold concentration exists
        n = float(r.get("n", 1.0)) or 1.0
        cands.append((float(value) * (th / (1.0 - th)) ** (1.0 / n), r))
    if not cands:
        return {"conc_M": None, "receptor": None, "placeholder": True,
                "occ_limit": float(th),
                "skipped_reason": "every vertebrate row for this compound is a placeholder"}
    real = [(c, r) for c, r in cands if not _is_placeholder(r)]
    pool = real or cands
    conc, row = min(pool, key=lambda t: t[0])
    return {
        "conc_M": float(conc),
        "receptor": row["receptor"],
        "param_type": row.get("param_type"),
        "param_value_M": float(row.get("param_value_M", row["ec50_M"])),
        "ec50_M": float(row.get("param_value_M", row["ec50_M"])),
        "placeholder": _is_placeholder(row),
        "occ_limit": float(th),
    }


# --------------------------------------------------------------------------
# circuit level
# --------------------------------------------------------------------------
def circuit_threshold_conc(
    assay: str = "subgraph",
    compound: str = "imidacloprid",
    readout: str = "auto",
    threshold_frac: float = 0.5,
    concs: Iterable[float] | None = None,
    graph: str | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Lowest concentration whose relative readout change reaches ``threshold_frac``.

    The relative change is ``|treated - vehicle| / |vehicle|`` with both arms
    run on the same graph at the same drive.  The crossing is interpolated
    linearly on ``log10(conc)`` between the last point below threshold and the
    first point at or above it (only when the curve is locally increasing;
    otherwise the ladder point itself is returned).

    Returns ``conc_M`` (``None`` if the threshold is never reached) and the
    whole ``curve``.
    """
    if readout == "auto":
        readout = DEFAULT_READOUT.get(assay, "mn9_hz")
    graph_name = graph if graph is not None else DEFAULT_GRAPH.get(assay, "named")
    xs = [float(c) for c in (concs if concs is not None else DEFAULT_CONCS)]
    warnings = list(BASE_WARNINGS)
    t0 = time.perf_counter()

    curve: list[dict[str, Any]] = []
    rel: list[float | None] = []
    vehicle_val: float | None = None
    for c in xs:
        res = drug_effect(assay, compound, c, readout, graph=graph_name, **kw)
        treated, vehicle = res["treated"], res["vehicle"]
        vehicle_val = vehicle if vehicle is not None else vehicle_val
        r = None
        if treated is not None and vehicle is not None and abs(vehicle) > 1e-12:
            r = abs(treated - vehicle) / abs(vehicle)
        rel.append(r)
        curve.append(
            {
                "conc_M": c,
                "treated": treated,
                "vehicle": vehicle,
                "effect": res["effect"],
                "rel_change": r,
            }
        )

    if vehicle_val is not None and abs(vehicle_val) <= 1e-12:
        warnings.append(
            f"vehicle {readout} is ~0 on graph {graph_name!r}; a relative "
            "threshold is undefined and no crossing is reported."
        )

    thr = float(threshold_frac)
    hit_conc: float | None = None
    hit_index: int | None = None
    for i, r in enumerate(rel):
        if r is not None and r >= thr:
            hit_index = i
            break
    finite = [(r, c) for r, c in zip(rel, xs) if r is not None]
    max_rel, max_rel_conc = max(finite, default=(None, None))
    if hit_index is None:
        warnings.append(
            f"{compound} never reaches a {thr:.0%} change of {readout} on "
            f"{graph_name} within {xs[0]:.1e}-{xs[-1]:.1e} M"
            + (f" (largest change {max_rel:.1%} at {max_rel_conc:.1e} M)." if max_rel is not None else ".")
        )
    elif hit_index == 0:
        hit_conc = xs[0]
        warnings.append(
            "threshold is already exceeded at the lowest tested concentration; "
            f"{hit_conc:.1e} M is an upper bound, not a crossing."
        )
    else:
        lo_r, hi_r = rel[hit_index - 1], rel[hit_index]
        lo_x, hi_x = math.log10(xs[hit_index - 1]), math.log10(xs[hit_index])
        if lo_r is not None and hi_r is not None and hi_r > lo_r:
            frac = (thr - lo_r) / (hi_r - lo_r)
            hit_conc = float(10.0 ** (lo_x + frac * (hi_x - lo_x)))
        else:
            hit_conc = xs[hit_index]
            warnings.append(
                "the readout is not locally monotone at the crossing; the "
                "ladder point is reported instead of an interpolation."
            )

    return {
        "assay": assay,
        "compound": compound,
        "readout": readout,
        "graph": graph_name,
        "threshold_frac": thr,
        "conc_M": hit_conc,
        "crossing_index": hit_index,
        "max_rel_change": max_rel,
        "max_rel_change_conc_M": max_rel_conc,
        "concs_M": xs,
        "curve": curve,
        "vehicle": vehicle_val,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }


def circuit_selectivity_index(
    compound: str = "imidacloprid",
    assay: str = "subgraph",
    readout: str = "auto",
    vert_occ_limit: float = 0.2,
    graph: str | None = None,
    concs: Iterable[float] | None = None,
    threshold_frac: float = 0.5,
    **kw: Any,
) -> dict[str, Any]:
    """log10 width of the window ``[C_circuit_effect, C_vertebrate_occupancy]``.

    Positive means the fly circuit is already perturbed at concentrations where
    the vertebrate receptor is still below ``vert_occ_limit`` occupancy.  The
    receptor-level SI (``log10(EC50_vert / EC50_insect)`` of the best pair) is
    reported next to it so the two can be compared directly.
    """
    thr = circuit_threshold_conc(
        assay, compound, readout=readout, threshold_frac=threshold_frac,
        concs=concs, graph=graph, **kw,
    )
    vert = vertebrate_threshold_conc(compound, vert_occ_limit)
    occ = _compound_rows(compound, 1e-6)
    pair_name, pair = _best_pair(occ)

    c_circ, c_vert = thr["conc_M"], vert["conc_M"]
    si = None
    if c_circ and c_vert:
        si = float(math.log10(c_vert) - math.log10(c_circ))

    warnings = list(thr["warnings"])
    if si is None:
        warnings.append(
            "circuit selectivity index is undefined: the circuit threshold or "
            "the vertebrate occupancy threshold could not be located."
        )
    if vert.get("placeholder"):
        warnings.append(
            "the vertebrate arm of this compound is a class placeholder; the "
            "circuit index is an artefact of that placeholder, not a margin."
        )
    warnings.append(
        "A positive circuit SI means only that this simulation moves before the "
        "teaching library's vertebrate receptor fills; it is not a therapeutic "
        "or environmental safety margin."
    )

    return {
        "compound": occ["key"],
        "name": occ["compound"],
        "class": occ.get("class"),
        "assay": assay,
        "readout": thr["readout"],
        "graph": thr["graph"],
        "threshold_frac": float(threshold_frac),
        "vert_occ_limit": float(vert_occ_limit),
        "circuit_threshold_M": c_circ,
        "max_rel_change": thr["max_rel_change"],
        "max_rel_change_conc_M": thr["max_rel_change_conc_M"],
        "vertebrate_threshold_M": c_vert,
        "vertebrate_receptor": vert.get("receptor"),
        "circuit_si_log10": si,
        "receptor_si_log10": (
            float(pair["log10_ec50_ratio_vert_over_insect"]) if pair else None
        ),
        "receptor_pair": pair_name,
        "receptor_pair_placeholder": bool(pair.get("placeholder")) if pair else True,
        "si_gap_circuit_minus_receptor": (
            float(si - pair["log10_ec50_ratio_vert_over_insect"])
            if (si is not None and pair and pair.get("log10_ec50_ratio_vert_over_insect") not in (None, float("-inf")))
            else None
        ),
        "curve": thr["curve"],
        "label": "model_derived",
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# landscape
# --------------------------------------------------------------------------
def graph_composition(graph_name: str = "named") -> dict[str, Any]:
    """E/I composition of a graph: node fractions and synapse-weighted fractions.

    The weighted fractions use the **presynaptic** transmitter of each edge
    times its synapse count, which is what the rate model actually scales.
    """
    from flylab.circuit.rate import load_graph

    g = load_graph(graph_name)
    nodes = g["nodes"]
    nt_of = {n["bodyId"]: (n.get("consensus_nt") or "unclear") for n in nodes}
    counts: dict[str, int] = {}
    for n in nodes:
        counts[n.get("consensus_nt") or "unclear"] = counts.get(n.get("consensus_nt") or "unclear", 0) + 1
    total_n = max(len(nodes), 1)

    w: dict[str, float] = {}
    total_w = 0.0
    for e in g["edges"]:
        nt = nt_of.get(e["pre"], "unclear")
        w[nt] = w.get(nt, 0.0) + float(e["weight"])
        total_w += float(e["weight"])
    total_w = total_w or 1.0

    def frac(d, tot, key):
        return float(d.get(key, 0) / tot)

    exc = frac(counts, total_n, "acetylcholine")
    inh = frac(counts, total_n, "gaba") + frac(counts, total_n, "glutamate")
    exc_w = frac(w, total_w, "acetylcholine")
    inh_w = frac(w, total_w, "gaba") + frac(w, total_w, "glutamate")
    return {
        "graph": graph_name,
        "n_nodes": g["n_nodes"],
        "n_edges": g["n_edges"],
        "node_fraction": {k: v / total_n for k, v in sorted(counts.items())},
        "synapse_fraction": {k: v / total_w for k, v in sorted(w.items())},
        "frac_ach_nodes": exc,
        "frac_gaba_nodes": frac(counts, total_n, "gaba"),
        "frac_glu_nodes": frac(counts, total_n, "glutamate"),
        "frac_ach_synapses": exc_w,
        "frac_gaba_synapses": frac(w, total_w, "gaba"),
        "frac_glu_synapses": frac(w, total_w, "glutamate"),
        "ei_ratio_nodes": (exc / inh) if inh > 0 else None,
        "ei_ratio_synapses": (exc_w / inh_w) if inh_w > 0 else None,
    }


def selectivity_landscape(
    assay: str = "subgraph",
    graphs: Sequence[str] = ("named", "taste_motor"),
    compounds: Sequence[str] | None = None,
    readout: str = "auto",
    vert_occ_limit: float = 0.2,
    threshold_frac: float = 0.5,
    concs: Iterable[float] | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Receptor SI and circuit SI for every compound on every graph.

    Each row also carries the graph's E/I composition (node and synapse
    weighted), which is the covariate for the question "does the neighborhood
    amplify or buffer receptor-level selectivity?".  Compounds whose insect
    rows are all class placeholders cannot move a fly circuit at all and are
    reported with ``skipped_reason`` instead of an index.
    """
    from flylab.pharm.occupancy import list_compounds

    t0 = time.perf_counter()
    names = list(compounds) if compounds is not None else list_compounds()
    table = receptor_selectivity_table(1e-6, compounds=names)
    by_compound = {r["compound"]: r for r in table["rows"]}
    comps = {gname: graph_composition(gname) for gname in graphs}

    rows: list[dict[str, Any]] = []
    warnings = list(BASE_WARNINGS)
    skipped: list[str] = []
    for name in names:
        rec = by_compound.get(name.lower().strip(), {})
        skip = None
        if rec.get("insect_all_placeholder"):
            skip = (
                "every insect receptor row for this compound is a class "
                "placeholder: it cannot move the fly circuit by construction"
            )
            skipped.append(name)
        for gname in graphs:
            comp = comps[gname]
            base = {
                "compound": name,
                "name": rec.get("name"),
                "class": rec.get("class"),
                "graph": gname,
                "assay": assay,
                "n_nodes": comp["n_nodes"],
                "n_edges": comp["n_edges"],
                "frac_ach_nodes": comp["frac_ach_nodes"],
                "frac_gaba_nodes": comp["frac_gaba_nodes"],
                "frac_glu_nodes": comp["frac_glu_nodes"],
                "frac_ach_synapses": comp["frac_ach_synapses"],
                "frac_gaba_synapses": comp["frac_gaba_synapses"],
                "frac_glu_synapses": comp["frac_glu_synapses"],
                "ei_ratio_synapses": comp["ei_ratio_synapses"],
                "receptor_si_log10": rec.get("receptor_si_log10"),
                "receptor_pair": rec.get("best_pair"),
                "max_insect_occupancy": rec.get("max_insect_occupancy"),
                "max_vertebrate_occupancy": rec.get("max_vertebrate_occupancy"),
            }
            if skip:
                rows.append(
                    {
                        **base,
                        "circuit_threshold_M": None,
                        "max_rel_change": None,
                        "max_rel_change_conc_M": None,
                        "vertebrate_threshold_M": None,
                        "circuit_si_log10": None,
                        "si_gap_circuit_minus_receptor": None,
                        "skipped_reason": skip,
                    }
                )
                continue
            si = circuit_selectivity_index(
                name, assay=assay, readout=readout, vert_occ_limit=vert_occ_limit,
                graph=gname, concs=concs, threshold_frac=threshold_frac, **kw,
            )
            rows.append(
                {
                    **base,
                    "readout": si["readout"],
                    "circuit_threshold_M": si["circuit_threshold_M"],
                    "max_rel_change": si["max_rel_change"],
                    "max_rel_change_conc_M": si["max_rel_change_conc_M"],
                    "vertebrate_threshold_M": si["vertebrate_threshold_M"],
                    "vertebrate_receptor": si["vertebrate_receptor"],
                    "circuit_si_log10": si["circuit_si_log10"],
                    "si_gap_circuit_minus_receptor": si["si_gap_circuit_minus_receptor"],
                    "circuit_beats_receptor": (
                        None
                        if (si["circuit_si_log10"] is None or si["receptor_si_log10"] is None)
                        else bool(si["circuit_si_log10"] > si["receptor_si_log10"])
                    ),
                    "skipped_reason": None,
                }
            )
    if skipped:
        warnings.append("skipped (insect rows all placeholder): " + ", ".join(sorted(set(skipped))))
    warnings.append(
        "Circuit SI and receptor SI are not interchangeable: the receptor index "
        "is a ratio of two teaching EC50s, the circuit index also carries the "
        "gain rules, the drive and the neighborhood's E/I composition."
    )
    return {
        "assay": assay,
        "graphs": list(graphs),
        "vert_occ_limit": float(vert_occ_limit),
        "threshold_frac": float(threshold_frac),
        "n_rows": len(rows),
        "composition": comps,
        "rows": rows,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }
