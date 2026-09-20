"""Transcriptome-weighted patches: a sensitivity analysis, not the primary model.

Data source: ``data/literature/receptor_expression_by_class.yaml``.  Read its
header before reading this module.  That file is an explicit **cross-atlas
inference**: none of the cited single-cell studies (Fly Cell Atlas 2022, Davie
2018, Allen 2020, Croset 2018, Konstantinides 2018, Enell 2007) is annotated
with MaleCNS ``superclass`` labels, most superclass rows are ``unknown``, and
motor-neuron receptor expression - the class MN9 belongs to - is a **confirmed
gap** with no source at all.

Consequently:

* The rest of FlyLab keeps **uniform** gains.  Nothing here changes any other
  assay's default.
* Every output of this module carries :data:`EXPRESSION_WARNING`.
* An ``absent`` call is never inferred from single-nucleus dropout; the dataset
  only ever records ``unknown``, and ``unknown`` maps to weight 1.0 by default,
  i.e. "assume the uniform model until measured".

What the layer does do is let one ask: *if* expression were as coarsely mapped
here, how much would the circuit readout move?  If the answer is "little", the
uniform-gain simplification is defensible; that is the whole claim.
"""

from __future__ import annotations

import copy
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import yaml

from flylab.pharm.mechanisms import default_gains, gains_from_occupancy

__all__ = [
    "EXPRESSION_PATH",
    "EXPRESSION_WARNING",
    "LEVEL_MAP",
    "RECEPTOR_GENES",
    "load_expression",
    "expression_table",
    "expression_weights",
    "weighted_gains",
    "run_weighted_subgraph_assay",
]

_DATA_DIRS = [
    Path(__file__).resolve().parents[2] / "data" / "literature",
    Path("data/literature"),
]
EXPRESSION_PATH = next(
    (
        d / "receptor_expression_by_class.yaml"
        for d in _DATA_DIRS
        if (d / "receptor_expression_by_class.yaml").exists()
    ),
    _DATA_DIRS[0] / "receptor_expression_by_class.yaml",
)

EXPRESSION_WARNING = (
    "COARSE CROSS-ATLAS INFERENCE - SENSITIVITY ANALYSIS ONLY. Per-node "
    "receptor weights come from class-level expression statements that were "
    "mapped onto MaleCNS superclasses by hand; no cited atlas uses MaleCNS "
    "labels, most superclasses are 'unknown', and receptor expression for "
    "identified adult motor neurons (the class MN9 belongs to) is a confirmed "
    "gap in the literature. The primary FlyLab model keeps UNIFORM gains; this "
    "layer exists to show how much that assumption could matter, and must not "
    "be reported as a per-cell measurement."
)

MOTOR_NEURON_WARNING = (
    "Motor-neuron gap: no cited source reports Rdl or nAChR expression for "
    "identified adult leg or labellar motor neurons, so MN9 and every other "
    "motor neuron carries the 'unknown' default weight. The readout of this "
    "assay is a motor neuron, so its own receptor complement is unmeasured."
)

#: Documented default ordinal -> multiplier map.
LEVEL_MAP: dict[str, float] = {
    "high": 1.0,
    "medium": 0.6,
    "low": 0.25,
    "absent": 0.0,
}

#: FlyLab receptor key -> gene-name prefixes used in the expression YAML.
RECEPTOR_GENES: dict[str, tuple[str, ...]] = {
    "insect_nAChR": ("nAChR",),
    "insect_RDL": ("Rdl",),
    "insect_GluCl": ("GluCl",),
    "insect_OctR": ("OctR", "Oct"),
    "insect_AChE": ("Ace",),
}

#: which gain each receptor's weight interpolates
RECEPTOR_GAIN = {
    "insect_nAChR": "g_ach",
    "insect_RDL": "g_gaba",
    "insect_GluCl": "g_glu",
    "insect_OctR": "g_oct",
}

AGGREGATIONS = ("mean", "max")


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
@lru_cache(maxsize=4)
def _read_expression(path_str: str) -> dict[str, Any]:
    p = Path(path_str)
    if not p.exists():  # pragma: no cover - only when the dataset is absent
        return {"by_superclass": [], "published_atlas_values": [], "_missing": str(p)}
    data = yaml.safe_load(p.read_text()) or {}
    if not isinstance(data, dict):  # pragma: no cover - defensive
        return {"by_superclass": []}
    data.setdefault("by_superclass", [])
    return data


def load_expression(path: Path | None = None) -> dict[str, Any]:
    """Deep copy of the expression dataset (defensive: never raises)."""
    return copy.deepcopy(_read_expression(str(path or EXPRESSION_PATH)))


def _genes_for(receptor: str) -> tuple[str, ...]:
    return RECEPTOR_GENES.get(receptor, (receptor,))


def _matching_genes(receptor: str, receptors_block: Mapping[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    prefixes = _genes_for(receptor)
    return [
        (gene, spec)
        for gene, spec in receptors_block.items()
        if isinstance(spec, dict) and str(gene).lower().startswith(tuple(p.lower() for p in prefixes))
    ]


def _level_weight(level: Any, level_map: Mapping[str, float], default_unknown: float) -> float | None:
    key = str(level or "unknown").strip().lower()
    if key in level_map:
        return float(level_map[key])
    return None  # unknown / unrecognised -> caller applies default_unknown


def _class_weight(
    receptor: str,
    receptors_block: Mapping[str, Any],
    level_map: Mapping[str, float],
    default_unknown: float,
    aggregation: str = "mean",
) -> tuple[float, str, list[dict[str, Any]]]:
    """``(weight, level_label, per-gene detail)`` for one superclass.

    Aggregation defaults to the **mean over genes with a known level**, not the
    max, because ``insect_nAChR`` is not one receptor: Croset et al. 2018 show
    alpha1/5/6/7 and beta1 in far more cells than alpha2/3/4, so a compound is
    not reaching a uniformly "high" target.  ``aggregation="max"`` gives the
    optimistic reading instead.
    """
    genes = _matching_genes(receptor, receptors_block)
    detail = []
    known: list[float] = []
    for gene, spec in genes:
        w = _level_weight(spec.get("level"), level_map, default_unknown)
        detail.append(
            {
                "gene": gene,
                "level": spec.get("level"),
                "weight": w,
                "fraction_cells_expressing": spec.get("fraction_cells_expressing"),
                "evidence_tier": spec.get("evidence_tier"),
                "source_key": spec.get("source_key"),
                "note": spec.get("note"),
            }
        )
        if w is not None:
            known.append(w)
    if not known:
        return float(default_unknown), "unknown", detail
    value = float(np.mean(known)) if aggregation == "mean" else float(np.max(known))
    label = min(LEVEL_MAP, key=lambda k: abs(LEVEL_MAP[k] - value))
    return value, label, detail


_INDEX_CACHE: dict[tuple, dict[str, dict[str, Any]]] = {}


def _superclass_index(
    receptor: str,
    level_map: Mapping[str, float] | None = None,
    default_unknown: float = 1.0,
    aggregation: str = "mean",
    path: Path | None = None,
) -> dict[str, dict[str, Any]]:
    lm = dict(level_map or LEVEL_MAP)
    cache_key = (
        receptor,
        tuple(sorted(lm.items())),
        float(default_unknown),
        aggregation,
        str(path or EXPRESSION_PATH),
    )
    if cache_key in _INDEX_CACHE:
        return _INDEX_CACHE[cache_key]
    data = load_expression(path)
    out: dict[str, dict[str, Any]] = {}
    for row in data.get("by_superclass") or []:
        if not isinstance(row, dict) or not row.get("superclass"):
            continue
        weight, label, detail = _class_weight(
            receptor, row.get("receptors") or {}, lm, default_unknown, aggregation
        )
        out[str(row["superclass"])] = {
            "weight": weight,
            "level": label,
            "known": any(d["weight"] is not None for d in detail),
            "genes": detail,
            "mapping": row.get("mapping"),
            "mapped_from": row.get("mapped_from"),
            "notes": row.get("notes"),
        }
    _INDEX_CACHE[cache_key] = out
    return out


# ---------------------------------------------------------------------------
# table + coverage
# ---------------------------------------------------------------------------
def _graph_superclasses(graph: str | None = None) -> list[str | None]:
    from flylab.circuit.rate import load_graph

    g = load_graph(graph)
    return [n.get("superclass") for n in g["nodes"]]


def coverage(
    receptor: str | None = None,
    graphs: Iterable[str] | None = None,
    default_unknown: float = 1.0,
    aggregation: str = "mean",
    path: Path | None = None,
) -> dict[str, Any]:
    """How many cells of the loaded graphs sit in classes with known expression."""
    from flylab.circuit.rate import GRAPHS

    names = list(graphs) if graphs is not None else list(GRAPHS)
    receptors = [receptor] if receptor else list(RECEPTOR_GENES)
    per_graph: dict[str, Any] = {}
    totals = {"n_cells": 0, "n_known": 0}
    for name in names:
        try:
            classes = _graph_superclasses(name)
        except Exception:  # pragma: no cover - graph not committed in this checkout
            continue
        block: dict[str, Any] = {"n_cells": len(classes)}
        for rec in receptors:
            index = _superclass_index(rec, None, default_unknown, aggregation, path)
            known = sum(1 for sc in classes if sc and index.get(sc, {}).get("known"))
            block[rec] = {
                "n_known": known,
                "n_unknown": len(classes) - known,
                "fraction_known": known / len(classes) if classes else 0.0,
            }
            totals["n_known"] += known
            totals["n_cells"] += len(classes)
        per_graph[name] = block
    overall = totals["n_known"] / totals["n_cells"] if totals["n_cells"] else 0.0
    return {
        "graphs": per_graph,
        "receptors": receptors,
        "fraction_known_overall": overall,
        "aggregation": aggregation,
        "default_unknown": default_unknown,
        "warning": EXPRESSION_WARNING,
    }


def expression_table(
    graphs: Iterable[str] | None = None,
    level_map: Mapping[str, float] | None = None,
    default_unknown: float = 1.0,
    aggregation: str = "mean",
    path: Path | None = None,
) -> dict[str, Any]:
    """JSON-ready rows per MaleCNS superclass x FlyLab receptor, plus coverage.

    Each row: ``superclass``, ``receptor``, ``level``, ``weight``,
    ``fraction_cells_expressing`` (published fraction or ``None`` - the dataset
    records almost none), ``evidence_tier``, ``source_key``, ``mapping`` and
    the per-gene detail behind the aggregate.
    """
    data = load_expression(path)
    lm = dict(level_map or LEVEL_MAP)
    rows: list[dict[str, Any]] = []
    for block in data.get("by_superclass") or []:
        if not isinstance(block, dict) or not block.get("superclass"):
            continue
        receptors_block = block.get("receptors") or {}
        for receptor in RECEPTOR_GENES:
            weight, label, detail = _class_weight(receptor, receptors_block, lm, default_unknown, aggregation)
            if not detail:
                continue
            fractions = [d["fraction_cells_expressing"] for d in detail if d["fraction_cells_expressing"] is not None]
            tiers = {d["evidence_tier"] for d in detail if d["evidence_tier"]}
            rows.append(
                {
                    "superclass": block["superclass"],
                    "receptor": receptor,
                    "level": label if any(d["weight"] is not None for d in detail) else "unknown",
                    "weight": weight,
                    "known": any(d["weight"] is not None for d in detail),
                    "fraction_cells_expressing": float(np.mean(fractions)) if fractions else None,
                    "evidence_tier": "class_placeholder" if "class_placeholder" in tiers or not tiers
                    else sorted(tiers)[0],
                    "source_key": next((d["source_key"] for d in detail if d["source_key"]), None),
                    "mapping": block.get("mapping"),
                    "mapped_from": block.get("mapped_from"),
                    "genes": detail,
                }
            )
    return {
        "rows": rows,
        "level_map": lm,
        "default_unknown": float(default_unknown),
        "aggregation": aggregation,
        "published_atlas_values": data.get("published_atlas_values", []),
        "motor_neuron_gap": data.get("motor_neuron_gap"),
        "coverage": coverage(None, graphs, default_unknown, aggregation, path),
        "not_sourced": data.get("not_sourced", []),
        "warnings": [EXPRESSION_WARNING, MOTOR_NEURON_WARNING],
    }


# ---------------------------------------------------------------------------
# per-node weights
# ---------------------------------------------------------------------------
def expression_weights(
    receptor: str,
    graph: str | None = None,
    default_unknown: float = 1.0,
    level_map: Mapping[str, float] | None = None,
    aggregation: str = "mean",
    path: Path | None = None,
) -> dict[int, float]:
    """Per-node receptor multiplier, keyed by ``bodyId``.

    The documented default ordinal map is ``high 1.0, medium 0.6, low 0.25,
    absent 0.0``; a class the dataset leaves ``unknown`` - which is most of them
    - takes ``default_unknown`` (1.0, i.e. the uniform model).  Nodes with no
    MaleCNS superclass annotation also take ``default_unknown``.
    """
    from flylab.circuit.rate import load_graph

    if default_unknown < 0:
        raise ValueError("default_unknown must be >= 0")
    index = _superclass_index(receptor, level_map, default_unknown, aggregation, path)
    g = load_graph(graph)
    out: dict[int, float] = {}
    for node in g["nodes"]:
        sc = node.get("superclass")
        entry = index.get(str(sc)) if sc else None
        out[int(node["bodyId"])] = float(entry["weight"]) if entry else float(default_unknown)
    return out


def weighted_gains(
    rows: Iterable[Mapping[str, Any]],
    graph: str | None = None,
    gains: Mapping[str, float] | None = None,
    default_unknown: float = 1.0,
    level_map: Mapping[str, float] | None = None,
    aggregation: str = "mean",
    path: Path | None = None,
) -> dict[int, dict[str, float]]:
    """Per-node gain vectors instead of one scalar per transmitter.

    For every node, each transmitter gain is interpolated between 1.0 (no
    receptor at all, so the drug cannot act there) and the full scalar gain,
    using that node's expression weight::

        g_node = 1 + (g_full - 1) * w_node

    ``g_ach`` uses the nAChR weight, ``g_gaba`` the Rdl weight, ``g_glu`` the
    GluCl weight.  ``ach_tone`` is folded into the cholinergic gain exactly as
    :meth:`flylab.circuit.rate.RateNetwork.scale_vector` does.  With every
    weight at 1.0 the result is the uniform gain for every node.
    """
    g = dict(default_gains())
    g.update(dict(gains) if gains is not None else gains_from_occupancy(list(rows)))
    full = {
        "g_ach": g["g_ach"] * g.get("ach_tone", 1.0),
        "g_gaba": g["g_gaba"],
        "g_glu": g["g_glu"],
    }
    weights = {
        "g_ach": expression_weights("insect_nAChR", graph, default_unknown, level_map, aggregation, path),
        "g_gaba": expression_weights("insect_RDL", graph, default_unknown, level_map, aggregation, path),
        "g_glu": expression_weights("insect_GluCl", graph, default_unknown, level_map, aggregation, path),
    }
    body_ids = list(weights["g_ach"])
    out: dict[int, dict[str, float]] = {}
    for body_id in body_ids:
        out[int(body_id)] = {
            key: 1.0 + (full[key] - 1.0) * weights[key].get(body_id, default_unknown)
            for key in ("g_ach", "g_gaba", "g_glu")
        }
    return out


# ---------------------------------------------------------------------------
# weighted assay
# ---------------------------------------------------------------------------
def _weighted_rates(
    net,
    drive_vec: np.ndarray,
    node_gains: Mapping[int, Mapping[str, float]],
    g_nav: float = 1.0,
    steps: int = 80,
    alpha: float = 0.3,
    r_max: float = 300.0,
) -> np.ndarray:
    """Mirror of :meth:`flylab.circuit.rate.RateNetwork.run` with per-node gains.

    ``flylab/circuit/rate.py`` is deliberately **not** edited.  Its update rule
    is reproduced here verbatim::

        r <- clip((1 - alpha) * r + alpha * g_nav * (drive + Weff @ r), 0, r_max)

    The only difference is ``Weff``.  ``RateNetwork.run`` scales each column
    (presynaptic node) by one transmitter-wide gain; here each entry is scaled
    by the gain the **postsynaptic** node has for that presynaptic
    transmitter - which is where a receptor actually sits::

        Weff[j, i] = W_signed[j, i] * g_node[j][gain_key(nt(i))]

    With every expression weight at 1.0 the two coincide exactly, and the test
    suite asserts that.
    """
    from flylab.circuit.rate import NT_GAIN_KEY

    n = len(net.body_ids)
    gain_matrix = np.ones((n, n), dtype=float)
    col_key = [NT_GAIN_KEY.get(net.nt_of_node[b]) for b in net.body_ids]
    per_node = np.array(
        [[node_gains.get(b, {}).get(k, 1.0) for k in ("g_ach", "g_gaba", "g_glu")] for b in net.body_ids],
        dtype=float,
    )
    key_col = {"g_ach": 0, "g_gaba": 1, "g_glu": 2}
    for i, key in enumerate(col_key):
        if key in key_col:
            gain_matrix[:, i] = per_node[:, key_col[key]]
        # transmitters FlyLab has no receptor weight for keep 1.0
    weff = net.W_signed * gain_matrix
    r = np.zeros(n, dtype=float)
    for _ in range(int(steps)):
        r = np.clip((1.0 - alpha) * r + alpha * (g_nav * (drive_vec + weff @ r)), 0.0, r_max)
    return r


def run_weighted_subgraph_assay(
    compound: str | None = None,
    conc_M: float = 0.0,
    graph: str | None = None,
    drive_hz: float = 40.0,
    steps: int = 80,
    default_unknown: float = 1.0,
    level_map: Mapping[str, float] | None = None,
    aggregation: str = "mean",
    library: dict[str, Any] | None = None,
    path: Path | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """The neighborhood rate assay with per-node, expression-weighted gains.

    Returns a notebook (``assay == "malecns_neighborhood_expression_weighted"``)
    carrying the uniform-gain result, the weighted result and
    ``delta_vs_uniform``, plus the coverage summary and the standing warning
    that this is a sensitivity analysis.
    """
    from flylab.assays.subgraph import named_mean, named_readout, run_subgraph_assay
    from flylab.circuit.rate import DEFAULT_GAINS, compute_gains, load_graph, rate_network, resolve_graph
    from flylab.notebook.schema import empty_notebook, set_map

    net = rate_network(graph)
    gains, occ = compute_gains(compound, conc_M, library=library)
    node_gains = weighted_gains(
        occ["receptors"] if occ else [],
        graph=graph,
        gains=gains,
        default_unknown=default_unknown,
        level_map=level_map,
        aggregation=aggregation,
        path=path,
    )
    drive = net.drive_vector(net.seed_drive(drive_hz))
    r_weighted = _weighted_rates(net, drive, node_gains, g_nav=gains["g_nav"], steps=steps)
    r_vehicle = net.run(drive, DEFAULT_GAINS, steps=steps)

    uniform_nb = run_subgraph_assay(compound, conc_M, drive_hz=drive_hz, steps=steps, graph=graph, library=library)
    uniform = uniform_nb["readouts"]

    named = named_readout(net, r_weighted)
    weighted = {
        "mn9_hz": named_mean(named, "MN9"),
        "dnp01_hz": named_mean(named, "DNp01"),
        "mean_hz": float(r_weighted.mean()),
        "max_hz": float(r_weighted.max()),
    }
    delta = {
        k: (None if weighted[k] is None or uniform.get(k) is None else float(weighted[k]) - float(uniform[k]))
        for k in ("mn9_hz", "dnp01_hz", "mean_hz", "max_hz")
    }
    values = [abs(v) for v in delta.values() if v is not None]

    path_obj = resolve_graph(graph)
    g = load_graph(path_obj)
    nb = empty_notebook("malecns_neighborhood_expression_weighted")
    set_map(nb, g["map"], f"{path_obj.stem}:expression_weighted", g["citation"])
    nb["compound"] = compound
    nb["concentration_M"] = conc_M
    nb["occupancy"] = occ["receptors"] if occ else []
    nb["gains"] = dict(gains)
    weight_values = {
        key: [node_gains[b][key] for b in net.body_ids if b in node_gains]
        for key in ("g_ach", "g_gaba", "g_glu")
    }
    nb["readouts"] = {
        "engine": "rate_expression_weighted",
        "n_nodes": g["n_nodes"],
        "n_edges": g["n_edges"],
        "drive_hz": float(drive_hz),
        "steps": int(steps),
        "named": named,
        **weighted,
        "uniform": {k: uniform.get(k) for k in ("mn9_hz", "dnp01_hz", "mean_hz", "max_hz")},
        "vehicle": {"mean_hz": float(r_vehicle.mean())},
        "delta_vs_uniform": delta,
        "max_abs_delta_vs_uniform": float(max(values)) if values else 0.0,
        "node_gain_summary": {
            key: {
                "min": float(np.min(v)) if v else None,
                "mean": float(np.mean(v)) if v else None,
                "max": float(np.max(v)) if v else None,
            }
            for key, v in weight_values.items()
        },
    }
    nb["expression"] = {
        "level_map": dict(level_map or LEVEL_MAP),
        "default_unknown": float(default_unknown),
        "aggregation": aggregation,
        "coverage": coverage(None, [graph] if graph else None, default_unknown, aggregation, path),
        "node_gains_sample": {
            str(b): node_gains[b] for b in list(node_gains)[:10]
        },
    }
    nb["warnings"] = [
        EXPRESSION_WARNING,
        MOTOR_NEURON_WARNING,
        "The uniform-gain result in readouts['uniform'] remains the primary "
        "FlyLab model; delta_vs_uniform is the sensitivity, not a correction.",
        "Hops-limited MaleCNS neighborhood, not the full 25M-edge CNS.",
        "Signs from predicted transmitters. ACh / RDL / GluCl gains are teaching patches.",
    ]
    return nb
