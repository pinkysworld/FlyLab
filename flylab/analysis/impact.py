"""Where a compound acts in the neighborhood graph.

Everything here is a *structural* readout of the model: an edge's effective
weight is ``sign(nt_pre) * synapse_count * gain(nt_pre)``.  Gains come from
``flylab.pharm.mechanisms.gains_from_occupancy`` when available (see
:func:`flylab.circuit.rate.compute_gains`).  These are teaching patches on a
hops-limited map, not measured synaptic changes.
"""
from __future__ import annotations

from typing import Any, Iterable, Sequence

import numpy as np

from flylab.circuit.rate import (
    DEFAULT_GAINS,
    NT_GAIN_KEY,
    SIGN,
    compute_gains,
    load_graph,
    rate_network,
    resolve_graph,
    seed_ids,
)

__all__ = [
    "nt_gain",
    "edge_impact",
    "node_impact",
    "path_impact",
    "grn_to_mn9_paths",
    "summarize_impact",
    "totals_by_transmitter",
]

BITTER_TYPES = ("LB1a", "LB1b", "LB1c", "LB1d")
SWEET_TYPES = ("LB3b", "LB3c")


def _gain_map(gains: dict[str, float] | None) -> dict[str, float]:
    g = dict(DEFAULT_GAINS)
    g.update({k: float(v) for k, v in (gains or {}).items() if k in DEFAULT_GAINS})
    return g


def nt_gain(nt: str | None, gains: dict[str, float] | None) -> float:
    """Synaptic multiplier applied to edges from a presynaptic cell of ``nt``."""
    g = _gain_map(gains)
    key = NT_GAIN_KEY.get(nt or "")
    if key is None:
        return 1.0
    val = g[key]
    if nt == "acetylcholine":
        val *= g["ach_tone"]
    return float(val)


def _node_table(graph: dict[str, Any]) -> dict[int, dict[str, Any]]:
    return {n["bodyId"]: n for n in graph["nodes"]}


# --------------------------------------------------------------------------
def edge_impact(graph: dict[str, Any], gains: dict[str, float] | None = None) -> list[dict[str, Any]]:
    """Per-edge effective weight before and after the drug's gains.

    Sorted by ``|delta|`` (largest change first).
    """
    nodes = _node_table(graph)
    rows: list[dict[str, Any]] = []
    for e in graph["edges"]:
        pre, post = e["pre"], e["post"]
        n_pre = nodes.get(pre, {})
        nt = n_pre.get("consensus_nt") or "unclear"
        sign = SIGN.get(nt, 0.0)
        w = float(e["weight"])
        before = sign * w
        gain = nt_gain(nt, gains)
        after = before * gain
        rows.append(
            {
                "pre": pre,
                "post": post,
                "pre_type": n_pre.get("type"),
                "post_type": nodes.get(post, {}).get("type"),
                "nt": nt,
                "weight": w,
                "gain": float(gain),
                "eff_weight_before": float(before),
                "eff_weight_after": float(after),
                "delta": float(after - before),
            }
        )
    rows.sort(key=lambda r: -abs(r["delta"]))
    return rows


def totals_by_transmitter(edges: Sequence[dict[str, Any]]) -> dict[str, dict[str, float]]:
    """Summed |effective weight| before/after, per presynaptic transmitter."""
    acc: dict[str, dict[str, float]] = {}
    for r in edges:
        a = acc.setdefault(r["nt"], {"n_edges": 0.0, "synapses": 0.0, "before": 0.0, "after": 0.0})
        a["n_edges"] += 1
        a["synapses"] += r["weight"]
        a["before"] += abs(r["eff_weight_before"])
        a["after"] += abs(r["eff_weight_after"])
    for a in acc.values():
        a["delta"] = a["after"] - a["before"]
        a["ratio"] = (a["after"] / a["before"]) if a["before"] else 1.0
    return dict(sorted(acc.items()))


def _as_rate_map(graph: dict[str, Any], rates) -> dict[int, float]:
    if isinstance(rates, dict):
        return {int(k): float(v) for k, v in rates.items()}
    arr = np.asarray(rates, dtype=float)
    return {n["bodyId"]: float(arr[i]) for i, n in enumerate(graph["nodes"])}


def node_impact(
    graph: dict[str, Any],
    rates_vehicle,
    rates_treated,
    top: int | None = None,
) -> list[dict[str, Any]]:
    """Per-node vehicle vs treated rate, sorted by ``|delta_hz|``."""
    veh = _as_rate_map(graph, rates_vehicle)
    tre = _as_rate_map(graph, rates_treated)
    rows = []
    for n in graph["nodes"]:
        b = n["bodyId"]
        v, t = veh.get(b, 0.0), tre.get(b, 0.0)
        rows.append(
            {
                "bodyId": b,
                "type": n.get("type"),
                "superclass": n.get("superclass"),
                "nt": n.get("consensus_nt") or "unclear",
                "rate_vehicle": v,
                "rate_treated": t,
                "delta_hz": t - v,
                "ratio": (t / v) if v > 1e-9 else None,
            }
        )
    rows.sort(key=lambda r: -abs(r["delta_hz"]))
    return rows[: int(top)] if top else rows


# --------------------------------------------------------------------------
def _incoming(graph: dict[str, Any]) -> dict[int, list[dict[str, Any]]]:
    inc: dict[int, list[dict[str, Any]]] = {}
    for e in graph["edges"]:
        inc.setdefault(e["post"], []).append(e)
    return inc


def _backward_paths(
    graph: dict[str, Any],
    targets: Iterable[int],
    gains: dict[str, float] | None,
    max_len: int = 2,
    beam: int = 20000,
    sources: set[int] | None = None,
) -> list[dict[str, Any]]:
    """Enumerate input paths of 1..``max_len`` edges ending on ``targets``.

    Walks backwards from each target, keeping at most ``beam`` partial paths per
    depth ranked by ``|signed product of effective weights|``.  Cycles are not
    followed (a body may appear once per path).
    """
    nodes = _node_table(graph)
    inc = _incoming(graph)
    targets = [t for t in targets if t in nodes]
    out: list[dict[str, Any]] = []

    # frontier entries: (path bodyIds src..target, before, after, min_weight)
    frontier = [([t], 1.0, 1.0, float("inf")) for t in targets]
    for _ in range(int(max_len)):
        nxt = []
        for path, before, after, minw in frontier:
            head = path[0]
            for e in inc.get(head, []):
                pre = e["pre"]
                if pre in path:
                    continue
                n_pre = nodes.get(pre)
                if n_pre is None:
                    continue
                nt = n_pre.get("consensus_nt") or "unclear"
                sign = SIGN.get(nt, 0.0)
                if sign == 0.0:
                    continue
                w = float(e["weight"])
                b2 = before * sign * w
                a2 = after * sign * w * nt_gain(nt, gains)
                m2 = min(minw, w)
                new = ([pre] + path, b2, a2, m2)
                nxt.append(new)
                if sources is None or pre in sources:
                    out.append(
                        {
                            "path": [pre] + path,
                            "types": [nodes.get(x, {}).get("type") for x in [pre] + path],
                            "nts": [nodes.get(x, {}).get("consensus_nt") or "unclear" for x in [pre] + path],
                            "length": len(path),
                            "min_weight": m2,
                            "product_before": float(b2),
                            "product_after": float(a2),
                            "delta": float(a2 - b2),
                        }
                    )
        nxt.sort(key=lambda r: -abs(r[2]))
        frontier = nxt[: int(beam)]
        if not frontier:
            break
    return out


def path_impact(
    graph: dict[str, Any],
    seeds: Iterable[int] | None = None,
    gains: dict[str, float] | None = None,
    max_len: int = 2,
    top: int = 10,
) -> list[dict[str, Any]]:
    """Top input paths into the seed cells (MN9 / DNp01 by default).

    Ranked by ``|product of effective weights|`` after the drug.
    """
    targets = list(seeds) if seeds is not None else seed_ids(graph)
    rows = _backward_paths(graph, targets, gains, max_len=max_len)
    rows.sort(key=lambda r: -abs(r["product_after"]))
    return rows[: int(top)]


def grn_to_mn9_paths(
    graph: dict[str, Any],
    max_len: int = 3,
    top: int = 20,
    gains: dict[str, float] | None = None,
    targets: Iterable[int] | None = None,
) -> dict[str, Any]:
    """Top labellar-GRN -> ... -> MN9 paths, split by modality hypothesis.

    ``sweet`` = LB3b/c, ``bitter`` = LB1a-d (Cell 2026 nomenclature; this
    assignment is a working hypothesis recorded in
    ``data/derived/malecns_gustatory_seeds.json``, **not** a MaleCNS
    annotation).  Returns signed products of effective weights before and after
    the drug's gains, plus the minimum synapse count along the path.
    """
    seeds = graph.get("seeds", {})
    tgt = list(targets) if targets is not None else list(seeds.get("MN9", []))
    sweet = {b for t in SWEET_TYPES for b in seeds.get(t, [])}
    bitter = {b for t in BITTER_TYPES for b in seeds.get(t, [])}
    sources = sweet | bitter
    if not tgt or not sources:
        return {
            "sweet": [],
            "bitter": [],
            "n_sweet": 0,
            "n_bitter": 0,
            "note": "graph has no labellar GRN seeds; use the taste_motor graph",
        }
    rows = _backward_paths(graph, tgt, gains, max_len=max_len, sources=sources)
    out_sweet = [r for r in rows if r["path"][0] in sweet]
    out_bitter = [r for r in rows if r["path"][0] in bitter]
    for r in out_sweet + out_bitter:
        r["modality"] = "sweet" if r["path"][0] in sweet else "bitter"
    out_sweet.sort(key=lambda r: -abs(r["product_after"]))
    out_bitter.sort(key=lambda r: -abs(r["product_after"]))
    return {
        "sweet": out_sweet[: int(top)],
        "bitter": out_bitter[: int(top)],
        "n_sweet": len(out_sweet),
        "n_bitter": len(out_bitter),
        "max_len": int(max_len),
        "modality_hypothesis": {"sweet": list(SWEET_TYPES), "bitter": list(BITTER_TYPES)},
        "note": (
            "Sweet/bitter assignment of labellar GRN types is a working "
            "hypothesis (Cell 2026 nomenclature), not a MaleCNS annotation."
        ),
    }


# --------------------------------------------------------------------------
def summarize_impact(
    compound: str | None = None,
    conc_M: float = 0.0,
    graph: str | None = None,
    graph_path_arg=None,
    drive_hz: float = 40.0,
    steps: int = 80,
    top_edges: int = 25,
    top_nodes: int = 25,
    top_paths: int = 10,
    max_len: int = 2,
) -> dict[str, Any]:
    """JSON-ready impact summary for one compound at one concentration."""
    path = resolve_graph(graph if graph is not None else graph_path_arg)
    g = load_graph(path)
    net = rate_network(path)
    gains, _ = compute_gains(compound, conc_M)
    drive = net.drive_vector(net.seed_drive(drive_hz))
    r_veh = net.run(drive, DEFAULT_GAINS, steps=steps)
    r_tre = net.run(drive, gains, steps=steps)

    edges = edge_impact(g, gains)
    return {
        "compound": compound,
        "concentration_M": conc_M,
        "graph": path.name,
        "map": g.get("map"),
        "gains": dict(gains),
        "n_edges": len(edges),
        "top_edges": edges[: int(top_edges)],
        "top_nodes": node_impact(g, r_veh, r_tre, top=top_nodes),
        "top_paths": path_impact(g, None, gains, max_len=max_len, top=top_paths),
        "totals_by_transmitter": totals_by_transmitter(edges),
        "label": "model_derived",
        "warnings": [
            "Hops-limited MaleCNS neighborhood, not the full CNS.",
            "Effective weights are teaching gain patches on predicted transmitters.",
        ],
    }
