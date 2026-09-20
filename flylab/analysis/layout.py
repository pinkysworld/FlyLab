"""Deterministic 2-D layout of a neighborhood graph for the bench viewer.

No force-directed step, no RNG that matters: the same graph always produces the
same positions, so the UI can diff two runs node-for-node.

Layering
--------
``named`` graph (MN9 / DNp01 seeds in the middle):

* column  0  seeds
* column -1  cells presynaptic to a seed
* column +1  cells postsynaptic to a seed
* column -2 / +2  anything further out

``taste_motor`` graph (labellar GRNs -> ... -> MN9):

* column  0  labellar GRN seeds (LB1a-d, LB3b/c) on the left
* column  N  MN9 / DNp01 on the right
* columns in between: relays, by hop distance from the GRNs

Within a column, cells are ordered by transmitter then by descending degree,
then by bodyId, and spread evenly on the y axis.
"""
from __future__ import annotations

from typing import Any, Iterable

from flylab.analysis.impact import BITTER_TYPES, SWEET_TYPES, nt_gain
from flylab.circuit.rate import SIGN, seed_ids

__all__ = ["graph_layout", "graph_for_viewer", "VIEWER_MIN_WEIGHT"]

#: edges below this synapse count are not sent to the viewer
VIEWER_MIN_WEIGHT = 5.0

COL_X = 320.0
ROW_Y = 26.0

NT_ORDER = {
    "acetylcholine": 0,
    "glutamate": 1,
    "gaba": 2,
    "octopamine": 3,
    "dopamine": 4,
    "serotonin": 5,
    "histamine": 6,
    "unclear": 7,
}


def _grn_seed_ids(graph: dict[str, Any]) -> list[int]:
    seeds = graph.get("seeds", {})
    out: list[int] = []
    for t in SWEET_TYPES + BITTER_TYPES:
        out.extend(seeds.get(t, []))
    return out


def _motor_seed_ids(graph: dict[str, Any]) -> list[int]:
    seeds = graph.get("seeds", {})
    return list(seeds.get("MN9", [])) + list(seeds.get("DNp01", []))


def _layers_named(graph: dict[str, Any]) -> dict[int, int]:
    seeds = set(seed_ids(graph))
    layer: dict[int, int] = {b: 0 for b in seeds}
    for e in graph["edges"]:
        pre, post = e["pre"], e["post"]
        if post in seeds and pre not in seeds:
            layer.setdefault(pre, -1)
    for e in graph["edges"]:
        pre, post = e["pre"], e["post"]
        if pre in seeds and post not in seeds and post not in layer:
            layer[post] = 1
    for n in graph["nodes"]:
        b = n["bodyId"]
        if b not in layer:
            layer[b] = -2
    return layer


def _layers_taste(graph: dict[str, Any]) -> dict[int, int]:
    """Forward hop distance from the GRN seeds; MN9/DNp01 pinned to the right."""
    grns = _grn_seed_ids(graph)
    motors = set(_motor_seed_ids(graph))
    out_adj: dict[int, list[int]] = {}
    for e in graph["edges"]:
        out_adj.setdefault(e["pre"], []).append(e["post"])
    layer: dict[int, int] = {b: 0 for b in grns}
    frontier = list(grns)
    depth = 0
    max_depth = 4
    while frontier and depth < max_depth:
        depth += 1
        nxt = []
        for b in frontier:
            for post in out_adj.get(b, []):
                if post not in layer:
                    layer[post] = depth
                    nxt.append(post)
        frontier = nxt
    right = max_depth
    for b in motors:
        layer[b] = right
    for n in graph["nodes"]:
        b = n["bodyId"]
        if b not in layer:
            layer[b] = right - 1
        elif layer[b] >= right and b not in motors:
            layer[b] = right - 1
    return layer


def _degree(graph: dict[str, Any]) -> dict[int, float]:
    deg: dict[int, float] = {}
    for e in graph["edges"]:
        w = float(e["weight"])
        deg[e["pre"]] = deg.get(e["pre"], 0.0) + w
        deg[e["post"]] = deg.get(e["post"], 0.0) + w
    return deg


def graph_layout(graph: dict[str, Any], seed: int = 0) -> dict[int, list[float]]:
    """Deterministic ``{bodyId: [x, y]}`` positions.

    ``seed`` is accepted for API symmetry; the layout has no random component,
    so the result does not depend on it.
    """
    has_grn = bool(_grn_seed_ids(graph))
    layer = _layers_taste(graph) if has_grn else _layers_named(graph)
    deg = _degree(graph)
    nt_of = {n["bodyId"]: (n.get("consensus_nt") or "unclear") for n in graph["nodes"]}

    cols: dict[int, list[int]] = {}
    for n in graph["nodes"]:
        b = n["bodyId"]
        cols.setdefault(layer.get(b, 0), []).append(b)

    pos: dict[int, list[float]] = {}
    for col, members in sorted(cols.items()):
        members.sort(key=lambda b: (NT_ORDER.get(nt_of.get(b, "unclear"), 9), -deg.get(b, 0.0), b))
        k = len(members)
        y0 = -(k - 1) * ROW_Y / 2.0
        for i, b in enumerate(members):
            pos[b] = [float(col * COL_X), float(y0 + i * ROW_Y)]
    return pos


def graph_for_viewer(
    graph: dict[str, Any],
    rates: Iterable[float] | dict[int, float] | None = None,
    gains: dict[str, float] | None = None,
    min_weight: float = VIEWER_MIN_WEIGHT,
    max_edges: int = 6000,
    seed: int = 0,
) -> dict[str, Any]:
    """Compact payload for the Circuit panel.

    ``rates`` may be a per-node array (graph node order) or a bodyId map.
    ``gains`` adds ``eff_weight`` per edge.  Edges below ``min_weight`` synapses
    are dropped and, on the dense taste_motor graph, only the ``max_edges``
    strongest survive so the payload stays small enough to draw.
    """
    pos = graph_layout(graph, seed=seed)
    seeds = set(seed_ids(graph))
    rate_map: dict[int, float] = {}
    if rates is not None:
        if isinstance(rates, dict):
            rate_map = {int(k): float(v) for k, v in rates.items()}
        else:
            vals = list(rates)
            rate_map = {n["bodyId"]: float(vals[i]) for i, n in enumerate(graph["nodes"]) if i < len(vals)}

    nodes = []
    for n in graph["nodes"]:
        b = n["bodyId"]
        x, y = pos.get(b, [0.0, 0.0])
        nodes.append(
            {
                "id": b,
                "type": n.get("type"),
                "superclass": n.get("superclass"),
                "nt": n.get("consensus_nt") or "unclear",
                "x": x,
                "y": y,
                "rate": rate_map.get(b),
                "is_seed": b in seeds,
            }
        )

    nt_of = {n["bodyId"]: (n.get("consensus_nt") or "unclear") for n in graph["nodes"]}
    edges = []
    for e in graph["edges"]:
        w = float(e["weight"])
        if w < float(min_weight):
            continue
        nt = nt_of.get(e["pre"], "unclear")
        base = SIGN.get(nt, 0.0) * w
        edges.append(
            {
                "source": e["pre"],
                "target": e["post"],
                "weight": w,
                "nt": nt,
                "eff_weight": float(base * nt_gain(nt, gains)),
            }
        )

    truncated = False
    if max_edges and len(edges) > int(max_edges):
        edges.sort(key=lambda r: -r["weight"])
        edges = edges[: int(max_edges)]
        truncated = True

    return {
        "truncated_edges": truncated,
        "map": graph.get("map"),
        "citation": graph.get("citation"),
        "n_nodes": len(nodes),
        "n_edges": len(edges),
        "min_weight": float(min_weight),
        "seeds": graph.get("seeds", {}),
        "gains": dict(gains) if gains else None,
        "nodes": nodes,
        "edges": edges,
    }
