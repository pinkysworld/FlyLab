"""Rate assay on the committed MaleCNS named-cell neighborhood.

Public API is unchanged; the network math now lives in
:mod:`flylab.circuit.rate` so the spiking assay, ensemble layer and analysis
package share one runtime.  ``run_subgraph_assay`` returns the historical
readouts bit-for-bit, plus ``gains``, ``readouts["by_superclass"]`` and
``readouts["top_changed"]``.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from flylab.circuit.rate import (
    CANDIDATES,
    DEFAULT_GAINS,
    SIGN,
    by_superclass,
    compute_gains,
    graph_path,
    load_graph,
    rate_network,
    top_changed,
)
from flylab.notebook.schema import empty_notebook

WARNINGS = [
    "Hops-limited MaleCNS neighborhood, not the full 25M-edge CNS.",
    "Signs from predicted transmitters. ACh gain and RDL gain are teaching patches.",
]


def _gains(compound, conc_M, library=None, rule_overrides=None):
    """Backward-compatible helper: ``(g_ach, g_gaba, occupancy_block)``."""
    gains, occ = compute_gains(compound, conc_M, library=library, rule_overrides=rule_overrides)
    return gains["g_ach"], gains["g_gaba"], occ


def named_readout(net, rates) -> dict[str, list[dict[str, Any]]]:
    named: dict[str, list[dict[str, Any]]] = {}
    for typ, ids in net.graph["seeds"].items():
        named[typ] = [
            {"bodyId": nid, "hz": float(rates[net.node_index[nid]]) if nid in net.node_index else None}
            for nid in ids
        ]
    return named


def named_mean(named: dict[str, list[dict[str, Any]]], typ: str) -> float | None:
    hz = [row["hz"] for row in named.get(typ, []) if row["hz"] is not None]
    return float(np.mean(hz)) if hz else None


def run_subgraph_assay(
    compound: str | None = None,
    conc_M: float = 0.0,
    graph_path_arg=None,
    drive_hz: float = 40.0,
    steps: int = 80,
    library: dict[str, Any] | None = None,
    rule_overrides: dict[str, float] | None = None,
    drive: dict[int, float] | None = None,
) -> dict[str, Any]:
    """Run the neighborhood rate model and return a notebook (schema 0.2/0.3).

    ``library`` / ``rule_overrides`` / ``drive`` are optional hooks used by the
    ensemble and sensitivity layers; the defaults reproduce the v0.4 numbers.
    """
    net = rate_network(graph_path_arg)
    g = net.graph
    gains, occ = compute_gains(compound, conc_M, library=library, rule_overrides=rule_overrides)
    drive_map = dict(drive) if drive else net.seed_drive(drive_hz)
    drive_vec = net.drive_vector(drive_map)

    r = net.run(drive_vec, gains, steps=steps)
    r_vehicle = net.run(drive_vec, DEFAULT_GAINS, steps=steps)

    named = named_readout(net, r)
    nb = empty_notebook("malecns_neighborhood")
    nb["map"] = {"name": g["map"], "version": "neighborhood", "citation": g["citation"]}
    nb["compound"] = compound
    nb["concentration_M"] = conc_M
    nb["occupancy"] = occ["receptors"] if occ else []
    nb.setdefault("gains", {})
    nb.setdefault("uncertainty", None)
    nb["gains"] = dict(gains)
    nb["readouts"] = {
        "n_nodes": g["n_nodes"],
        "n_edges": g["n_edges"],
        "named": named,
        "mean_hz": float(r.mean()),
        "max_hz": float(r.max()),
        "g_ach": float(gains["g_ach"]),
        "g_gaba": float(gains["g_gaba"]),
        "mn9_hz": named_mean(named, "MN9"),
        "dnp01_hz": named_mean(named, "DNp01"),
        "drive_hz": float(drive_hz),
        "steps": int(steps),
        "by_superclass": by_superclass(net, r),
        "top_changed": top_changed(net, r_vehicle, r, k=15),
        "vehicle": {
            "mean_hz": float(r_vehicle.mean()),
            "mn9_hz": named_mean(named_readout(net, r_vehicle), "MN9"),
            "dnp01_hz": named_mean(named_readout(net, r_vehicle), "DNp01"),
        },
    }
    nb["warnings"] = list(WARNINGS)
    return nb


def dose_response_subgraph(compound: str, concs=None):
    concs = concs or [10**e for e in range(-10, -4)]
    points = []
    for c in concs:
        nb = run_subgraph_assay(compound, c)
        mn9 = nb["readouts"]["named"].get("MN9", [])
        hz = [row["hz"] for row in mn9 if row["hz"] is not None]
        points.append(
            {
                "conc_M": c,
                "mean_hz": nb["readouts"]["mean_hz"],
                "mn9_hz": float(np.mean(hz)) if hz else None,
                "g_ach": nb["readouts"]["g_ach"],
                "g_gaba": nb["readouts"]["g_gaba"],
            }
        )
    return points


__all__ = [
    "SIGN",
    "CANDIDATES",
    "graph_path",
    "load_graph",
    "run_subgraph_assay",
    "dose_response_subgraph",
    "named_readout",
    "named_mean",
]
