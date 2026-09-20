"""LIF spiking assay on a committed MaleCNS neighborhood graph."""
from __future__ import annotations

from typing import Any

import numpy as np

from flylab.circuit.lif import (
    BACKGROUND_FRACTION_DEFAULT,
    DRIVE_FANIN_DEFAULT,
    DRIVE_W_DEFAULT,
    LIFNetwork,
    W_SCALE_DEFAULT,
    lif_network,
)
from flylab.circuit.rate import DEFAULT_GAINS, compute_gains, load_graph, resolve_graph
from flylab.notebook.schema import empty_notebook

try:  # notebook schema 0.3
    from flylab.notebook.schema import set_map
except ImportError:  # pragma: no cover - schema 0.2 fallback

    def set_map(nb, name, version, citation):
        nb["map"] = {"name": name, "version": version, "citation": citation}
        return nb


ASSAY = "malecns_neighborhood_lif"

WARNINGS = [
    "Hops-limited MaleCNS neighborhood, not the full 25M-edge CNS.",
    "Signs from predicted transmitters. ACh / RDL / GluCl gains are teaching patches.",
    "LIF parameters are Shiu-style literature defaults; the weight scale is a "
    "single hand calibration for a plausible vehicle MN9 rate. Nothing here is "
    "fitted to animal data and no spike rate is a measurement from a living fly.",
    "Every cell receives a uniform background Poisson drive standing in for the "
    "inputs lost when the neighborhood was cut out of the CNS.",
]

#: keep the notebook JSON small
MAX_SEED_RASTER_NODES = 80
MAX_EXTRA_RASTER_NODES = 30
MAX_SPIKES_PER_NODE = 800
MAX_RASTER_SPIKES = 15000


def _default_drive(graph: dict[str, Any], drive_hz: float, background_hz: float) -> dict[int, float]:
    drive = {n["bodyId"]: float(background_hz) for n in graph["nodes"]}
    for ids in graph.get("seeds", {}).values():
        for body_id in ids:
            drive[body_id] = float(drive_hz)
    return drive


def _named(net: LIFNetwork, rates: np.ndarray) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    for typ, ids in net.graph.get("seeds", {}).items():
        out[typ] = [
            {"bodyId": b, "hz": float(rates[net.node_index[b]]) if b in net.node_index else None}
            for b in ids
        ]
    return out


def _named_mean(named: dict[str, list[dict[str, Any]]], typ: str) -> float | None:
    hz = [r["hz"] for r in named.get(typ, []) if r["hz"] is not None]
    return float(np.mean(hz)) if hz else None


def _by_superclass(net: LIFNetwork, rates: np.ndarray) -> dict[str, float]:
    acc: dict[str, list[float]] = {}
    for i, sc in enumerate(net.superclasses):
        acc.setdefault(sc or "unknown", []).append(float(rates[i]))
    return {k: float(np.mean(v)) for k, v in sorted(acc.items())}


def _spikes_by_index(result) -> dict[int, list[float]]:
    out: dict[int, list[float]] = {}
    for t, i in result.spikes:
        out.setdefault(i, []).append(round(float(t), 1))
    return out


def _raster(net: LIFNetwork, result, rates: np.ndarray) -> list[dict[str, Any]]:
    by_idx = _spikes_by_index(result)
    seed_idx = [net.node_index[b] for b in net.seed_ids if b in net.node_index]
    order = [int(i) for i in np.argsort(-rates)]
    picks: list[int] = []
    for i in seed_idx[:MAX_SEED_RASTER_NODES]:
        if i not in picks:
            picks.append(i)
    limit = len(picks) + MAX_EXTRA_RASTER_NODES
    for i in order:
        if len(picks) >= limit:
            break
        if i not in picks and rates[i] > 0:
            picks.append(i)
    rows: list[dict[str, Any]] = []
    budget = MAX_RASTER_SPIKES
    for i in picks:
        st = by_idx.get(i, [])[:MAX_SPIKES_PER_NODE]
        if len(st) > budget:
            st = st[: max(budget, 0)]
        budget -= len(st)
        body_id = net.body_ids[i]
        rows.append(
            {
                "bodyId": body_id,
                "type": net.types[i],
                "nt": net.nt_of_node[body_id],
                "is_seed": body_id in net.seed_ids,
                "hz": float(rates[i]),
                "spikes_ms": st,
            }
        )
        if budget <= 0:
            break
    return rows


def _psth(net: LIFNetwork, result, bin_ms: float = 10.0) -> dict[str, Any]:
    n_bins = max(int(round(result.t_ms / bin_ms)), 1)
    by_idx = _spikes_by_index(result)
    cells: dict[str, list[dict[str, Any]]] = {}
    for typ, ids in net.graph.get("seeds", {}).items():
        rows = []
        for b in ids:
            i = net.node_index.get(b)
            if i is None:
                continue
            counts = [0] * n_bins
            for t in by_idx.get(i, []):
                k = min(int(t // bin_ms), n_bins - 1)
                counts[k] += 1
            rows.append({"bodyId": b, "counts": counts})
        cells[typ] = rows
    return {"bin_ms": bin_ms, "n_bins": n_bins, "cells": cells}


def _notebook(assay: str, seed: int):
    try:
        return empty_notebook(assay, seed=seed)
    except TypeError:  # pragma: no cover - schema 0.2
        return empty_notebook(assay)


def run_spiking_assay(
    compound: str | None = None,
    conc_M: float = 0.0,
    graph_path_arg=None,
    drive_hz: float = 40.0,
    t_ms: float = 500.0,
    seed: int = 0,
    drive: dict[int, float] | None = None,
    graph: str | None = None,
    background_hz: float | None = None,
    library: dict[str, Any] | None = None,
    rule_overrides: dict[str, float] | None = None,
    w_scale: float = W_SCALE_DEFAULT,
    drive_w: float = DRIVE_W_DEFAULT,
    drive_fanin: int = DRIVE_FANIN_DEFAULT,
    graph_obj: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run the Shiu-style LIF network and return a notebook.

    ``drive`` overrides the default drive map (seeds at ``drive_hz``, every
    other cell at ``background_hz``).  Deterministic for a given ``seed``.
    """
    if graph_obj is not None:
        path = resolve_graph(graph if graph is not None else graph_path_arg)
        g = graph_obj
        net = LIFNetwork(g, seed=seed, w_scale=w_scale, drive_w=drive_w, drive_fanin=drive_fanin)
    else:
        path = resolve_graph(graph if graph is not None else graph_path_arg)
        g = load_graph(path)
        net = lif_network(g, seed=seed, w_scale=w_scale, drive_w=drive_w, drive_fanin=drive_fanin)
    if background_hz is None:
        background_hz = BACKGROUND_FRACTION_DEFAULT * float(drive_hz)
    drive_map = dict(drive) if drive else _default_drive(g, drive_hz, background_hz)

    gains, occ = compute_gains(compound, conc_M, library=library, rule_overrides=rule_overrides)
    result = net.run(drive_map, gains, t_ms=t_ms)
    vehicle = net.run(drive_map, DEFAULT_GAINS, t_ms=t_ms)

    rates = result.rates_hz
    named = _named(net, rates)
    named_veh = _named(net, vehicle.rates_hz)

    nb = _notebook(ASSAY, seed)
    set_map(nb, g["map"], f"{path.stem}:lif", g["citation"])
    nb["compound"] = compound
    nb["concentration_M"] = conc_M
    nb["occupancy"] = occ["receptors"] if occ else []
    nb.setdefault("gains", {})
    nb.setdefault("uncertainty", None)
    nb["gains"] = dict(gains)
    nb["readouts"] = {
        "engine": "lif",
        "n_nodes": len(net),
        "n_edges": g["n_edges"],
        "t_ms": float(result.t_ms),
        "analysis_window_ms": list(result.window_ms),
        "drive_hz": float(drive_hz),
        "background_hz": float(background_hz),
        "seed": int(seed),
        "named": named,
        "mn9_hz": _named_mean(named, "MN9"),
        "dnp01_hz": _named_mean(named, "DNp01"),
        "mean_hz": float(rates.mean()),
        "max_hz": float(rates.max()),
        "frac_active": float((rates > 0).mean()),
        "n_spikes": len(result.spikes),
        "g_ach": float(gains["g_ach"]),
        "g_gaba": float(gains["g_gaba"]),
        "by_superclass": _by_superclass(net, rates),
        "raster": _raster(net, result, rates),
        "psth": _psth(net, result),
        "vehicle": {
            "mn9_hz": _named_mean(named_veh, "MN9"),
            "dnp01_hz": _named_mean(named_veh, "DNp01"),
            "mean_hz": float(vehicle.rates_hz.mean()),
        },
        "lif_params": {
            "tau_m_ms": 20.0,
            "v_rest_mV": -52.0,
            "v_th_mV": -45.0,
            "v_reset_mV": -52.0,
            "t_ref_ms": 2.2,
            "tau_s_ms": 5.0,
            "dt_ms": net.dt_ms,
            "w_scale_mV_per_synapse": net.w_scale,
            "drive_w_mV": net.drive_w,
            "drive_fanin": net.drive_fanin,
        },
    }
    nb["warnings"] = list(WARNINGS)
    return nb


def dose_response_spiking(compound: str, concs=None, **kw):
    concs = concs or [10**e for e in range(-10, -4)]
    rows = []
    for c in concs:
        nb = run_spiking_assay(compound, c, **kw)
        r = nb["readouts"]
        rows.append(
            {
                "conc_M": c,
                "mn9_hz": r["mn9_hz"],
                "dnp01_hz": r["dnp01_hz"],
                "mean_hz": r["mean_hz"],
                "g_ach": r["g_ach"],
                "g_gaba": r["g_gaba"],
            }
        )
    return rows


__all__ = ["run_spiking_assay", "dose_response_spiking", "ASSAY"]
