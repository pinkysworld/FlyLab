"""Map-extracted labellar-GRN -> MN9 assay on the taste_motor neighborhood.

This drives the *real* MaleCNS labellar gustatory receptor neurons (LB1a-d,
LB3b/c body IDs) and reads MN9, instead of driving MN9 directly.  It is the map
counterpart of ``flylab/assays/taste.py`` (``reduced_taste_v0``).

Two things must stay explicit in every notebook this writes:

1. The sweet/bitter assignment of the labellar types is a **working hypothesis**
   (Cell 2026 nomenclature, recorded in
   ``data/derived/malecns_gustatory_seeds.json``).  MaleCNS does not annotate
   modality.
2. ``reduced_taste_v0`` remains the directional control for the bitter veto.
   This assay reports whatever the extracted path does; it does not assert that
   the map path vetoes MN9.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from flylab.circuit.lif import LIFNetwork, lif_network
from flylab.circuit.rate import (
    DEFAULT_GAINS,
    RateNetwork,
    compute_gains,
    load_graph,
    rate_network,
    resolve_graph,
)
from flylab.notebook.schema import empty_notebook

try:  # notebook schema 0.3
    from flylab.notebook.schema import set_map
except ImportError:  # pragma: no cover - schema 0.2 fallback

    def set_map(nb, name, version, citation):
        nb["map"] = {"name": name, "version": version, "citation": citation}
        return nb


ASSAY = "malecns_taste_map"
SWEET_TYPES = ("LB3b", "LB3c")
BITTER_TYPES = ("LB1a", "LB1b", "LB1c", "LB1d")
MOTOR_TYPES = ("MN9", "DNp01")

WARNINGS = [
    "Sweet/bitter assignment of the labellar GRN types (LB3b/c sweet, LB1a-d "
    "bitter) is a working hypothesis from the Cell 2026 nomenclature, not a "
    "MaleCNS annotation.",
    "reduced_taste_v0 remains the directional control for the bitter veto; this "
    "assay reports what the extracted path does and asserts nothing about it.",
    "Hops-limited MaleCNS neighborhood with a weight floor, not the full "
    "25M-edge CNS; missing inhibition can change the sign of any readout.",
    "ACh / RDL / GluCl gains are teaching patches and no rate here is a "
    "measurement from a living fly.",
]


def _type_ids(graph: dict[str, Any], types) -> list[int]:
    seeds = graph.get("seeds", {})
    out: list[int] = []
    for t in types:
        out.extend(seeds.get(t, []))
    return out


def _drive_map(graph: dict[str, Any], sugar_hz: float, bitter_hz: float) -> dict[int, float]:
    """Drive sweet and bitter GRNs only - never MN9 / DNp01."""
    drive: dict[int, float] = {}
    for b in _type_ids(graph, SWEET_TYPES):
        drive[b] = float(sugar_hz)
    for b in _type_ids(graph, BITTER_TYPES):
        drive[b] = float(bitter_hz)
    return drive


def _mn9(graph: dict[str, Any], index: dict[int, int], rates) -> float | None:
    hz = [float(rates[index[b]]) for b in graph.get("seeds", {}).get("MN9", []) if b in index]
    return float(np.mean(hz)) if hz else None


def run_taste_map_assay(
    compound: str | None = None,
    conc_M: float = 0.0,
    sugar_hz: float = 150.0,
    bitter_hz: float = 0.0,
    engine: str = "rate",
    seed: int = 0,
    graph: str | None = "taste_motor",
    steps: int = 80,
    t_ms: float = 400.0,
    library: dict[str, Any] | None = None,
    rule_overrides: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Sugar / sugar+bitter / vehicle / no-drive MN9 rates on the map path.

    ``engine`` is ``"rate"`` (leaky rate network) or ``"lif"`` (Shiu-style
    spiking).  ``bitter_hz`` of 0 means the bitter probe condition reuses
    ``sugar_hz`` as the bitter drive, matching ``assays/taste.py``.
    """
    path = resolve_graph(graph)
    g = load_graph(path)
    gains, occ = compute_gains(compound, conc_M, library=library, rule_overrides=rule_overrides)
    bitter_probe_hz = float(bitter_hz) if bitter_hz and bitter_hz > 0 else float(sugar_hz)

    if engine not in ("rate", "lif"):
        raise ValueError("engine must be 'rate' or 'lif'")

    if engine == "rate":
        net: RateNetwork | LIFNetwork = rate_network(path)
        index = net.node_index

        def run(drive, gain):
            return net.run(net.drive_vector(drive), gain, steps=steps)

    else:
        net = lif_network(g, seed=seed)
        net.seed = int(seed)
        index = net.node_index

        def run(drive, gain):
            return net.run(drive, gain, t_ms=t_ms).rates_hz

    d_sugar = _drive_map(g, sugar_hz, 0.0)
    d_both = _drive_map(g, sugar_hz, bitter_probe_hz)
    d_none = _drive_map(g, 0.0, 0.0)

    r_sugar = run(d_sugar, gains)
    r_both = run(d_both, gains)
    r_veh_sugar = run(d_sugar, DEFAULT_GAINS)
    r_none = run(d_none, DEFAULT_GAINS)

    mn9_sugar = _mn9(g, index, r_sugar)
    mn9_both = _mn9(g, index, r_both)
    mn9_veh = _mn9(g, index, r_veh_sugar)
    mn9_none = _mn9(g, index, r_none)
    veto = (mn9_both / mn9_sugar) if (mn9_sugar and mn9_sugar > 1e-6) else None

    seed_set = {b for ids in g.get("seeds", {}).values() for b in ids}
    body_ids = net.body_ids
    delta = np.asarray(r_sugar, dtype=float) - np.asarray(r_none, dtype=float)
    relay_order = [int(i) for i in np.argsort(-np.abs(delta)) if body_ids[int(i)] not in seed_set]
    relays = []
    for i in relay_order[:15]:
        b = body_ids[i]
        relays.append(
            {
                "bodyId": b,
                "type": net.types[i],
                "superclass": net.superclasses[i],
                "nt": net.nt_of_node[b],
                "hz_sugar": float(r_sugar[i]),
                "hz_sugar_bitter": float(r_both[i]),
                "hz_no_drive": float(r_none[i]),
                "delta_hz": float(delta[i]),
            }
        )

    nb = empty_notebook(ASSAY, seed=seed) if engine == "lif" else empty_notebook(ASSAY)
    set_map(nb, g["map"], "taste_motor_neighborhood", g["citation"])
    nb["compound"] = compound
    nb["concentration_M"] = conc_M
    nb["occupancy"] = occ["receptors"] if occ else []
    nb.setdefault("gains", {})
    nb.setdefault("uncertainty", None)
    nb["gains"] = dict(gains)
    nb["readouts"] = {
        "engine": engine,
        "n_nodes": g["n_nodes"],
        "n_edges": g["n_edges"],
        "sugar_hz": float(sugar_hz),
        "bitter_hz": float(bitter_probe_hz),
        "seed": int(seed),
        "sweet_types": list(SWEET_TYPES),
        "bitter_types": list(BITTER_TYPES),
        "n_sweet_grn": len(_type_ids(g, SWEET_TYPES)),
        "n_bitter_grn": len(_type_ids(g, BITTER_TYPES)),
        "mn9_hz": mn9_sugar,
        "mn9_sugar_hz": mn9_sugar,
        "mn9_sugar_bitter_hz": mn9_both,
        "mn9_vehicle_sugar_hz": mn9_veh,
        "mn9_no_drive_hz": mn9_none,
        "bitter_veto_ratio": veto,
        "g_ach": float(gains["g_ach"]),
        "g_gaba": float(gains["g_gaba"]),
        "mean_hz": float(np.mean(r_sugar)),
        "top_relays": relays,
    }
    nb["warnings"] = list(WARNINGS)
    return nb


__all__ = ["run_taste_map_assay", "ASSAY", "SWEET_TYPES", "BITTER_TYPES"]
