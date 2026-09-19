from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from flylab.notebook.schema import empty_notebook
from flylab.pharm.occupancy import compare_compound

SIGN = {
    "acetylcholine": 1.0,
    "glutamate": -0.4,
    "gaba": -1.0,
    "histamine": -0.5,
    "dopamine": 0.2,
    "serotonin": 0.2,
    "octopamine": 0.2,
}
CANDIDATES = [
    Path("data/derived/malecns_named_neighborhood.json"),
    Path(__file__).resolve().parents[2] / "data/derived/malecns_named_neighborhood.json",
    Path.home() / ".flylab/malecns_named_neighborhood.json",
]


def graph_path(path=None) -> Path:
    if path:
        p = Path(path)
        if p.exists():
            return p
    for p in CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError("neighborhood JSON missing. Run extract-malecns-subgraph Action.")


def load_graph(path=None):
    return json.loads(graph_path(path).read_text())


def _gains(compound, conc_M):
    if not compound:
        return 1.0, 1.0, None
    occ = compare_compound(compound, conc_M)
    g_ach, g_gaba = 1.0, 1.0
    for row in occ["receptors"]:
        th, d = row["occupancy"], row["direction"]
        if row["receptor"] == "insect_nAChR":
            if d == "agonist":
                g_ach = max(0.05, 1.0 + 0.4 * th - 1.6 * th * th)
            elif d == "antagonist":
                g_ach = max(0.05, 1.0 - th)
        if row["receptor"] == "insect_RDL":
            if d == "antagonist":
                g_gaba = max(0.05, 1.0 - th)
            elif d == "agonist":
                g_gaba = max(0.05, 1.0 + 0.4 * th)
    return g_ach, g_gaba, occ


def run_subgraph_assay(compound=None, conc_M=0.0, graph_path_arg=None, drive_hz=40.0, steps=80):
    g = load_graph(graph_path_arg)
    nodes = g["nodes"]
    index = {n["bodyId"]: i for i, n in enumerate(nodes)}
    n = len(nodes)
    W = np.zeros((n, n), dtype=float)
    for e in g["edges"]:
        i, j = index.get(e["pre"]), index.get(e["post"])
        if i is None or j is None:
            continue
        nt = nodes[i].get("consensus_nt") or "unclear"
        W[j, i] += SIGN.get(nt, 0.0) * float(e["weight"])
    denom = np.maximum(np.abs(W).sum(axis=1, keepdims=True), 1.0)
    W = W / denom
    g_ach, g_gaba, occ = _gains(compound, conc_M)
    scale = []
    for nd in nodes:
        nt = nd.get("consensus_nt")
        if nt == "acetylcholine":
            scale.append(g_ach)
        elif nt == "gaba":
            scale.append(g_gaba)
        else:
            scale.append(1.0)
    scale = np.array(scale)
    r = np.zeros(n)
    seed_ids = {i for ids in g["seeds"].values() for i in ids}
    drive = np.array([drive_hz if nd["bodyId"] in seed_ids else 0.0 for nd in nodes])
    for _ in range(steps):
        r = np.clip(0.7 * r + 0.3 * (drive + (W * scale) @ r), 0.0, 300.0)
    named = {}
    for typ, ids in g["seeds"].items():
        named[typ] = [{"bodyId": nid, "hz": float(r[index[nid]]) if nid in index else None} for nid in ids]
    nb = empty_notebook("malecns_neighborhood")
    nb["map"] = {"name": g["map"], "version": "neighborhood", "citation": g["citation"]}
    nb["compound"] = compound
    nb["concentration_M"] = conc_M
    nb["occupancy"] = occ["receptors"] if occ else []
    nb["readouts"] = {
        "n_nodes": g["n_nodes"],
        "n_edges": g["n_edges"],
        "named": named,
        "mean_hz": float(r.mean()),
        "max_hz": float(r.max()),
        "g_ach": float(g_ach),
        "g_gaba": float(g_gaba),
    }
    nb["warnings"] = [
        "Hops-limited MaleCNS neighborhood, not the full 25M-edge CNS.",
        "Signs from predicted transmitters. ACh gain and RDL gain are teaching patches.",
    ]
    return nb


def dose_response_subgraph(compound: str, concs=None):
    concs = concs or [10**e for e in range(-10, -4)]
    points = []
    for c in concs:
        nb = run_subgraph_assay(compound, c)
        mn9 = nb["readouts"]["named"].get("MN9", [])
        hz = [row["hz"] for row in mn9 if row["hz"] is not None]
        points.append({
            "conc_M": c,
            "mean_hz": nb["readouts"]["mean_hz"],
            "mn9_hz": float(np.mean(hz)) if hz else None,
            "g_ach": nb["readouts"]["g_ach"],
            "g_gaba": nb["readouts"]["g_gaba"],
        })
    return points
