from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from flylab.notebook.schema import empty_notebook
from flylab.pharm.occupancy import compare_compound

SIGN = {"acetylcholine": 1.0, "glutamate": -0.4, "gaba": -1.0, "histamine": -0.5, "dopamine": 0.2, "serotonin": 0.2, "octopamine": 0.2}
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

def _g_ach(compound, conc_M):
    if not compound:
        return 1.0, None
    occ = compare_compound(compound, conc_M)
    row = next(r for r in occ["receptors"] if r["receptor"] == "insect_nAChR")
    th, d = row["occupancy"], row["direction"]
    if d == "agonist":
        gain = max(0.05, 1.0 + 0.4 * th - 1.6 * th * th)
    elif d == "antagonist":
        gain = max(0.05, 1.0 - th)
    else:
        gain = 1.0
    return gain, occ

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
    g_ach, occ = _g_ach(compound, conc_M)
    scale = np.array([g_ach if nd.get("consensus_nt") == "acetylcholine" else 1.0 for nd in nodes])
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
    nb["readouts"] = {"n_nodes": g["n_nodes"], "n_edges": g["n_edges"], "named": named, "mean_hz": float(r.mean()), "max_hz": float(r.max())}
    nb["warnings"] = ["Hops-limited MaleCNS neighborhood, not the full 25M-edge CNS.", "Signs from predicted transmitters."]
    return nb
