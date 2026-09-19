"""Cut named-cell neighborhoods from MaleCNS weights via Arrow batches."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import pandas as pd
import pyarrow as pa

from flylab.maps.malecns import FILES, MAP_CITATION, MAP_ID, default_dir

DEFAULT_TYPES = ("MN9", "DNp01")
GUSTATORY_TYPES = ("LB1a", "LB1b", "LB1c", "LB1d", "LB3b", "LB3c")


def seed_ids(ann: pd.DataFrame, types: Iterable[str]) -> dict[str, list[int]]:
    """Exact type match, then prefix match (LB1a matches LB1a_L)."""
    col = ann["type"].fillna("").astype(str)
    out: dict[str, list[int]] = {}
    for t in types:
        exact = ann[col == t]
        if len(exact):
            hits = exact
        else:
            hits = ann[col.str.startswith(t)]
        out[t] = sorted({int(x) for x in hits["bodyId"].tolist()})
    return out


def _touching_edges(weights_path: Path, ids: set[int], min_weight: int):
    src = pa.memory_map(str(weights_path), "r")
    reader = pa.ipc.open_file(src)
    rows = []
    idset = ids
    for i in range(reader.num_record_batches):
        b = reader.get_batch(i)
        pre = b.column("body_pre").to_pylist()
        post = b.column("body_post").to_pylist()
        w = b.column("weight").to_pylist()
        for a, c, wt in zip(pre, post, w):
            if wt < min_weight:
                continue
            if a in idset or c in idset:
                rows.append((int(a), int(c), int(wt)))
    return rows


def extract_neighborhood(dest=None, types=DEFAULT_TYPES, hops=1, min_weight=5, out=None):
    dest = dest or default_dir()
    weights_path = dest / FILES["weights"]
    if not weights_path.exists():
        raise FileNotFoundError(f"{weights_path} missing")
    ann = pd.read_feather(dest / FILES["annotations"])
    nt = pd.read_feather(dest / FILES["neurotransmitters"])
    traced = ann[ann["status"] == "Traced"]
    seeds = seed_ids(traced, types)
    keep = {i for ids in seeds.values() for i in ids}
    if not keep:
        raise ValueError(f"no seeds for types={list(types)}")
    partner = set(keep)
    all_edges = []
    frontier = set(keep)
    for _ in range(max(hops, 0)):
        if not frontier:
            break
        rows = _touching_edges(weights_path, frontier, min_weight)
        all_edges.extend(rows)
        nxt = set()
        for a, c, wt in rows:
            if a not in partner:
                nxt.add(a)
            if c not in partner:
                nxt.add(c)
            partner.add(a)
            partner.add(c)
        frontier = nxt
    seen = set()
    edges = []
    for a, c, wt in all_edges:
        if a in partner and c in partner and (a, c) not in seen:
            seen.add((a, c))
            edges.append({"pre": a, "post": c, "weight": float(wt)})
    meta = traced.merge(nt, left_on="bodyId", right_on="body", how="left").set_index("bodyId")
    nodes = []
    for body in sorted(partner):
        if body in meta.index:
            r = meta.loc[body]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            nodes.append({
                "bodyId": int(body),
                "type": None if pd.isna(r.get("type")) else str(r["type"]),
                "superclass": None if pd.isna(r.get("superclass")) else str(r["superclass"]),
                "consensus_nt": None if pd.isna(r.get("consensus_nt")) else str(r["consensus_nt"]),
            })
        else:
            nodes.append({"bodyId": int(body), "type": None, "superclass": None, "consensus_nt": None})
    payload = {
        "map": MAP_ID,
        "citation": MAP_CITATION,
        "types": list(types),
        "hops": hops,
        "min_weight": min_weight,
        "n_nodes": len(nodes),
        "n_edges": len(edges),
        "seeds": seeds,
        "nodes": nodes,
        "edges": edges,
    }
    dest_out = out or Path("data/derived/malecns_named_neighborhood.json")
    dest_out.parent.mkdir(parents=True, exist_ok=True)
    dest_out.write_text(json.dumps(payload))
    payload["path"] = str(dest_out)
    return payload
