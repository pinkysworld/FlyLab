"""Cut named-cell neighborhoods from MaleCNS weights via Arrow batches.

Two products live here:

* :func:`extract_neighborhood` -- the original hops-limited cut that produced
  the committed ``named`` and ``taste_motor`` graphs.  Unchanged.
* :func:`extract_scaled_cut` / :func:`extract_ladder` -- a **ladder** of
  nested cuts at prescribed node budgets, built so that the only thing that
  differs between rungs is how many cells were kept.  See
  :data:`LADDER_RECIPE` for the rule and why it is that rule.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
import pyarrow as pa

from flylab.maps.ladder import (
    COMMIT_BYTE_BUDGET,
    DEFAULT_TYPES,
    GUSTATORY_TYPES,
    LADDER_MIN_WEIGHT,
    LADDER_RECIPE,
    LADDER_SEED_TYPES,
    LADDER_SIZES,
    TASTE_MOTOR_TYPES,
    ladder_filename,
)
from flylab.maps.malecns import (
    FILES,
    MAP_CITATION,
    MAP_ID,
    default_dir,
    edge_arrays,
    traced_body_ids,
)


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


def _induced_edges(weights_path: Path, ids: set[int], min_weight: int):
    """Edges with BOTH ends inside ``ids`` (induced subgraph), via Arrow compute."""
    import pyarrow.compute as pc

    keep = pa.array(sorted(ids), type=pa.int64())
    reader = pa.ipc.open_file(pa.memory_map(str(weights_path), "r"))
    rows = []
    for i in range(reader.num_record_batches):
        b = reader.get_batch(i)
        pre = b.column("body_pre").cast(pa.int64())
        post = b.column("body_post").cast(pa.int64())
        mask = pc.and_(pc.is_in(pre, keep), pc.is_in(post, keep))
        mask = pc.and_(mask, pc.greater_equal(b.column("weight"), min_weight))
        f = b.filter(mask)
        for a, c, wt in zip(f.column("body_pre").to_pylist(), f.column("body_post").to_pylist(), f.column("weight").to_pylist()):
            rows.append((int(a), int(c), int(wt)))
    return rows


def extract_neighborhood(dest=None, types=DEFAULT_TYPES, hops=1, min_weight=5, out=None, closure_min_weight=None):
    """Cut the ``hops``-hop neighborhood of the seed ``types``.

    ``closure_min_weight``: if set, add every edge of at least that weight whose two
    ends are both already in the node set (induced closure). No new nodes are added,
    but multi-synapse paths between seed groups (e.g. GRN -> relay -> MN9) become
    visible. Off by default so the original MN9/DNp01 product is unchanged.
    """
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
    if closure_min_weight is not None:
        all_edges.extend(_induced_edges(weights_path, partner, int(closure_min_weight)))
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
        "closure_min_weight": closure_min_weight,
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


# ==========================================================================
# the scale ladder
# ==========================================================================
def _growth_order(
    pre: "np.ndarray",
    post: "np.ndarray",
    w: "np.ndarray",
    bodies: "np.ndarray",
    seed_bodies: "np.ndarray",
    budget: int,
) -> tuple["np.ndarray", list[int]]:
    """Indices into ``bodies``, in ladder order, plus the size of each rank.

    ``bodies`` must be sorted; ``pre``/``post`` are body ids.  See
    :data:`LADDER_RECIPE`.
    """
    n = int(bodies.size)
    ip = np.searchsorted(bodies, pre)
    iq = np.searchsorted(bodies, post)
    wf = np.asarray(w, dtype=np.float64)
    sel = np.zeros(n, dtype=bool)
    pos = np.searchsorted(bodies, seed_bodies)
    pos = np.clip(pos, 0, max(n - 1, 0))
    seed_idx = pos[bodies[pos] == seed_bodies]
    sel[seed_idx] = True
    order: list[int] = sorted(int(i) for i in seed_idx)
    ranks = [len(order)]
    while len(order) < int(budget):
        touch = sel[ip] | sel[iq]
        if not touch.any():
            break
        score = np.bincount(ip[touch], weights=wf[touch], minlength=n)
        score += np.bincount(iq[touch], weights=wf[touch], minlength=n)
        score[sel] = -1.0
        cand = np.nonzero(score > 0)[0]
        if cand.size == 0:
            break
        # descending score, ascending bodyId: lexsort's last key is primary
        ring = cand[np.lexsort((bodies[cand], -score[cand]))]
        take = ring[: int(budget) - len(order)]
        sel[take] = True
        order.extend(int(i) for i in take)
        ranks.append(int(take.size))
    return np.asarray(order, dtype=np.int64), ranks


def _node_records(
    body_ids: "np.ndarray", ann: pd.DataFrame, nt: pd.DataFrame
) -> list[dict[str, Any]]:
    """``bodyId``/``type``/``superclass``/``consensus_nt`` rows, sorted by bodyId."""
    want = pd.DataFrame({"bodyId": np.asarray(body_ids, dtype="int64")})
    a = ann[["bodyId", "type", "superclass"]].drop_duplicates("bodyId")
    t = nt[["body", "consensus_nt"]].drop_duplicates("body")
    m = want.merge(a, on="bodyId", how="left").merge(
        t, left_on="bodyId", right_on="body", how="left"
    )
    out: list[dict[str, Any]] = []
    for r in m.itertuples():
        out.append(
            {
                "bodyId": int(r.bodyId),
                "type": None if pd.isna(r.type) else str(r.type),
                "superclass": None if pd.isna(r.superclass) else str(r.superclass),
                "consensus_nt": None if pd.isna(r.consensus_nt) else str(r.consensus_nt),
            }
        )
    return out


def _composition(nodes: list[dict[str, Any]], edges: list[dict[str, Any]]) -> dict[str, Any]:
    """Transmitter composition of a cut, by cell count and by outgoing weight.

    The weighted shares are the ones the drug actually acts through: the gain
    patch multiplies a presynaptic cell's outgoing edges, so what matters is
    each transmitter's share of total *outgoing synaptic weight*, not its share
    of cells.  Referee finding B2 turns on exactly this distinction.
    """
    nt_of = {nd["bodyId"]: (nd.get("consensus_nt") or "unclear") for nd in nodes}
    counts: dict[str, int] = {}
    for v in nt_of.values():
        counts[v] = counts.get(v, 0) + 1
    out_w: dict[str, float] = {}
    total = 0.0
    for e in edges:
        k = nt_of.get(e["pre"], "unclear")
        wt = float(e["weight"])
        out_w[k] = out_w.get(k, 0.0) + wt
        total += wt
    n = max(len(nodes), 1)
    return {
        "n_nodes": len(nodes),
        "cell_counts": dict(sorted(counts.items())),
        "cell_shares": {k: v / n for k, v in sorted(counts.items())},
        "out_weight": dict(sorted(out_w.items())),
        "out_weight_shares": (
            {k: v / total for k, v in sorted(out_w.items())} if total > 0 else {}
        ),
        "total_out_weight": total,
        "note": (
            "out_weight_shares is the composition the gain patch acts through "
            "(gains scale a presynaptic cell's outgoing edges); cell_shares is "
            "the histogram a label permutation preserves."
        ),
    }


def _census_of_payload(path: Path, seed_types: Sequence[str]) -> dict[str, Any] | None:
    """:func:`flylab.analysis.dependence.cut_census` on a just-written cut.

    Imported lazily (the analysis package pulls in the whole engine) and the
    memo caches are dropped afterwards so the file can be rewritten with the
    census inside it without leaving a stale graph cached in-process.
    """
    try:
        from flylab.analysis.dependence import cut_census
    except Exception:  # pragma: no cover - analysis layer optional at extract time
        return None
    try:
        census = cut_census(str(path), seed_types=list(seed_types))
    except Exception:  # pragma: no cover
        return None
    finally:
        _forget_cached_graph(path)
    census.pop("graph", None)
    return census


def _forget_cached_graph(path: Path) -> None:
    from flylab.circuit.rate import forget_graph

    forget_graph(path)
    key = str(Path(path).resolve())
    try:
        from flylab.analysis import nullmodels as nm

        nm._STATE_CACHE.pop(key, None)
    except Exception:  # pragma: no cover
        pass


def extract_scaled_cut(
    n_target: int,
    dest: Path | None = None,
    types: Sequence[str] = LADDER_SEED_TYPES,
    min_weight: int = LADDER_MIN_WEIGHT,
    out: Path | str | None = None,
    census: bool = True,
    _cache: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """One rung of the scale ladder: the induced cut on the first ``n_target`` cells.

    The rule and its justification are in :data:`LADDER_RECIPE`, and the
    payload records it verbatim under ``recipe`` together with the cut's
    :func:`~flylab.analysis.dependence.cut_census` statistics and its
    transmitter composition, so a reader of the scaling table can see *how*
    the cut was made and *what changed structurally* beside the verdict.

    ``_cache`` is an opaque dict reused by :func:`extract_ladder` so the 1.0 GB
    weight matrix is streamed once for the whole ladder instead of once a rung.
    """
    dest = Path(dest) if dest else default_dir()
    cache = _cache if _cache is not None else {}
    if "edges" not in cache:
        ann_all = pd.read_feather(dest / FILES["annotations"])
        cache["ann"] = ann_all[ann_all["status"] == "Traced"]
        cache["nt"] = pd.read_feather(dest / FILES["neurotransmitters"])
        traced = traced_body_ids(dest)
        cache["edges"] = edge_arrays(dest, min_weight=int(min_weight), keep=traced)
        cache["min_weight"] = int(min_weight)
    if cache.get("min_weight") != int(min_weight):
        raise ValueError(
            f"cached edges were filtered at min_weight={cache.get('min_weight')}, "
            f"not {int(min_weight)}"
        )
    ann, nt = cache["ann"], cache["nt"]
    pre, post, w = cache["edges"]

    seeds = seed_ids(ann, types)
    seed_bodies = np.array(
        sorted({int(i) for ids in seeds.values() for i in ids}), dtype=np.int64
    )
    if seed_bodies.size == 0:
        raise ValueError(f"no traced seeds for types={list(types)}")
    if "bodies" not in cache:
        cache["bodies"] = np.union1d(np.union1d(pre, post), seed_bodies)
    bodies = cache["bodies"]

    order, ranks = _growth_order(pre, post, w, bodies, seed_bodies, int(n_target))
    keep_bodies = np.sort(bodies[order])
    from flylab.maps.malecns import _in_sorted

    mask = _in_sorted(pre, keep_bodies) & _in_sorted(post, keep_bodies)
    ep, eq, ew = pre[mask], post[mask], w[mask]
    edges = [
        {"pre": int(a), "post": int(b), "weight": int(round(float(c)))}
        for a, b, c in zip(ep.tolist(), eq.tolist(), ew.tolist())
    ]
    nodes = _node_records(keep_bodies, ann, nt)
    kept_seeds = {}
    for t, ids in seeds.items():
        arr = np.asarray(sorted(int(i) for i in ids), dtype=np.int64)
        kept_seeds[t] = (
            [int(x) for x in arr[_in_sorted(arr, keep_bodies)]] if arr.size else []
        )
    payload: dict[str, Any] = {
        "map": MAP_ID,
        "citation": MAP_CITATION,
        "cut": "scale_ladder",
        "n_target": int(n_target),
        "types": list(types),
        "min_weight": int(min_weight),
        "hops": None,
        "closure_min_weight": int(min_weight),
        "growth": "ranked_bfs_induced",
        "recipe": LADDER_RECIPE,
        "rank_sizes": ranks,
        "n_ranks": len(ranks),
        "n_nodes": len(nodes),
        "n_edges": len(edges),
        "seeds": kept_seeds,
        "nodes": nodes,
        "edges": edges,
    }
    payload["composition"] = _composition(nodes, edges)
    out_path = Path(out) if out is not None else Path("data/derived") / ladder_filename(n_target)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload))
    if census:
        c = _census_of_payload(out_path, list(kept_seeds))
        if c is not None:
            payload["census"] = c
            out_path.write_text(json.dumps(payload))
    payload["path"] = str(out_path)
    payload["bytes"] = out_path.stat().st_size
    return payload


def ladder_summary(payload: dict[str, Any]) -> dict[str, Any]:
    """The payload without ``nodes``/``edges``: what a manifest or a log wants."""
    return {k: v for k, v in payload.items() if k not in ("nodes", "edges", "recipe")}


def extract_ladder(
    sizes: Sequence[int] = LADDER_SIZES,
    dest: Path | None = None,
    types: Sequence[str] = LADDER_SEED_TYPES,
    min_weight: int = LADDER_MIN_WEIGHT,
    out_dir: Path | str = Path("data/derived"),
    census: bool = True,
) -> dict[str, Any]:
    """Build the whole ladder, streaming the weight matrix once.

    Returns ``{"cuts": [summary, ...], "recipe": ..., "sizes": ...}``; each cut
    is written to ``out_dir / ladder_filename(size)``.  Only the rungs small
    enough to belong in git should be committed -- see ``bytes`` in each
    summary and the ``COMMIT_BYTE_BUDGET`` the workflow applies.
    """
    cache: dict[str, Any] = {}
    out_dir = Path(out_dir)
    cuts = []
    for k in sorted(int(s) for s in sizes):
        p = extract_scaled_cut(
            k,
            dest=dest,
            types=types,
            min_weight=min_weight,
            out=out_dir / ladder_filename(k),
            census=census,
            _cache=cache,
        )
        cuts.append(ladder_summary(p))
    return {
        "sizes": [int(s) for s in sizes],
        "seed_types": list(types),
        "min_weight": int(min_weight),
        "recipe": LADDER_RECIPE,
        "map": MAP_ID,
        "citation": MAP_CITATION,
        "cuts": cuts,
        "nested": True,
        "label": "model_derived",
    }


__all__ = [
    "DEFAULT_TYPES",
    "GUSTATORY_TYPES",
    "TASTE_MOTOR_TYPES",
    "seed_ids",
    "extract_neighborhood",
    "LADDER_SIZES",
    "LADDER_SEED_TYPES",
    "LADDER_MIN_WEIGHT",
    "LADDER_RECIPE",
    "COMMIT_BYTE_BUDGET",
    "ladder_filename",
    "ladder_summary",
    "extract_scaled_cut",
    "extract_ladder",
]
