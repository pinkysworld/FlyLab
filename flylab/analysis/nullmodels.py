"""Connectome null models: "does the MaleCNS wiring actually matter?"

Every circuit readout in FlyLab is produced by pushing drive through a *real*
MaleCNS neighborhood.  That only means something if a **destroyed** version of
the same neighborhood gives a different answer.  This module builds four
degradations of the graph and measures the drug effect on each of them:

==============================  ==========================================
mode                            what is destroyed / what is kept
==============================  ==========================================
``sign_permute``                consensus transmitter labels are permuted
                                across nodes.  The NT histogram (and hence
                                the global E/I balance) is preserved exactly;
                                *which* cell is cholinergic is destroyed.
``weight_permute``              synapse counts are permuted across edges.
                                Topology and transmitters are kept; the
                                weight-topology pairing is destroyed.
``rewire_degree_preserving``    double-edge swaps (``~10 x E`` accepted
                                swaps).  Every node keeps its out-degree
                                (pre) and in-degree (post) and its
                                transmitter; the wiring pattern is destroyed.
``erdos_renyi``                 same N and E, endpoints drawn uniformly,
                                weights drawn from the empirical weight list.
                                Only the size of the graph survives.
==============================  ==========================================

Invariants shared by every mode: the node set, the node order, the ``seeds``
block, ``n_nodes`` and ``n_edges`` are unchanged, so **seeds stay seeds** and
the same named cells are driven and read in the shuffled graph as in the real
one.  Self-loops and duplicated (pre, post) pairs are never created.

Nothing in here is fitted to animal data.  A z-score against a shuffled
connectome is a statement about this model, not about a fly, so every result
carries ``label = "model_derived"`` and a ``warnings`` list.
"""
from __future__ import annotations

import contextlib
import json
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

__all__ = [
    "MODES",
    "ASSAYS",
    "DEFAULT_READOUT",
    "shuffle_graph",
    "null_distribution",
    "null_panel",
    "connectome_information_score",
]

#: the four degradations, in increasing order of destruction
MODES: tuple[str, ...] = (
    "sign_permute",
    "weight_permute",
    "rewire_degree_preserving",
    "erdos_renyi",
)

#: assays this module can drive
ASSAYS: tuple[str, ...] = ("subgraph", "spiking", "taste_map")

#: readout that actually carries the dose-response, per assay.  On the
#: seed-driven ``subgraph`` rate assay MN9 is clamped by its own drive, so the
#: network mean is the informative readout (see
#: ``flylab.assays.ensemble.DEFAULT_READOUT``).
DEFAULT_READOUT: dict[str, str] = {
    "subgraph": "mean_hz",
    "spiking": "mn9_hz",
    "taste_map": "mn9_hz",
}

#: default graph per assay
DEFAULT_GRAPH: dict[str, str] = {
    "subgraph": "named",
    "spiking": "named",
    "taste_map": "taste_motor",
}

BASE_WARNINGS = [
    "Null models destroy the MaleCNS wiring on purpose: they are controls, not "
    "biology. Nothing here is a measurement from a living fly.",
    "The p-value is an empirical two-sided permutation p, (k+1)/(n+1); with n "
    "shuffles its floor is 1/(n+1).",
    "model_derived: z-scores describe this simulation and its teaching EC50 "
    "library, not a measured drug effect.",
]

MODE_NOTES: dict[str, str] = {
    "sign_permute": (
        "sign_permute keeps the transmitter histogram exactly, so a small |z| "
        "means the drug effect follows global E/I balance rather than the "
        "identity of the cholinergic cells."
    ),
    "weight_permute": (
        "weight_permute keeps topology and transmitters, so it isolates the "
        "pairing between synapse counts and wiring."
    ),
    "rewire_degree_preserving": (
        "rewire_degree_preserving keeps every node's in/out degree and "
        "transmitter, so it isolates the wiring pattern itself."
    ),
    "erdos_renyi": (
        "erdos_renyi keeps only N, E and the weight histogram; it is the "
        "weakest null and should be the easiest to beat."
    ),
}


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _rng(rng: Any) -> np.random.Generator:
    if isinstance(rng, np.random.Generator):
        return rng
    return np.random.default_rng(rng)


def _f(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


# --------------------------------------------------------------------------
# graph shuffles
# --------------------------------------------------------------------------
def _double_edge_swap(
    edges: list[dict[str, Any]],
    rng: np.random.Generator,
    n_swaps: int,
    max_tries_per_swap: int = 20,
) -> tuple[list[dict[str, Any]], int, int]:
    """Degree-preserving double-edge swaps.

    ``(a->b, c->d)`` becomes ``(a->d, c->b)``: node ``a`` and ``c`` keep their
    out-degree, ``b`` and ``d`` keep their in-degree.  The weight travels with
    its presynaptic node.  Swaps that would make a self-loop or a duplicate
    edge are rejected.
    """
    pre = [int(e["pre"]) for e in edges]
    post = [int(e["post"]) for e in edges]
    present = {(p, q) for p, q in zip(pre, post)}
    m = len(edges)
    done = 0
    tries = 0
    if m < 2:
        return [dict(e) for e in edges], 0, 0
    budget = int(n_swaps) * int(max_tries_per_swap)
    chunk = max(int(n_swaps), 1024)
    while done < n_swaps and tries < budget:
        picks = rng.integers(0, m, size=(min(chunk, budget - tries), 2))
        for i, j in picks:
            tries += 1
            i, j = int(i), int(j)
            if i == j:
                continue
            a, b = pre[i], post[i]
            c, d = pre[j], post[j]
            if a == d or c == b:
                continue
            if (a, d) in present or (c, b) in present:
                continue
            present.discard((a, b))
            present.discard((c, d))
            present.add((a, d))
            present.add((c, b))
            post[i], post[j] = d, b
            done += 1
            if done >= n_swaps:
                break
    out = []
    for e, p, q in zip(edges, pre, post):
        row = dict(e)
        row["pre"], row["post"] = p, q
        out.append(row)
    return out, done, tries


def _erdos_renyi(
    body_ids: Sequence[int],
    edges: Sequence[dict[str, Any]],
    rng: np.random.Generator,
) -> list[dict[str, Any]]:
    """Same N and E; endpoints uniform, weights drawn from the empirical list."""
    n = len(body_ids)
    m = len(edges)
    weights = np.array([float(e["weight"]) for e in edges], dtype=float)
    ids = np.asarray(body_ids)
    taken: set[tuple[int, int]] = set()
    out: list[dict[str, Any]] = []
    drawn = rng.choice(weights, size=m, replace=True) if m else np.array([])
    tries = 0
    while len(out) < m and tries < 200 * max(m, 1):
        need = m - len(out)
        aa = rng.integers(0, n, size=need)
        bb = rng.integers(0, n, size=need)
        for a, b in zip(aa, bb):
            tries += 1
            if a == b:
                continue
            key = (int(ids[a]), int(ids[b]))
            if key in taken:
                continue
            taken.add(key)
            out.append({"pre": key[0], "post": key[1], "weight": float(drawn[len(out)])})
            if len(out) >= m:
                break
    return out


def shuffle_graph(
    graph: dict[str, Any],
    mode: str = "sign_permute",
    rng: Any = 0,
) -> dict[str, Any]:
    """Return a degraded copy of ``graph`` for one null ``mode``.

    ``rng`` is a ``numpy.random.Generator`` or any seed accepted by
    ``numpy.random.default_rng``; the result is deterministic given the seed.
    The input graph is never mutated (the module-level graph cache is shared).

    See the module docstring for what each mode keeps and destroys.  The node
    list, its order, the ``seeds`` block, ``n_nodes`` and ``n_edges`` are
    preserved by every mode, so seeds stay seeds.
    """
    if mode not in MODES:
        raise ValueError(f"unknown null mode {mode!r}; expected one of {MODES}")
    r = _rng(rng)
    nodes = [dict(n) for n in graph["nodes"]]
    edges = [dict(e) for e in graph["edges"]]
    meta: dict[str, Any] = {"mode": mode}

    if mode == "sign_permute":
        nts = [n.get("consensus_nt") for n in nodes]
        perm = r.permutation(len(nodes))
        for i, node in enumerate(nodes):
            node["consensus_nt"] = nts[int(perm[i])]
        meta["permuted_nodes"] = len(nodes)
    elif mode == "weight_permute":
        w = [float(e["weight"]) for e in edges]
        perm = r.permutation(len(edges))
        for i, e in enumerate(edges):
            e["weight"] = w[int(perm[i])]
        meta["permuted_edges"] = len(edges)
    elif mode == "rewire_degree_preserving":
        target = 10 * len(edges)
        edges, done, tries = _double_edge_swap(edges, r, n_swaps=target)
        meta.update({"swaps_requested": target, "swaps_done": done, "attempts": tries})
    else:  # erdos_renyi
        edges = _erdos_renyi([n["bodyId"] for n in nodes], edges, r)
        meta["resampled_edges"] = len(edges)

    out = dict(graph)
    out["nodes"] = nodes
    out["edges"] = edges
    out["n_nodes"] = len(nodes)
    out["n_edges"] = len(edges)
    out["seeds"] = {k: list(v) for k, v in (graph.get("seeds") or {}).items()}
    out["null_model"] = meta
    return out


# --------------------------------------------------------------------------
# running an assay on a shuffled graph
# --------------------------------------------------------------------------
def _forget_graph(path: Path, graph_obj: dict[str, Any] | None = None) -> None:
    """Drop a temp graph from the rate / LIF memo caches."""
    from flylab.circuit import lif as lif_mod
    from flylab.circuit import rate as rate_mod

    key = str(Path(path).resolve())
    cached = rate_mod._GRAPH_CACHE.pop(key, None)
    rate_mod._NET_CACHE.pop(key, None)
    for g in (cached, graph_obj):
        if g is None:
            continue
        for k in [k for k in lif_mod._LIF_CACHE if k and k[0] == id(g)]:
            lif_mod._LIF_CACHE.pop(k, None)


@contextlib.contextmanager
def _graph_source(assay: str, graph_obj: dict[str, Any] | None, graph: str | None):
    """Yield the assay kwargs that make the assay read ``graph_obj``.

    ``subgraph`` and ``spiking`` accept ``graph_obj=`` directly.
    ``run_taste_map_assay`` does not, and this module must not edit it, so the
    shuffled graph is written to a temp JSON and passed as a path (which
    ``flylab.circuit.rate.resolve_graph`` accepts).  The temp file and its
    memoised network are removed afterwards.
    """
    name = graph if graph is not None else DEFAULT_GRAPH[assay]
    if graph_obj is None:
        yield {"graph": name}
        return
    if assay in ("subgraph", "spiking"):
        yield {"graph": name, "graph_obj": graph_obj}
        return
    tmpdir = Path(tempfile.mkdtemp(prefix="flylab_null_"))
    path = tmpdir / "shuffled_graph.json"
    path.write_text(json.dumps(graph_obj))
    try:
        yield {"graph": str(path)}
    finally:
        _forget_graph(path, graph_obj)
        shutil.rmtree(tmpdir, ignore_errors=True)


def _run_assay(
    assay: str,
    compound: str | None,
    conc_M: float,
    source_kw: dict[str, Any],
    seed: int = 0,
    **kw: Any,
) -> dict[str, Any]:
    if assay == "subgraph":
        from flylab.assays.subgraph import run_subgraph_assay

        return run_subgraph_assay(compound, conc_M, **source_kw, **kw)
    if assay == "spiking":
        from flylab.assays.spiking import run_spiking_assay

        return run_spiking_assay(compound, conc_M, seed=seed, **source_kw, **kw)
    if assay == "taste_map":
        from flylab.assays.taste_map import run_taste_map_assay

        return run_taste_map_assay(compound, conc_M, seed=seed, **source_kw, **kw)
    raise ValueError(f"unknown assay {assay!r}; expected one of {ASSAYS}")


def readout_value(nb: dict[str, Any], readout: str) -> float | None:
    """One scalar readout out of a notebook (``readouts`` then ``gains``)."""
    r = nb.get("readouts") or {}
    if readout in r:
        return _f(r[readout])
    return _f((nb.get("gains") or {}).get(readout))


def vehicle_value(assay: str, nb: dict[str, Any], readout: str) -> float | None:
    """Vehicle value already embedded in a treated notebook, if there is one.

    ``subgraph`` and ``spiking`` carry ``readouts["vehicle"]``; ``taste_map``
    carries ``mn9_vehicle_sugar_hz``.  Anything else needs a second run.
    """
    r = nb.get("readouts") or {}
    veh = r.get("vehicle")
    if isinstance(veh, dict) and readout in veh:
        return _f(veh[readout])
    if assay == "taste_map" and readout in ("mn9_hz", "mn9_sugar_hz"):
        return _f(r.get("mn9_vehicle_sugar_hz"))
    return None


def drug_effect(
    assay: str,
    compound: str | None,
    conc_M: float,
    readout: str,
    graph_obj: dict[str, Any] | None = None,
    graph: str | None = None,
    seed: int = 0,
    **kw: Any,
) -> dict[str, Any]:
    """``treated - vehicle`` for one readout, both on the *same* graph."""
    with _graph_source(assay, graph_obj, graph) as source_kw:
        nb = _run_assay(assay, compound, conc_M, source_kw, seed=seed, **kw)
        treated = readout_value(nb, readout)
        vehicle = vehicle_value(assay, nb, readout)
        if vehicle is None:
            veh_nb = _run_assay(assay, None, 0.0, source_kw, seed=seed, **kw)
            vehicle = readout_value(veh_nb, readout)
    effect = None if (treated is None or vehicle is None) else float(treated - vehicle)
    return {"treated": treated, "vehicle": vehicle, "effect": effect}


# --------------------------------------------------------------------------
# null distribution
# --------------------------------------------------------------------------
def null_distribution(
    assay: str = "subgraph",
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    mode: str = "sign_permute",
    n: int = 50,
    seed: int = 0,
    readout: str = "auto",
    graph: str | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Drug effect on the real graph vs ``n`` shuffled graphs.

    The effect is always ``treated - vehicle`` for ``readout``, with both arms
    run on the *same* graph, so a shuffle can only change the effect through
    the wiring.  ``readout="auto"`` uses :data:`DEFAULT_READOUT`.

    Returns ``real_effect``, the ``null_effects`` list, ``null_mean``,
    ``null_sd``, ``z``, the empirical two-sided ``p_two_sided`` ``(k+1)/(n+1)``,
    plus ``mode``, ``readout``, ``n`` and ``runtime_s``.
    """
    from flylab.circuit.rate import load_graph

    if assay not in ASSAYS:
        raise ValueError(f"unknown assay {assay!r}; expected one of {ASSAYS}")
    if mode not in MODES:
        raise ValueError(f"unknown null mode {mode!r}; expected one of {MODES}")
    if readout == "auto":
        readout = DEFAULT_READOUT[assay]
    t0 = time.perf_counter()
    graph_name = graph if graph is not None else DEFAULT_GRAPH[assay]
    base = load_graph(graph_name)

    warnings = list(BASE_WARNINGS) + [MODE_NOTES[mode]]
    real = drug_effect(assay, compound, conc_M, readout, graph=graph_name, **kw)

    nulls: list[float | None] = []
    for i in range(int(n)):
        shuffled = shuffle_graph(base, mode, np.random.default_rng([int(seed), i]))
        res = drug_effect(
            assay, compound, conc_M, readout, graph_obj=shuffled, graph=graph_name, **kw
        )
        nulls.append(res["effect"])

    ok = np.array([v for v in nulls if v is not None], dtype=float)
    n_ok = int(ok.size)
    null_mean = float(ok.mean()) if n_ok else None
    null_sd = float(ok.std(ddof=1)) if n_ok > 1 else 0.0
    real_effect = real["effect"]

    z: float | None = None
    p: float | None = None
    if real_effect is None:
        warnings.append(
            f"readout {readout!r} is undefined on the real graph for this "
            "compound; no z-score could be formed."
        )
    elif n_ok == 0:
        warnings.append("every shuffled graph gave an undefined readout; no null.")
    else:
        if null_sd > 0:
            z = float((real_effect - null_mean) / null_sd)
        else:
            warnings.append(
                "the null distribution has zero spread; z is undefined (the "
                "shuffles all gave the same effect)."
            )
        k = int(np.sum(np.abs(ok - null_mean) >= abs(real_effect - null_mean)))
        p = float((k + 1) / (n_ok + 1))
    if n_ok < int(n):
        warnings.append(f"{int(n) - n_ok} of {int(n)} shuffles gave an undefined readout.")

    return {
        "assay": assay,
        "compound": compound,
        "conc_M": float(conc_M),
        "graph": graph_name,
        "mode": mode,
        "readout": readout,
        "n": int(n),
        "n_ok": n_ok,
        "seed": int(seed),
        "real_effect": real_effect,
        "real_treated": real["treated"],
        "real_vehicle": real["vehicle"],
        "null_effects": nulls,
        "null_mean": null_mean,
        "null_sd": null_sd,
        "z": z,
        "p_two_sided": p,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# panel over modes
# --------------------------------------------------------------------------
def null_panel(
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    modes: Iterable[str] = MODES,
    n: int = 30,
    seed: int = 0,
    n_taste: int | None = None,
    include_taste_map: bool = True,
    **kw: Any,
) -> dict[str, Any]:
    """All null modes for the ``subgraph`` mean rate and the taste-map veto.

    The taste-map arm is the bitter veto ratio (``mn9_sugar_bitter / mn9_sugar``)
    on the 1841-node ``taste_motor`` graph, which needs two runs per shuffle and
    a temp JSON, so it defaults to ``n_taste = min(n, 10)`` shuffles.
    """
    t0 = time.perf_counter()
    modes = tuple(modes)
    n_taste = int(n_taste) if n_taste is not None else min(int(n), 10)
    rows: list[dict[str, Any]] = []
    warnings = list(BASE_WARNINGS)

    arms: list[tuple[str, str, int]] = [("subgraph", "mean_hz", int(n))]
    if include_taste_map:
        arms.append(("taste_map", "bitter_veto_ratio", n_taste))

    for assay, readout, n_i in arms:
        for mode in modes:
            res = null_distribution(
                assay, compound, conc_M, mode, n=n_i, seed=seed, readout=readout, **kw
            )
            rows.append(
                {
                    "assay": assay,
                    "readout": readout,
                    "mode": mode,
                    "n": res["n"],
                    "n_ok": res["n_ok"],
                    "real_effect": res["real_effect"],
                    "null_mean": res["null_mean"],
                    "null_sd": res["null_sd"],
                    "z": res["z"],
                    "p_two_sided": res["p_two_sided"],
                    "runtime_s": res["runtime_s"],
                }
            )
            for w in res["warnings"]:
                if w not in warnings:
                    warnings.append(w)

    return {
        "compound": compound,
        "conc_M": float(conc_M),
        "modes": list(modes),
        "n": int(n),
        "n_taste": n_taste,
        "seed": int(seed),
        "rows": rows,
        "summary": connectome_information_score(panel=rows),
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }


def connectome_information_score(
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    panel: dict[str, Any] | list[dict[str, Any]] | None = None,
    assay: str | None = None,
    modes: Iterable[str] = MODES,
    n: int = 30,
    seed: int = 0,
    z_threshold: float = 2.0,
    **kw: Any,
) -> dict[str, Any]:
    """Fraction of null modes the real graph beats at ``|z| > z_threshold``.

    ``panel`` may be a :func:`null_panel` result or its ``rows``; when it is
    ``None`` the panel is computed first.  ``assay`` filters the rows to one
    assay arm (default: all rows).  1.0 means the connectome mattered against
    every null; 0.0 means a shuffled graph reproduced the drug effect.
    """
    if panel is None:
        panel = null_panel(compound, conc_M, modes=modes, n=n, seed=seed, **kw)
    rows = panel["rows"] if isinstance(panel, dict) else list(panel)
    if assay is not None:
        rows = [r for r in rows if r.get("assay") == assay]
    scored = [r for r in rows if r.get("z") is not None]
    hits = [r for r in scored if abs(float(r["z"])) > float(z_threshold)]
    warnings = list(BASE_WARNINGS) + [
        "The score counts modes, not biology: it asks how many degradations of "
        "the MaleCNS wiring fail to reproduce this model's drug effect."
    ]
    if len(scored) < len(rows):
        warnings.append(f"{len(rows) - len(scored)} of {len(rows)} rows had an undefined z.")
    return {
        "score": (len(hits) / len(scored)) if scored else None,
        "n_rows": len(rows),
        "n_scored": len(scored),
        "n_beaten": len(hits),
        "z_threshold": float(z_threshold),
        "per_row": [
            {
                "assay": r.get("assay"),
                "readout": r.get("readout"),
                "mode": r.get("mode"),
                "z": r.get("z"),
                "beaten": (r.get("z") is not None and abs(float(r["z"])) > float(z_threshold)),
            }
            for r in rows
        ],
        "label": "model_derived",
        "warnings": warnings,
    }
