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

Statistics
----------
The headline number is the **empirical two-sided permutation probability**

    p = (k + 1) / (n_ok + 1),   k = #{ |null_i - null_mean| >= |real - null_mean| }

reported together with its **resolution** ``1 / (n_ok + 1)``: with 50 shuffles
no effect can be more extreme than p = 0.02 whatever the z-score, which is why
the principal analyses are run at n = 500-1000.  ``z`` is reported second, as a
standardised distance, and is *not* a probability: on degraded graphs the null
is often far from normal and a |z| of 40 only means the null spread collapsed.

Every ``null_distribution`` also carries a **convergence check**: the statistic
is recomputed on the first n/4 and n/2 shuffles (the draws are i.i.d., so a
prefix is itself a valid smaller permutation sample) and the result is flagged
``stabilised`` when p moved by less than ``tol`` (default 0.02, i.e. two
percentage points) between every checkpoint and the full sample, and the
verdict at ``alpha`` never changed.  ``n_stabilised`` is the smallest
checkpoint from which that holds, which is the answer to "how many permutations
are needed before the topology conclusion stops moving?".

Performance
-----------
High-permutation runs go through a small in-module rate engine
(:class:`_GraphState`) instead of the notebook assay wrappers:

* the graph is parsed into index/weight arrays **once** per run and each
  shuffle only permutes those arrays (no dict copies, no JSON round-trip, no
  temp files, no re-parsing);
* the row-normalised signed matrix is kept in edge-list form and the recurrent
  step is a ``bincount`` scatter instead of a dense ``N x N`` mat-vec;
* the pharmacological gains are concentration-dependent only, so they are
  computed once and reused by every shuffle;
* only the runs the requested readout actually needs are executed (two for a
  rate readout, four for the bitter-veto ratio) instead of a full notebook.

The engine reproduces ``flylab.circuit.rate.RateNetwork`` to ~1e-13 Hz (it is
the same arithmetic in a different summation order); ``engine="assay"`` forces
the original notebook path, and unsupported combinations (the LIF ``spiking``
assay, an unusual readout) fall back to it automatically.  ``n_jobs > 1``
spreads the shuffles over processes; because shuffle *i* is seeded by
``(seed, i)`` alone, the parallel result is **identical** to the serial one.

Nothing in here is fitted to animal data.  A p-value against a shuffled
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
    "DEFAULT_GRAPH",
    "FAST_READOUTS",
    "DEFAULT_TOL",
    "DEFAULT_ALPHA",
    "DEFAULT_CHECKPOINTS",
    "shuffle_graph",
    "null_distribution",
    "null_panel",
    "connectome_information_score",
    "drug_effect",
    "fast_drug_effect",
    "readout_value",
    "vehicle_value",
    "empirical_p",
    "p_resolution",
    "convergence_report",
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

#: readouts the in-module rate engine can produce, per assay.  Anything else
#: (and the LIF ``spiking`` assay) falls back to the notebook assay path.
FAST_READOUTS: dict[str, tuple[str, ...]] = {
    "subgraph": ("mean_hz", "max_hz", "mn9_hz", "dnp01_hz"),
    "taste_map": (
        "mn9_hz",
        "mn9_sugar_hz",
        "mn9_sugar_bitter_hz",
        "bitter_veto_ratio",
        "mean_hz",
    ),
}

#: assay keyword arguments the fast engine understands; anything else forces
#: the notebook path (``engine="lif"`` in particular).
_FAST_KW: dict[str, tuple[str, ...]] = {
    "subgraph": ("drive_hz", "steps", "drive", "library", "rule_overrides"),
    "taste_map": (
        "sugar_hz",
        "bitter_hz",
        "steps",
        "library",
        "rule_overrides",
        "seed",
        "engine",
    ),
}

#: default convergence tolerance on the empirical p (two percentage points)
DEFAULT_TOL = 0.02
#: default significance threshold used for the "verdict stable" check
DEFAULT_ALPHA = 0.05
#: fractions of ``n`` at which the statistic is recomputed
DEFAULT_CHECKPOINTS: tuple[float, ...] = (0.25, 0.5, 1.0)

BASE_WARNINGS = [
    "Null models destroy the MaleCNS wiring on purpose: they are controls, not "
    "biology. Nothing here is a measurement from a living fly.",
    "The headline statistic is the empirical two-sided permutation p, "
    "(k+1)/(n+1); with n shuffles its resolution (and floor) is 1/(n+1), so a "
    "p at the floor means 'not resolvable with this n', not 'p = 0'.",
    "z is reported second and is not a probability: the null distribution on a "
    "degraded graph is usually not normal, and a large |z| can simply mean the "
    "shuffles collapsed onto one value.",
    "model_derived: these numbers describe this simulation and its teaching "
    "EC50 library, not a measured drug effect.",
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


def p_resolution(n_ok: int) -> float:
    """Resolution (and floor) of an empirical permutation p with ``n_ok`` draws."""
    return 1.0 / (int(n_ok) + 1)


def empirical_p(
    real: float | None,
    nulls: Sequence[float],
    centre: float | None = None,
) -> dict[str, Any]:
    """Two-sided empirical permutation probability ``(k + 1) / (n + 1)``.

    ``k`` counts null draws at least as far from ``centre`` (the null mean by
    default) as the real effect is.  The ``+1`` in both places counts the
    observed statistic itself, which is what keeps the estimate valid (and its
    floor honest) for any ``n``.
    """
    ok = np.asarray([v for v in nulls if v is not None], dtype=float)
    n_ok = int(ok.size)
    if real is None or n_ok == 0:
        return {"p": None, "k": None, "n": n_ok, "resolution": p_resolution(n_ok)}
    c = float(ok.mean()) if centre is None else float(centre)
    k = int(np.sum(np.abs(ok - c) >= abs(float(real) - c)))
    return {
        "p": float((k + 1) / (n_ok + 1)),
        "k": k,
        "n": n_ok,
        "resolution": p_resolution(n_ok),
    }


def convergence_report(
    real: float | None,
    nulls: Sequence[float | None],
    tol: float = DEFAULT_TOL,
    alpha: float = DEFAULT_ALPHA,
    checkpoints: Sequence[float] = DEFAULT_CHECKPOINTS,
) -> dict[str, Any]:
    """Has the permutation p stopped moving?

    The statistic is recomputed on the first ``f * n`` shuffles for each
    fraction ``f`` in ``checkpoints`` (default n/4, n/2, n).  The draws are
    i.i.d. given the seed, so a prefix is itself a valid smaller permutation
    sample; no extra circuit runs are needed.

    ``stabilised`` is True when, for every checkpoint, ``|p_k - p_final| <=
    tol`` **and** the verdict ``p <= alpha`` is the same as at the full sample.
    ``n_stabilised`` is the smallest checkpoint n (never the full sample
    itself) from which that holds for every later checkpoint - the number of
    permutations this conclusion actually needed.  ``None`` means not even the
    half-sample agreed, so the answer still moves at this n.
    """
    draws = list(nulls)
    total = len(draws)
    rows: list[dict[str, Any]] = []
    for frac in checkpoints:
        n_k = int(round(float(frac) * total))
        if n_k < 1:
            continue
        prefix = draws[:n_k]
        ok = np.asarray([v for v in prefix if v is not None], dtype=float)
        if ok.size == 0 or real is None:
            rows.append(
                {
                    "fraction": float(frac),
                    "n": n_k,
                    "n_ok": int(ok.size),
                    "p": None,
                    "p_resolution": p_resolution(int(ok.size)),
                    "z": None,
                    "beats_null": None,
                }
            )
            continue
        mean = float(ok.mean())
        sd = float(ok.std(ddof=1)) if ok.size > 1 else 0.0
        pp = empirical_p(real, ok, centre=mean)
        rows.append(
            {
                "fraction": float(frac),
                "n": n_k,
                "n_ok": int(ok.size),
                "p": pp["p"],
                "p_resolution": pp["resolution"],
                "z": (float((real - mean) / sd) if sd > 0 else None),
                "beats_null": (pp["p"] is not None and pp["p"] <= float(alpha)),
            }
        )
    if not rows or rows[-1]["p"] is None:
        return {
            "tolerance": float(tol),
            "alpha": float(alpha),
            "checkpoints": rows,
            "p_max_shift": None,
            "stabilised": False,
            "n_stabilised": None,
            "verdict_stable": None,
            "resolution_limited": None,
            "note": "no usable null draws; convergence could not be assessed.",
        }
    final = rows[-1]
    shifts = [abs(r["p"] - final["p"]) for r in rows[:-1] if r["p"] is not None]
    verdicts = [r["beats_null"] for r in rows if r["beats_null"] is not None]
    verdict_stable = len(set(verdicts)) <= 1
    stabilised = bool(all(s <= float(tol) for s in shifts) and verdict_stable)
    # The full sample trivially agrees with itself, so it is not a candidate:
    # "stabilised at n" may only be claimed when an *earlier* checkpoint and
    # everything after it already agreed with the final answer.
    n_stab: int | None = None
    for i, row in enumerate(rows[:-1]):
        if row["p"] is None:
            continue
        if all(
            r["p"] is not None
            and abs(r["p"] - final["p"]) <= float(tol)
            and (r["beats_null"] == final["beats_null"])
            for r in rows[i:]
        ):
            n_stab = int(row["n"])
            break
    # A p pinned at its own floor has not converged for a different reason:
    # it is limited by the permutation resolution, not by sampling noise, and
    # the cure is more shuffles rather than a wider tolerance.
    resolution_limited = bool(final["p"] <= final["p_resolution"] + 1e-12)
    note = (
        "p recomputed on nested prefixes of the same shuffle stream; "
        f"stabilised means every checkpoint sat within {float(tol)} of the "
        "full-sample p and never changed the verdict at "
        f"alpha = {float(alpha)}."
    )
    if resolution_limited:
        note += (
            " p is at its resolution floor 1/(n+1) at the full sample, so it "
            "is still falling with n: read it as 'p <= floor' and raise n to "
            "resolve it."
        )
    return {
        "tolerance": float(tol),
        "alpha": float(alpha),
        "checkpoints": rows,
        "p_max_shift": (max(shifts) if shifts else 0.0),
        "stabilised": stabilised,
        "n_stabilised": n_stab,
        "verdict_stable": verdict_stable,
        "resolution_limited": resolution_limited,
        "note": note,
    }


# --------------------------------------------------------------------------
# array-level graph state and rate engine
# --------------------------------------------------------------------------
class _GraphState:
    """Index/weight arrays for one graph, plus the rate runtime that uses them.

    This is the null-model fast path: parse the connectome once, then permute
    arrays instead of rebuilding dictionaries and dense matrices per shuffle.
    The arithmetic is the same as :class:`flylab.circuit.rate.RateNetwork` -
    signed presynaptic weights, rows normalised by ``max(sum |w|, 1)``, per-
    presynaptic-node gain scaling, leaky update, clipped at ``[0, r_max]`` -
    written as an edge-list scatter rather than a dense mat-vec.
    """

    __slots__ = (
        "n",
        "body_ids",
        "index",
        "nt",
        "pre",
        "post",
        "w",
        "seeds",
        "meta",
        "_sw",
        "_denom",
    )

    def __init__(
        self,
        n: int,
        body_ids: np.ndarray,
        index: dict[int, int],
        nt: np.ndarray,
        pre: np.ndarray,
        post: np.ndarray,
        w: np.ndarray,
        seeds: dict[str, list[int]],
        meta: dict[str, Any] | None = None,
    ) -> None:
        self.n = int(n)
        self.body_ids = body_ids
        self.index = index
        self.nt = nt
        self.pre = pre
        self.post = post
        self.w = w
        self.seeds = seeds
        self.meta = dict(meta or {})
        self._sw: np.ndarray | None = None
        self._denom: np.ndarray | None = None

    # -- construction -----------------------------------------------------
    @classmethod
    def from_graph(cls, graph: dict[str, Any]) -> "_GraphState":
        nodes = graph["nodes"]
        body_ids = np.array([int(nd["bodyId"]) for nd in nodes], dtype=np.int64)
        index = {int(b): i for i, b in enumerate(body_ids.tolist())}
        # the raw label is kept (``None`` included, as the graph JSON has it):
        # ``SIGN``/``NT_GAIN_KEY`` miss it exactly as they miss "unclear", and
        # ``shuffle_graph`` must hand back the histogram it was given.
        nt = np.empty(len(nodes), dtype=object)
        for i, nd in enumerate(nodes):
            nt[i] = nd.get("consensus_nt")
        pre: list[int] = []
        post: list[int] = []
        w: list[float] = []
        for e in graph["edges"]:
            i = index.get(int(e["pre"]))
            j = index.get(int(e["post"]))
            if i is None or j is None:
                continue
            pre.append(i)
            post.append(j)
            w.append(float(e["weight"]))
        seeds = {k: [int(v) for v in vals] for k, vals in (graph.get("seeds") or {}).items()}
        return cls(
            n=len(nodes),
            body_ids=body_ids,
            index=index,
            nt=nt,
            pre=np.asarray(pre, dtype=np.int64),
            post=np.asarray(post, dtype=np.int64),
            w=np.asarray(w, dtype=float),
            seeds=seeds,
        )

    def replace(
        self,
        nt: np.ndarray | None = None,
        pre: np.ndarray | None = None,
        post: np.ndarray | None = None,
        w: np.ndarray | None = None,
        meta: dict[str, Any] | None = None,
    ) -> "_GraphState":
        return _GraphState(
            n=self.n,
            body_ids=self.body_ids,
            index=self.index,
            nt=self.nt if nt is None else nt,
            pre=self.pre if pre is None else pre,
            post=self.post if post is None else post,
            w=self.w if w is None else w,
            seeds=self.seeds,
            meta=meta,
        )

    # -- rate runtime -----------------------------------------------------
    def _signed(self) -> tuple[np.ndarray, np.ndarray]:
        """``(signed edge weights, per-postsynaptic normaliser)``, memoised."""
        if self._sw is None or self._denom is None:
            from flylab.circuit.rate import SIGN

            sign = np.array([SIGN.get(x, 0.0) for x in self.nt.tolist()], dtype=float)
            sw = sign[self.pre] * self.w
            denom = np.maximum(
                np.bincount(self.post, weights=np.abs(sw), minlength=self.n), 1.0
            )
            self._sw, self._denom = sw, denom
        return self._sw, self._denom

    def scale_vector(self, gains: dict[str, float]) -> tuple[np.ndarray, dict[str, float]]:
        from flylab.circuit.rate import DEFAULT_GAINS, NT_GAIN_KEY

        g = dict(DEFAULT_GAINS)
        g.update({k: float(v) for k, v in (gains or {}).items() if k in DEFAULT_GAINS})
        out = np.ones(self.n, dtype=float)
        for i, nt in enumerate(self.nt.tolist()):
            key = NT_GAIN_KEY.get(nt)
            if key is None:
                continue
            val = g[key]
            if nt == "acetylcholine":
                val *= g["ach_tone"]
            out[i] = val
        return out, g

    def drive_vector(self, mapping: dict[int, float] | None) -> np.ndarray:
        d = np.zeros(self.n, dtype=float)
        for body_id, hz in (mapping or {}).items():
            i = self.index.get(int(body_id))
            if i is not None:
                d[i] = float(hz)
        return d

    def seed_indices(self, types: Iterable[str]) -> list[int]:
        out: list[int] = []
        for t in types:
            for b in self.seeds.get(t, []):
                i = self.index.get(int(b))
                if i is not None and i not in out:
                    out.append(i)
        return out

    def run(
        self,
        drive: np.ndarray,
        gains: dict[str, float] | None,
        steps: int = 80,
        alpha: float = 0.3,
        r_max: float = 300.0,
    ) -> np.ndarray:
        sw, denom = self._signed()
        scale, g = self.scale_vector(gains or {})
        weff = sw * scale[self.pre] / denom[self.post]
        g_nav = g["g_nav"]
        pre, post, n = self.pre, self.post, self.n
        r = np.zeros(n, dtype=float)
        for _ in range(int(steps)):
            rec = np.bincount(post, weights=weff * r[pre], minlength=n)
            r = np.clip((1.0 - alpha) * r + alpha * (g_nav * (drive + rec)), 0.0, r_max)
        return r


_STATE_CACHE: dict[str, _GraphState] = {}


def _state_for(graph_name: str | None) -> _GraphState:
    """Memoised :class:`_GraphState` for a committed graph (per process)."""
    from flylab.circuit.rate import load_graph, resolve_graph

    key = str(resolve_graph(graph_name).resolve())
    if key not in _STATE_CACHE:
        _STATE_CACHE[key] = _GraphState.from_graph(load_graph(key))
    return _STATE_CACHE[key]


# --------------------------------------------------------------------------
# graph shuffles (array level; the dict API below is a thin wrapper)
# --------------------------------------------------------------------------
def _double_edge_swap_arrays(
    pre: np.ndarray,
    post: np.ndarray,
    rng: np.random.Generator,
    n_swaps: int,
    max_tries_per_swap: int = 20,
) -> tuple[np.ndarray, int, int]:
    """Degree-preserving double-edge swaps on index arrays.

    ``(a->b, c->d)`` becomes ``(a->d, c->b)``: node ``a`` and ``c`` keep their
    out-degree, ``b`` and ``d`` keep their in-degree.  The weight travels with
    its presynaptic node.  Swaps that would make a self-loop or a duplicate
    edge are rejected.  Endpoint pairs are packed into one int64 key, which is
    what makes this cheap enough for a thousand shuffles.
    """
    m = int(pre.size)
    if m < 2:
        return post.copy(), 0, 0
    prel = pre.tolist()
    postl = post.tolist()
    shift = 1 << 32
    present = {a * shift + b for a, b in zip(prel, postl)}
    done = 0
    tries = 0
    budget = int(n_swaps) * int(max_tries_per_swap)
    chunk = max(int(n_swaps), 1024)
    while done < n_swaps and tries < budget:
        picks = rng.integers(0, m, size=(min(chunk, budget - tries), 2))
        for i, j in picks.tolist():
            tries += 1
            if i == j:
                continue
            a, b = prel[i], postl[i]
            c, d = prel[j], postl[j]
            if a == d or c == b:
                continue
            k1 = a * shift + d
            k2 = c * shift + b
            if k1 in present or k2 in present:
                continue
            present.discard(a * shift + b)
            present.discard(c * shift + d)
            present.add(k1)
            present.add(k2)
            postl[i], postl[j] = d, b
            done += 1
            if done >= n_swaps:
                break
    return np.asarray(postl, dtype=np.int64), done, tries


def _erdos_renyi_arrays(
    n_nodes: int,
    weights: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Same N and E; endpoints uniform, weights drawn from the empirical list."""
    m = int(weights.size)
    drawn = rng.choice(weights, size=m, replace=True) if m else np.asarray([], dtype=float)
    taken: set[int] = set()
    pre: list[int] = []
    post: list[int] = []
    w: list[float] = []
    shift = 1 << 32
    tries = 0
    while len(pre) < m and tries < 200 * max(m, 1):
        need = m - len(pre)
        aa = rng.integers(0, n_nodes, size=need)
        bb = rng.integers(0, n_nodes, size=need)
        for a, b in zip(aa.tolist(), bb.tolist()):
            tries += 1
            if a == b:
                continue
            key = a * shift + b
            if key in taken:
                continue
            taken.add(key)
            w.append(float(drawn[len(pre)]))
            pre.append(a)
            post.append(b)
            if len(pre) >= m:
                break
    return (
        np.asarray(pre, dtype=np.int64),
        np.asarray(post, dtype=np.int64),
        np.asarray(w, dtype=float),
    )


def _shuffle_state(state: _GraphState, mode: str, rng: Any = 0) -> _GraphState:
    """One degraded copy of ``state``; see :func:`shuffle_graph` for the modes."""
    if mode not in MODES:
        raise ValueError(f"unknown null mode {mode!r}; expected one of {MODES}")
    r = _rng(rng)
    meta: dict[str, Any] = {"mode": mode}
    if mode == "sign_permute":
        perm = r.permutation(state.n)
        meta["permuted_nodes"] = state.n
        return state.replace(nt=state.nt[perm], meta=meta)
    if mode == "weight_permute":
        perm = r.permutation(int(state.w.size))
        meta["permuted_edges"] = int(state.w.size)
        return state.replace(w=state.w[perm], meta=meta)
    if mode == "rewire_degree_preserving":
        target = 10 * int(state.pre.size)
        post, done, tries = _double_edge_swap_arrays(state.pre, state.post, r, n_swaps=target)
        meta.update({"swaps_requested": target, "swaps_done": done, "attempts": tries})
        return state.replace(post=post, meta=meta)
    pre, post, w = _erdos_renyi_arrays(state.n, state.w, r)
    meta["resampled_edges"] = int(pre.size)
    return state.replace(pre=pre, post=post, w=w, meta=meta)


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

    This is the dict-level API, kept for the server, the browser bridge and the
    tests.  The null loops run the same shuffles on arrays
    (:func:`_shuffle_state`) and never materialise these dictionaries; the two
    paths consume the RNG identically and give the same graph for a given seed.
    """
    if mode not in MODES:
        raise ValueError(f"unknown null mode {mode!r}; expected one of {MODES}")
    base = _GraphState.from_graph(graph)
    st = _shuffle_state(base, mode, rng)

    nodes = [dict(nd) for nd in graph["nodes"]]
    if mode == "sign_permute":
        for nd, nt in zip(nodes, st.nt.tolist()):
            nd["consensus_nt"] = nt

    body_ids = base.body_ids.tolist()
    if mode == "erdos_renyi":
        edges = [
            {"pre": body_ids[p], "post": body_ids[q], "weight": float(wt)}
            for p, q, wt in zip(st.pre.tolist(), st.post.tolist(), st.w.tolist())
        ]
    else:
        edges = []
        for e, p, q, wt in zip(
            graph["edges"], st.pre.tolist(), st.post.tolist(), st.w.tolist()
        ):
            row = dict(e)
            row["pre"], row["post"], row["weight"] = body_ids[p], body_ids[q], float(wt)
            edges.append(row)

    out = dict(graph)
    out["nodes"] = nodes
    out["edges"] = edges
    out["n_nodes"] = len(nodes)
    out["n_edges"] = len(edges)
    out["seeds"] = {k: list(v) for k, v in (graph.get("seeds") or {}).items()}
    out["null_model"] = dict(st.meta)
    return out


# --------------------------------------------------------------------------
# running an assay on a shuffled graph (notebook path)
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
    """``treated - vehicle`` for one readout, both on the *same* graph.

    This is the notebook path: it runs the real assay wrapper and is what the
    selectivity layer and the ``engine="assay"`` null runs use.  For the null
    loops :func:`fast_drug_effect` computes the same contrast directly on the
    rate engine.
    """
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
# the fast (engine) path
# --------------------------------------------------------------------------
def _fast_supported(assay: str, readout: str, kw: dict[str, Any]) -> tuple[bool, str]:
    if assay not in FAST_READOUTS:
        return False, f"assay {assay!r} has no rate-engine path (LIF is run through the assay)"
    if readout not in FAST_READOUTS[assay]:
        return False, f"readout {readout!r} is not one the rate engine produces"
    allowed = _FAST_KW[assay]
    extra = [k for k in kw if k not in allowed]
    if extra:
        return False, f"assay keyword(s) {sorted(extra)} are not supported by the rate engine"
    if assay == "taste_map" and str(kw.get("engine", "rate")) != "rate":
        return False, "taste_map engine='lif' is run through the assay"
    return True, ""


def _plan(state: _GraphState, assay: str, readout: str, gains: dict[str, float], kw: dict[str, Any]):
    """Drive vectors and the run list a readout needs (graph-invariant)."""
    from flylab.circuit.rate import DEFAULT_GAINS

    steps = int(kw.get("steps", 80))
    if assay == "subgraph":
        drive_map = kw.get("drive")
        if drive_map:
            drive = state.drive_vector(dict(drive_map))
        else:
            hz = float(kw.get("drive_hz", 40.0))
            drive = np.zeros(state.n, dtype=float)
            for i in state.seed_indices(list(state.seeds)):
                drive[i] = hz
        return {
            "assay": assay,
            "readout": readout,
            "steps": steps,
            "gains": dict(gains),
            "vehicle_gains": dict(DEFAULT_GAINS),
            "drive": drive,
            "mn9": state.seed_indices(["MN9"]),
            "dnp01": state.seed_indices(["DNp01"]),
        }
    # taste_map
    from flylab.assays.taste_map import BITTER_TYPES, SWEET_TYPES

    sugar_hz = float(kw.get("sugar_hz", 150.0))
    bitter_hz = float(kw.get("bitter_hz", 0.0) or 0.0)
    bitter_probe = bitter_hz if bitter_hz > 0 else sugar_hz
    sweet = state.seed_indices(SWEET_TYPES)
    bitter = state.seed_indices(BITTER_TYPES)
    d_sugar = np.zeros(state.n, dtype=float)
    d_sugar[sweet] = sugar_hz
    d_both = d_sugar.copy()
    d_both[bitter] = bitter_probe
    return {
        "assay": assay,
        "readout": readout,
        "steps": steps,
        "gains": dict(gains),
        "vehicle_gains": dict(DEFAULT_GAINS),
        "d_sugar": d_sugar,
        "d_both": d_both,
        "mn9": state.seed_indices(["MN9"]),
        "needs_both": readout in ("mn9_sugar_bitter_hz", "bitter_veto_ratio"),
    }


def _mean_at(rates: np.ndarray, idx: Sequence[int]) -> float | None:
    if not len(idx):
        return None
    return float(np.mean(rates[list(idx)]))


def _fast_effect(state: _GraphState, plan: dict[str, Any]) -> dict[str, Any]:
    """``treated - vehicle`` for one readout, straight off the rate engine."""
    steps = plan["steps"]
    readout = plan["readout"]
    if plan["assay"] == "subgraph":
        r_t = state.run(plan["drive"], plan["gains"], steps=steps)
        r_v = state.run(plan["drive"], plan["vehicle_gains"], steps=steps)

        def pick(r: np.ndarray) -> float | None:
            if readout == "mean_hz":
                return float(r.mean())
            if readout == "max_hz":
                return float(r.max())
            if readout == "mn9_hz":
                return _mean_at(r, plan["mn9"])
            return _mean_at(r, plan["dnp01"])

        treated, vehicle = pick(r_t), pick(r_v)
    else:
        r_ts = state.run(plan["d_sugar"], plan["gains"], steps=steps)
        r_vs = state.run(plan["d_sugar"], plan["vehicle_gains"], steps=steps)
        r_tb = r_vb = None
        if plan["needs_both"]:
            r_tb = state.run(plan["d_both"], plan["gains"], steps=steps)
            r_vb = state.run(plan["d_both"], plan["vehicle_gains"], steps=steps)

        def pick2(r_s: np.ndarray, r_b: np.ndarray | None) -> float | None:
            if readout == "mean_hz":
                return float(np.mean(r_s))
            if readout in ("mn9_hz", "mn9_sugar_hz"):
                return _mean_at(r_s, plan["mn9"])
            both = _mean_at(r_b, plan["mn9"]) if r_b is not None else None
            if readout == "mn9_sugar_bitter_hz":
                return both
            sugar = _mean_at(r_s, plan["mn9"])
            if sugar is None or both is None or not (sugar > 1e-6):
                return None
            return float(both / sugar)

        treated, vehicle = pick2(r_ts, r_tb), pick2(r_vs, r_vb)
    effect = None if (treated is None or vehicle is None) else float(treated - vehicle)
    return {"treated": treated, "vehicle": vehicle, "effect": effect}


def fast_drug_effect(
    assay: str,
    compound: str | None,
    conc_M: float,
    readout: str,
    graph: str | None = None,
    state: _GraphState | None = None,
    gains: dict[str, float] | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """:func:`drug_effect` computed on the in-module rate engine.

    Same contrast, same arithmetic, ~5-100x faster and with no temp files; it
    agrees with the notebook path to ~1e-13 Hz (different summation order).
    ``state`` lets a caller reuse an already-shuffled graph.
    """
    from flylab.circuit.rate import compute_gains

    st = state if state is not None else _state_for(graph)
    if gains is None:
        gains, _ = compute_gains(
            compound,
            conc_M,
            library=kw.get("library"),
            rule_overrides=kw.get("rule_overrides"),
        )
    plan = _plan(st, assay, readout, gains, kw)
    return _fast_effect(st, plan)


# --------------------------------------------------------------------------
# the shuffle loop, serial and parallel
# --------------------------------------------------------------------------
def _effects_for_indices(spec: dict[str, Any], indices: Sequence[int]) -> list[float | None]:
    """Drug effect on shuffle ``i`` of ``spec`` for each ``i`` in ``indices``.

    Shuffle ``i`` depends only on ``(seed, i)``, so this is independent of how
    the indices are partitioned - which is what makes the parallel path give
    bit-identical results to the serial one.
    """
    mode = spec["mode"]
    seed = int(spec["seed"])
    kw = dict(spec.get("kw") or {})
    out: list[float | None] = []
    if spec["engine"] == "fast":
        base = _state_for(spec["graph"])
        plan = _plan(base, spec["assay"], spec["readout"], spec["gains"], kw)
        for i in indices:
            st = _shuffle_state(base, mode, np.random.default_rng([seed, int(i)]))
            out.append(_fast_effect(st, plan)["effect"])
        return out
    from flylab.circuit.rate import load_graph

    graph = load_graph(spec["graph"])
    for i in indices:
        shuffled = shuffle_graph(graph, mode, np.random.default_rng([seed, int(i)]))
        res = drug_effect(
            spec["assay"],
            spec["compound"],
            spec["conc_M"],
            spec["readout"],
            graph_obj=shuffled,
            graph=spec["graph"],
            **kw,
        )
        out.append(res["effect"])
    return out


def _shuffle_block(payload: tuple[dict[str, Any], list[int]]) -> list[float | None]:
    """Module-level worker (must be importable for ProcessPoolExecutor)."""
    spec, indices = payload
    return _effects_for_indices(spec, indices)


def _run_shuffles(
    spec: dict[str, Any],
    n: int,
    n_jobs: int,
    warnings: list[str],
) -> list[float | None]:
    indices = list(range(int(n)))
    if int(n_jobs) <= 1 or len(indices) < 2:
        return _effects_for_indices(spec, indices)
    jobs = min(int(n_jobs), len(indices))
    size = (len(indices) + jobs - 1) // jobs
    blocks = [indices[i : i + size] for i in range(0, len(indices), size)]
    try:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=jobs) as pool:
            chunks = list(pool.map(_shuffle_block, [(spec, b) for b in blocks]))
        out: list[float | None] = []
        for c in chunks:
            out.extend(c)
        warnings.append(
            f"shuffles were run on {jobs} processes; shuffle i is seeded by "
            "(seed, i) alone, so the draws are identical to a serial run."
        )
        return out
    except Exception as exc:  # pragma: no cover - platform dependent
        warnings.append(f"parallel execution failed ({exc!r}); fell back to serial.")
        return _effects_for_indices(spec, indices)


# --------------------------------------------------------------------------
# null distribution
# --------------------------------------------------------------------------
def null_distribution(
    assay: str = "subgraph",
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    mode: str = "sign_permute",
    n: int = 200,
    seed: int = 0,
    readout: str = "auto",
    graph: str | None = None,
    n_jobs: int = 1,
    engine: str = "auto",
    tol: float = DEFAULT_TOL,
    alpha: float = DEFAULT_ALPHA,
    checkpoints: Sequence[float] = DEFAULT_CHECKPOINTS,
    **kw: Any,
) -> dict[str, Any]:
    """Drug effect on the real graph vs ``n`` shuffled graphs.

    The effect is always ``treated - vehicle`` for ``readout``, with both arms
    run on the *same* graph, so a shuffle can only change the effect through
    the wiring.  ``readout="auto"`` uses :data:`DEFAULT_READOUT`.

    Args:
        n: number of shuffles.  The empirical p cannot be smaller than
            ``1/(n+1)``: 50 shuffles floor it at 0.02, so headline claims are
            run at 500-1000.
        n_jobs: 1 (default) runs in-process and is deterministic; >1 spreads
            the shuffles over that many processes and returns **identical**
            numbers for the same seed.
        engine: ``"auto"`` (default) uses the in-module rate engine when the
            assay/readout/kwargs allow it and the notebook assay otherwise;
            ``"fast"`` demands the engine, ``"assay"`` forces the notebook path.
        tol, alpha, checkpoints: convergence check (see
            :func:`convergence_report`).

    Returns:
        dict with the empirical two-sided ``p_two_sided`` **first**, its
        ``p_resolution`` ``1/(n_ok+1)``, the ``stabilised`` flag and the
        ``convergence`` trace, then ``z`` and the null moments, the raw
        ``null_effects``, ``runtime_s``, ``label`` and ``warnings``.
    """
    if assay not in ASSAYS:
        raise ValueError(f"unknown assay {assay!r}; expected one of {ASSAYS}")
    if mode not in MODES:
        raise ValueError(f"unknown null mode {mode!r}; expected one of {MODES}")
    if engine not in ("auto", "fast", "assay"):
        raise ValueError(f"unknown engine {engine!r}; expected auto, fast or assay")
    if readout == "auto":
        readout = DEFAULT_READOUT[assay]
    t0 = time.perf_counter()
    graph_name = graph if graph is not None else DEFAULT_GRAPH[assay]
    warnings = list(BASE_WARNINGS) + [MODE_NOTES[mode]]

    ok_fast, why = _fast_supported(assay, readout, kw)
    if engine == "fast" and not ok_fast:
        raise ValueError(f"engine='fast' is not available here: {why}")
    use_fast = ok_fast if engine == "auto" else (engine == "fast")
    if engine == "auto" and not ok_fast:
        warnings.append(f"rate engine unavailable ({why}); ran the notebook assay path.")

    if use_fast:
        from flylab.circuit.rate import compute_gains

        gains, _occ = compute_gains(
            compound,
            conc_M,
            library=kw.get("library"),
            rule_overrides=kw.get("rule_overrides"),
        )
        base = _state_for(graph_name)
        plan = _plan(base, assay, readout, gains, kw)
        real = _fast_effect(base, plan)
    else:
        gains = None
        real = drug_effect(assay, compound, conc_M, readout, graph=graph_name, **kw)

    spec = {
        "assay": assay,
        "compound": compound,
        "conc_M": float(conc_M),
        "readout": readout,
        "graph": graph_name,
        "mode": mode,
        "seed": int(seed),
        "engine": "fast" if use_fast else "assay",
        "gains": dict(gains) if gains else None,
        "kw": dict(kw),
    }
    nulls = _run_shuffles(spec, int(n), int(n_jobs), warnings)

    ok = np.array([v for v in nulls if v is not None], dtype=float)
    n_ok = int(ok.size)
    null_mean = float(ok.mean()) if n_ok else None
    null_median = float(np.median(ok)) if n_ok else None
    null_sd = float(ok.std(ddof=1)) if n_ok > 1 else 0.0
    real_effect = real["effect"]

    z: float | None = None
    p: float | None = None
    k: int | None = None
    if real_effect is None:
        warnings.append(
            f"readout {readout!r} is undefined on the real graph for this "
            "compound; no permutation p and no z-score could be formed."
        )
    elif n_ok == 0:
        warnings.append("every shuffled graph gave an undefined readout; no null.")
    else:
        pp = empirical_p(real_effect, ok, centre=null_mean)
        p, k = pp["p"], pp["k"]
        if null_sd > 0:
            z = float((real_effect - null_mean) / null_sd)
        else:
            warnings.append(
                "the null distribution has zero spread; z is undefined (the "
                "shuffles all gave the same effect)."
            )
        if p is not None and p <= p_resolution(n_ok) + 1e-12:
            warnings.append(
                f"p is at its resolution floor {p_resolution(n_ok):.4f}: with "
                f"{n_ok} shuffles this is the smallest p obtainable, so report "
                "it as 'p <= floor', not as an exact value."
            )
        if null_sd > 10.0 * max(abs(real_effect), 1e-12):
            warnings.append(
                "the null distribution is heavy-tailed (sd >> |real effect|): "
                "on degraded graphs a ratio readout can blow up, so z and p "
                "are not interpretable here."
            )
    if n_ok < int(n):
        warnings.append(f"{int(n) - n_ok} of {int(n)} shuffles gave an undefined readout.")

    conv = convergence_report(real_effect, nulls, tol=tol, alpha=alpha, checkpoints=checkpoints)
    if not conv["stabilised"] and n_ok and p is not None:
        shift = conv["p_max_shift"]
        moved = f"{shift:.3f}" if shift is not None else "an unknown amount"
        verdict = (
            "the verdict at alpha never changed, so only the precision of p is "
            "at stake"
            if conv.get("verdict_stable")
            else "the verdict at alpha also changed between checkpoints"
        )
        warnings.append(
            f"the permutation p had not stabilised at this n (it moved by {moved} "
            f"between checkpoints, tolerance {float(tol):.3f}); {verdict}. "
            "Increase n before quoting the value."
        )

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
        "engine": spec["engine"],
        "n_jobs": int(n_jobs),
        "real_effect": real_effect,
        # empirical permutation probability first, z second (reviewer request)
        "p_two_sided": p,
        "p_resolution": p_resolution(n_ok),
        "p_k": k,
        "alpha": float(alpha),
        "beats_null": (p is not None and p <= float(alpha)),
        "stabilised": bool(conv["stabilised"]),
        "n_stabilised": conv["n_stabilised"],
        "resolution_limited": bool(conv.get("resolution_limited")),
        "verdict_stable": conv.get("verdict_stable"),
        "convergence": conv,
        "z": z,
        "null_mean": null_mean,
        "null_median": null_median,
        "null_sd": null_sd,
        "real_treated": real["treated"],
        "real_vehicle": real["vehicle"],
        "null_effects": nulls,
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
    n: int = 200,
    seed: int = 0,
    n_taste: int | None = None,
    include_taste_map: bool = True,
    n_jobs: int = 1,
    **kw: Any,
) -> dict[str, Any]:
    """All null modes for the ``subgraph`` mean rate and the taste-map veto.

    The taste-map arm is the bitter veto ratio (``mn9_sugar_bitter / mn9_sugar``)
    on the 1841-node ``taste_motor`` graph, which needs four rate runs per
    shuffle, so it defaults to ``n_taste = min(n, 300)``.

    Every row carries the empirical ``p_two_sided`` and its ``p_resolution``
    first, then ``z``, the null moments, the ``stabilised`` flag **and the raw
    ``null_effects``** - the paper's F6 histogram is rebuilt from the panel and
    could not be before.
    """
    t0 = time.perf_counter()
    modes = tuple(modes)
    n_taste = int(n_taste) if n_taste is not None else min(int(n), 300)
    rows: list[dict[str, Any]] = []
    warnings = list(BASE_WARNINGS)

    arms: list[tuple[str, str, int]] = [("subgraph", "mean_hz", int(n))]
    if include_taste_map:
        arms.append(("taste_map", "bitter_veto_ratio", n_taste))

    for assay, readout, n_i in arms:
        for mode in modes:
            res = null_distribution(
                assay,
                compound,
                conc_M,
                mode,
                n=n_i,
                seed=seed,
                readout=readout,
                n_jobs=n_jobs,
                **kw,
            )
            rows.append(
                {
                    "assay": assay,
                    "readout": readout,
                    "mode": mode,
                    "n": res["n"],
                    "n_ok": res["n_ok"],
                    "real_effect": res["real_effect"],
                    "p_two_sided": res["p_two_sided"],
                    "p_resolution": res["p_resolution"],
                    "beats_null": res["beats_null"],
                    "stabilised": res["stabilised"],
                    "n_stabilised": res["n_stabilised"],
                    "z": res["z"],
                    "null_mean": res["null_mean"],
                    "null_median": res["null_median"],
                    "null_sd": res["null_sd"],
                    "null_effects": res["null_effects"],
                    "engine": res["engine"],
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
    n: int = 200,
    seed: int = 0,
    z_threshold: float = 2.0,
    alpha: float = DEFAULT_ALPHA,
    **kw: Any,
) -> dict[str, Any]:
    """Fraction of null modes the real graph beats.

    ``score`` counts modes with ``|z| > z_threshold`` (the v0.5 definition, kept
    so the paper's history stays comparable); ``score_p`` counts modes with
    empirical ``p <= alpha`` and is the one to lead with, because it does not
    assume the null is normal.

    ``panel`` may be a :func:`null_panel` result or its ``rows``; when it is
    ``None`` the panel is computed first.  ``assay`` filters the rows to one
    assay arm (default: all rows).  1.0 means the connectome mattered against
    every null; 0.0 means a shuffled graph reproduced the drug effect.
    

    .. deprecated:: 0.6
        Superseded by :func:`flylab.analysis.dependence.dependence_profile`
        and :func:`~flylab.analysis.dependence.necessary_information_level`,
        which report an empirical permutation p per mode and name the weakest
        graph model that already reproduces the effect. This single score
        collapses that into one number and is no longer used by the paper.
    """
    if panel is None:
        panel = null_panel(compound, conc_M, modes=modes, n=n, seed=seed, **kw)
    rows = panel["rows"] if isinstance(panel, dict) else list(panel)
    if assay is not None:
        rows = [r for r in rows if r.get("assay") == assay]
    scored = [r for r in rows if r.get("z") is not None]
    hits = [r for r in scored if abs(float(r["z"])) > float(z_threshold)]
    p_scored = [r for r in rows if r.get("p_two_sided") is not None]
    p_hits = [r for r in p_scored if float(r["p_two_sided"]) <= float(alpha)]
    warnings = list(BASE_WARNINGS) + [
        "The score counts modes, not biology: it asks how many degradations of "
        "the MaleCNS wiring fail to reproduce this model's drug effect."
    ]
    if len(scored) < len(rows):
        warnings.append(f"{len(rows) - len(scored)} of {len(rows)} rows had an undefined z.")
    unstable = [r for r in rows if r.get("stabilised") is False]
    if unstable:
        warnings.append(
            f"{len(unstable)} of {len(rows)} rows had not stabilised at their n; "
            "their contribution to the score is provisional."
        )
    return {
        "score_p": (len(p_hits) / len(p_scored)) if p_scored else None,
        "n_beaten_p": len(p_hits),
        "n_scored_p": len(p_scored),
        "alpha": float(alpha),
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
                "p_two_sided": r.get("p_two_sided"),
                "z": r.get("z"),
                "beaten_p": (
                    r.get("p_two_sided") is not None
                    and float(r["p_two_sided"]) <= float(alpha)
                ),
                "beaten": (r.get("z") is not None and abs(float(r["z"])) > float(z_threshold)),
            }
            for r in rows
        ],
        "label": "model_derived",
        "warnings": warnings,
    }
