"""Rate network on the MaleCNS named-cell neighborhood.

Extracted verbatim (numerically) from the original ``flylab/assays/subgraph.py``
inline loop so that the spiking assay, the ensemble/sensitivity layer and the
analysis package all share one runtime.

Model
-----
``W_signed[j, i]`` is the signed synapse count from presynaptic node ``i`` onto
postsynaptic node ``j``.  The sign comes from the *presynaptic* consensus
transmitter (see :data:`SIGN`).  Rows are then normalised by
``max(sum_i |W[j, i]|, 1)`` so no cell can receive more than unit total drive.

The leaky update, for ``alpha = 0.3``::

    r <- clip((1 - alpha) * r + alpha * g_nav * (drive + (W_signed * scale) @ r), 0, 300)

``scale`` is a per-*presynaptic*-node vector built from the pharmacological
gains:

==================  ====================================
presynaptic nt      multiplier
==================  ====================================
acetylcholine       ``g_ach * ach_tone``
gaba                ``g_gaba``
glutamate           ``g_glu``
octopamine          ``g_oct``
anything else       ``1.0``
==================  ====================================

``g_nav`` is *not* a synapse scaling: it multiplies the total (recurrent +
external) input as a global excitability factor, which is how a Nav modulator
is represented at rate level.  ``ach_tone`` (AChE inhibition) multiplies ACh
edge strength on top of the receptor gain.

With all gains at 1.0 this reproduces the pre-refactor numbers exactly.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

SIGN: dict[str, float] = {
    "acetylcholine": 1.0,
    "glutamate": -0.4,
    "gaba": -1.0,
    "histamine": -0.5,
    "dopamine": 0.2,
    "serotonin": 0.2,
    "octopamine": 0.2,
}

#: transmitter -> gain key used to scale that transmitter's outgoing edges
NT_GAIN_KEY: dict[str, str] = {
    "acetylcholine": "g_ach",
    "gaba": "g_gaba",
    "glutamate": "g_glu",
    "octopamine": "g_oct",
}

DEFAULT_GAINS: dict[str, float] = {
    "g_ach": 1.0,
    "g_gaba": 1.0,
    "g_glu": 1.0,
    "g_oct": 1.0,
    "g_nav": 1.0,
    "ach_tone": 1.0,
}

CANDIDATES = [
    Path("data/derived/malecns_named_neighborhood.json"),
    Path(__file__).resolve().parents[2] / "data/derived/malecns_named_neighborhood.json",
    Path.home() / ".flylab/malecns_named_neighborhood.json",
]


def graph_path(path: str | Path | None = None) -> Path:
    if path:
        p = Path(path)
        if p.exists():
            return p
    for p in CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError("neighborhood JSON missing. Run extract-malecns-subgraph Action.")


_GRAPH_CACHE: dict[str, dict[str, Any]] = {}


def load_graph(path: str | Path | None = None) -> dict[str, Any]:
    """Load (and memoise) the committed neighborhood JSON.

    The returned dict is shared; callers must treat it as read-only.
    """
    key = str(graph_path(path).resolve())
    if key not in _GRAPH_CACHE:
        _GRAPH_CACHE[key] = json.loads(Path(key).read_text())
    return _GRAPH_CACHE[key]


def seed_ids(graph: dict[str, Any]) -> list[int]:
    out: list[int] = []
    for ids in graph.get("seeds", {}).values():
        for i in ids:
            if i not in out:
                out.append(i)
    return out


# --------------------------------------------------------------------------
# gains
# --------------------------------------------------------------------------
def fallback_gains(
    rows: list[dict[str, Any]],
    agonist_a: float = 0.4,
    agonist_b: float = 1.6,
) -> dict[str, float]:
    """Inline copy of the historical patch rules (used until / unless
    ``flylab.pharm.mechanisms.gains_from_occupancy`` is importable).

    * insect nAChR agonist     -> ``g_ach  = max(0.05, 1 + a*th - b*th^2)``
    * insect nAChR antagonist  -> ``g_ach  = max(0.05, 1 - th)``
    * insect RDL antagonist    -> ``g_gaba = max(0.05, 1 - th)``
    * insect RDL agonist       -> ``g_gaba = max(0.05, 1 + 0.4*th)``

    ``agonist_a`` / ``agonist_b`` default to the frozen 0.4 / 1.6 coefficients.
    They exist only so :func:`flylab.assays.ensemble.sensitivity` can probe the
    coefficient without editing the rule itself.
    """
    gains = dict(DEFAULT_GAINS)
    for row in rows or []:
        th = float(row.get("occupancy", 0.0))
        d = row.get("direction")
        rec = row.get("receptor")
        if rec == "insect_nAChR":
            if d == "agonist":
                gains["g_ach"] = max(0.05, 1.0 + agonist_a * th - agonist_b * th * th)
            elif d == "antagonist":
                gains["g_ach"] = max(0.05, 1.0 - th)
        elif rec == "insect_RDL":
            if d == "antagonist":
                gains["g_gaba"] = max(0.05, 1.0 - th)
            elif d == "agonist":
                gains["g_gaba"] = max(0.05, 1.0 + 0.4 * th)
    return gains


def compute_gains(
    compound: str | None,
    conc_M: float = 0.0,
    library: dict[str, Any] | None = None,
    rule_overrides: dict[str, float] | None = None,
) -> tuple[dict[str, float], dict[str, Any] | None]:
    """``(gains, occupancy_block)`` for one compound at one concentration.

    Prefers ``flylab.pharm.mechanisms.gains_from_occupancy`` (owned by the pharm
    agent, single source of truth).  Falls back to :func:`fallback_gains` while
    that module does not exist yet.  ``rule_overrides`` forces the fallback rule
    with modified coefficients and is only used by the sensitivity analysis.
    """
    if not compound:
        return dict(DEFAULT_GAINS), None
    from flylab.pharm.occupancy import compare_compound

    occ = compare_compound(compound, conc_M, library=library) if library else compare_compound(compound, conc_M)
    rows = occ["receptors"]
    if rule_overrides:
        return fallback_gains(rows, **rule_overrides), occ
    try:  # pragma: no cover - depends on the concurrently built pharm module
        from flylab.pharm.mechanisms import gains_from_occupancy

        raw = gains_from_occupancy(rows) or {}
        gains = {k: float(raw.get(k, v)) for k, v in DEFAULT_GAINS.items()}
    except Exception:
        gains = fallback_gains(rows)
    return gains, occ


# --------------------------------------------------------------------------
# network
# --------------------------------------------------------------------------
class RateNetwork:
    """Row-normalised signed rate network over the neighborhood graph."""

    def __init__(self, graph: dict[str, Any], sign: dict[str, float] | None = None):
        self.graph = graph
        self.sign = dict(sign or SIGN)
        self.nodes: list[dict[str, Any]] = graph["nodes"]
        self.body_ids: list[int] = [n["bodyId"] for n in self.nodes]
        self.node_index: dict[int, int] = {b: i for i, b in enumerate(self.body_ids)}
        self.nt_of_node: dict[int, str] = {
            n["bodyId"]: (n.get("consensus_nt") or "unclear") for n in self.nodes
        }
        self.types: list[str | None] = [n.get("type") for n in self.nodes]
        self.superclasses: list[str | None] = [n.get("superclass") for n in self.nodes]
        self.seed_ids: list[int] = seed_ids(graph)

        n = len(self.nodes)
        W = np.zeros((n, n), dtype=float)
        for e in graph["edges"]:
            i = self.node_index.get(e["pre"])
            j = self.node_index.get(e["post"])
            if i is None or j is None:
                continue
            nt = self.nt_of_node[self.body_ids[i]]
            W[j, i] += self.sign.get(nt, 0.0) * float(e["weight"])
        denom = np.maximum(np.abs(W).sum(axis=1, keepdims=True), 1.0)
        self.W_raw = W
        self.W_signed = W / denom
        self._nt_vec = np.array([self.nt_of_node[b] for b in self.body_ids], dtype=object)

    # -- helpers ----------------------------------------------------------
    def __len__(self) -> int:
        return len(self.nodes)

    def scale_vector(self, gains: dict[str, float] | None) -> np.ndarray:
        """Per-presynaptic-node synaptic multiplier from the gain dict."""
        g = dict(DEFAULT_GAINS)
        g.update({k: float(v) for k, v in (gains or {}).items() if k in DEFAULT_GAINS})
        out = np.ones(len(self.nodes), dtype=float)
        for idx, nt in enumerate(self._nt_vec):
            key = NT_GAIN_KEY.get(nt)
            if key is None:
                continue
            val = g[key]
            if nt == "acetylcholine":
                val *= g["ach_tone"]
            out[idx] = val
        return out

    def drive_vector(self, drive: dict[int, float] | None, default_hz: float = 0.0) -> np.ndarray:
        d = np.full(len(self.nodes), float(default_hz))
        if drive:
            for body_id, hz in drive.items():
                i = self.node_index.get(int(body_id))
                if i is not None:
                    d[i] = float(hz)
        return d

    def seed_drive(self, drive_hz: float) -> dict[int, float]:
        return {b: float(drive_hz) for b in self.seed_ids}

    # -- run --------------------------------------------------------------
    def run(
        self,
        drive: dict[int, float] | np.ndarray,
        gains: dict[str, float] | None = None,
        steps: int = 80,
        alpha: float = 0.3,
        r_max: float = 300.0,
    ) -> np.ndarray:
        """Iterate the leaky rate update and return per-node rates (Hz)."""
        g = dict(DEFAULT_GAINS)
        g.update({k: float(v) for k, v in (gains or {}).items() if k in DEFAULT_GAINS})
        d = drive if isinstance(drive, np.ndarray) else self.drive_vector(drive)
        Weff = self.W_signed * self.scale_vector(g)
        g_nav = g["g_nav"]
        r = np.zeros(len(self.nodes), dtype=float)
        for _ in range(int(steps)):
            r = np.clip((1.0 - alpha) * r + alpha * (g_nav * (d + Weff @ r)), 0.0, r_max)
        return r

    def rates_by_body(self, rates: np.ndarray) -> dict[int, float]:
        return {b: float(rates[i]) for i, b in enumerate(self.body_ids)}


_NET_CACHE: dict[str, RateNetwork] = {}


def rate_network(path: str | Path | None = None) -> RateNetwork:
    """Memoised :class:`RateNetwork` for the committed neighborhood graph."""
    key = str(graph_path(path).resolve())
    if key not in _NET_CACHE:
        _NET_CACHE[key] = RateNetwork(load_graph(key))
    return _NET_CACHE[key]


def by_superclass(net: RateNetwork, rates: np.ndarray) -> dict[str, float]:
    """Mean rate per MaleCNS superclass (``unknown`` for unannotated cells)."""
    acc: dict[str, list[float]] = {}
    for i, sc in enumerate(net.superclasses):
        acc.setdefault(sc or "unknown", []).append(float(rates[i]))
    return {k: float(np.mean(v)) for k, v in sorted(acc.items())}


def top_changed(
    net: RateNetwork,
    rates_vehicle: np.ndarray,
    rates_treated: np.ndarray,
    k: int = 15,
) -> list[dict[str, Any]]:
    """Top-``k`` nodes by ``|rate_treated - rate_vehicle|``."""
    delta = np.asarray(rates_treated, dtype=float) - np.asarray(rates_vehicle, dtype=float)
    order = np.argsort(-np.abs(delta))[: int(k)]
    out = []
    for i in order:
        i = int(i)
        out.append(
            {
                "bodyId": net.body_ids[i],
                "type": net.types[i],
                "superclass": net.superclasses[i],
                "nt": net.nt_of_node[net.body_ids[i]],
                "rate_vehicle": float(rates_vehicle[i]),
                "rate_treated": float(rates_treated[i]),
                "delta_hz": float(delta[i]),
            }
        )
    return out


__all__ = [
    "SIGN",
    "NT_GAIN_KEY",
    "DEFAULT_GAINS",
    "RateNetwork",
    "rate_network",
    "load_graph",
    "graph_path",
    "seed_ids",
    "compute_gains",
    "fallback_gains",
    "by_superclass",
    "top_changed",
]
