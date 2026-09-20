"""Leaky integrate-and-fire network on the MaleCNS named-cell neighborhood.

Parameters follow the leaky integrate-and-fire formulation used by
Shiu, Sterne, Spiller et al., *Nature* 634 (2024) for the *Drosophila* brain
("A leaky integrate-and-fire computational model based on the connectome of the
entire adult Drosophila brain reveals insights into sensorimotor processing").
They are literature-order defaults, **not** a fit to any FlyLab measurement:

===================  ==========  ==========================================
parameter            value       note
===================  ==========  ==========================================
``tau_m``            20 ms       membrane time constant
``v_rest``           -52 mV      resting potential
``v_th``             -45 mV      spike threshold (7 mV above rest)
``v_reset``          -52 mV      reset to rest
``t_ref``            2.2 ms      absolute refractory period
``tau_s``            5 ms        exponential synaptic current
``dt``               0.1 ms      forward-Euler step
===================  ==========  ==========================================

Synapse model: current-based exponential synapses.  A presynaptic spike on
edge ``i -> j`` injects a charge of ``W_eff[j, i]`` mV into ``j``, i.e.

    dI/dt = -I / tau_s,   I[j] += W_eff[j, i] / tau_s  on a spike
    dV/dt = (v_rest - V) / tau_m + g_nav * I

so ``W_eff`` is in millivolts of total depolarisation per presynaptic spike
(peak PSP is a little smaller because the membrane leaks while the synapse
decays).  ``W_eff[j, i] = w_scale * synapse_count * SIGN[nt_i] * gain(nt_i)``
with the same signs and gain keys as :mod:`flylab.circuit.rate`.

External drive is an independent Poisson process per driven cell at
``drive_hz``; each external event injects ``drive_w`` mV.

Calibration (``W_SCALE_DEFAULT`` / ``DRIVE_W_DEFAULT``): chosen once so that
the *vehicle* MN9 cells fire in a plausible 5-60 Hz band under the standard
40 Hz seed drive, and so the rest of the neighborhood is neither silent nor
saturated.  Nothing here is fitted to animal data - see
``docs/RESEARCH_MAP.md`` parameter policy.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from flylab.circuit.rate import DEFAULT_GAINS, NT_GAIN_KEY, SIGN, seed_ids

# --- Shiu-style cell parameters -------------------------------------------
TAU_M_MS = 20.0
V_REST_MV = -52.0
V_TH_MV = -45.0
V_RESET_MV = -52.0
T_REF_MS = 2.2
TAU_S_MS = 5.0
DT_MS = 0.1

# --- FlyLab calibration constants (see module docstring) ------------------
#: mV of injected charge per presynaptic *synapse* (per spike)
W_SCALE_DEFAULT = 0.0125
#: mV of injected charge per external Poisson drive event
DRIVE_W_DEFAULT = 0.75
#: number of independent external synapses per driven cell (drive "bundle")
DRIVE_FANIN_DEFAULT = 20


@dataclass
class SpikeResult:
    """Spikes and rates from one :meth:`LIFNetwork.run`."""

    spikes: list[tuple[float, int]]
    rates_hz: np.ndarray
    t_ms: float
    window_ms: tuple[float, float] = (0.0, 0.0)
    body_ids: list[int] = field(default_factory=list)

    def spike_times(self, index: int) -> list[float]:
        return [t for t, i in self.spikes if i == index]

    def rate_of(self, body_id: int) -> float | None:
        try:
            return float(self.rates_hz[self.body_ids.index(body_id)])
        except ValueError:
            return None


class LIFNetwork:
    """Dense-matrix LIF over the ~1126-node neighborhood graph.

    Dense is fine at this size (1126 x 1126 float32 = 5 MB); per step only the
    columns of cells that actually spiked are summed, so a 500 ms run is well
    under a second.
    """

    def __init__(
        self,
        graph: dict[str, Any],
        seed: int = 0,
        w_scale: float = W_SCALE_DEFAULT,
        dt_ms: float = DT_MS,
        drive_w: float = DRIVE_W_DEFAULT,
        drive_fanin: int = DRIVE_FANIN_DEFAULT,
        sign: dict[str, float] | None = None,
    ):
        self.graph = graph
        self.seed = int(seed)
        self.w_scale = float(w_scale)
        self.dt_ms = float(dt_ms)
        self.drive_w = float(drive_w)
        self.drive_fanin = int(drive_fanin)
        self.sign = dict(sign or SIGN)

        self.nodes: list[dict[str, Any]] = graph["nodes"]
        self.body_ids: list[int] = [n["bodyId"] for n in self.nodes]
        self.node_index: dict[int, int] = {b: i for i, b in enumerate(self.body_ids)}
        self.nt_of_node: dict[int, str] = {
            n["bodyId"]: (n.get("consensus_nt") or "unclear") for n in self.nodes
        }
        self.types = [n.get("type") for n in self.nodes]
        self.superclasses = [n.get("superclass") for n in self.nodes]
        self.seed_ids = seed_ids(graph)

        n = len(self.nodes)
        C = np.zeros((n, n), dtype=np.float32)  # signed synapse counts
        for e in graph["edges"]:
            i = self.node_index.get(e["pre"])
            j = self.node_index.get(e["post"])
            if i is None or j is None:
                continue
            C[j, i] += self.sign.get(self.nt_of_node[self.body_ids[i]], 0.0) * float(e["weight"])
        self.C_signed = C
        self._nt_vec = [self.nt_of_node[b] for b in self.body_ids]

    def __len__(self) -> int:
        return len(self.nodes)

    # -- gains -------------------------------------------------------------
    def gain_vector(self, gains: dict[str, float] | None) -> np.ndarray:
        g = dict(DEFAULT_GAINS)
        g.update({k: float(v) for k, v in (gains or {}).items() if k in DEFAULT_GAINS})
        out = np.ones(len(self.nodes), dtype=np.float32)
        for idx, nt in enumerate(self._nt_vec):
            key = NT_GAIN_KEY.get(nt)
            if key is None:
                continue
            val = g[key]
            if nt == "acetylcholine":
                val *= g["ach_tone"]
            out[idx] = val
        return out

    def effective_weights(self, gains: dict[str, float] | None) -> np.ndarray:
        return self.C_signed * (self.w_scale * self.gain_vector(gains))

    def seed_drive(self, drive_hz: float) -> dict[int, float]:
        return {b: float(drive_hz) for b in self.seed_ids}

    # -- simulation --------------------------------------------------------
    def run(
        self,
        drive_hz: dict[int, float],
        gains: dict[str, float] | None = None,
        t_ms: float = 500.0,
        settle_ms: float | None = None,
    ) -> SpikeResult:
        """Simulate ``t_ms`` of network activity.

        ``drive_hz`` maps bodyId -> external Poisson rate (Hz); any node may be
        driven, not just the seeds.  Deterministic for a given ``seed``.
        Rates are counted over ``[settle_ms, t_ms]`` (default: first 20% of the
        run, capped at 100 ms, is discarded as transient).
        """
        g = dict(DEFAULT_GAINS)
        g.update({k: float(v) for k, v in (gains or {}).items() if k in DEFAULT_GAINS})
        n = len(self.nodes)
        dt = self.dt_ms
        n_steps = max(int(round(float(t_ms) / dt)), 1)
        if settle_ms is None:
            settle_ms = min(100.0, 0.2 * float(t_ms))
        settle_step = int(round(float(settle_ms) / dt))

        W = self.effective_weights(g)
        g_nav = np.float32(g["g_nav"])

        drive_idx = []
        drive_lam = []
        for body_id, hz in (drive_hz or {}).items():
            i = self.node_index.get(int(body_id))
            if i is not None and float(hz) > 0:
                drive_idx.append(i)
                drive_lam.append(float(hz) * self.drive_fanin * dt / 1000.0)
        drive_idx_arr = np.array(drive_idx, dtype=np.int64)
        drive_lam_arr = np.array(drive_lam, dtype=np.float64)
        drive_w = np.float32(self.drive_w / TAU_S_MS)

        rng = np.random.default_rng(self.seed)
        V = np.full(n, V_REST_MV, dtype=np.float32)
        I = np.zeros(n, dtype=np.float32)
        ref = np.zeros(n, dtype=np.float32)  # remaining refractory time (ms)
        counts = np.zeros(n, dtype=np.int64)
        spikes: list[tuple[float, int]] = []

        decay = np.float32(np.exp(-dt / TAU_S_MS))
        leak = np.float32(dt / TAU_M_MS)
        v_rest = np.float32(V_REST_MV)
        w_over_tau = np.float32(1.0 / TAU_S_MS)

        for step in range(n_steps):
            t = (step + 1) * dt
            # external Poisson drive
            if drive_idx_arr.size:
                ext = rng.poisson(drive_lam_arr)
                nz = ext > 0
                if nz.any():
                    I[drive_idx_arr[nz]] += drive_w * ext[nz].astype(np.float32)
            # membrane update
            free = ref <= 0.0
            V[free] += leak * (v_rest - V[free]) + np.float32(dt) * g_nav * I[free]
            ref[~free] -= dt
            V[~free] = V_RESET_MV
            # spikes
            fired = np.flatnonzero(V >= V_TH_MV)
            if fired.size:
                V[fired] = V_RESET_MV
                ref[fired] = T_REF_MS
                I += W[:, fired].sum(axis=1) * w_over_tau
                if step >= settle_step:
                    counts[fired] += 1
                for i in fired:
                    spikes.append((float(t), int(i)))
            I *= decay

        window_ms = (float(settle_step * dt), float(n_steps * dt))
        dur_s = max((window_ms[1] - window_ms[0]) / 1000.0, 1e-9)
        rates = counts.astype(float) / dur_s
        return SpikeResult(
            spikes=spikes,
            rates_hz=rates,
            t_ms=float(n_steps * dt),
            window_ms=window_ms,
            body_ids=list(self.body_ids),
        )


_LIF_CACHE: dict[tuple, LIFNetwork] = {}


def lif_network(graph: dict[str, Any], seed: int = 0, **kw) -> LIFNetwork:
    """Memoised LIF network (the signed synapse matrix is seed-independent)."""
    key = (id(graph), float(kw.get("w_scale", W_SCALE_DEFAULT)), float(kw.get("drive_w", DRIVE_W_DEFAULT)))
    net = _LIF_CACHE.get(key)
    if net is None:
        net = LIFNetwork(graph, seed=seed, **kw)
        _LIF_CACHE[key] = net
    net.seed = int(seed)
    return net


__all__ = [
    "LIFNetwork",
    "SpikeResult",
    "lif_network",
    "W_SCALE_DEFAULT",
    "DRIVE_W_DEFAULT",
    "TAU_M_MS",
    "V_REST_MV",
    "V_TH_MV",
    "V_RESET_MV",
    "T_REF_MS",
    "TAU_S_MS",
    "DT_MS",
]
