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
``drive_fanin * drive_hz``; each external event injects ``drive_w`` mV.  Cells
that are not seeds still get a uniform background Poisson drive
(``BACKGROUND_FRACTION_DEFAULT * drive_hz``) standing in for the inputs that
were cut away when the neighborhood was excised from the full CNS - without it
the 942 input-only cells are silent and the drug cannot act on them at all.

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
#: mV of injected charge per presynaptic *synapse* (per spike).
#: 0.05 mV/synapse: a 140-synapse connection then carries ~7 mV, i.e. one
#: presynaptic spike is about threshold-sized, and the strongest edge in the
#: committed neighborhood (739 synapses) is reliably suprathreshold.
W_SCALE_DEFAULT = 0.05
#: mV of injected charge per external Poisson drive event
DRIVE_W_DEFAULT = 0.5
#: number of independent external synapses per driven cell (drive "bundle").
#: fanin x drive_w = 10 mV.s, so a cell driven at f Hz sits at
#: v_rest + 0.2*f mV; 35 Hz is the rheobase of an isolated cell.
DRIVE_FANIN_DEFAULT = 20
#: background Poisson rate, as a fraction of ``drive_hz``, applied to every
#: cell to stand in for the inputs lost when the neighborhood was cut out of
#: the 25M-edge CNS.  0.65 x 40 Hz = 26 Hz keeps ~60% of the named
#: neighborhood spiking and vehicle MN9 near 30 Hz without saturating.
BACKGROUND_FRACTION_DEFAULT = 0.65

# --- WebAssembly / browser defaults ---------------------------------------
#: Simulated window used by the static (Pyodide) build.  The integration step
#: stays at :data:`DT_MS` = 0.1 ms -- a coarser step would be a *different*
#: model, not the same one running faster -- so the browser shortens the
#: window instead.  200 ms still clears the 100 ms settle cap with a 100 ms
#: counting window.
BROWSER_T_MS = 200.0
#: The browser build must NOT relax the step: 0.1 ms here is the same constant
#: the native build integrates with.
BROWSER_DT_MS = DT_MS
#: The browser build runs the sparse kernel (same spikes, less work per step).
BROWSER_SPARSE = True


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
        sparse: bool = False,
    ):
        self.graph = graph
        self.seed = int(seed)
        self.w_scale = float(w_scale)
        self.dt_ms = float(dt_ms)
        self.drive_w = float(drive_w)
        self.drive_fanin = int(drive_fanin)
        self.sign = dict(sign or SIGN)
        #: default kernel for :meth:`run` (``run(sparse=...)`` overrides it)
        self.sparse = bool(sparse)

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
        #: CSR-style arrays over the *columns* of ``C_signed`` (built lazily)
        self._csr: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None

    def __len__(self) -> int:
        return len(self.nodes)

    # -- sparse representation --------------------------------------------
    def csr(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(indptr, indices, data)`` of the signed synapse matrix, by column.

        One entry per non-zero ``C_signed[j, i]``; column ``i`` (the
        *presynaptic* cell, the one that spikes) owns the slice
        ``indices[indptr[i]:indptr[i + 1]]`` of postsynaptic rows ``j`` with
        weights ``data[...]``.  Column-major is the useful orientation here
        because a step touches only the columns of cells that fired.

        Built once from the dense matrix and cached, so the dense and sparse
        kernels are two readings of exactly the same numbers.  On the committed
        1126-node neighborhood this is 1230 non-zeros against 1.27M dense
        cells, i.e. ~0.1% fill.
        """
        if self._csr is None:
            C = self.C_signed
            rows, cols = np.nonzero(C)
            order = np.argsort(cols, kind="stable")  # group by presynaptic column
            rows, cols = rows[order], cols[order]
            data = C[rows, cols].astype(np.float32, copy=True)
            indptr = np.zeros(C.shape[1] + 1, dtype=np.int64)
            np.add.at(indptr, cols + 1, 1)
            indptr = np.cumsum(indptr)
            self._csr = (indptr, rows.astype(np.int64), data)
        return self._csr

    @property
    def nnz(self) -> int:
        """Number of non-zero signed connections."""
        return int(self.csr()[2].size)

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
        sparse: bool | None = None,
    ) -> SpikeResult:
        """Simulate ``t_ms`` of network activity.

        ``drive_hz`` maps bodyId -> external Poisson rate (Hz); any node may be
        driven, not just the seeds.  Deterministic for a given ``seed``.
        Rates are counted over ``[settle_ms, t_ms]`` (default: first 20% of the
        run, capped at 100 ms, is discarded as transient).

        ``sparse`` selects the synaptic kernel: ``False`` sums the columns of
        the dense effective-weight matrix (the default, bit-identical to every
        previous release), ``True`` walks only the non-zero entries of the
        columns that fired.  Both reduce each column in the same ascending row
        order, so the two kernels produce identical spike trains for the same
        seed; ``None`` means "whatever this network was built with".
        """
        g = dict(DEFAULT_GAINS)
        g.update({k: float(v) for k, v in (gains or {}).items() if k in DEFAULT_GAINS})
        n = len(self.nodes)
        dt = self.dt_ms
        n_steps = max(int(round(float(t_ms) / dt)), 1)
        if settle_ms is None:
            settle_ms = min(100.0, 0.2 * float(t_ms))
        settle_step = int(round(float(settle_ms) / dt))

        use_sparse = self.sparse if sparse is None else bool(sparse)
        if use_sparse:
            indptr, indices, data = self.csr()
            col_scale = (self.w_scale * self.gain_vector(g)).astype(np.float32)
            # per-entry effective weight: data is C_signed, scaled by its own
            # presynaptic column's gain -- exactly ``effective_weights()``
            w_sparse = data * np.repeat(col_scale, np.diff(indptr))
            contrib = np.zeros(n, dtype=np.float32)
            W = None
        else:
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
                if use_sparse:
                    contrib[:] = 0.0
                    for i in fired:
                        lo, hi = indptr[i], indptr[i + 1]
                        if hi > lo:
                            contrib[indices[lo:hi]] += w_sparse[lo:hi]
                    I += contrib * w_over_tau
                else:
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
    sparse = bool(kw.pop("sparse", False))
    key = (id(graph), float(kw.get("w_scale", W_SCALE_DEFAULT)), float(kw.get("drive_w", DRIVE_W_DEFAULT)))
    net = _LIF_CACHE.get(key)
    if net is None:
        net = LIFNetwork(graph, seed=seed, **kw)
        _LIF_CACHE[key] = net
    net.seed = int(seed)
    net.sparse = sparse
    return net


__all__ = [
    "LIFNetwork",
    "SpikeResult",
    "lif_network",
    "W_SCALE_DEFAULT",
    "DRIVE_W_DEFAULT",
    "DRIVE_FANIN_DEFAULT",
    "BACKGROUND_FRACTION_DEFAULT",
    "TAU_M_MS",
    "V_REST_MV",
    "V_TH_MV",
    "V_RESET_MV",
    "T_REF_MS",
    "TAU_S_MS",
    "DT_MS",
    "BROWSER_T_MS",
    "BROWSER_DT_MS",
    "BROWSER_SPARSE",
]
