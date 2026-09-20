"""Variance-based global uncertainty attribution (Sobol'-style) for FlyLab.

What this answers
-----------------
The one-at-a-time tornado in :func:`flylab.assays.ensemble.sensitivity` varies
one parameter at a time around one operating point.  That is not an uncertainty
budget: it cannot see interactions, it depends on the base point, and it cannot
say what *share* of the output's variance each assumption owns.  This module
does the global version: it samples every named assumption jointly and
decomposes the variance of a circuit readout into first-order and total-order
Sobol' indices, with bootstrap confidence intervals and a convergence
diagnostic, and prints the result as an **uncertainty budget**.

The nine factors
----------------
======================  ============================================================
``potency``             log10 shift of every sourced potency/affinity value of the
                        compound (+/- ``potency_log10``, default +/-0.6 ~ x/4)
``hill_n``              multiplier on every Hill coefficient
                        (``10**+/-hill_log10``, default x/2)
``gain_transform``      **categorical**: which member of
                        :data:`flylab.analysis.robustness.SHAPES` maps engagement
                        onto gain
``gain_coef``           multiplier on every gain-rule coefficient (0.5 - 1.5x)
``drive``               external drive rate (0.5 - 1.5x of the assay default)
``expression``          receptor-expression assumption: lambda in [0, 1]
                        interpolating each node between the uniform model
                        (lambda = 0, what FlyLab ships) and the cross-atlas
                        expression weighting of :mod:`flylab.pharm.expression`
``weight_threshold``    minimum synapse count kept in the graph (base .. 4x base)
``transmitter``         probability that a predicted transmitter is wrong
                        (0 .. ``transmitter_p_max``, default 0.20); that fraction
                        of nodes is relabelled from the graph's own transmitter
                        distribution
``lif_seed``            the LIF network's RNG stream
======================  ============================================================

Design choice: Jansen estimators on Latin-hypercube base matrices
-----------------------------------------------------------------
FlyLab's science core may not import ``scipy`` (browser build), so
``scipy.stats.qmc.Sobol`` is unavailable and a Joe-Kuo direction-number table
would have to be vendored.  Instead the standard Saltelli *cross-sampling*
design is built on two independent **Latin-hypercube** base matrices ``A`` and
``B`` (``n_base x k`` each), plus the ``k`` matrices ``AB^j`` = ``A`` with
column ``j`` taken from ``B``: ``n_base * (k + 2)`` model evaluations.  LHS
stratifies each factor's marginal, which at the small sample sizes this model
can afford is closer to a low-discrepancy sequence than plain Monte-Carlo, and
the cross-sampling identity that the estimators rest on does not require the
base samples to be independent *across rows*, only that ``A`` and ``B`` are
independent of each other.

The estimators are Jansen's (1999), which are the recommended pair for small
first-order and large total-order indices:

.. math::

    S_j = \\frac{V(Y) - \\frac{1}{2N}\\sum (f(B) - f(A_B^j))^2}{V(Y)}
    \\qquad
    T_j = \\frac{\\frac{1}{2N}\\sum (f(A) - f(A_B^j))^2}{V(Y)}

Both are unbiased only in the limit; at finite ``N`` an index that is truly 0
is estimated as a small positive **or negative** number.  Estimates are
therefore *not* clipped -- a negative ``S_j`` is reported as-is and is itself a
convergence diagnostic.  The ``lif_seed`` factor has no effect at all on the
deterministic rate engine, so its indices are a built-in **null-factor
control**: they measure the estimator's own noise floor.

The categorical factor
----------------------
``gain_transform`` is unordered.  It is indexed by
``SHAPE_ORDER[floor(u * len(SHAPE_ORDER))]`` with :data:`SHAPE_ORDER` fixed and
documented, so the sampler visits the five transformations with equal
probability.  Its Sobol' indices must be read as **group indices**: "how much
of the output variance is owned by the choice of transformation class", not as
a derivative.  The numerical value of the index does not depend on the order
(any permutation of an equiprobable categorical gives the same variance
decomposition), but the *stratification* LHS applies to that column does, which
is why the order is frozen here rather than sampled.

Runtime
-------
``sobol_analysis()`` defaults to ``n_base=128``, i.e. ``128 * 11 = 1408``
evaluations of the rate readout at roughly 15-25 ms each: **30-60 s**.
``engine="lif"`` costs about 0.5 s per evaluation, so use ``n_base<=32`` there
(~5 min) -- it is the only way ``lif_seed`` gets a non-zero index.

One practical note: each evaluation is 80 small dense matrix-vector products,
which are memory-bound, so a multithreaded BLAS spends more time synchronising
than computing and a loaded machine can slow a run down by an order of
magnitude.  Setting ``OPENBLAS_NUM_THREADS=1`` (or the OpenMP/MKL equivalent)
before importing numpy makes a long run several times faster and does not
change a single digit of the result.
"""
from __future__ import annotations

import math
import time
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from flylab.analysis.robustness import SHAPES, make_spec, spec_gains

__all__ = [
    "FACTORS",
    "RESOLUTION_STATES",
    "RESOLVED_LABEL",
    "UNRESOLVED_LABEL",
    "NULL_CONTROL_LABEL",
    "NULL_FACTOR",
    "resolution_summary",
    "FACTOR_NAMES",
    "SHAPE_ORDER",
    "RateSurrogate",
    "baseline_unit_vector",
    "lhs",
    "saltelli_matrices",
    "jansen_indices",
    "sobol_analysis",
    "uncertainty_budget",
    "to_markdown",
    "ishigami",
    "ishigami_unit",
    "ISHIGAMI_REFERENCE",
]

#: frozen order of the categorical ``gain_transform`` factor
SHAPE_ORDER: tuple[str, ...] = (
    "linear",
    "saturating",
    "weak_biphasic",
    "flylab_biphasic",
    "strong_biphasic",
)
assert set(SHAPE_ORDER) == set(SHAPES)

#: index of the shipped FlyLab transformation inside :data:`SHAPE_ORDER`
DEFAULT_SHAPE_INDEX = SHAPE_ORDER.index("flylab_biphasic")

#: gain keys that scale a presynaptic transmitter's outgoing edges, in the
#: column order the surrogate uses internally
_GAIN_COLUMNS: tuple[str, ...] = ("g_ach", "g_gaba", "g_glu", "g_oct")
_GAIN_COLUMN_INDEX = {k: i for i, k in enumerate(_GAIN_COLUMNS)}

#: how many dense thresholded weight matrices to keep (each is ~10 MB)
_MATRIX_CACHE_SIZE = 8

FACTORS: tuple[dict[str, Any], ...] = (
    {
        "name": "potency",
        "kind": "continuous",
        "what": "sourced potency / affinity value (EC50, IC50, Kd, Ki) of the compound",
        "range": "+/- potency_log10 decades (default 0.6, about a factor of 4 either way)",
        "baseline_u": 0.5,
    },
    {
        "name": "hill_n",
        "kind": "continuous",
        "what": "Hill coefficient of every sourced row",
        "range": "x 10**+/-hill_log10 (default 0.3, i.e. x/2)",
        "baseline_u": 0.5,
    },
    {
        "name": "gain_transform",
        "kind": "categorical",
        "what": "which transformation maps receptor engagement onto synaptic gain",
        "range": f"uniform over {len(SHAPE_ORDER)} prespecified shapes: {', '.join(SHAPE_ORDER)}",
        "baseline_u": (DEFAULT_SHAPE_INDEX + 0.5) / len(SHAPE_ORDER),
    },
    {
        "name": "gain_coef",
        "kind": "continuous",
        "what": "multiplier on every gain-rule coefficient",
        "range": "0.5 - 1.5x the shipped coefficients",
        "baseline_u": 0.5,
    },
    {
        "name": "drive",
        "kind": "continuous",
        "what": "external drive applied to the seed cells",
        "range": "0.5 - 1.5x the assay default",
        "baseline_u": 0.5,
    },
    {
        "name": "expression",
        "kind": "continuous",
        "what": "receptor-expression assumption (uniform gains vs cross-atlas weighting)",
        "range": "lambda 0 (uniform, shipped) - 1 (fully expression weighted)",
        "baseline_u": 0.0,
    },
    {
        "name": "weight_threshold",
        "kind": "continuous",
        "what": "minimum synapse count an edge must carry to be kept",
        "range": "graph min_weight .. threshold_max_mult x min_weight (default 4x)",
        "baseline_u": 0.0,
    },
    {
        "name": "transmitter",
        "kind": "continuous",
        "what": "probability that a predicted consensus transmitter is wrong",
        "range": "0 - transmitter_p_max (default 0.20) of nodes relabelled",
        "baseline_u": 0.0,
    },
    {
        "name": "lif_seed",
        "kind": "seed",
        "what": "LIF stochasticity (Poisson drive and initial state)",
        "range": "RNG stream; NULL FACTOR on the deterministic rate engine",
        "baseline_u": 0.0,
    },
)

FACTOR_NAMES: tuple[str, ...] = tuple(f["name"] for f in FACTORS)

BASE_WARNINGS = [
    "model_derived: these are shares of the variance of a SIMULATED readout "
    "under stated input ranges. They are not shares of biological variability, "
    "and a factor's index is only as meaningful as the range assumed for it.",
    "Input ranges are modelling assumptions declared in FACTORS; widening a "
    "factor's range raises its index almost mechanically.",
    "Sobol indices for the categorical gain_transform factor are group indices "
    "over the prespecified transformation family, not derivatives.",
    "Jansen estimates are unclipped: a small negative first-order index means "
    "the index is indistinguishable from zero at this sample size.",
]


# --------------------------------------------------------------------------
# sampling designs
# --------------------------------------------------------------------------
def lhs(n: int, k: int, rng: np.random.Generator) -> np.ndarray:
    """Latin-hypercube sample of ``n`` points in the unit ``k``-cube."""
    out = np.empty((n, k), dtype=float)
    cut = np.arange(n, dtype=float)
    for j in range(k):
        out[:, j] = (cut + rng.random(n)) / n
        rng.shuffle(out[:, j])
    return out


def saltelli_matrices(n_base: int, k: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """``(A, B, [AB^0 .. AB^{k-1}])`` cross-sampling design.

    ``n_base * (k + 2)`` rows in total.
    """
    rng = np.random.default_rng(seed)
    A = lhs(n_base, k, rng)
    B = lhs(n_base, k, rng)
    AB = []
    for j in range(k):
        M = A.copy()
        M[:, j] = B[:, j]
        AB.append(M)
    return A, B, AB


def jansen_indices(
    yA: np.ndarray, yB: np.ndarray, yAB: Sequence[np.ndarray]
) -> tuple[np.ndarray, np.ndarray, float]:
    """Jansen (1999) first-order and total-order indices, unclipped.

    Returns ``(S, T, var)`` where ``var`` is the pooled output variance of
    ``[yA, yB]``.
    """
    yA = np.asarray(yA, dtype=float)
    yB = np.asarray(yB, dtype=float)
    var = float(np.var(np.concatenate([yA, yB]), ddof=1))
    k = len(yAB)
    S = np.zeros(k)
    T = np.zeros(k)
    if var <= 0 or not math.isfinite(var):
        return S, T, var
    n = yA.size
    for j, yj in enumerate(yAB):
        yj = np.asarray(yj, dtype=float)
        S[j] = (var - np.sum((yB - yj) ** 2) / (2.0 * n)) / var
        T[j] = (np.sum((yA - yj) ** 2) / (2.0 * n)) / var
    return S, T, var


# --------------------------------------------------------------------------
# the FlyLab model as a function of the unit cube
# --------------------------------------------------------------------------
class RateSurrogate:
    """The FlyLab rate assay, re-expressed as ``f: [0,1]^9 -> readout``.

    This is not an approximation: it is the same signed, row-normalised,
    80-step leaky rate update that :class:`flylab.circuit.rate.RateNetwork`
    runs, assembled once from the committed graph so that a Sobol' design of a
    few thousand evaluations finishes in under a minute.  It reproduces
    :func:`flylab.assays.subgraph.run_subgraph_assay` bit-for-bit at
    :func:`baseline_unit_vector` -- the test suite asserts it.

    ``engine="lif"`` runs the real spiking assay instead (0.5 s per evaluation)
    and is the only way ``lif_seed`` can carry variance.
    """

    def __init__(
        self,
        compound: str = "imidacloprid",
        conc_M: float = 1e-6,
        readout: str = "mean_hz",
        graph: str = "named",
        engine: str = "rate",
        steps: int = 80,
        alpha: float = 0.3,
        r_max: float = 300.0,
        drive_hz: float = 40.0,
        potency_log10: float = 0.6,
        hill_log10: float = 0.3,
        gain_coef_range: tuple[float, float] = (0.5, 1.5),
        drive_range: tuple[float, float] = (0.5, 1.5),
        threshold_max_mult: float = 4.0,
        transmitter_p_max: float = 0.20,
        seed: int = 0,
    ):
        from flylab.circuit.rate import SIGN, NT_GAIN_KEY, load_graph

        if engine not in ("rate", "lif"):
            raise ValueError("engine must be 'rate' or 'lif'")
        self.compound = compound
        self.conc_M = float(conc_M)
        self.readout = readout
        self.graph_name = graph
        self.engine = engine
        self.steps = int(steps)
        self.alpha = float(alpha)
        self.r_max = float(r_max)
        self.drive_hz = float(drive_hz)
        self.potency_log10 = float(potency_log10)
        self.hill_log10 = float(hill_log10)
        self.gain_coef_range = tuple(gain_coef_range)
        self.drive_range = tuple(drive_range)
        self.threshold_max_mult = float(threshold_max_mult)
        self.transmitter_p_max = float(transmitter_p_max)
        self.seed = int(seed)
        self.n_calls = 0

        g = load_graph(graph)
        self.graph = g
        self._sign = dict(SIGN)
        self._nt_gain_key = dict(NT_GAIN_KEY)
        nodes = g["nodes"]
        self.body_ids = [int(n["bodyId"]) for n in nodes]
        self.index = {b: i for i, b in enumerate(self.body_ids)}
        self.n = len(nodes)
        self.nt0 = np.array(
            [str(n.get("consensus_nt") or "unclear") for n in nodes], dtype=object
        )

        pre, post, w = [], [], []
        for e in g["edges"]:
            i = self.index.get(int(e["pre"]))
            j = self.index.get(int(e["post"]))
            if i is None or j is None:
                continue
            pre.append(i)
            post.append(j)
            w.append(float(e["weight"]))
        self.pre = np.asarray(pre, dtype=np.intp)
        self.post = np.asarray(post, dtype=np.intp)
        self.w = np.asarray(w, dtype=float)
        self.base_min_weight = float(g.get("min_weight") or (self.w.min() if self.w.size else 1.0))

        # seed drive
        seeds = {int(b) for ids in g.get("seeds", {}).values() for b in ids}
        self.seed_mask = np.zeros(self.n, dtype=float)
        for b in seeds:
            i = self.index.get(b)
            if i is not None:
                self.seed_mask[i] = 1.0
        self.mn9_idx = [self.index[b] for b in g.get("seeds", {}).get("MN9", []) if b in self.index]
        self.dnp01_idx = [self.index[b] for b in g.get("seeds", {}).get("DNp01", []) if b in self.index]

        # fixed transmitter-relabelling stream: node order and replacement
        # labels are drawn once, so the ``transmitter`` factor is the *rate*
        # p and p = 0 is exactly the shipped graph.
        rng = np.random.default_rng([self.seed, 0xA11CE])
        self._relabel_order = rng.permutation(self.n)
        labels, counts = np.unique(self.nt0, return_counts=True)
        probs = counts / counts.sum()
        self._relabel_to = rng.choice(labels, size=self.n, p=probs)

        # expression weights (cached; lambda = 0 gives the uniform model)
        self._expr = self._expression_weights()

        # base occupancy rows
        self._rows0 = self._base_rows()
        self._cacheW: dict[float, np.ndarray] = {}
        self._cache_nt: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        self._unique_w = np.unique(self.w) if self.w.size else np.asarray([1.0])
        self._Wbuf = np.zeros((self.n, self.n), dtype=float)

    # -- inputs ----------------------------------------------------------
    def _expression_weights(self) -> dict[str, np.ndarray]:
        out = {k: np.ones(self.n, dtype=float) for k in ("g_ach", "g_gaba", "g_glu")}
        try:
            from flylab.pharm.expression import expression_weights

            for receptor, key in (
                ("insect_nAChR", "g_ach"),
                ("insect_RDL", "g_gaba"),
                ("insect_GluCl", "g_glu"),
            ):
                wmap = expression_weights(receptor, self.graph_name)
                out[key] = np.array([float(wmap.get(b, 1.0)) for b in self.body_ids], dtype=float)
        except Exception:  # pragma: no cover - dataset absent
            pass
        return out

    def _base_rows(self) -> list[dict[str, Any]]:
        from flylab.pharm.occupancy import compare_compound

        if not self.compound:
            return []
        return [dict(r) for r in compare_compound(self.compound, self.conc_M)["receptors"]]

    def _rows(self, potency_shift: float, hill_mult: float) -> list[dict[str, Any]]:
        """Base rows with the potency and Hill factors applied.

        Not-modelled rows (``engagement is None``) stay ``None``: perturbing a
        row with no sourced value would invent evidence.
        """
        from flylab.pharm.occupancy import hill_occupancy

        out = []
        for row in self._rows0:
            r = dict(row)
            eng = r.get("engagement", r.get("occupancy"))
            value = r.get("param_value_M", r.get("value_M", r.get("ec50_M")))
            if eng is None or value is None:
                r["engagement"] = None
                r["occupancy"] = None
                out.append(r)
                continue
            v = float(value) * (10.0 ** potency_shift)
            n = float(r.get("n", 1.0) or 1.0) * hill_mult
            e = hill_occupancy(self.conc_M, v, n)
            r["engagement"] = e
            r["occupancy"] = e
            out.append(r)
        return out

    def decode(self, u: Sequence[float]) -> dict[str, Any]:
        """Unit vector -> the named settings it stands for (for reporting)."""
        u = list(u)
        idx = min(int(u[2] * len(SHAPE_ORDER)), len(SHAPE_ORDER) - 1)
        lo, hi = self.gain_coef_range
        dlo, dhi = self.drive_range
        return {
            "potency_shift_log10": self.potency_log10 * (2.0 * u[0] - 1.0),
            "hill_mult": 10.0 ** (self.hill_log10 * (2.0 * u[1] - 1.0)),
            "gain_transform": SHAPE_ORDER[idx],
            "gain_coef": lo + u[3] * (hi - lo),
            "drive_hz": self.drive_hz * (dlo + u[4] * (dhi - dlo)),
            "expression_lambda": float(u[5]),
            "min_weight": self.base_min_weight * (1.0 + u[6] * (self.threshold_max_mult - 1.0)),
            "transmitter_p": self.transmitter_p_max * float(u[7]),
            "lif_seed": int(float(u[8]) * (2 ** 31 - 1)),
        }

    # -- network ---------------------------------------------------------
    def _sign_vector(self, transmitter_p: float) -> tuple[np.ndarray, np.ndarray]:
        """``(sign per node, transmitter label per node)`` at error rate ``p``.

        The relabelling order and the replacement labels were drawn once in
        ``__init__``, so the sets are nested in ``p`` and ``p = 0`` is exactly
        the shipped graph.
        """
        return self._relabel(int(math.floor(float(transmitter_p) * self.n)))

    def _relabel(self, n_flip: int) -> tuple[np.ndarray, np.ndarray]:
        """``(sign per node, transmitter label per node)`` with exactly
        ``n_flip`` nodes relabelled, taken from the fixed stream."""
        nt = self.nt0
        if n_flip > 0:
            nt = nt.copy()
            sel = self._relabel_order[: int(n_flip)]
            nt[sel] = self._relabel_to[sel]
        sign = np.array([self._sign.get(x, 0.0) for x in nt], dtype=float)
        return sign, nt

    def _node_arrays(self, transmitter_p: float) -> tuple[np.ndarray, np.ndarray]:
        """Cached ``(sign, gain-key index)`` per node for one error rate.

        The gain-key index points into :data:`_GAIN_COLUMNS` (with the last
        slot meaning "this transmitter has no gain", multiplier 1.0).
        """
        n_flip = int(math.floor(float(transmitter_p) * self.n))
        hit = self._cache_nt.get(n_flip)
        if hit is None:
            sign, nt = self._relabel(n_flip)
            idx = np.array(
                [_GAIN_COLUMN_INDEX.get(self._nt_gain_key.get(x, ""), len(_GAIN_COLUMNS)) for x in nt],
                dtype=np.intp,
            )
            if len(self._cache_nt) > 256:  # pragma: no cover - bounded cache
                self._cache_nt.clear()
            hit = (sign, idx)
            self._cache_nt[n_flip] = hit
        return hit

    def _effective_threshold(self, min_weight: float) -> float:
        """Snap a threshold to the smallest edge weight it keeps.

        Only the *set* of surviving edges matters, and edge weights are
        integers, so this is exact and makes the matrix cache hit.
        """
        i = int(np.searchsorted(self._unique_w, float(min_weight), side="left"))
        return float(self._unique_w[i]) if i < self._unique_w.size else float("inf")

    def _unsigned_matrix(self, min_weight: float) -> np.ndarray:
        key = self._effective_threshold(min_weight)
        M = self._cacheW.get(key)
        if M is None:
            mask = self.w >= key
            M = np.zeros((self.n, self.n), dtype=float)
            np.add.at(M, (self.post[mask], self.pre[mask]), self.w[mask])
            # a dense n x n matrix is ~10 MB here, so the cache is deliberately
            # small: holding dozens of them costs more in allocation pressure
            # than rebuilding one (a few ms).
            if len(self._cacheW) >= _MATRIX_CACHE_SIZE:
                self._cacheW.pop(next(iter(self._cacheW)))
            self._cacheW[key] = M
        return M

    def _rates(self, settings: Mapping[str, Any], gains: Mapping[str, float]) -> np.ndarray:
        """The rate update of :meth:`flylab.circuit.rate.RateNetwork.run`,
        verbatim, with the sampled graph and gains."""
        sign, key_idx = self._node_arrays(settings["transmitter_p"])
        M = self._unsigned_matrix(settings["min_weight"])
        W = self._Wbuf
        np.multiply(M, sign[None, :], out=W)
        denom = np.maximum(np.abs(W).sum(axis=1, keepdims=True), 1.0)
        W /= denom

        lam = float(settings["expression_lambda"])
        vals = np.array(
            [
                gains["g_ach"] * gains.get("ach_tone", 1.0),
                gains["g_gaba"],
                gains["g_glu"],
                gains["g_oct"],
                1.0,  # transmitters FlyLab has no receptor gain for
            ],
            dtype=float,
        )
        if lam <= 0.0:
            W *= vals[key_idx][None, :]
        else:
            # per-node (postsynaptic) receptor weight, blended towards the
            # cross-atlas expression model; lambda = 0 reproduces the uniform
            # column scaling above exactly.
            per_key = np.ones((self.n, len(vals)), dtype=float)
            for col, key in enumerate(_GAIN_COLUMNS):
                w_eff = 1.0 + lam * (self._expr[key] - 1.0) if key in self._expr else 1.0
                per_key[:, col] = 1.0 + (vals[col] - 1.0) * w_eff
            W *= per_key[:, key_idx]

        d = self.seed_mask * float(settings["drive_hz"])
        g_nav = float(gains["g_nav"])
        r = np.zeros(self.n, dtype=float)
        a = self.alpha
        for _ in range(self.steps):
            r = np.clip((1.0 - a) * r + a * (g_nav * (d + W @ r)), 0.0, self.r_max)
        return r

    # -- evaluation ------------------------------------------------------
    def gains_for(self, settings: Mapping[str, Any]) -> dict[str, float]:
        rows = self._rows(settings["potency_shift_log10"], settings["hill_mult"])
        spec = make_spec(settings["gain_transform"], float(settings["gain_coef"]))
        return spec_gains(rows, spec)

    def __call__(self, u: Sequence[float]) -> float:
        self.n_calls += 1
        settings = self.decode(u)
        if self.engine == "lif":
            return self._lif(settings)
        gains = self.gains_for(settings)
        r = self._rates(settings, gains)
        if self.readout == "mean_hz":
            return float(r.mean())
        if self.readout == "mn9_hz":
            return float(np.mean(r[self.mn9_idx])) if self.mn9_idx else float("nan")
        if self.readout == "dnp01_hz":
            return float(np.mean(r[self.dnp01_idx])) if self.dnp01_idx else float("nan")
        if self.readout == "max_hz":
            return float(r.max())
        raise ValueError(f"unknown readout {self.readout!r}")

    def _perturbed_library(self, potency_shift: float, hill_mult: float) -> dict[str, Any] | None:
        """A library copy with this compound's sourced values shifted.

        Only rows that carry a value move; a not-modelled row stays
        not-modelled, because perturbing it would invent evidence.
        """
        if potency_shift == 0.0 and hill_mult == 1.0:
            return None
        import copy

        from flylab.pharm.occupancy import load_library

        lib = copy.deepcopy(load_library())
        entry = lib.get("compounds", {}).get(str(self.compound).lower().strip())
        if entry is None:
            return None
        for spec in entry.get("receptors", {}).values():
            value = spec.get("value_M", spec.get("ec50_M"))
            if value is None:
                continue
            shifted = float(value) * (10.0 ** potency_shift)
            spec["value_M"] = shifted
            spec["ec50_M"] = shifted  # deprecated alias, kept in step
            spec["n"] = float(spec.get("n", 1.0) or 1.0) * hill_mult
        return lib

    def _lif(self, settings: Mapping[str, Any]) -> float:
        """The real Shiu-style LIF assay under the same settings.

        Two factors do not reach this path and their indices on the LIF engine
        must be read as zero by construction, not as a finding: ``expression``
        (the spiking runtime has no per-node gain hook and
        :mod:`flylab.pharm.expression` is a rate-model sensitivity layer), and
        nothing else -- potency, Hill, transformation, coefficient, drive,
        threshold, transmitter and the seed all apply.
        """
        from flylab.analysis.robustness import mechanism_spec
        from flylab.assays.spiking import run_spiking_assay

        spec = make_spec(settings["gain_transform"], float(settings["gain_coef"]))
        graph_obj = self._threshold_graph(settings)
        library = self._perturbed_library(
            float(settings["potency_shift_log10"]), float(settings["hill_mult"])
        )
        with mechanism_spec(spec):
            nb = run_spiking_assay(
                self.compound, self.conc_M,
                drive_hz=float(settings["drive_hz"]),
                seed=int(settings["lif_seed"]) % (2 ** 31 - 1),
                graph=self.graph_name,
                graph_obj=graph_obj,
                library=library,
            )
        v = (nb.get("readouts") or {}).get(self.readout)
        return float(v) if v is not None else float("nan")

    def _threshold_graph(self, settings: Mapping[str, Any]) -> dict[str, Any] | None:
        thr = float(settings["min_weight"])
        p = float(settings["transmitter_p"])
        if thr <= self.base_min_weight and p <= 0:
            return None
        g = dict(self.graph)
        g["edges"] = [e for e in self.graph["edges"] if float(e["weight"]) >= thr]
        g["n_edges"] = len(g["edges"])
        if p > 0:
            _, nt = self._sign_vector(p)
            g["nodes"] = [
                {**node, "consensus_nt": str(nt[i])} for i, node in enumerate(self.graph["nodes"])
            ]
        return g


def baseline_unit_vector() -> list[float]:
    """The point of the unit cube that reproduces the shipped FlyLab model.

    Potency, Hill, gain coefficient and drive sit at the centre of their
    ranges; the categorical factor sits in the ``flylab_biphasic`` cell; the
    expression, weight-threshold, transmitter-error and LIF-seed factors sit at
    the lower edge, which is "the assumption FlyLab actually makes".
    """
    return [float(f["baseline_u"]) for f in FACTORS]


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------
def _bootstrap(
    yA: np.ndarray, yB: np.ndarray, yAB: Sequence[np.ndarray], n_boot: int, seed: int, ci: tuple[float, float]
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    n = yA.size
    k = len(yAB)
    S_boot = np.empty((n_boot, k))
    T_boot = np.empty((n_boot, k))
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        s, t, _ = jansen_indices(yA[idx], yB[idx], [y[idx] for y in yAB])
        S_boot[b] = s
        T_boot[b] = t
    return {
        "first_ci": np.percentile(S_boot, ci, axis=0).T.tolist(),
        "total_ci": np.percentile(T_boot, ci, axis=0).T.tolist(),
    }


def sobol_analysis(
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    readout: str = "mean_hz",
    n_base: int = 128,
    seed: int = 0,
    n_boot: int = 200,
    ci: tuple[float, float] = (2.5, 97.5),
    engine: str = "rate",
    graph: str = "named",
    model: Callable[[Sequence[float]], float] | None = None,
    factor_names: Sequence[str] | None = None,
    progress: Callable[[int, int], None] | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Variance-based attribution over :data:`FACTORS` for one compound/readout.

    ``model`` and ``factor_names`` override the FlyLab model with any
    ``f: [0,1]^k -> float`` (used by the Ishigami validation test).

    Returns first-order and total-order indices with bootstrap CIs, a
    convergence diagnostic (indices recomputed at ``n_base/4`` and
    ``n_base/2``), the uncertainty budget and the readable warnings.
    """
    t0 = time.perf_counter()
    if model is None:
        surrogate = RateSurrogate(
            compound=compound, conc_M=conc_M, readout=readout, engine=engine, graph=graph, seed=seed, **kw
        )
        f = surrogate
        names = list(FACTOR_NAMES)
    else:
        surrogate = None
        f = model
        names = list(factor_names if factor_names is not None else [f"x{i}" for i in range(3)])

    k = len(names)
    A, B, AB = saltelli_matrices(int(n_base), k, seed=seed)
    total_rows = int(n_base) * (k + 2)
    done = 0

    def evaluate(M: np.ndarray) -> np.ndarray:
        nonlocal done
        out = np.empty(M.shape[0], dtype=float)
        for i in range(M.shape[0]):
            out[i] = float(f(M[i]))
            done += 1
            if progress and done % 200 == 0:
                progress(done, total_rows)
        return out

    yA = evaluate(A)
    yB = evaluate(B)
    yAB = [evaluate(M) for M in AB]

    S, T, var = jansen_indices(yA, yB, yAB)
    boot = _bootstrap(yA, yB, yAB, int(n_boot), seed + 1, ci) if n_boot else {"first_ci": None, "total_ci": None}

    # convergence: the same estimator on the first quarter and half of the rows
    convergence = []
    prev_S = prev_T = None
    for frac in (0.25, 0.5, 1.0):
        m = max(4, int(round(frac * int(n_base))))
        s, t, _ = jansen_indices(yA[:m], yB[:m], [y[:m] for y in yAB])
        row = {
            "n_base": m,
            "first_order": s.tolist(),
            "total_order": t.tolist(),
            "max_abs_delta_first": None if prev_S is None else float(np.max(np.abs(s - prev_S))),
            "max_abs_delta_total": None if prev_T is None else float(np.max(np.abs(t - prev_T))),
        }
        convergence.append(row)
        prev_S, prev_T = s, t

    rows: list[dict[str, Any]] = []
    for j, name in enumerate(names):
        meta = next((dict(fa) for fa in FACTORS if fa["name"] == name), {"name": name})
        rows.append(
            {
                "factor": name,
                "kind": meta.get("kind", "continuous"),
                "what": meta.get("what"),
                "range": meta.get("range"),
                "first_order": float(S[j]),
                "first_order_ci": (boot["first_ci"][j] if boot["first_ci"] else None),
                "total_order": float(T[j]),
                "total_order_ci": (boot["total_ci"][j] if boot["total_ci"] else None),
                "interaction": float(T[j] - S[j]),
            }
        )
        ci = rows[-1]["first_order_ci"]
        excludes_zero = (
            None if ci is None else bool(float(ci[0]) > 0.0 or float(ci[1]) < 0.0)
        )
        rows[-1]["first_order_ci_excludes_zero"] = excludes_zero
    rows.sort(key=lambda r: -r["total_order"])

    sum_first = float(sum(r["first_order"] for r in rows))
    warnings = list(BASE_WARNINGS)
    resolution = resolution_summary(rows, engine=engine)
    for row in rows:
        # the three-state view is the single source of truth for a row's state
        row["resolution_status"] = resolution["by_factor"][row["factor"]]
        row["resolved"] = bool(row["resolution_status"] == RESOLVED_LABEL)
    warnings.append(resolution["statement"])
    noise = next((r for r in rows if r["factor"] == "lif_seed"), None)
    if noise is not None and engine == "rate":
        warnings.append(
            "lif_seed has no effect on the deterministic rate engine, so its "
            f"indices (S={noise['first_order']:+.3f}, T={noise['total_order']:+.3f}) "
            "are a null-factor control: they show the SIGN and rough SIZE of the "
            "estimator's error on a factor that is exactly zero. They are one "
            "draw of that error, not a symmetric tolerance band, and an index "
            "is NOT resolved merely by exceeding them. Resolution is decided per "
            "factor by whether its own first-order confidence interval excludes "
            "zero. Run with engine='lif' to give LIF stochasticity a real index."
        )
    if engine == "lif":
        warnings.append(
            "engine='lif': the expression factor is not applied on the spiking "
            "runtime (it is a rate-model sensitivity layer), so its indices are "
            "zero by construction there, not by measurement."
        )
    last = convergence[-1]
    if last["max_abs_delta_first"] is not None and last["max_abs_delta_first"] > 0.10:
        warnings.append(
            f"not converged: doubling n_base moved a first-order index by "
            f"{last['max_abs_delta_first']:.2f}; raise n_base before quoting these shares."
        )
    if var <= 0:
        warnings.append("output variance is zero: the readout does not move over the sampled ranges.")
    if resolution["noise_floor"] > 0:
        warnings.append(
            f"estimator noise floor {resolution['noise_floor']:+.3f}, measured as "
            f"the most negative first-order estimate ({resolution['noise_floor_factor']}); "
            "a true first-order index cannot be negative, so that magnitude is a "
            "lower bound on the estimator's error at this sample size. It is NOT "
            "a symmetric tolerance band and exceeding it does not resolve a "
            "factor; only a confidence interval excluding zero does."
        )
    if float(1.0 - sum_first) != float(1.0 - sum(max(0.0, r["first_order"]) for r in rows)):
        warnings.append(
            f"interaction share is {float(1.0 - sum_first):.3f} on the raw "
            "first-order estimates and "
            f"{float(1.0 - sum(max(0.0, r['first_order']) for r in rows)):.3f} with "
            "negative estimates clipped to zero; the difference is estimator "
            "noise being counted as interaction."
        )

    return {
        "compound": compound,
        "conc_M": float(conc_M),
        "readout": readout,
        "engine": engine,
        "graph": graph,
        "n_base": int(n_base),
        "n_factors": k,
        "n_evaluations": total_rows,
        "design": "Saltelli cross-sampling on two Latin-hypercube base matrices; Jansen (1999) estimators",
        "estimator": "jansen1999",
        "seed": int(seed),
        "output_mean": float(np.mean(np.concatenate([yA, yB]))),
        "output_variance": float(var),
        "output_sd": float(math.sqrt(var)) if var > 0 else 0.0,
        "baseline_value": float(f(baseline_unit_vector())) if model is None else None,
        "sum_first_order": sum_first,
        "sum_first_order_clipped": float(sum(max(0.0, r["first_order"]) for r in rows)),
        "interaction_share": float(1.0 - sum_first),
        # negative first-order estimates are estimator noise, and leaving them
        # in inflates the interaction residual; both are reported
        "interaction_share_clipped": float(
            1.0 - sum(max(0.0, r["first_order"]) for r in rows)
        ),
        "rows": rows,
        "resolution": resolution,
        "summary": resolution["statement"],
        "convergence": convergence,
        "budget": uncertainty_budget(
            {
                "rows": rows,
                "sum_first_order": sum_first,
                "output_variance": var,
                "engine": engine,
                "resolution": resolution,
            }
        ),
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }


#: the three states a factor's first-order index can be in.  They are about
#: what this SAMPLE could resolve, not about how important the factor is.
RESOLVED_LABEL = "resolved"
UNRESOLVED_LABEL = "unresolved at this sample size"
NULL_CONTROL_LABEL = "null control"

RESOLUTION_STATES: tuple[str, ...] = (
    RESOLVED_LABEL,
    UNRESOLVED_LABEL,
    NULL_CONTROL_LABEL,
)

#: the factor that is zero by construction on the deterministic rate engine
NULL_FACTOR = "lif_seed"


def resolution_summary(
    rows: Sequence[Mapping[str, Any]],
    null_factor: str = NULL_FACTOR,
    engine: str = "rate",
) -> dict[str, Any]:
    """Which factors this sample could resolve, as states rather than numbers.

    A column of small signed first-order indices invites exactly the two
    misreadings this module has already been bitten by: reading a negative
    estimate as a negative effect, and reading "bigger than the null factor"
    as "resolved".  Each factor is therefore given one explicit state:

    ``resolved``
        its bootstrap first-order confidence interval excludes zero.  This is
        the *only* state that licenses quoting a share of the variance.
    ``null control``
        the deliberately null factor itself (the spiking RNG seed, which the
        deterministic rate engine ignores and whose true index is exactly 0),
        plus any factor whose estimate is no larger in magnitude than it and
        whose interval includes zero.  Reported separately because such a
        factor is not merely unresolved: it is indistinguishable from a
        quantity known to be zero.
    ``unresolved at this sample size``
        everything else: the interval includes zero, so the sample cannot tell
        the factor's contribution apart from nothing.

    ``null control`` is a refinement *inside* the unresolved set, never an
    escape from it: a factor in that state is also listed in ``unresolved``.
    The null factor's own index is one draw of the estimator's error, not a
    symmetric tolerance band, which is why nothing here is phrased as a
    +-floor.
    """
    rows = list(rows)
    by_name = {r["factor"]: r for r in rows}
    null_row = by_name.get(null_factor) if engine == "rate" else None
    # The estimator's floor is NOT the null factor's own index: a first-order
    # index that is truly >= 0 can only be estimated below zero by estimator
    # error, so the largest magnitude among the NEGATIVE estimates is a
    # measured lower bound on that error, and it can be much larger than the
    # null factor's own draw (on the shipped run: gain_coef at -0.046 against
    # lif_seed at -0.007, 6.8x). The floor is the larger of the two.
    negatives = [abs(float(r["first_order"])) for r in rows if float(r["first_order"]) < 0]
    noise_floor = max(negatives) if negatives else 0.0
    noise_floor_factor = None
    if negatives:
        noise_floor_factor = min(
            (r for r in rows if float(r["first_order"]) < 0),
            key=lambda r: float(r["first_order"]),
        )["factor"]
    null_mag = None
    if null_row is not None or negatives:
        null_mag = max(
            noise_floor,
            0.0 if null_row is None else abs(float(null_row["first_order"])),
        )

    out: list[dict[str, Any]] = []
    for r in rows:
        ci = r.get("first_order_ci")
        excludes_zero = (
            None if ci is None else bool(float(ci[0]) > 0.0 or float(ci[1]) < 0.0)
        )
        if null_row is not None and r["factor"] == null_factor:
            state = NULL_CONTROL_LABEL
        elif excludes_zero:
            state = RESOLVED_LABEL
        elif null_mag is not None and abs(float(r["first_order"])) <= null_mag:
            state = NULL_CONTROL_LABEL
        else:
            state = UNRESOLVED_LABEL
        out.append(
            {
                "factor": r["factor"],
                "state": state,
                "resolved": bool(state == RESOLVED_LABEL),
                "first_order": float(r["first_order"]),
                "first_order_ci": ([float(ci[0]), float(ci[1])] if ci else None),
                "first_order_ci_excludes_zero": excludes_zero,
                "total_order": float(r["total_order"]),
                "comparable_to_null_control": (
                    None
                    if null_mag is None
                    else bool(abs(float(r["first_order"])) <= null_mag)
                ),
            }
        )
    resolved = [r["factor"] for r in out if r["state"] == RESOLVED_LABEL]
    nulls = [r["factor"] for r in out if r["state"] == NULL_CONTROL_LABEL]
    unresolved = [r["factor"] for r in out if r["state"] == UNRESOLVED_LABEL]
    ranked = sorted(
        (r for r in out if r["state"] == RESOLVED_LABEL),
        key=lambda r: -r["first_order"],
    )
    if resolved:
        statement = (
            "first-order confidence intervals exclude zero for "
            + ", ".join(
                f"{r['factor']} (S1 {r['first_order']:+.3f}, CI "
                f"[{r['first_order_ci'][0]:+.3f}, {r['first_order_ci'][1]:+.3f}])"
                for r in ranked
            )
            + "; every other effect is unresolved at this sample size"
            + (
                " (" + ", ".join(sorted(unresolved + nulls)) + ")"
                if (unresolved or nulls)
                else ""
            )
            + ". 'Unresolved' means the sample cannot separate the effect from "
            "zero; it is not a claim that the effect is zero, and a negative "
            "point estimate is estimator noise, not a negative contribution."
        )
    else:
        statement = (
            "no factor's first-order confidence interval excludes zero: every "
            "effect is unresolved at this sample size, and no share of the "
            "variance should be quoted. Raise n_base."
        )
    if nulls:
        statement += (
            " Indistinguishable from the estimator's own noise (floor "
            f"{noise_floor:+.3f}, set by {noise_floor_factor}"
            + (
                f"; the deliberately null factor {null_factor} estimates "
                f"{float(null_row['first_order']):+.3f}"
                if null_row is not None
                else ""
            )
            + "): "
            + ", ".join(sorted(nulls))
            + "."
        )
    return {
        "states": out,
        "by_factor": {r["factor"]: r["state"] for r in out},
        "noise_floor": float(noise_floor),
        "noise_floor_factor": noise_floor_factor,
        "noise_floor_definition": (
            "the largest magnitude among the negative first-order estimates, "
            "which is a measured lower bound on the estimator's error because "
            "a true first-order index cannot be negative; it is compared "
            "against the null factor's own index and the larger is used"
        ),
        "resolved": resolved,
        "unresolved": sorted(unresolved + nulls),
        "unresolved_strict": sorted(unresolved),
        "null_control": sorted(nulls),
        "null_factor": (null_factor if null_row is not None else None),
        "null_factor_magnitude": null_mag,
        "statement": statement,
    }


def uncertainty_budget(result: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Share of output variance per source, plus the interaction residual.

    ``share_of_variance`` is the first-order index: the fraction of the
    readout's variance that resolving that factor *alone* would remove.  The
    shares do not sum to 1; the remainder is carried by interactions between
    factors and is reported as its own row.
    """
    var = float(result.get("output_variance", 0.0))
    rows = list(result["rows"])
    states = (
        result.get("resolution")
        or resolution_summary(rows, engine=str(result.get("engine", "rate")))
    )["by_factor"]
    out = [
        {
            "source": r["factor"],
            "kind": r.get("kind"),
            "state": states.get(r["factor"]),
            "resolved": bool(states.get(r["factor"]) == RESOLVED_LABEL),
            "share_of_variance": r["first_order"],
            "share_with_interactions": r["total_order"],
            "variance_removed": r["first_order"] * var,
            "ci": r.get("first_order_ci"),
        }
        for r in sorted(rows, key=lambda r: -r["first_order"])
    ]
    resid = 1.0 - float(result.get("sum_first_order", sum(r["first_order"] for r in rows)))
    out.append(
        {
            "source": "interactions (higher order)",
            "kind": "residual",
            "state": None,
            "resolved": None,
            "share_of_variance": resid,
            "share_with_interactions": None,
            "variance_removed": resid * var,
            "ci": None,
        }
    )
    return out


def to_markdown(result: Mapping[str, Any]) -> str:
    """The uncertainty budget as a Markdown table."""
    lines = [
        f"### Uncertainty budget - {result.get('compound')} "
        f"{result.get('conc_M'):g} M, readout {result.get('readout')} "
        f"({result.get('engine')} engine, n_base={result.get('n_base')}, "
        f"{result.get('n_evaluations')} evaluations)",
        "",
        f"output mean {result.get('output_mean'):.4g}, variance {result.get('output_variance'):.4g}",
        "",
        "| source | state | share of variance (S1) | 95% CI | with interactions (ST) |",
        "|---|---|---|---|---|",
    ]
    by_name = {r["factor"]: r for r in result["rows"]}
    for row in result["budget"]:
        ci = row.get("ci")
        ci_s = f"[{ci[0]:+.3f}, {ci[1]:+.3f}]" if ci else "-"
        st = by_name.get(row["source"], {}).get("total_order")
        lines.append(
            f"| {row['source']} | {row.get('state') or '-'} | "
            f"{row['share_of_variance']:+.3f} | {ci_s} | "
            + (f"{st:+.3f} |" if st is not None else "- |")
        )
    if result.get("summary"):
        lines += ["", f"_{result['summary']}_"]
    return "\n".join(lines)


# --------------------------------------------------------------------------
# estimator validation: the Ishigami function
# --------------------------------------------------------------------------
ISHIGAMI_A = 7.0
ISHIGAMI_B = 0.1


def ishigami(x: Sequence[float], a: float = ISHIGAMI_A, b: float = ISHIGAMI_B) -> float:
    """``sin(x1) + a sin^2(x2) + b x3^4 sin(x1)`` on ``x in [-pi, pi]^3``.

    The standard analytic benchmark for variance-based sensitivity: ``x3`` has
    a first-order index of exactly 0 but a non-zero total-order index, so it
    catches an estimator that confuses "no main effect" with "no effect".
    """
    x1, x2, x3 = float(x[0]), float(x[1]), float(x[2])
    return math.sin(x1) + a * math.sin(x2) ** 2 + b * (x3 ** 4) * math.sin(x1)


def ishigami_unit(u: Sequence[float], a: float = ISHIGAMI_A, b: float = ISHIGAMI_B) -> float:
    """:func:`ishigami` reparameterised onto the unit cube."""
    x = [(-math.pi + 2.0 * math.pi * float(ui)) for ui in u[:3]]
    return ishigami(x, a, b)


def _ishigami_reference(a: float = ISHIGAMI_A, b: float = ISHIGAMI_B) -> dict[str, list[float]]:
    p4, p8 = math.pi ** 4, math.pi ** 8
    v1 = 0.5 * (1.0 + b * p4 / 5.0) ** 2
    v2 = a * a / 8.0
    v13 = b * b * p8 * (1.0 / 18.0 - 1.0 / 50.0)
    vt = v1 + v2 + v13
    return {
        "first_order": [v1 / vt, v2 / vt, 0.0],
        "total_order": [(v1 + v13) / vt, v2 / vt, v13 / vt],
        "variance": vt,
    }


#: analytic Sobol' indices of :func:`ishigami` for ``a=7, b=0.1``
ISHIGAMI_REFERENCE: dict[str, list[float]] = _ishigami_reference()
