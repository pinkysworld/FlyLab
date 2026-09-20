r"""Rate network on the MaleCNS named-cell neighborhood.

Extracted verbatim (numerically) from the original ``flylab/assays/subgraph.py``
inline loop so that the spiking assay, the ensemble/sensitivity layer and the
analysis package all share one runtime.

Model
-----
``W_signed[j, i]`` is the signed synapse count from presynaptic node ``i`` onto
postsynaptic node ``j``.  The sign comes from the *presynaptic* consensus
transmitter (see :data:`SIGN` and, with a rationale per entry,
:data:`SIGN_TABLE`).  Rows are then normalised -- by default by
``max(sum_i |W[j, i]|, 1)`` -- so no cell can receive more than unit total
drive; see `Row normalisation`_ below, which is a substantive modelling choice
and not a detail.

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

.. _Row normalisation:

Row normalisation, and what it makes the readout measure
--------------------------------------------------------
**This step was undocumented through v0.6 and is not in the manuscript's §2.2,
§S3 or Table T2, which give the update rule without it.**  It is stated here in
full because it changes the interpretation of every rate readout in the paper.

After the signed matrix is built, each *row* (each postsynaptic cell) is
divided by its own total absolute input::

    denom[j]      = max(sum_i |W[j, i]|, 1.0)     # normalise="row_abs", default
    W_signed[j,i] = W[j, i] / denom[j]

The pharmacological gains are applied **afterwards**, as a presynaptic *column*
scale (:meth:`RateNetwork.scale_vector`).  Cell ``j``'s recurrent input is
therefore

.. math::

    \\sum_i \\frac{|W_{ji}|}{\\sum_k |W_{jk}|}\; \\mathrm{sign}(nt_i)\\,
    g(nt_i)\\, r_i ,

that is, **a composition-weighted average of its presynaptic gains**, the
weights being the share of ``j``'s own total input that each presynaptic
transmitter class supplies.  Every cell is, in effect, a unit-gain integrator
of its own input mix.

Four consequences, all of which a reader of the results needs:

1. *Absolute connection strength is removed from a cell's operating point.*  A
   cell receiving 10 000 synapses and a cell receiving 10 behave identically if
   their input *composition* is identical.  What survives the normalisation is
   which cell contacts which, and in what proportion, not how strongly.
2. *The drug enters each cell only through a weighted mean of transmitter
   gains, with weights that are a property of the cut and not of the compound.*
   Two compounds with the same gain vector are indistinguishable at every cell,
   and the network effect is close to a smooth function of
   ``(g_ach, g_gaba, g_glu, g_oct, g_nav)`` alone -- which is precisely what the
   composition-only ablation level (``B_composition_only``) computes from the
   global census.  The normalisation therefore makes "composition-dominated"
   close to a built-in property of the readout, and the second research
   question's composition-dominated majority must be read as a statement about
   *this* readout under *this* normalisation.  It is not, by itself, evidence
   that the wiring does not matter.
3. *It is what makes the iteration stable and the statistics affordable.*  Row
   sums of ``|W_signed|`` are at most 1, so the recurrent operator is a
   contraction (spectral radius 0.51 on the ``named`` cut) and 80 steps settle.
   The same matrix unnormalised has spectral radius ~150; the iteration is then
   supercritical and the clip at ``r_max`` does the work that the normalisation
   used to do (at 1 uM fipronil, 184 of 1126 cells sit at ``r_max``).
4. *It has no measurement behind it.*  It is an asserted modelling choice of
   the same kind as the transmitter signs, and it is currently outside the
   prespecified specification-robustness family.

Alternatives are therefore exposed by name rather than hard-wired
(:data:`NORMALISATIONS`, ``RateNetwork(..., normalise=...)``):

==============  ==========================================================
``normalise``   ``denom[j]``, and what it preserves
==============  ==========================================================
``"row_abs"``   ``max(sum_i |W[j, i]|, 1.0)``.  The shipped default and the
                one every published FlyLab number was computed with.  Input
                *composition* only; total input strength is divided out.
``"none"``      ``1.0``.  Raw signed synapse counts.  Strongly connected
                cells dominate, the operator is supercritical, and the
                readout is shaped by the ``r_max`` clip.
``"degree"``    ``max(k_in[j], 1.0)`` with ``k_in`` the number of
                presynaptic partners.  A degree-corrected middle case: a
                cell with strong inputs stays strongly driven (relative
                weights are preserved, unlike ``row_abs``) but hubs are
                discounted by their fan-in.
==============  ==========================================================

:func:`readout_under_normalisations` runs one compound under all three and is
the function the manuscript should quote when it states whether a result
depends on the normalisation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

#: Per-transmitter signs used to build the signed weight matrix, with a
#: rationale and an evidence status for each entry.  Five of these magnitudes
#: (glutamate, histamine, dopamine, serotonin, octopamine) appear nowhere in
#: the manuscript; they are published here so that they can be quoted.
#:
#: ``evidence`` is one of:
#:
#: * ``"convention"`` -- a sign convention that is standard for insect circuits
#:   and is stated as such in the paper;
#: * ``"asserted"`` -- a magnitude chosen by this project with no measurement
#:   and no citation behind it.
#:
#: **None of these coefficients is currently a member of the prespecified
#: specification-robustness family (§2.4), which varies only the
#: engagement-to-gain transformation.**  The sign-permutation null model bounds
#: how much the *assignment* of transmitters to cells matters; it does not
#: bound these numbers, because it permutes labels and leaves the table itself
#: untouched.
SIGN_TABLE: dict[str, dict[str, Any]] = {
    "acetylcholine": {
        "sign": 1.0,
        "role": "fast excitatory",
        "evidence": "convention",
        "rationale": (
            "Acetylcholine is the principal fast excitatory transmitter of the "
            "insect CNS; +1 is the reference against which every other entry is "
            "scaled, so its magnitude is a normalisation, not a claim."
        ),
    },
    "gaba": {
        "sign": -1.0,
        "role": "fast inhibitory",
        "evidence": "convention",
        "rationale": (
            "GABA acts on RDL, a GABA-gated chloride channel; inhibitory at "
            "full weight, i.e. taken to be as strong per synapse as ACh."
        ),
    },
    "glutamate": {
        "sign": -0.4,
        "role": "inhibitory (GluCl convention), attenuated",
        "evidence": "asserted",
        "rationale": (
            "Insect glutamate is signed inhibitory at GluCl-bearing targets "
            "(the convention the manuscript states in §4). The magnitude 0.4 "
            "rather than 1.0 is an asserted hedge against the cells whose "
            "glutamatergic output is excitatory, for which the global sign is "
            "simply wrong. No measurement supports the value 0.4."
        ),
    },
    "histamine": {
        "sign": -0.5,
        "role": "fast inhibitory (photoreceptor/mechanosensory)",
        "evidence": "asserted",
        "rationale": (
            "Histamine gates chloride channels (ort/HisCl) and is inhibitory; "
            "the magnitude 0.5 is asserted, chosen below GABA because "
            "histaminergic drive in this cut is sparse and its postsynaptic "
            "channel complement is not known per cell."
        ),
    },
    "dopamine": {
        "sign": 0.2,
        "role": "modulatory, net excitatory in this model",
        "evidence": "asserted",
        "rationale": (
            "Dopamine is modulatory and its sign at a given synapse depends on "
            "the receptor subtype expressed there, which this model does not "
            "represent. +0.2 is an asserted small positive stand-in; a rate "
            "model has no slow-modulation machinery to put it in instead."
        ),
    },
    "serotonin": {
        "sign": 0.2,
        "role": "modulatory, net excitatory in this model",
        "evidence": "asserted",
        "rationale": (
            "As for dopamine: receptor-subtype dependent in reality, collapsed "
            "here to an asserted small positive weight."
        ),
    },
    "octopamine": {
        "sign": 0.2,
        "role": "modulatory, net excitatory in this model",
        "evidence": "asserted",
        "rationale": (
            "Octopamine is the insect arousal amine and is treated as weakly "
            "excitatory; +0.2 is asserted and is the sign on which the "
            "formamidine (``g_oct``) mechanism acts."
        ),
    },
}

#: transmitter -> signed weight multiplier (derived from :data:`SIGN_TABLE` so
#: the two cannot drift).  A transmitter absent from this table contributes 0.
SIGN: dict[str, float] = {k: float(v["sign"]) for k, v in SIGN_TABLE.items()}


def sign_table_rows() -> list[dict[str, Any]]:
    """:data:`SIGN_TABLE` as a list of records, for the paper and the UI.

    Every row carries ``transmitter``, ``sign``, ``role``, ``evidence`` and
    ``rationale``.  The five rows with ``evidence == "asserted"`` are the
    coefficients that no manuscript table currently prints and that no
    robustness family currently varies.
    """
    return [
        {"transmitter": nt, "sign": float(row["sign"]), "role": row["role"],
         "evidence": row["evidence"], "rationale": row["rationale"],
         "in_specification_family": False}
        for nt, row in SIGN_TABLE.items()
    ]


#: available row-normalisation modes (see the module docstring)
NORMALISATIONS: tuple[str, ...] = ("row_abs", "none", "degree")
DEFAULT_NORMALISATION = "row_abs"

#: one-line description of each mode, for notebooks and tables
NORMALISATION_NOTES: dict[str, str] = {
    "row_abs": (
        "denom[j] = max(sum_i |W[j, i]|, 1): each cell's recurrent input is a "
        "composition-weighted average of its presynaptic gains. Shipped "
        "default; every published FlyLab number uses it."
    ),
    "none": (
        "denom[j] = 1: raw signed synapse counts. The recurrent operator is "
        "supercritical (spectral radius ~150 on the named cut) and the r_max "
        "clip, not the wiring, bounds the rates."
    ),
    "degree": (
        "denom[j] = max(in-degree[j], 1): degree-corrected. Relative input "
        "strengths survive (unlike row_abs) but hubs are discounted by fan-in."
    ),
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

#: named graph -> committed JSON filename
GRAPHS: dict[str, str] = {
    "named": "malecns_named_neighborhood.json",
    "taste_motor": "malecns_taste_motor_neighborhood.json",
}
DEFAULT_GRAPH = "named"

#: directories searched for a committed graph JSON, in order
DATA_DIRS = [
    Path("data/derived"),
    Path(__file__).resolve().parents[2] / "data/derived",
    Path.home() / ".flylab",
]

CANDIDATES = [d / GRAPHS["named"] for d in DATA_DIRS]


def resolve_graph(name_or_path: str | Path | None = None) -> Path:
    """Resolve ``"named"`` / ``"taste_motor"`` / an explicit path to a JSON file.

    ``None`` resolves to the default (``"named"``) graph, which is the one the
    v0.4 regression numbers were taken on.
    """
    key = name_or_path if name_or_path is not None else DEFAULT_GRAPH
    if isinstance(key, str) and key in GRAPHS:
        for d in DATA_DIRS:
            p = d / GRAPHS[key]
            if p.exists():
                return p
        raise FileNotFoundError(
            f"graph {key!r} ({GRAPHS[key]}) missing. Run the extract-malecns-subgraph Action."
        )
    p = Path(key)
    if p.exists():
        return p
    # unknown string that is not a path: fall back to the default graph
    for d in DATA_DIRS:
        q = d / GRAPHS[DEFAULT_GRAPH]
        if q.exists():
            return q
    raise FileNotFoundError("neighborhood JSON missing. Run extract-malecns-subgraph Action.")


def graph_path(path: str | Path | None = None) -> Path:
    return resolve_graph(path)


_GRAPH_CACHE: dict[str, dict[str, Any]] = {}


def load_graph(path: str | Path | None = None) -> dict[str, Any]:
    """Load (and memoise) the committed neighborhood JSON.

    The returned dict is shared; callers must treat it as read-only.
    """
    key = str(resolve_graph(path).resolve())
    if key not in _GRAPH_CACHE:
        _GRAPH_CACHE[key] = json.loads(Path(key).read_text())
    return _GRAPH_CACHE[key]


def forget_graph(path: str | Path) -> None:
    """Drop one graph (and every normalisation of it) from the memo caches.

    Needed when a file is written, read back and then rewritten in the same
    process -- the extractor does this to fold a cut's census into its own
    metadata.
    """
    key = str(Path(path).resolve())
    _GRAPH_CACHE.pop(key, None)
    for k in [k for k in _NET_CACHE if k and k[0] == key]:
        _NET_CACHE.pop(k, None)


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
        raw = row.get("engagement", row.get("occupancy"))
        if raw is None:
            # schema v3: no sourced value -> receptor not modelled, gain unchanged
            continue
        th = float(raw)
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
    r"""Row-normalised signed rate network over the neighborhood graph.

    ``W_signed[j, i]`` is ``sign(nt_i) * synapses(i->j) / denom[j]``, and the
    pharmacological gains multiply the *columns* afterwards.  Cell ``j``'s
    recurrent input is therefore a **composition-weighted average of its
    presynaptic gains**: with the default ``normalise="row_abs"`` the weights
    are the share of ``j``'s total absolute input carried by each presynaptic
    transmitter class, so the drug reaches ``j`` only through that weighted
    mean and not through how strongly ``j`` is wired.

    This is the step the manuscript does not state (§2.2, §S3, T2 give the
    update rule without it), and it is the reason a "composition-dominated"
    verdict is close to built into the readout rather than discovered in it:
    two compounds with the same gain vector are indistinguishable at every
    cell, and the network response is nearly a smooth function of the gain
    vector alone -- the same object the composition-only ablation level
    computes from the global census.  See the module docstring for the full
    statement and for what each alternative preserves.

    Args:
        graph: a committed neighborhood graph dict (nodes/edges/seeds).
        sign: transmitter -> signed multiplier; defaults to :data:`SIGN`,
            whose per-entry rationale is in :data:`SIGN_TABLE`.
        normalise: one of :data:`NORMALISATIONS`.  ``"row_abs"`` (default)
            divides each row by ``max(sum_i |W[j, i]|, 1)`` and is what every
            published FlyLab number was computed with; ``"none"`` leaves the
            raw signed synapse counts (supercritical: the ``r_max`` clip, not
            the wiring, bounds the rates); ``"degree"`` divides by the number
            of presynaptic partners, keeping relative input strengths while
            discounting hubs by their fan-in.

    Attributes:
        W_raw: the signed, *un*-normalised matrix.
        W_signed: the matrix actually iterated (``W_raw / denom``).
        row_denominator: the per-row denominator that was applied.
        normalise: the mode name.
        in_degree: number of presynaptic partners per node.
    """

    def __init__(
        self,
        graph: dict[str, Any],
        sign: dict[str, float] | None = None,
        normalise: str = DEFAULT_NORMALISATION,
    ):
        if normalise not in NORMALISATIONS:
            raise ValueError(
                f"unknown normalise={normalise!r}; expected one of {NORMALISATIONS}"
            )
        self.graph = graph
        self.normalise = normalise
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
        self.W_raw = W
        self.in_degree = (np.abs(W) > 0).sum(axis=1)
        denom = self._row_denominator(W)
        self.row_denominator = denom.ravel()
        self.W_signed = W / denom
        self._nt_vec = np.array([self.nt_of_node[b] for b in self.body_ids], dtype=object)

    # -- normalisation ----------------------------------------------------
    def _row_denominator(self, W: np.ndarray) -> np.ndarray:
        """Per-row divisor for the selected mode, as an ``(n, 1)`` column.

        ``row_abs`` (default) returns ``max(sum_i |W[j, i]|, 1)``: the readout
        then measures each cell's input *composition* and not its input
        strength.  ``none`` returns 1 (raw synapse counts).  ``degree`` returns
        ``max(in-degree, 1)``, which keeps relative input strengths but
        discounts fan-in.
        """
        if self.normalise == "none":
            return np.ones((W.shape[0], 1), dtype=float)
        if self.normalise == "degree":
            return np.maximum(self.in_degree.reshape(-1, 1).astype(float), 1.0)
        return np.maximum(np.abs(W).sum(axis=1, keepdims=True), 1.0)

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


_NET_CACHE: dict[tuple[str, str], RateNetwork] = {}


def rate_network(
    path: str | Path | None = None, normalise: str = DEFAULT_NORMALISATION
) -> RateNetwork:
    """Memoised :class:`RateNetwork` for a committed neighborhood graph.

    ``path`` may be ``"named"``, ``"taste_motor"``, an explicit path, or None.
    ``normalise`` selects the row normalisation (:data:`NORMALISATIONS`); each
    mode is cached separately and the default is the shipped ``"row_abs"``.
    """
    key = (str(resolve_graph(path).resolve()), str(normalise))
    if key not in _NET_CACHE:
        _NET_CACHE[key] = RateNetwork(load_graph(key[0]), normalise=normalise)
    return _NET_CACHE[key]


#: readouts :func:`readout_under_normalisations` reports, all in Hz
NORMALISATION_READOUTS: tuple[str, ...] = ("mean_hz", "mn9_hz", "dnp01_hz", "max_hz")


def _subgraph_readouts(net: RateNetwork, rates: np.ndarray) -> dict[str, float | None]:
    """``mean_hz`` / ``mn9_hz`` / ``dnp01_hz`` / ``max_hz`` for one run."""
    seeds = net.graph.get("seeds", {})

    def mean_of(typ: str) -> float | None:
        idx = [net.node_index[b] for b in seeds.get(typ, []) if b in net.node_index]
        return float(np.mean(rates[idx])) if idx else None

    return {
        "mean_hz": float(rates.mean()),
        "mn9_hz": mean_of("MN9"),
        "dnp01_hz": mean_of("DNp01"),
        "max_hz": float(rates.max()),
    }


def effect_under_normalisation(
    gains: dict[str, float] | None,
    graph: str | Path | None = None,
    normalise: str = DEFAULT_NORMALISATION,
    drive_hz: float = 40.0,
    steps: int = 80,
    readout: str = "mean_hz",
) -> float | None:
    """``treated - vehicle`` for an arbitrary gain dict under one normalisation.

    The same contrast the assays report (both arms on the same graph, seeds
    driven at ``drive_hz``), exposed so that the analysis layer can ask what a
    result looks like when the row normalisation is changed.  With the defaults
    it reproduces ``run_subgraph_assay``'s ``mean_hz`` contrast exactly.
    """
    net = rate_network(graph, normalise=normalise)
    drive = net.drive_vector(net.seed_drive(drive_hz))
    treated = _subgraph_readouts(net, net.run(drive, gains, steps=steps))[readout]
    vehicle = _vehicle_readouts(net, drive_hz, steps)[readout]
    if treated is None or vehicle is None:
        return None
    return float(treated - vehicle)


_VEHICLE_CACHE: dict[tuple[int, str, float, int], dict[str, float | None]] = {}


def _vehicle_readouts(
    net: RateNetwork, drive_hz: float = 40.0, steps: int = 80
) -> dict[str, float | None]:
    """Memoised vehicle (all gains 1.0) readouts for one network and drive.

    The vehicle arm is identical for every compound, so the ablation and
    normalisation sweeps would otherwise re-run it thousands of times.
    """
    key = (id(net), net.normalise, float(drive_hz), int(steps))
    if key not in _VEHICLE_CACHE:
        drive = net.drive_vector(net.seed_drive(drive_hz))
        _VEHICLE_CACHE[key] = _subgraph_readouts(net, net.run(drive, DEFAULT_GAINS, steps=steps))
    return _VEHICLE_CACHE[key]


def readout_under_normalisations(
    compound: str | None,
    conc_M: float = 1e-6,
    graph: str | Path | None = None,
    drive_hz: float = 40.0,
    steps: int = 80,
    normalisations: tuple[str, ...] | list[str] = NORMALISATIONS,
    library: dict[str, Any] | None = None,
    gains: dict[str, float] | None = None,
) -> dict[str, Any]:
    """One compound's rate readouts under every row normalisation.

    The row normalisation (module docstring) is an undocumented modelling
    choice that decides what the rate readout measures: under the shipped
    ``"row_abs"`` mode a cell's recurrent input is a composition-weighted
    average of its presynaptic gains, which makes a composition-dominated
    verdict close to a property of the readout.  This function runs the same
    compound, the same concentration and the same drug contrast
    (``treated - vehicle``) under each mode so a claim can be checked against
    the alternatives rather than asserted under one of them.

    Returns a JSON-serialisable block with, per mode, the vehicle and treated
    readouts, the effect, how many cells sit at ``r_max`` (the clip does the
    stabilising that ``row_abs`` otherwise does), and a ``statements`` list
    saying whether the *direction* of each readout survives.

    ``gains`` overrides the compound's own gains (used by the ablation layer).
    """
    modes = [m for m in normalisations]
    for m in modes:
        if m not in NORMALISATIONS:
            raise ValueError(f"unknown normalisation {m!r}; expected one of {NORMALISATIONS}")
    if gains is None:
        gains, _occ = compute_gains(compound, conc_M, library=library)
    out: dict[str, Any] = {}
    for mode in modes:
        net = rate_network(graph, normalise=mode)
        drive = net.drive_vector(net.seed_drive(drive_hz))
        r_t = net.run(drive, gains, steps=steps)
        r_v = net.run(drive, DEFAULT_GAINS, steps=steps)
        treated = _subgraph_readouts(net, r_t)
        vehicle = _subgraph_readouts(net, r_v)
        out[mode] = {
            "normalise": mode,
            "note": NORMALISATION_NOTES[mode],
            "vehicle": vehicle,
            "treated": treated,
            "effect": {
                k: (None if treated[k] is None or vehicle[k] is None else float(treated[k] - vehicle[k]))
                for k in NORMALISATION_READOUTS
            },
            "n_at_r_max_treated": int((r_t >= 299.999).sum()),
            "n_at_r_max_vehicle": int((r_v >= 299.999).sum()),
            "row_denominator_median": float(np.median(net.row_denominator)),
        }

    ref = DEFAULT_NORMALISATION if DEFAULT_NORMALISATION in out else modes[0]
    statements: list[str] = []
    agree: dict[str, bool] = {}
    for key in NORMALISATION_READOUTS:
        vals = {m: out[m]["effect"][key] for m in modes if out[m]["effect"][key] is not None}
        if len(vals) < 2:
            continue
        signs = {np.sign(round(v, 9)) for v in vals.values()}
        same = len(signs) == 1
        agree[key] = bool(same)
        rendered = ", ".join(f"{m} {vals[m]:+.3f}" for m in modes if m in vals)
        statements.append(
            f"{compound or 'vehicle'} at {conc_M:.1e} M, {key}: "
            + ("the sign of the effect is the same under every normalisation"
               if same else
               "THE SIGN OF THE EFFECT DEPENDS ON THE NORMALISATION")
            + f" ({rendered} Hz)."
        )
    statements.append(
        "Evidence against a result being normalisation-independent would be a "
        "sign flip between modes, or a mode in which the effect is within "
        "numerical noise of zero while it is not under 'row_abs'."
    )
    return {
        "compound": compound,
        "conc_M": float(conc_M),
        "graph": str(resolve_graph(graph).name),
        "gains": dict(gains),
        "reference_normalisation": ref,
        "normalisations": out,
        "effect_sign_agrees": agree,
        "statements": statements,
        "label": "model_derived",
        "warnings": [
            "The row normalisation is an asserted modelling choice, is absent "
            "from the manuscript's 2.2/S3/T2, and is not a member of the "
            "specification-robustness family.",
            "Under 'none' and 'degree' the recurrent operator is supercritical "
            "on this cut, so the r_max clip shapes the readout; the levels are "
            "not comparable across modes, only the signs and the orderings are.",
        ],
    }


# --------------------------------------------------------------------------
# sparse (edge-list) runtime, for cuts too large for a dense N x N matrix
# --------------------------------------------------------------------------
#: Default number of edges processed per scatter chunk.  The recurrent step
#: builds two temporaries of this length; 4M edges is ~64 MB of float64, which
#: keeps the whole-CNS run inside a few GB instead of tens.
SPARSE_CHUNK = 4_000_000


def sparse_memory_gb(n_nodes: int, n_edges: int, chunk: int = SPARSE_CHUNK) -> float:
    """Resident bytes (GB) :class:`SparseRateNetwork` needs, as a formula.

    Per edge, 32 B: ``pre`` and ``post`` as int32 (8 B), the raw weight, the
    signed weight and the gain-scaled ``weff`` as float64 (24 B).  Per node,
    40 B: five float64 vectors (rates, drive, scale, denominator, scratch).
    Plus two ``chunk``-length float64 temporaries for the scatter.

    This is the engine's own footprint.  A whole-CNS *load* peaks higher
    (measured 1.9 GB against 0.9 GB here) because the streaming reader holds
    int64 body-id arrays and the annotation frames before it hands them over;
    :func:`flylab.assays.fullcns.run_fullcns_assay` reports the measured
    ``ru_maxrss`` beside this prediction rather than instead of it.
    """
    return float(32 * int(n_edges) + 40 * int(n_nodes) + 16 * int(chunk)) / 1e9


class SparseRateNetwork:
    r"""The same rate model as :class:`RateNetwork`, on an edge list.

    :class:`RateNetwork` materialises a dense ``N x N`` matrix, which is fine
    for the committed 1-2k-node cuts (10 MB) and impossible at whole-CNS scale
    (165 122 nodes would be 218 TB).  This class holds the connectome as the
    three arrays the null-model fast path already uses -- ``pre``, ``post``,
    ``weight`` -- and writes the recurrent step as a ``bincount`` scatter::

        rec[j] = sum_{e: post[e] == j} weff[e] * r[pre[e]]
        r     <- clip((1 - alpha) * r + alpha * g_nav * (drive + rec), 0, r_max)

    with ``weff[e] = sign(nt[pre[e]]) * w[e] * scale[pre[e]] / denom[post[e]]``:
    identical arithmetic to the dense engine, differing only in summation
    order (agreement is ~1e-12 Hz; ``tests/test_circuit_rate.py`` pins it).

    **This is an addition, not a change.**  Nothing in the shipped pipeline
    routes through it; :class:`RateNetwork` is untouched and remains the
    default everywhere, so the frozen numerical regression is unaffected.

    Args:
        n: number of nodes.
        pre / post: edge endpoint *indices* (not body ids).
        weight: synapse counts, unsigned.
        nt: per-node consensus transmitter labels (object or str array).
        body_ids: optional per-node body ids, for ``rates_by_body``.
        seeds: optional ``{type: [bodyId, ...]}`` as the graph JSON has it.
        normalise: one of :data:`NORMALISATIONS`; the same three modes and the
            same denominators as the dense engine.
        chunk: edges per scatter chunk (bounds the temporaries).
    """

    def __init__(
        self,
        n: int,
        pre: np.ndarray,
        post: np.ndarray,
        weight: np.ndarray,
        nt: np.ndarray | list[Any],
        body_ids: np.ndarray | list[int] | None = None,
        seeds: dict[str, list[int]] | None = None,
        sign: dict[str, float] | None = None,
        normalise: str = DEFAULT_NORMALISATION,
        chunk: int = SPARSE_CHUNK,
    ) -> None:
        if normalise not in NORMALISATIONS:
            raise ValueError(
                f"unknown normalise={normalise!r}; expected one of {NORMALISATIONS}"
            )
        self.n = int(n)
        self.pre = np.asarray(pre, dtype=np.int32)
        self.post = np.asarray(post, dtype=np.int32)
        self.w = np.asarray(weight, dtype=np.float64)
        self.nt = np.asarray(list(nt), dtype=object)
        self.body_ids = (
            np.arange(self.n, dtype=np.int64)
            if body_ids is None
            else np.asarray(body_ids, dtype=np.int64)
        )
        self.node_index = {int(b): i for i, b in enumerate(self.body_ids.tolist())}
        self.seeds = {k: [int(v) for v in vals] for k, vals in (seeds or {}).items()}
        self.sign = dict(sign or SIGN)
        self.normalise = normalise
        self.chunk = int(chunk)
        sign_vec = np.array(
            [self.sign.get(x, 0.0) for x in self.nt.tolist()], dtype=np.float64
        )
        self.sw = sign_vec[self.pre] * self.w
        # the dense engine counts a presynaptic *partner* only when its signed
        # contribution is non-zero (a transmitter absent from SIGN contributes
        # 0), so the degree normalisation must count the same thing
        self.in_degree = np.bincount(
            self.post[self.sw != 0.0], minlength=self.n
        ) if self.sw.size else np.zeros(self.n, dtype=np.int64)
        self.row_denominator = self._denominator()

    # -- construction -----------------------------------------------------
    @classmethod
    def from_graph(
        cls,
        graph: dict[str, Any],
        normalise: str = DEFAULT_NORMALISATION,
        chunk: int = SPARSE_CHUNK,
    ) -> "SparseRateNetwork":
        """Build from a committed neighborhood graph dict."""
        nodes = graph["nodes"]
        body_ids = np.array([int(nd["bodyId"]) for nd in nodes], dtype=np.int64)
        index = {int(b): i for i, b in enumerate(body_ids.tolist())}
        nt = [nd.get("consensus_nt") or "unclear" for nd in nodes]
        pre, post, w = [], [], []
        for e in graph["edges"]:
            i, j = index.get(int(e["pre"])), index.get(int(e["post"]))
            if i is None or j is None:
                continue
            pre.append(i)
            post.append(j)
            w.append(float(e["weight"]))
        return cls(
            n=len(nodes),
            pre=np.asarray(pre, dtype=np.int32),
            post=np.asarray(post, dtype=np.int32),
            weight=np.asarray(w, dtype=np.float64),
            nt=nt,
            body_ids=body_ids,
            seeds=graph.get("seeds") or {},
            normalise=normalise,
            chunk=chunk,
        )

    # -- helpers ----------------------------------------------------------
    def __len__(self) -> int:
        return self.n

    def _denominator(self) -> np.ndarray:
        if self.normalise == "none":
            return np.ones(self.n, dtype=np.float64)
        if self.normalise == "degree":
            return np.maximum(self.in_degree.astype(np.float64), 1.0)
        return np.maximum(
            np.bincount(self.post, weights=np.abs(self.sw), minlength=self.n), 1.0
        )

    def scale_vector(self, gains: dict[str, float] | None) -> np.ndarray:
        g = dict(DEFAULT_GAINS)
        g.update({k: float(v) for k, v in (gains or {}).items() if k in DEFAULT_GAINS})
        out = np.ones(self.n, dtype=np.float64)
        for idx, nt in enumerate(self.nt.tolist()):
            key = NT_GAIN_KEY.get(nt)
            if key is None:
                continue
            val = g[key]
            if nt == "acetylcholine":
                val *= g["ach_tone"]
            out[idx] = val
        return out

    def drive_vector(
        self, drive: dict[int, float] | None, default_hz: float = 0.0
    ) -> np.ndarray:
        d = np.full(self.n, float(default_hz), dtype=np.float64)
        for body_id, hz in (drive or {}).items():
            i = self.node_index.get(int(body_id))
            if i is not None:
                d[i] = float(hz)
        return d

    def seed_indices(self, types: Any = None) -> list[int]:
        names = list(self.seeds) if types is None else list(types)
        out: list[int] = []
        for t in names:
            for b in self.seeds.get(t, []):
                i = self.node_index.get(int(b))
                if i is not None and i not in out:
                    out.append(i)
        return out

    def seed_drive(self, drive_hz: float, types: Any = None) -> dict[int, float]:
        idx = self.seed_indices(types)
        return {int(self.body_ids[i]): float(drive_hz) for i in idx}

    def memory_gb(self) -> float:
        return sparse_memory_gb(self.n, int(self.pre.size), self.chunk)

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
        scale = self.scale_vector(g)
        weff = self.sw * scale[self.pre] / self.row_denominator[self.post]
        g_nav = g["g_nav"]
        pre, post, n = self.pre, self.post, self.n
        m = int(pre.size)
        chunk = max(int(self.chunk), 1)
        r = np.zeros(n, dtype=np.float64)
        for _ in range(int(steps)):
            if m <= chunk:
                rec = np.bincount(post, weights=weff * r[pre], minlength=n)
            else:
                rec = np.zeros(n, dtype=np.float64)
                for a in range(0, m, chunk):
                    b = min(a + chunk, m)
                    rec += np.bincount(
                        post[a:b], weights=weff[a:b] * r[pre[a:b]], minlength=n
                    )
            r = np.clip((1.0 - alpha) * r + alpha * (g_nav * (d + rec)), 0.0, r_max)
        return r

    def rates_by_body(self, rates: np.ndarray) -> dict[int, float]:
        return {int(b): float(rates[i]) for i, b in enumerate(self.body_ids.tolist())}


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
    "SIGN_TABLE",
    "sign_table_rows",
    "NORMALISATIONS",
    "DEFAULT_NORMALISATION",
    "NORMALISATION_NOTES",
    "NORMALISATION_READOUTS",
    "effect_under_normalisation",
    "readout_under_normalisations",
    "NT_GAIN_KEY",
    "GRAPHS",
    "DEFAULT_GRAPH",
    "resolve_graph",
    "DEFAULT_GAINS",
    "RateNetwork",
    "SparseRateNetwork",
    "SPARSE_CHUNK",
    "sparse_memory_gb",
    "forget_graph",
    "rate_network",
    "load_graph",
    "graph_path",
    "seed_ids",
    "compute_gains",
    "fallback_gains",
    "by_superclass",
    "top_changed",
]
