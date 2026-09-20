r"""Whole-CNS rate assay: one run over the complete traced MaleCNS connectome.

What this is
------------
Every other FlyLab circuit assay runs on a *cut* -- a few thousand cells chosen
around a seed set.  This one runs the same signed rate model over all 165 122
traced cells and every synapse between them, by streaming the 1.0 GB weight
matrix batch by batch into the edge arrays
:class:`flylab.circuit.rate.SparseRateNetwork` iterates.  The matrix is never
loaded through pandas (``pd.read_feather`` on it costs ~3.6 GB of int64 columns
before anything useful happens) and no dense ``N x N`` matrix is ever formed:
at this scale that would be 218 GB of float64 for a matrix that is 99.9 % zeros.

What this is **not**
--------------------
It is a *single forward run*, and that is the whole of it:

* **No dependence testing.**  A 1000-permutation connectome-dependence profile
  at this scale is days of compute (see
  :func:`flylab.analysis.scale.feasibility_frontier`), so none of the verdicts
  the rest of FlyLab reports -- distinguishable / equivalent within tolerance /
  indeterminate -- exists for this readout.
* **No null models.**  Nothing here has been compared against a rewired,
  weight-permuted, transmitter-permuted or Erdos-Renyi control, so there is no
  evidence that any number below depends on the wiring rather than on the
  transmitter census.
* **No specification robustness.**  The mechanism-rule family
  (:mod:`flylab.analysis.robustness`) is not run, so the sign and the size of
  the effect under the other admissible specifications are unknown here.
* **No uncertainty budget, no ablation ladder, no value-of-information.**

The gain patch is applied **uniformly** to every cell carrying a given
transmitter.  Receptor expression is unmapped for about four fifths of these
cells: of the 165 122 traced cells, a MaleCNS superclass with *any* published
nAChR statement covers 27.4 %, RDL 25.1 %, GluCl 8.0 %, and OctR and AChE 0 %
(:mod:`flylab.pharm.expression`, which is itself a coarse cross-atlas
inference).  The uniform patch is therefore not a per-cell model of drug
action; it is one global assumption applied 165 122 times.

The result is consequently **a prediction that cannot be interrogated by the
instruments the rest of FlyLab relies on.**  Use it to see what the model does
when the cut boundary is removed, not as evidence about the connectome.
"""
from __future__ import annotations

import resource
import time
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from flylab.circuit.rate import (
    DEFAULT_GAINS,
    DEFAULT_NORMALISATION,
    SPARSE_CHUNK,
    SparseRateNetwork,
    compute_gains,
    sparse_memory_gb,
)
from flylab.maps.malecns import (
    FILES,
    MAP_CITATION,
    MAP_ID,
    N_EDGES_V1,
    NAMED_TYPES,
    default_dir,
    edge_arrays,
    estimate_edge_memory_gb,
    traced_body_ids,
    weights_path,
)
from flylab.notebook.schema import empty_notebook, set_map

#: Refuse rather than thrash: the default ceiling on the edge arrays alone.
#: ``min_weight=1`` (the entire matrix) needs ~1.8 GB of arrays and ~4.5 GB
#: resident once the signed and gain-scaled copies exist, so 8 GB is a real
#: ceiling on a 16 GB machine and not a formality.
DEFAULT_MEMORY_CEILING_GB = 8.0

#: Default synapse-count floor.  ``1`` means *every* traced-to-traced edge in
#: v1.0; the committed cuts use 5.  This is the one knob that trades fidelity
#: for memory, and it is recorded in the notebook.
DEFAULT_MIN_WEIGHT = 1

FULLCNS_WARNINGS: list[str] = [
    "WHOLE-CNS SINGLE RUN: no dependence testing, no null models and no "
    "specification-robustness analysis are available at this scale. Nothing "
    "in this notebook has been compared against a rewired, weight-permuted, "
    "transmitter-permuted or Erdos-Renyi control, so there is no evidence "
    "here that any readout depends on the wiring rather than on the "
    "transmitter census.",
    "The gain patch is applied UNIFORMLY across cells. Receptor expression is "
    "unmapped for about four fifths of them: a MaleCNS superclass with any "
    "published statement covers 27.4 % of traced cells for nAChR, 25.1 % for "
    "RDL, 8.0 % for GluCl and 0 % for OctR and AChE. The patch is one global "
    "assumption applied 165 122 times, not a per-cell model of drug action.",
    "This result is therefore A PREDICTION THAT CANNOT BE INTERROGATED by the "
    "instruments the rest of FlyLab relies on (dependence profile, null "
    "ladder, ablation ladder, specification family, uncertainty budget, "
    "value of information). Do not quote it beside a cut-level result as if "
    "the two carried the same evidential weight.",
    "Signs come from predicted consensus transmitters, and five of the seven "
    "sign magnitudes (glutamate, histamine, dopamine, serotonin, octopamine) "
    "are asserted with no measurement behind them (rate.SIGN_TABLE).",
    "The engine row-normalises: each cell's recurrent input is a "
    "composition-weighted average of its presynaptic gains, so absolute "
    "connection strength is divided out. See flylab/circuit/rate.py.",
    "Rates are model output, not measured firing rates. Only the sign and the "
    "ordering of a change mean anything; the level does not.",
]


def estimate_fullcns_cost(
    min_weight: int = DEFAULT_MIN_WEIGHT,
    steps: int = 80,
    chunk: int = SPARSE_CHUNK,
    n_nodes: int = 165_122,
) -> dict[str, Any]:
    """Predicted arrays, resident memory and runtime for one whole-CNS run.

    The runtime constant is measured on this engine
    (:data:`SECONDS_PER_EDGE_STEP`), not extrapolated from the dense path --
    the dense :class:`~flylab.circuit.rate.RateNetwork` costs ~8 us per edge,
    the sparse scatter ~5.5e-3 us per edge per step.
    """
    from flylab.maps.malecns import TRACED_EDGES_BY_MIN_WEIGHT

    n_edges = TRACED_EDGES_BY_MIN_WEIGHT.get(int(min_weight), N_EDGES_V1)
    arrays_gb = estimate_edge_memory_gb(min_weight, traced_only=True)
    resident_gb = sparse_memory_gb(n_nodes, n_edges, chunk)
    return {
        "min_weight": int(min_weight),
        "n_nodes": int(n_nodes),
        "n_edges_traced": int(n_edges),
        "edge_arrays_gb": arrays_gb,
        "predicted_resident_gb": resident_gb,
        "steps": int(steps),
        "predicted_runtime_s": float(SECONDS_PER_EDGE_STEP * n_edges * int(steps)),
        "note": (
            "n_edges_traced counts rows with BOTH ends among the 165 122 traced "
            "bodies -- 25 563 197 at weight 1, against 151 856 684 rows in the "
            "file. Memory is for the arrays, not for the pandas frame the assay "
            "never builds; measured peak RSS for a weight-1 run is 1.9 GB, "
            "above predicted_resident_gb because the annotation frames and the "
            "transient int64 body-id arrays are counted too."
        ),
    }


#: Measured cost of one scatter step, per edge, at whole-CNS size (s).  On the
#: small cuts the same kernel costs ~3.4e-9 s; the gather ``r[pre]`` misses
#: cache once the rate vector no longer fits in L2, so the large-scale constant
#: is the one to extrapolate with.  Measured: 30.3 s for 80 steps over
#: 25 563 197 edges.
SECONDS_PER_EDGE_STEP = 1.48e-8


def _peak_rss_gb() -> float:
    """Peak resident set size of this process, in GB (Linux ru_maxrss is KB)."""
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1e6


def load_full_cns(
    dest: Path | None = None,
    min_weight: int = DEFAULT_MIN_WEIGHT,
    memory_ceiling_gb: float = DEFAULT_MEMORY_CEILING_GB,
    normalise: str = DEFAULT_NORMALISATION,
    chunk: int = SPARSE_CHUNK,
    seed_types: Sequence[str] = NAMED_TYPES,
) -> tuple[SparseRateNetwork, dict[str, Any]]:
    """Stream the whole traced connectome into a :class:`SparseRateNetwork`.

    Raises ``FileNotFoundError`` naming ``flylab download-malecns --full`` if
    the weight matrix is absent (it is ~1.0 GB and is never committed), and
    ``MemoryError`` with the numbers if the edge arrays would exceed
    ``memory_ceiling_gb``.
    """
    import pandas as pd

    dest = Path(dest) if dest else default_dir()
    weights_path(dest)  # explicit, early failure with the download instruction
    t0 = time.perf_counter()

    traced = traced_body_ids(dest)
    ann = pd.read_feather(
        dest / FILES["annotations"], columns=["bodyId", "type", "superclass", "status"]
    )
    ann = ann[ann["status"] == "Traced"].drop_duplicates("bodyId")
    nt = pd.read_feather(
        dest / FILES["neurotransmitters"], columns=["body", "consensus_nt"]
    ).drop_duplicates("body")
    meta = (
        pd.DataFrame({"bodyId": traced})
        .merge(ann[["bodyId", "type", "superclass"]], on="bodyId", how="left")
        .merge(nt, left_on="bodyId", right_on="body", how="left")
    )
    t_meta = time.perf_counter() - t0

    t1 = time.perf_counter()
    pre_b, post_b, w = edge_arrays(
        dest, min_weight=int(min_weight), keep=traced, memory_ceiling_gb=memory_ceiling_gb
    )
    t_stream = time.perf_counter() - t1

    pre = np.searchsorted(traced, pre_b).astype(np.int32)
    post = np.searchsorted(traced, post_b).astype(np.int32)
    del pre_b, post_b

    nt_labels = meta["consensus_nt"].fillna("unclear").astype(str).tolist()
    types = [None if pd.isna(x) else str(x) for x in meta["type"].tolist()]
    superclasses = [None if pd.isna(x) else str(x) for x in meta["superclass"].tolist()]

    seeds: dict[str, list[int]] = {}
    tser = meta["type"].fillna("").astype(str)
    for t in seed_types:
        hit = meta.loc[tser == t, "bodyId"]
        if hit.empty:
            hit = meta.loc[tser.str.startswith(t), "bodyId"]
        seeds[t] = sorted(int(x) for x in hit.tolist())

    net = SparseRateNetwork(
        n=int(traced.size),
        pre=pre,
        post=post,
        weight=w,
        nt=nt_labels,
        body_ids=traced,
        seeds=seeds,
        normalise=normalise,
        chunk=chunk,
    )
    info = {
        "n_nodes": int(traced.size),
        "n_edges": int(net.pre.size),
        "min_weight": int(min_weight),
        "normalise": normalise,
        "chunk": int(chunk),
        "types": types,
        "superclasses": superclasses,
        "seeds": {k: len(v) for k, v in seeds.items()},
        "load_metadata_s": float(t_meta),
        "load_weights_s": float(t_stream),
        "load_total_s": float(time.perf_counter() - t0),
        "edge_arrays_gb": float(net.pre.nbytes + net.post.nbytes + net.w.nbytes) / 1e9,
        "predicted_resident_gb": net.memory_gb(),
        "peak_rss_gb_after_load": _peak_rss_gb(),
        "weights_file": str(dest / FILES["weights"]),
    }
    return net, info


def _means_by(labels: Sequence[Any], rates: np.ndarray) -> dict[str, float]:
    acc: dict[str, list[float]] = {}
    for i, lab in enumerate(labels):
        key = "unknown" if lab is None or lab != lab or lab == "" else str(lab)
        acc.setdefault(key, []).append(float(rates[i]))
    return {k: float(np.mean(v)) for k, v in sorted(acc.items())}


def _named(net: SparseRateNetwork, rates: np.ndarray, typ: str) -> float | None:
    idx = net.seed_indices([typ])
    return float(np.mean(rates[idx])) if idx else None


def run_fullcns_assay(
    compound: str | None = None,
    conc_M: float = 0.0,
    dest: Path | None = None,
    min_weight: int = DEFAULT_MIN_WEIGHT,
    drive_hz: float = 40.0,
    drive: dict[int, float] | None = None,
    steps: int = 80,
    memory_ceiling_gb: float = DEFAULT_MEMORY_CEILING_GB,
    normalise: str = DEFAULT_NORMALISATION,
    chunk: int = SPARSE_CHUNK,
    library: dict[str, Any] | None = None,
    net: SparseRateNetwork | None = None,
    net_info: dict[str, Any] | None = None,
    top_k: int = 25,
) -> dict[str, Any]:
    """One treated + one vehicle run over the complete traced connectome.

    The experiment is deliberately the *same* one the cut assays run -- the
    named seed cells driven at ``drive_hz``, ``treated - vehicle`` on the same
    graph -- with the cut boundary removed, so the two are comparable.

    Returns a schema-0.3 notebook carrying the usual provenance, the measured
    runtime and peak RSS, and :data:`FULLCNS_WARNINGS` stating plainly what is
    *not* available at this scale.  Pass ``net``/``net_info`` from
    :func:`load_full_cns` to run several compounds on one load.
    """
    t_start = time.perf_counter()
    if net is None:
        net, net_info = load_full_cns(
            dest=dest,
            min_weight=min_weight,
            memory_ceiling_gb=memory_ceiling_gb,
            normalise=normalise,
            chunk=chunk,
        )
    info = dict(net_info or {})

    gains, occ = compute_gains(compound, conc_M, library=library)
    drive_map = dict(drive) if drive else net.seed_drive(drive_hz)
    d = net.drive_vector(drive_map)

    t0 = time.perf_counter()
    r_t = net.run(d, gains, steps=steps)
    t_treated = time.perf_counter() - t0
    t0 = time.perf_counter()
    r_v = net.run(d, DEFAULT_GAINS, steps=steps)
    t_vehicle = time.perf_counter() - t0

    delta = r_t - r_v
    order = np.argsort(-np.abs(delta))[: int(top_k)]
    types = info.get("types") or [None] * net.n
    superclasses = info.get("superclasses") or [None] * net.n
    top_changed = [
        {
            "bodyId": int(net.body_ids[i]),
            "type": types[i],
            "superclass": superclasses[i],
            "nt": str(net.nt[i]),
            "rate_vehicle": float(r_v[i]),
            "rate_treated": float(r_t[i]),
            "delta_hz": float(delta[i]),
        }
        for i in order.tolist()
    ]

    nb = empty_notebook("malecns_full_cns")
    set_map(nb, MAP_ID, "v1.0-full", MAP_CITATION)
    nb["compound"] = compound
    nb["concentration_M"] = conc_M
    nb["occupancy"] = occ["receptors"] if occ else []
    nb["gains"] = dict(gains)
    runtime = {
        "load_metadata_s": info.get("load_metadata_s"),
        "load_weights_s": info.get("load_weights_s"),
        "load_total_s": info.get("load_total_s"),
        "run_treated_s": float(t_treated),
        "run_vehicle_s": float(t_vehicle),
        "total_s": float(time.perf_counter() - t_start),
        "seconds_per_edge_step": (
            float(t_vehicle) / (max(int(net.pre.size), 1) * max(int(steps), 1))
        ),
    }
    memory = {
        "edge_arrays_gb": info.get("edge_arrays_gb"),
        "predicted_resident_gb": info.get("predicted_resident_gb"),
        "peak_rss_gb": _peak_rss_gb(),
        "memory_ceiling_gb": float(memory_ceiling_gb),
        "chunk_edges": int(chunk),
        "note": (
            "peak_rss_gb is the whole process (ru_maxrss), so it includes the "
            "pandas annotation frames and the transient int64 body-id arrays "
            "the streaming loader builds and frees, not just the engine."
        ),
    }
    nb["readouts"] = {
        "n_nodes": int(net.n),
        "n_edges": int(net.pre.size),
        "min_weight": int(min_weight),
        "normalise": normalise,
        "drive_hz": float(drive_hz),
        "n_driven_cells": int(len(drive_map)),
        "steps": int(steps),
        "mean_hz": float(r_t.mean()),
        "max_hz": float(r_t.max()),
        "mn9_hz": _named(net, r_t, "MN9"),
        "dnp01_hz": _named(net, r_t, "DNp01"),
        "n_active_treated": int((r_t > 0.1).sum()),
        "n_active_vehicle": int((r_v > 0.1).sum()),
        "n_at_r_max_treated": int((r_t >= 299.999).sum()),
        "n_at_r_max_vehicle": int((r_v >= 299.999).sum()),
        "vehicle": {
            "mean_hz": float(r_v.mean()),
            "max_hz": float(r_v.max()),
            "mn9_hz": _named(net, r_v, "MN9"),
            "dnp01_hz": _named(net, r_v, "DNp01"),
        },
        "delta_mean_hz": float(r_t.mean() - r_v.mean()),
        "by_superclass": _means_by(superclasses, r_t),
        "by_superclass_vehicle": _means_by(superclasses, r_v),
        "top_changed": top_changed,
        "seed_counts": info.get("seeds"),
        "runtime": runtime,
        "memory": memory,
        "receptor_expression_coverage": RECEPTOR_COVERAGE,
    }
    nb["warnings"] = list(FULLCNS_WARNINGS)
    reach_v = float((r_v > 0.1).sum()) / max(net.n, 1)
    nb["readouts"]["reach_vehicle"] = reach_v
    nb["readouts"]["reach_treated"] = float((r_t > 0.1).sum()) / max(net.n, 1)
    if reach_v < 0.05:
        nb["warnings"].append(
            f"ACTIVITY DOES NOT SPREAD: driving {len(drive_map)} seed cells at "
            f"{float(drive_hz):.0f} Hz leaves only {int((r_v > 0.1).sum())} of "
            f"{net.n} cells above 0.1 Hz in the vehicle arm ({reach_v:.2%}). The "
            "row-normalised operator is a contraction, so the drive decays "
            "within a few hops and the whole-CNS network mean is dominated by "
            "the silent majority. The network-mean readout is therefore NOT "
            "comparable with the same readout on a cut, and a 'whole-CNS "
            "effect' computed this way is an effect on a few hundred cells "
            "divided by 165 122."
        )
    if occ:
        nb["warnings"].append(occ["disclaimer"])
    nb["label"] = "model_derived"
    return nb


#: Share of the 165 122 traced cells whose MaleCNS superclass carries *any*
#: published expression statement for each FlyLab receptor, computed from
#: ``data/literature/receptor_expression_by_class.yaml`` against the v1.0
#: superclass census.  These are the numbers behind the "unmapped for about
#: four fifths of them" warning.
RECEPTOR_COVERAGE: dict[str, float] = {
    "insect_nAChR": 0.2744,
    "insect_RDL": 0.2505,
    "insect_GluCl": 0.0796,
    "insect_OctR": 0.0,
    "insect_AChE": 0.0,
}


def receptor_coverage_now(dest: Path | None = None) -> dict[str, Any]:
    """Recompute :data:`RECEPTOR_COVERAGE` from the atlas census (slow path)."""
    from flylab.maps.malecns import load_census
    from flylab.pharm.expression import RECEPTOR_GENES, _superclass_index

    census = load_census(dest)
    counts = census["superclass_counts"]
    total = sum(counts.values()) or 1
    out = {}
    for rec in RECEPTOR_GENES:
        index = _superclass_index(rec, None, 1.0, "mean", None)
        known = sum(v for k, v in counts.items() if index.get(k, {}).get("known"))
        out[rec] = known / total
    return {"n_traced": total, "fraction_class_mapped": out}


__all__ = [
    "run_fullcns_assay",
    "load_full_cns",
    "estimate_fullcns_cost",
    "receptor_coverage_now",
    "FULLCNS_WARNINGS",
    "RECEPTOR_COVERAGE",
    "DEFAULT_MEMORY_CEILING_GB",
    "DEFAULT_MIN_WEIGHT",
    "SECONDS_PER_EDGE_STEP",
]
