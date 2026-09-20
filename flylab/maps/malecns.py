from __future__ import annotations
import os
import urllib.request
from pathlib import Path
from typing import Any

BASE = "https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome"
FILES = {
    "annotations": "body-annotations-male-cns-v1.0-minconf-0.5.feather",
    "neurotransmitters": "body-neurotransmitters-male-cns-v1.0.feather",
    "weights": "connectome-weights-male-cns-v1.0-minconf-0.5.feather",
}
MAP_ID = "male-cns:v1.0"
MAP_CITATION = "MaleCNS v1.0, HHMI Janelia FlyEM / Cambridge / MRC LMB / Google Research. CC-BY 4.0."

def default_dir() -> Path:
    env = os.environ.get("FLYLAB_MALECNS")
    if env:
        return Path(env)
    artifacts = Path("/home/workdir/artifacts/malecns_data")
    if (artifacts / FILES["annotations"]).exists():
        return artifacts
    return Path.home() / ".flylab" / "malecns_v1"

def download(kind: str = "atlas", dest: Path | None = None) -> Path:
    dest = dest or default_dir()
    dest.mkdir(parents=True, exist_ok=True)
    names = ["annotations", "neurotransmitters"] + (["weights"] if kind == "full" else [])
    for key in names:
        path = dest / FILES[key]
        if path.exists() and path.stat().st_size > 1_000_000:
            continue
        urllib.request.urlretrieve(f"{BASE}/{FILES[key]}", path)
    return dest

CENSUS_FALLBACK = Path(__file__).resolve().parents[2] / "data" / "derived" / "malecns_census_v1.json"
GUSTATORY_SEEDS = Path(__file__).resolve().parents[2] / "data" / "derived" / "malecns_gustatory_seeds.json"
NAMED_TYPES = ("MN9", "DNp01", "LB1a", "LB1b", "LB1c", "LB1d", "LB3b", "LB3c")


def load_census_fallback(path: Path | None = None) -> dict[str, Any]:
    """Committed census built once from the public atlas (no feather download needed)."""
    import json
    p = Path(path) if path else CENSUS_FALLBACK
    if not p.exists():
        raise FileNotFoundError(f"census fallback missing at {p}")
    return json.loads(p.read_text())


def load_gustatory_seeds(path: Path | None = None) -> dict[str, Any]:
    """Labellar GRN body IDs by MaleCNS type, plus the sweet/bitter working hypothesis."""
    import json
    p = Path(path) if path else GUSTATORY_SEEDS
    return json.loads(p.read_text())


def load_census(dest: Path | None = None, allow_fallback: bool = True) -> dict[str, Any]:
    dest = dest or default_dir()
    ann_path, nt_path = dest / FILES["annotations"], dest / FILES["neurotransmitters"]
    if not ann_path.exists() or not nt_path.exists():
        if allow_fallback and CENSUS_FALLBACK.exists():
            return load_census_fallback()
        raise FileNotFoundError(f"MaleCNS atlas not in {dest}. Run flylab download-malecns")
    import pandas as pd
    ann = pd.read_feather(ann_path)
    nt = pd.read_feather(nt_path)
    traced = ann[ann["status"] == "Traced"].copy()
    merged = traced.merge(nt, left_on="bodyId", right_on="body", how="left")
    named = {}
    for key in NAMED_TYPES:
        hits = merged[merged["type"].fillna("") == key]
        named[key] = [{"bodyId": int(r.bodyId), "superclass": r.superclass, "consensus_nt": r.consensus_nt} for r in hits.itertuples()]
    return {
        "map": MAP_ID, "citation": MAP_CITATION, "source": "local_atlas", "n_traced": int(len(merged)),
        "n_annotated_rows": int(len(ann)),
        "neurotransmitter_counts": merged["consensus_nt"].fillna("unclear").value_counts().astype(int).to_dict(),
        "superclass_counts": merged["superclass"].fillna("unknown").value_counts().astype(int).to_dict(),
        "class_counts": merged["class"].fillna("unknown").value_counts().astype(int).to_dict(),
        "named_cells": named, "weights_present": (dest / FILES["weights"]).exists(),
    }


# --------------------------------------------------------------------------
# streaming access to the full weight matrix
# --------------------------------------------------------------------------
#: The v1.0 flat connectome, as shipped: 151 856 684 rows of
#: ``(body_pre, body_post, weight)`` at ``minconf 0.5``, in 2 318 Arrow record
#: batches of 65 536 rows.  Read through :func:`stream_weight_batches` the file
#: is never materialised as a pandas frame: ``pd.read_feather`` on it costs
#: about 3.6 GB of int64 columns plus the frame's own copy, where the arrays
#: the rate engine actually wants (``int32`` indices, ``float32`` weights) cost
#: 1.8 GB for the whole matrix and far less for any threshold above 1.
N_EDGES_V1 = 151_856_684
N_TRACED_V1 = 165_122

#: Rows surviving each synapse-count floor, measured on v1.0 over **all**
#: bodies.  Most of them are not traced cells: two thirds of the weight-1 rows
#: have at least one end outside the 165 122 traced bodies.
EDGES_BY_MIN_WEIGHT: dict[int, int] = {
    1: 151_856_684,
    2: 57_670_765,
    3: 23_014_406,
    5: 7_622_864,
    10: 2_799_910,
}

#: The count that actually matters for a whole-CNS run: rows with **both** ends
#: among the traced bodies, measured on v1.0.  The engine's committed cuts use
#: ``min_weight=5``; a whole-CNS run defaults to 1, i.e. 25 563 197 edges.
TRACED_EDGES_BY_MIN_WEIGHT: dict[int, int] = {
    1: 25_563_197,
    2: 15_270_273,
    3: 10_511_038,
    5: 6_235_682,
    10: 2_749_407,
}

#: Bytes per edge held in memory by :func:`edge_arrays` (int32 pre, int32 post,
#: float32 weight).  Used by the memory estimators so a caller can refuse a run
#: instead of thrashing.
BYTES_PER_EDGE = 12


def weights_path(dest: Path | None = None) -> Path:
    """Path of the full weight matrix, or a clear error saying how to get it.

    The matrix is ~1.0 GB and is **never committed**; it is downloaded by
    ``flylab download-malecns --full`` (or by the extract Action).
    """
    dest = Path(dest) if dest else default_dir()
    path = dest / FILES["weights"]
    if not path.exists():
        raise FileNotFoundError(
            f"MaleCNS weight matrix not found at {path}. It is ~1.0 GB and is "
            "never committed to the repository. Run `flylab download-malecns "
            "--full` (or set FLYLAB_MALECNS to a directory that already has "
            f"{FILES['weights']})."
        )
    return path


def estimate_edge_memory_gb(min_weight: int = 1, traced_only: bool = False) -> float:
    """GB of RAM :func:`edge_arrays` needs at this synapse-count floor.

    ``traced_only`` uses :data:`TRACED_EDGES_BY_MIN_WEIGHT`, i.e. what a run
    restricted to traced cells actually holds.
    """
    table = TRACED_EDGES_BY_MIN_WEIGHT if traced_only else EDGES_BY_MIN_WEIGHT
    n = table.get(int(min_weight))
    if n is None:
        # linear-in-log interpolation is not worth it; be conservative
        n = table[max(k for k in table if k <= int(min_weight))]
    return float(n * BYTES_PER_EDGE) / 1e9


def stream_weight_batches(
    dest: Path | None = None,
    min_weight: int = 1,
    keep: Any = None,
):
    """Yield ``(pre, post, weight)`` numpy arrays, one Arrow record batch at a time.

    The file is memory-mapped and read batch by batch, so peak memory is one
    batch (65 536 rows, ~1.5 MB) plus whatever the caller keeps.  Nothing goes
    through pandas.

    ``min_weight`` drops rows below that synapse count inside the loop.
    ``keep``, if given, is a **sorted** int64 array of body ids; only rows with
    both ends in it survive (membership by ``searchsorted``, not ``isin``,
    so it is O(E log |keep|) rather than quadratic).
    """
    import numpy as np
    import pyarrow as pa

    path = weights_path(dest)
    reader = pa.ipc.open_file(pa.memory_map(str(path), "r"))
    keep_arr = None if keep is None else np.asarray(keep, dtype=np.int64)
    for i in range(reader.num_record_batches):
        b = reader.get_batch(i)
        pre = b.column("body_pre").to_numpy(zero_copy_only=False)
        post = b.column("body_post").to_numpy(zero_copy_only=False)
        w = b.column("weight").to_numpy(zero_copy_only=False)
        if min_weight > 1:
            m = w >= int(min_weight)
            pre, post, w = pre[m], post[m], w[m]
        if keep_arr is not None and pre.size:
            m = _in_sorted(pre, keep_arr) & _in_sorted(post, keep_arr)
            pre, post, w = pre[m], post[m], w[m]
        if pre.size:
            yield pre, post, w


def _in_sorted(x: Any, sorted_keys: Any) -> Any:
    """``np.isin`` for a sorted key array, in O(len(x) log len(keys))."""
    import numpy as np

    x = np.asarray(x)
    if sorted_keys.size == 0 or x.size == 0:
        return np.zeros(x.shape, dtype=bool)
    idx = np.searchsorted(sorted_keys, x)
    idx[idx >= sorted_keys.size] = 0
    return sorted_keys[idx] == x


def traced_body_ids(dest: Path | None = None) -> Any:
    """Sorted int64 array of ``status == "Traced"`` body ids (165 122 in v1.0)."""
    import numpy as np
    import pandas as pd

    dest = Path(dest) if dest else default_dir()
    ann_path = dest / FILES["annotations"]
    if not ann_path.exists():
        raise FileNotFoundError(
            f"MaleCNS annotations not found at {ann_path}. Run `flylab download-malecns`."
        )
    ann = pd.read_feather(ann_path, columns=["bodyId", "status"])
    return np.sort(ann.loc[ann["status"] == "Traced", "bodyId"].to_numpy().astype("int64"))


def edge_arrays(
    dest: Path | None = None,
    min_weight: int = 1,
    keep: Any = None,
    memory_ceiling_gb: float = 8.0,
) -> tuple[Any, Any, Any]:
    """The whole (filtered) weight matrix as three concatenated numpy arrays.

    Returns ``(body_pre, body_post, weight)`` with body ids as ``int64`` and
    weights as ``float32``.  Built by :func:`stream_weight_batches`, so the
    file is never loaded through pandas.

    ``memory_ceiling_gb`` is a refusal, not a target: if the *estimated* array
    footprint exceeds it the call raises ``MemoryError`` with the numbers,
    rather than swapping the machine.
    """
    import numpy as np

    need = estimate_edge_memory_gb(min_weight)
    if need > float(memory_ceiling_gb):
        raise MemoryError(
            f"the MaleCNS weight matrix at min_weight={min_weight} needs about "
            f"{need:.2f} GB of edge arrays (all bodies), over the "
            f"{float(memory_ceiling_gb):.2f} GB "
            "ceiling. Raise memory_ceiling_gb deliberately, or raise min_weight "
            f"(rows by floor: {EDGES_BY_MIN_WEIGHT})."
        )
    pres, posts, ws = [], [], []
    for pre, post, w in stream_weight_batches(dest, min_weight=min_weight, keep=keep):
        pres.append(pre.astype(np.int64, copy=False))
        posts.append(post.astype(np.int64, copy=False))
        ws.append(w.astype(np.float32, copy=False))
    if not pres:
        empty_i = np.empty(0, dtype=np.int64)
        return empty_i, empty_i.copy(), np.empty(0, dtype=np.float32)
    return np.concatenate(pres), np.concatenate(posts), np.concatenate(ws)
