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
