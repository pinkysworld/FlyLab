"""Batch designs: compounds x concentrations x replicates -> table + CSV."""
from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any, Literal, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator

from flylab.assays.ensemble import ASSAYS, READOUT_KEYS, readouts_of, run_assay, sample_library

__all__ = ["ExperimentDesign", "run_experiment", "design_from_yaml", "rows_to_csv"]


class ExperimentDesign(BaseModel):
    """Validated batch design (see ``docs/DESIGN_v0.5.md``)."""

    assay: Literal["subgraph", "spiking", "taste", "taste_map"] = "subgraph"
    compounds: list[str] = Field(default_factory=lambda: ["imidacloprid"])
    concs_M: list[float] = Field(default_factory=lambda: [1e-9, 1e-8, 1e-7, 1e-6, 1e-5])
    replicates: int = 1
    seed: int = 0
    readouts: list[str] = Field(default_factory=lambda: list(READOUT_KEYS))
    jitter_log10: float = 0.0
    drive_hz: float | None = None
    include_vehicle: bool = True
    graph: str | None = None
    keep_notebooks: bool = False
    options: dict[str, Any] = Field(default_factory=dict)

    @field_validator("compounds")
    @classmethod
    def _compounds_nonempty(cls, v):
        v = [c for c in v if c]
        if not v:
            raise ValueError("design needs at least one compound")
        return v

    @field_validator("concs_M")
    @classmethod
    def _concs_ok(cls, v):
        if not v:
            raise ValueError("design needs at least one concentration")
        if any(float(c) < 0 for c in v):
            raise ValueError("concentrations must be >= 0")
        return [float(c) for c in v]

    @field_validator("replicates")
    @classmethod
    def _reps_ok(cls, v):
        if int(v) < 1:
            raise ValueError("replicates must be >= 1")
        return int(v)

    @field_validator("assay")
    @classmethod
    def _assay_ok(cls, v):
        if v not in ASSAYS:
            raise ValueError(f"assay must be one of {ASSAYS}")
        return v


def rows_to_csv(rows: Sequence[dict[str, Any]], fields: Sequence[str] | None = None) -> str:
    if not rows:
        return ""
    keys = list(fields) if fields else list(rows[0].keys())
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=keys, extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in keys})
    return buf.getvalue()


def _run_cell(design: ExperimentDesign, compound: str | None, conc_M: float, rep: int):
    """One replicate: seed offset, plus library jitter when ``jitter_log10 > 0``."""
    seed = int(design.seed) + int(rep)
    lib = None
    if compound and design.jitter_log10 > 0:
        rng = np.random.default_rng([int(design.seed), rep, int(abs(hash(compound)) % 10_000)])
        lib = sample_library(rng, float(design.jitter_log10))
    kw = dict(design.options)
    if design.graph is not None and design.assay in ("subgraph", "spiking", "taste_map"):
        kw["graph"] = design.graph
    return run_assay(
        design.assay,
        compound,
        conc_M,
        library=lib,
        drive_hz=design.drive_hz,
        seed=seed,
        **kw,
    )


def run_experiment(design: dict[str, Any] | ExperimentDesign) -> dict[str, Any]:
    """Run a batch design.

    Returns ``{design, rows, vehicle_rows, summary, notebooks, csv, warnings}``.
    ``rows`` has exactly ``len(compounds) * len(concs_M) * replicates`` entries;
    the vehicle (``compound = None``) replicate set is kept separately in
    ``vehicle_rows`` so that row count stays a clean product, and is included in
    the CSV and the summary.
    """
    d = design if isinstance(design, ExperimentDesign) else ExperimentDesign(**dict(design))
    readouts = [r for r in d.readouts] or list(READOUT_KEYS)

    rows: list[dict[str, Any]] = []
    vehicle_rows: list[dict[str, Any]] = []
    notebooks: list[dict[str, Any]] = []

    for compound in d.compounds:
        for conc in d.concs_M:
            for rep in range(d.replicates):
                nb = _run_cell(d, compound, conc, rep)
                vals = readouts_of(d.assay, nb)
                rows.append(
                    {
                        "compound": compound,
                        "conc_M": conc,
                        "replicate": rep,
                        **{k: vals.get(k) for k in readouts},
                    }
                )
                if d.keep_notebooks:
                    notebooks.append(nb)

    if d.include_vehicle:
        for rep in range(d.replicates):
            nb = _run_cell(d, None, 0.0, rep)
            vals = readouts_of(d.assay, nb)
            vehicle_rows.append(
                {
                    "compound": None,
                    "conc_M": 0.0,
                    "replicate": rep,
                    **{k: vals.get(k) for k in readouts},
                }
            )
            if d.keep_notebooks:
                notebooks.append(nb)

    summary: list[dict[str, Any]] = []
    groups: dict[tuple, list[dict[str, Any]]] = {}
    for r in rows + vehicle_rows:
        groups.setdefault((r["compound"], r["conc_M"]), []).append(r)
    for (compound, conc), grp in groups.items():
        entry: dict[str, Any] = {"compound": compound, "conc_M": conc, "n": len(grp)}
        for k in readouts:
            v = np.array(
                [g[k] for g in grp if g.get(k) is not None and np.isfinite(g[k])], dtype=float
            )
            entry[f"{k}_mean"] = float(v.mean()) if v.size else None
            entry[f"{k}_sd"] = float(v.std(ddof=1)) if v.size > 1 else 0.0
        summary.append(entry)
    summary.sort(key=lambda e: (str(e["compound"]), e["conc_M"]))

    fields = ["compound", "conc_M", "replicate"] + readouts
    return {
        "design": d.model_dump(),
        "rows": rows,
        "vehicle_rows": vehicle_rows,
        "summary": summary,
        "notebooks": notebooks,
        "csv": rows_to_csv(rows + vehicle_rows, fields),
        "n_rows": len(rows),
        "label": "model_derived",
        "warnings": [
            "Every value is simulated: teaching EC50 library, hops-limited map, "
            "gain patches. No row is a measurement from a living fly.",
        ],
    }


def design_from_yaml(path: str | Path) -> dict[str, Any]:
    """Load and validate a design YAML/JSON file; returns the design dict."""
    import yaml

    data = yaml.safe_load(Path(path).read_text()) or {}
    if "design" in data and isinstance(data["design"], dict):
        data = data["design"]
    return ExperimentDesign(**data).model_dump()
