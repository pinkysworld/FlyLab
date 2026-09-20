"""Batch designs: compounds x concentrations x replicates -> table + CSV.

``ExperimentDesign`` is a pydantic model when pydantic is importable and an
equivalent dataclass otherwise.  The dataclass path exists so the whole science
core runs where pydantic is not available (WebAssembly / Pyodide, see
``flylab/browser/bridge.py``); both paths have the same field names, defaults,
coercion, ``ValueError`` messages and ``.model_dump()`` output.
"""
from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any, Literal, Sequence

import numpy as np

from flylab.assays.ensemble import ASSAYS, READOUT_KEYS, readouts_of, run_assay, sample_library

__all__ = [
    "ExperimentDesign",
    "run_experiment",
    "design_from_yaml",
    "rows_to_csv",
    "HAVE_PYDANTIC",
]


# --------------------------------------------------------------------------
# field rules, shared by the pydantic model and the dataclass fallback
# --------------------------------------------------------------------------
def _check_compounds(v: Any) -> list[str]:
    v = [c for c in (v or []) if c]
    if not v:
        raise ValueError("design needs at least one compound")
    return [str(c) for c in v]


def _check_concs(v: Any) -> list[float]:
    if not v:
        raise ValueError("design needs at least one concentration")
    if any(float(c) < 0 for c in v):
        raise ValueError("concentrations must be >= 0")
    return [float(c) for c in v]


def _check_replicates(v: Any) -> int:
    if int(v) < 1:
        raise ValueError("replicates must be >= 1")
    return int(v)


def _check_assay(v: Any) -> str:
    if v not in ASSAYS:
        raise ValueError(f"assay must be one of {ASSAYS}")
    return str(v)


def _default_readouts() -> list[str]:
    return list(READOUT_KEYS)


def _default_compounds() -> list[str]:
    return ["imidacloprid"]


def _default_concs() -> list[float]:
    return [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]


try:  # pragma: no cover - exercised by both paths in tests/test_browser_bridge.py
    from pydantic import BaseModel, Field, field_validator

    HAVE_PYDANTIC = True
except Exception:  # pragma: no cover - selected on Pyodide / minimal installs
    HAVE_PYDANTIC = False


if HAVE_PYDANTIC:

    class ExperimentDesign(BaseModel):  # type: ignore[no-redef]
        """Validated batch design (see ``docs/DESIGN_v0.5.md``)."""

        assay: Literal["subgraph", "spiking", "taste", "taste_map"] = "subgraph"
        compounds: list[str] = Field(default_factory=_default_compounds)
        concs_M: list[float] = Field(default_factory=_default_concs)
        replicates: int = 1
        seed: int = 0
        readouts: list[str] = Field(default_factory=_default_readouts)
        jitter_log10: float = 0.0
        drive_hz: float | None = None
        include_vehicle: bool = True
        graph: str | None = None
        keep_notebooks: bool = False
        options: dict[str, Any] = Field(default_factory=dict)

        @field_validator("compounds")
        @classmethod
        def _compounds_nonempty(cls, v):
            return _check_compounds(v)

        @field_validator("concs_M")
        @classmethod
        def _concs_ok(cls, v):
            return _check_concs(v)

        @field_validator("replicates")
        @classmethod
        def _reps_ok(cls, v):
            return _check_replicates(v)

        @field_validator("assay")
        @classmethod
        def _assay_ok(cls, v):
            return _check_assay(v)

else:  # pragma: no cover - selected only when pydantic is absent

    from dataclasses import asdict, dataclass, field

    def _as_bool(v: Any) -> bool:
        if isinstance(v, str):
            low = v.strip().lower()
            if low in ("true", "yes", "on", "1"):
                return True
            if low in ("false", "no", "off", "0", ""):
                return False
            raise ValueError(f"expected a boolean, got {v!r}")
        return bool(v)

    def _as_opt_float(v: Any) -> float | None:
        return None if v is None else float(v)

    def _as_opt_str(v: Any) -> str | None:
        return None if v is None else str(v)

    @dataclass
    class ExperimentDesign:  # type: ignore[no-redef]
        """Validated batch design (pydantic-free twin of the model above).

        Same field names, defaults and error messages; extra keyword arguments
        are ignored, as pydantic's default ``extra="ignore"`` does.
        """

        assay: str = "subgraph"
        compounds: list[str] = field(default_factory=_default_compounds)
        concs_M: list[float] = field(default_factory=_default_concs)
        replicates: int = 1
        seed: int = 0
        readouts: list[str] = field(default_factory=_default_readouts)
        jitter_log10: float = 0.0
        drive_hz: float | None = None
        include_vehicle: bool = True
        graph: str | None = None
        keep_notebooks: bool = False
        options: dict[str, Any] = field(default_factory=dict)

        #: field -> coercion/validation callable, applied in declaration order
        _RULES = {
            "assay": _check_assay,
            "compounds": _check_compounds,
            "concs_M": _check_concs,
            "replicates": _check_replicates,
            "seed": int,
            "readouts": lambda v: [str(x) for x in (v or [])],
            "jitter_log10": float,
            "drive_hz": _as_opt_float,
            "include_vehicle": _as_bool,
            "graph": _as_opt_str,
            "keep_notebooks": _as_bool,
            "options": lambda v: dict(v or {}),
        }

        def __init__(self, **kw: Any) -> None:
            defaults = {
                "assay": "subgraph",
                "compounds": _default_compounds(),
                "concs_M": _default_concs(),
                "replicates": 1,
                "seed": 0,
                "readouts": _default_readouts(),
                "jitter_log10": 0.0,
                "drive_hz": None,
                "include_vehicle": True,
                "graph": None,
                "keep_notebooks": False,
                "options": {},
            }
            for name, default in defaults.items():
                if name in kw and kw[name] is not None:
                    value = ExperimentDesign._RULES[name](kw[name])
                elif name in kw and kw[name] is None and name in ("drive_hz", "graph"):
                    value = None
                else:
                    value = default
                object.__setattr__(self, name, value)

        def model_dump(self, **_: Any) -> dict[str, Any]:
            """pydantic-compatible dict of the validated design."""
            return asdict(self)

        def dict(self, **kw: Any) -> dict[str, Any]:  # pragma: no cover - alias
            return self.model_dump(**kw)


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
