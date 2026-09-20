"""FlyLab bench HTTP API (FastAPI).

Every route is a thin, validated wrapper around a library function: the server
owns request shapes, bounds and error mapping, never pharmacology or circuit
maths.  Endpoints whose backing module is still being written by another agent
import lazily inside the handler and answer ``501 {"detail": "module not
available yet"}`` instead of failing at import time, so the bench keeps running
while the v0.5 analysis layer lands.

Error contract: ``KeyError`` / ``FileNotFoundError`` -> 404, ``ValueError``
(including pydantic validation raised inside a handler) -> 400, both carrying
the library's own message.
"""

from __future__ import annotations

import csv
import io
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Literal

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from flylab.analysis.impact import grn_to_mn9_paths, summarize_impact
from flylab.browser import bridge
from flylab.analysis.layout import graph_for_viewer
from flylab.assays.ensemble import circuit_ic50, run_ensemble, sensitivity
from flylab.assays.experiment import ExperimentDesign, rows_to_csv, run_experiment
from flylab.assays.spiking import run_spiking_assay
from flylab.assays.subgraph import dose_response_subgraph, run_subgraph_assay
from flylab.assays.taste import dose_response, library_keys, run_taste_assay
from flylab.assays.taste_map import run_taste_map_assay
from flylab.assays.wholens import run_wholens_assay
from flylab.circuit.rate import DEFAULT_GAINS, GRAPHS, compute_gains, load_graph, rate_network
from flylab.maps.malecns import MAP_CITATION, MAP_ID, load_census, load_gustatory_seeds
from flylab.pharm.exposure import ROUTE_DEFAULTS, exposure_profile
from flylab.pharm.mechanisms import GAIN_KEYS, mechanism_table_rows
from flylab.pharm.occupancy import (
    compare_compound,
    library_sha256,
    list_compounds,
    load_library,
    occupancy_curve,
    receptor_table,
    selectivity_pairs,
)
from flylab.pharm.uncertainty import occupancy_ci

VERSION = "0.5.0"
STATIC = Path(__file__).parent / "static"

#: bounds shared with the pydantic models and documented in /api/meta
MAX_N_REP = 64
MAX_T_MS = 3000.0
MAX_EXPERIMENT_ROWS = 2000
MAX_CONC_M = 1.0

NOT_READY = "module not available yet"

app = FastAPI(
    title="FlyLab",
    version=VERSION,
    description="Virtual Drosophila pharmacology bench on MaleCNS v1.0.",
)

# The bench is a local tool: a researcher may open the static page from a
# different localhost port (or a notebook) than the API. Allow that, nothing else.
app.add_middleware(
    CORSMiddleware,
    allow_origin_regex=r"^https?://(localhost|127\.0\.0\.1|\[::1\])(:\d+)?$",
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.mount("/static", StaticFiles(directory=STATIC), name="static")


# --------------------------------------------------------------------------
# error mapping
# --------------------------------------------------------------------------
def _message(exc: Exception) -> str:
    if isinstance(exc, KeyError) and exc.args:
        return str(exc.args[0])
    return str(exc) or exc.__class__.__name__


@app.exception_handler(KeyError)
def _key_error(_: Request, exc: KeyError):
    return JSONResponse(status_code=404, content={"detail": _message(exc)})


@app.exception_handler(FileNotFoundError)
def _missing_file(_: Request, exc: FileNotFoundError):
    return JSONResponse(status_code=404, content={"detail": _message(exc)})


@app.exception_handler(ValueError)
def _value_error(_: Request, exc: ValueError):
    return JSONResponse(status_code=400, content={"detail": _message(exc)})


def _lazy(dotted: str, name: str) -> Callable[..., Any]:
    """Import ``name`` from ``dotted`` or raise 501.

    Used by every forward-compatible endpoint so a module another agent has not
    landed yet degrades to a documented 501 instead of breaking the bench.
    """
    try:
        module = __import__(dotted, fromlist=[name])
        return getattr(module, name)
    except (ImportError, AttributeError):
        raise HTTPException(status_code=501, detail=NOT_READY)


def _call(fn: Callable[..., Any], /, **kw: Any) -> Any:
    """Call a forward-compatible library function with keyword arguments only.

    The analysis modules land separately, so their parameter *order* is not
    something this server may assume. A keyword this build does not know yet is
    a contract mismatch, not a user error, so it answers 501 rather than 500.
    """
    try:
        return fn(**kw)
    except TypeError as exc:
        if "unexpected keyword" in str(exc) or "required positional" in str(exc):
            raise HTTPException(
                status_code=501, detail=f"module signature not compatible yet: {exc}"
            )
        raise


def _with_genotype(fn: Callable[..., Any], genotype: str | None, /, *args: Any, **kw: Any) -> Any:
    """Call ``fn`` forwarding ``genotype=`` only when one was requested.

    ``flylab.pharm.genotype`` is still being written; until the assays grow the
    keyword, asking for a genotype answers 501 rather than silently ignoring it.
    """
    if genotype is None:
        return fn(*args, **kw)
    try:
        return fn(*args, genotype=genotype, **kw)
    except TypeError as exc:
        if "genotype" in str(exc):
            raise HTTPException(status_code=501, detail="genotype not supported yet")
        raise


def _timed(fn: Callable[[], Any]) -> tuple[Any, float]:
    t0 = time.perf_counter()
    out = fn()
    return out, (time.perf_counter() - t0) * 1000.0


def _stamp(nb: Any, runtime_ms: float) -> Any:
    """Attach the server-side runtime to a notebook without touching its schema."""
    if isinstance(nb, dict):
        nb.setdefault("readouts", {})
        if isinstance(nb["readouts"], dict):
            nb["readouts"]["runtime_ms"] = round(runtime_ms, 2)
    return nb


# --------------------------------------------------------------------------
# request models
# --------------------------------------------------------------------------
GraphName = Literal["named", "taste_motor"]
AssayName = Literal["subgraph", "spiking", "taste", "taste_map"]


class AssayRequest(BaseModel):
    """Shared shape of the reduced-taste / whole-CNS / subgraph assays."""

    compound: str | None = None
    conc_M: float = Field(0.0, ge=0, le=MAX_CONC_M)
    sugar_hz: float = Field(150.0, ge=0, le=1000)
    bitter_hz: float = Field(0.0, ge=0, le=1000)
    drive_hz: float = Field(40.0, ge=0, le=1000)
    steps: int = Field(80, ge=1, le=1000)
    graph: GraphName | None = None
    genotype: str | None = None


class SpikingRequest(BaseModel):
    compound: str | None = None
    conc_M: float = Field(0.0, ge=0, le=MAX_CONC_M)
    drive_hz: float = Field(40.0, ge=0, le=1000)
    t_ms: float = Field(500.0, gt=0, le=MAX_T_MS)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    graph: GraphName | None = None
    genotype: str | None = None


class TasteMapRequest(BaseModel):
    compound: str | None = None
    conc_M: float = Field(0.0, ge=0, le=MAX_CONC_M)
    sugar_hz: float = Field(150.0, ge=0, le=1000)
    bitter_hz: float = Field(0.0, ge=0, le=1000)
    engine: Literal["rate", "lif"] = "rate"
    seed: int = Field(0, ge=0, le=2**31 - 1)
    graph: GraphName | None = "taste_motor"
    genotype: str | None = None


class EnsembleRequest(BaseModel):
    assay: AssayName = "subgraph"
    compound: str | None = None
    conc_M: float = Field(0.0, ge=0, le=MAX_CONC_M)
    n_rep: int = Field(8, ge=2, le=MAX_N_REP)
    ec50_sd_log10: float = Field(0.3, ge=0, le=3)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    graph: GraphName | None = None
    genotype: str | None = None


class SensitivityRequest(BaseModel):
    assay: AssayName = "subgraph"
    compound: str = "imidacloprid"
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    readout: str = "mean_hz"
    factor: float = Field(2.0, gt=1, le=100)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    graph: GraphName | None = None
    genotype: str | None = None


class Ic50Request(BaseModel):
    assay: AssayName = "subgraph"
    compound: str = "imidacloprid"
    readout: str = "mean_hz"
    concs_M: list[float] | None = None
    n_boot: int = Field(200, ge=0, le=2000)
    n_rep: int = Field(4, ge=1, le=MAX_N_REP)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    graph: GraphName | None = None
    genotype: str | None = None


class ImpactRequest(BaseModel):
    compound: str | None = None
    conc_M: float = Field(0.0, ge=0, le=MAX_CONC_M)
    graph: GraphName | None = "named"
    drive_hz: float = Field(40.0, ge=0, le=1000)
    top_edges: int = Field(25, ge=1, le=200)
    top_nodes: int = Field(25, ge=1, le=200)
    top_paths: int = Field(10, ge=1, le=100)
    max_len: int = Field(2, ge=1, le=4)


class ExposureRequest(BaseModel):
    compound: str = "imidacloprid"
    dose: float = Field(1.0, ge=0, le=1e6)
    route: Literal["feeding", "topical", "bath"] = "feeding"
    t_h: float = Field(24.0, gt=0, le=240)
    dt_h: float = Field(0.05, gt=0, le=24)


class ExperimentRequest(BaseModel):
    """Design JSON, validated again by :class:`ExperimentDesign` in the handler."""

    assay: AssayName = "subgraph"
    compounds: list[str] = Field(default_factory=lambda: ["imidacloprid"])
    concs_M: list[float] = Field(default_factory=lambda: [1e-9, 1e-8, 1e-7, 1e-6, 1e-5])
    replicates: int = Field(1, ge=1, le=MAX_N_REP)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    readouts: list[str] | None = None
    jitter_log10: float = Field(0.0, ge=0, le=3)
    drive_hz: float | None = Field(None, ge=0, le=1000)
    include_vehicle: bool = True
    graph: GraphName | None = None
    keep_notebooks: bool = False
    genotype: str | None = None


class LiveLabRequest(BaseModel):
    notebook: dict[str, Any]
    csv_text: str
    assay_name: str = "live_lab"
    source: str = "imported CSV"


class NullRequest(BaseModel):
    assay: AssayName = "subgraph"
    compound: str = "imidacloprid"
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    mode: str = "sign_permute"
    n: int = Field(20, ge=1, le=500)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    graph: GraphName | None = None


class NullPanelRequest(BaseModel):
    compound: str = "imidacloprid"
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    n: int = Field(8, ge=1, le=500)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    #: optional narrowing, because a full panel is minutes of compute
    modes: list[str] | None = None
    include_taste_map: bool | None = None


class LandscapeRequest(BaseModel):
    assay: AssayName = "subgraph"
    graphs: list[str] | None = None
    compounds: list[str] | None = None


class CircuitSiRequest(BaseModel):
    compound: str = "imidacloprid"
    assay: AssayName = "subgraph"
    graph: GraphName | None = "named"


class MixtureComponent(BaseModel):
    compound: str
    conc_M: float = Field(0.0, ge=0, le=MAX_CONC_M)


class MixtureRequest(BaseModel):
    components: list[MixtureComponent] = Field(default_factory=list)
    assay: AssayName = "subgraph"
    graph: GraphName | None = None
    model: Literal["bliss", "loewe"] = "bliss"


# --------------------------------------------------------------------------
# page
# --------------------------------------------------------------------------
@app.get("/", include_in_schema=False)
def index():
    return FileResponse(STATIC / "index.html")


@app.get("/api/health")
def health():
    return {"ok": True, "version": VERSION}


# --------------------------------------------------------------------------
# meta / library
# --------------------------------------------------------------------------
def _graph_meta() -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key in GRAPHS:
        try:
            g = load_graph(key)
        except FileNotFoundError as exc:
            out[key] = {"available": False, "error": str(exc)}
            continue
        out[key] = {
            "available": True,
            "file": GRAPHS[key],
            "map": g.get("map"),
            "citation": g.get("citation"),
            "n_nodes": g.get("n_nodes"),
            "n_edges": g.get("n_edges"),
            "min_weight": g.get("min_weight"),
            "hops": g.get("hops"),
            "seeds": {t: len(ids) for t, ids in (g.get("seeds") or {}).items()},
            "n_seed_cells": sum(len(ids) for ids in (g.get("seeds") or {}).values()),
        }
    return out


def _compound_meta(lib: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for key in list_compounds(lib):
        spec = lib["compounds"][key]
        targets = [
            {
                "receptor": receptor,
                # schema v3: param_type says what the number is; ec50_M is the
                # deprecated alias and is None on a not-modelled row.
                "param_type": r.get("param_type", "EC50"),
                "param_value_M": float(r["value_M"]),
                "ec50_M": float(r["value_M"]),
                "relation": r.get("relation"),
                "species": r.get("species"),
                "direction": r.get("direction", "none"),
                "evidence_tier": r.get("evidence_tier", "class_placeholder"),
                "source": r.get("source", ""),
            }
            for receptor, r in (spec.get("receptors") or {}).items()
            if receptor.startswith("insect_")
            and r.get("direction") not in (None, "none")
            and r.get("value_M", r.get("ec50_M")) is not None
        ]
        targets.sort(key=lambda r: r["param_value_M"])
        rows.append(
            {
                "key": key,
                "name": spec.get("name", key),
                "class": spec.get("class"),
                "cas": spec.get("cas"),
                "insect_targets": targets,
                "n_receptors": len(spec.get("receptors") or {}),
                "evidence_tier": (
                    "literature_order"
                    if any(t["evidence_tier"] == "literature_order" for t in targets)
                    else "class_placeholder"
                ),
            }
        )
    return rows


@app.get("/api/meta")
def meta():
    """Everything the UI needs to draw its controls before the first run."""
    lib = load_library()
    try:
        seeds = load_gustatory_seeds()
    except (FileNotFoundError, OSError):
        seeds = {"note": "gustatory seed file not present"}
    return {
        "version": VERSION,
        "notebook_version": "0.3",
        "map": {"id": MAP_ID, "citation": MAP_CITATION},
        "graphs": _graph_meta(),
        "library": {
            "version": lib.get("library_version"),
            "schema_version": lib.get("schema_version"),
            "sha256": library_sha256(),
            "n_compounds": len(lib.get("compounds") or {}),
        },
        "compounds": _compound_meta(lib),
        "receptors": receptor_table(lib),
        "selectivity_pairs": selectivity_pairs(lib),
        "mechanisms": mechanism_table_rows(),
        "gain_keys": list(GAIN_KEYS),
        "assays": ["subgraph", "spiking", "taste", "taste_map"],
        "engines": ["rate", "lif"],
        "routes": ROUTE_DEFAULTS,
        "gustatory_hypothesis": seeds,
        "limits": {
            "max_n_rep": MAX_N_REP,
            "max_t_ms": MAX_T_MS,
            "max_experiment_rows": MAX_EXPERIMENT_ROWS,
            "max_conc_M": MAX_CONC_M,
        },
    }


@app.get("/api/census")
def census():
    """MaleCNS traced-cell census (local atlas when present, else the committed fallback)."""
    return load_census()


@app.get("/api/drugs")
def drugs():
    return {"compounds": library_keys()}


@app.get("/api/drugs/{key}")
def drug(key: str):
    lib = load_library()
    k = key.lower().strip()
    if k not in lib["compounds"]:
        raise HTTPException(status_code=404, detail=f"unknown compound {key!r}")
    return {"key": k, **lib["compounds"][k]}


@app.get("/api/occupancy")
def occupancy(compound: str, conc_M: float = 0.0):
    if conc_M < 0:
        raise HTTPException(status_code=400, detail="conc_M must be >= 0")
    return compare_compound(compound, conc_M)


@app.get("/api/occupancy/curve")
def occupancy_curve_route(compound: str):
    return {
        "compound": compound,
        "points": occupancy_curve(compound),
        "receptors": receptor_table(),
    }


@app.get("/api/occupancy/ci")
def occupancy_ci_route(
    compound: str,
    conc_M: float = 0.0,
    sd_log10: float = 0.3,
    n: int = 500,
    seed: int = 0,
):
    if conc_M < 0:
        raise HTTPException(status_code=400, detail="conc_M must be >= 0")
    if not 2 <= n <= 5000:
        raise HTTPException(status_code=400, detail="n must be between 2 and 5000")
    return occupancy_ci(compound, conc_M, sd_log10=sd_log10, n=n, seed=seed)


# --------------------------------------------------------------------------
# assays
# --------------------------------------------------------------------------
@app.post("/api/assay/taste")
def assay_taste(req: AssayRequest):
    nb, ms = _timed(
        lambda: _with_genotype(
            run_taste_assay, req.genotype, req.compound, req.conc_M, req.sugar_hz, req.bitter_hz
        )
    )
    return _stamp(nb, ms)


@app.get("/api/assay/taste/dose-response")
def assay_taste_dose(compound: str):
    return {"compound": compound, "points": dose_response(compound, [10**e for e in range(-10, -4)])}


@app.post("/api/assay/cns")
def assay_cns(req: AssayRequest):
    nb, ms = _timed(
        lambda: _with_genotype(run_wholens_assay, req.genotype, req.compound, req.conc_M)
    )
    return _stamp(nb, ms)


@app.post("/api/assay/subgraph")
def assay_subgraph(req: AssayRequest):
    nb, ms = _timed(
        lambda: _with_genotype(
            run_subgraph_assay,
            req.genotype,
            req.compound,
            req.conc_M,
            drive_hz=req.drive_hz,
            steps=req.steps,
            graph=req.graph,
        )
    )
    return _stamp(nb, ms)


@app.get("/api/assay/subgraph/dose-response")
def assay_subgraph_dose(compound: str):
    return {"compound": compound, "points": dose_response_subgraph(compound)}


@app.post("/api/assay/spiking")
def assay_spiking(req: SpikingRequest):
    nb, ms = _timed(
        lambda: _with_genotype(
            run_spiking_assay,
            req.genotype,
            req.compound,
            req.conc_M,
            drive_hz=req.drive_hz,
            t_ms=req.t_ms,
            seed=req.seed,
            graph=req.graph,
        )
    )
    return _stamp(nb, ms)


@app.post("/api/assay/taste-map")
def assay_taste_map(req: TasteMapRequest):
    nb, ms = _timed(
        lambda: _with_genotype(
            run_taste_map_assay,
            req.genotype,
            req.compound,
            req.conc_M,
            sugar_hz=req.sugar_hz,
            bitter_hz=req.bitter_hz,
            engine=req.engine,
            seed=req.seed,
            graph=req.graph or "taste_motor",
        )
    )
    return _stamp(nb, ms)


@app.post("/api/assay/ensemble")
def assay_ensemble(req: EnsembleRequest):
    kw: dict[str, Any] = {}
    if req.graph and req.assay in ("subgraph", "spiking", "taste_map"):
        kw["graph"] = req.graph
    nb, ms = _timed(
        lambda: _with_genotype(
            run_ensemble,
            req.genotype,
            req.assay,
            req.compound,
            req.conc_M,
            n_rep=req.n_rep,
            ec50_sd_log10=req.ec50_sd_log10,
            seed=req.seed,
            **kw,
        )
    )
    return _stamp(nb, ms)


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------
@app.post("/api/analysis/sensitivity")
def analysis_sensitivity(req: SensitivityRequest):
    kw: dict[str, Any] = {}
    if req.graph and req.assay in ("subgraph", "spiking", "taste_map"):
        kw["graph"] = req.graph
    rows, ms = _timed(
        lambda: _with_genotype(
            sensitivity,
            req.genotype,
            req.assay,
            req.compound,
            req.conc_M,
            readout=req.readout,
            factor=req.factor,
            seed=req.seed,
            **kw,
        )
    )
    return {
        "assay": req.assay,
        "compound": req.compound,
        "concentration_M": req.conc_M,
        "readout": req.readout,
        "factor": req.factor,
        "rows": rows,
        "runtime_ms": round(ms, 2),
        "label": "model_derived",
        "warnings": [
            "One-at-a-time tornado over teaching parameters; it reports model "
            "sensitivity, not biological variability.",
        ],
    }


@app.post("/api/analysis/ic50")
def analysis_ic50(req: Ic50Request):
    if req.concs_M is not None:
        if not req.concs_M:
            raise HTTPException(status_code=400, detail="concs_M must not be empty")
        if len(req.concs_M) > 64:
            raise HTTPException(status_code=400, detail="concs_M is limited to 64 points")
        if any(c < 0 for c in req.concs_M):
            raise HTTPException(status_code=400, detail="concentrations must be >= 0")
    kw: dict[str, Any] = {}
    if req.graph and req.assay in ("subgraph", "spiking", "taste_map"):
        kw["graph"] = req.graph
    out, ms = _timed(
        lambda: _with_genotype(
            circuit_ic50,
            req.genotype,
            req.assay,
            req.compound,
            readout=req.readout,
            concs=req.concs_M,
            n_boot=req.n_boot,
            n_rep=req.n_rep,
            seed=req.seed,
            **kw,
        )
    )
    out["runtime_ms"] = round(ms, 2)
    return out


# --------------------------------------------------------------------------
# graph viewer
# --------------------------------------------------------------------------
@app.get("/api/graph")
def graph(
    graph: GraphName = "named",
    min_weight: float = 5.0,
    max_edges: int = 2000,
    compound: str | None = None,
    conc_M: float = 0.0,
    drive_hz: float = 40.0,
    steps: int = 80,
):
    """Viewer payload: positions, transmitter, vehicle and treated rates.

    When ``compound`` is given the nodes carry both rates and the edges carry
    both effective weights, so the UI can toggle vehicle/treated with no second
    request.
    """
    if min_weight < 0:
        raise HTTPException(status_code=400, detail="min_weight must be >= 0")
    if not 1 <= max_edges <= 20000:
        raise HTTPException(status_code=400, detail="max_edges must be between 1 and 20000")
    if conc_M < 0:
        raise HTTPException(status_code=400, detail="conc_M must be >= 0")

    t0 = time.perf_counter()
    g = load_graph(graph)
    net = rate_network(graph)
    drive = net.drive_vector(net.seed_drive(drive_hz))
    r_veh = net.run(drive, DEFAULT_GAINS, steps=steps)
    gains, _ = compute_gains(compound, conc_M) if compound else (dict(DEFAULT_GAINS), None)
    r_tre = net.run(drive, gains, steps=steps) if compound else r_veh

    payload = graph_for_viewer(
        g, rates=r_tre, gains=gains, min_weight=min_weight, max_edges=max_edges
    )
    vehicle = graph_for_viewer(
        g, rates=r_veh, gains=None, min_weight=min_weight, max_edges=max_edges
    )
    veh_rate = {n["id"]: n["rate"] for n in vehicle["nodes"]}
    for node in payload["nodes"]:
        node["rate_vehicle"] = veh_rate.get(node["id"])
    veh_eff = {(e["source"], e["target"]): e["eff_weight"] for e in vehicle["edges"]}
    for edge in payload["edges"]:
        base = veh_eff.get((edge["source"], edge["target"]), edge["eff_weight"])
        edge["eff_weight_vehicle"] = base
        edge["eff_weight_delta"] = edge["eff_weight"] - base

    payload["graph"] = graph
    payload["compound"] = compound
    payload["concentration_M"] = conc_M
    payload["drive_hz"] = drive_hz
    payload["runtime_ms"] = round((time.perf_counter() - t0) * 1000.0, 2)
    payload["warnings"] = [
        "Hops-limited MaleCNS neighborhood with a weight floor, not the full CNS.",
        "Node rates come from the rate model; they are not measured spike rates.",
    ]
    if payload.get("truncated_edges"):
        payload["warnings"].append(
            f"Only the {max_edges} heaviest edges above {min_weight} synapses are drawn."
        )
    return payload


@app.post("/api/graph/impact")
def graph_impact(req: ImpactRequest):
    name = req.graph or "named"
    out, ms = _timed(
        lambda: summarize_impact(
            req.compound,
            req.conc_M,
            graph=name,
            drive_hz=req.drive_hz,
            top_edges=req.top_edges,
            top_nodes=req.top_nodes,
            top_paths=req.top_paths,
            max_len=req.max_len,
        )
    )
    if name == "taste_motor":
        g = load_graph(name)
        gains, _ = compute_gains(req.compound, req.conc_M)
        out["grn_to_mn9_paths"] = grn_to_mn9_paths(
            g, max_len=3, top=req.top_paths, gains=gains
        )
    out["runtime_ms"] = round(ms, 2)
    return out


# --------------------------------------------------------------------------
# exposure
# --------------------------------------------------------------------------
@app.post("/api/exposure")
def exposure(req: ExposureRequest):
    if req.dt_h > req.t_h:
        raise HTTPException(status_code=400, detail="dt_h must be <= t_h")
    out, ms = _timed(
        lambda: exposure_profile(req.compound, req.dose, req.route, t_h=req.t_h, dt_h=req.dt_h)
    )
    out["runtime_ms"] = round(ms, 2)
    return out


# --------------------------------------------------------------------------
# experiment
# --------------------------------------------------------------------------
def _design(req: ExperimentRequest) -> ExperimentDesign:
    n_rows = len(req.compounds) * len(req.concs_M) * req.replicates
    if n_rows > MAX_EXPERIMENT_ROWS:
        raise HTTPException(
            status_code=400,
            detail=(
                f"design would produce {n_rows} rows; the limit is "
                f"{MAX_EXPERIMENT_ROWS}. Reduce compounds, concentrations or replicates."
            ),
        )
    if any(c < 0 for c in req.concs_M):
        raise HTTPException(status_code=400, detail="concentrations must be >= 0")
    data = req.model_dump(exclude_none=True, exclude={"genotype", "readouts"})
    if req.readouts:
        data["readouts"] = req.readouts
    return ExperimentDesign(**data)


def _run_design(req: ExperimentRequest) -> dict[str, Any]:
    design = _design(req)
    out, ms = _timed(lambda: _with_genotype(run_experiment, req.genotype, design))
    out["runtime_ms"] = round(ms, 2)
    return out


@app.post("/api/experiment")
def experiment(req: ExperimentRequest):
    out = _run_design(req)
    if not req.keep_notebooks:
        out.pop("notebooks", None)
    return out


@app.post("/api/experiment/csv")
def experiment_csv(req: ExperimentRequest):
    out = _run_design(req)
    text = out.get("csv") or rows_to_csv(out.get("rows", []) + out.get("vehicle_rows", []))
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return PlainTextResponse(
        text,
        media_type="text/csv; charset=utf-8",
        headers={"Content-Disposition": f'attachment; filename="flylab_experiment_{stamp}.csv"'},
    )


# --------------------------------------------------------------------------
# notebook: live lab import
# --------------------------------------------------------------------------
def _coerce(value: str) -> Any:
    v = (value or "").strip()
    if v == "":
        return None
    try:
        return int(v)
    except ValueError:
        pass
    try:
        f = float(v)
    except ValueError:
        return v
    return f


@app.post("/api/notebook/live-lab")
def notebook_live_lab(req: LiveLabRequest):
    """Attach a **real** measured table to a notebook's ``live_lab`` slot.

    This endpoint only parses and stores what the user pasted. It never
    simulates, fills or extrapolates a value: ``live_lab`` stays ``null``
    unless a human supplies a table.
    """
    text = req.csv_text or ""
    if not text.strip():
        raise HTTPException(status_code=400, detail="csv_text is empty")
    try:
        dialect = csv.Sniffer().sniff(text[:4096], delimiters=",;\t|")
    except csv.Error:
        dialect = csv.excel
    reader = csv.DictReader(io.StringIO(text), dialect=dialect)
    if not reader.fieldnames:
        raise HTTPException(status_code=400, detail="csv_text has no header row")
    columns = [c.strip() for c in reader.fieldnames if c is not None]
    rows: list[dict[str, Any]] = []
    for raw in reader:
        rows.append(
            {
                (k.strip() if k else ""): _coerce(v if isinstance(v, str) else "")
                for k, v in raw.items()
                if k is not None
            }
        )
    if not rows:
        raise HTTPException(status_code=400, detail="csv_text has a header but no data rows")

    notebook = dict(req.notebook)
    notebook["live_lab"] = {
        "assay_name": req.assay_name,
        "source": req.source,
        "imported_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "columns": columns,
        "n_rows": len(rows),
        "table": rows,
    }
    warnings = list(notebook.get("warnings") or [])
    note = (
        "live_lab holds an imported table supplied by the user. FlyLab never "
        "generates live-animal values; provenance of these rows is the user's."
    )
    if note not in warnings:
        warnings.append(note)
    notebook["warnings"] = warnings
    return notebook


# --------------------------------------------------------------------------
# forward-compatible endpoints (v0.5 analysis layer, landing separately)
# --------------------------------------------------------------------------
@app.post("/api/analysis/null")
def analysis_null(req: NullRequest):
    fn = _lazy("flylab.analysis.nullmodels", "null_distribution")
    kw: dict[str, Any] = {}
    if req.graph:
        kw["graph"] = req.graph
    return _call(
        fn,
        assay=req.assay,
        compound=req.compound,
        conc_M=req.conc_M,
        mode=req.mode,
        n=req.n,
        seed=req.seed,
        **kw,
    )


@app.post("/api/analysis/null-panel")
def analysis_null_panel(req: NullPanelRequest):
    fn = _lazy("flylab.analysis.nullmodels", "null_panel")
    kw: dict[str, Any] = {}
    if req.modes:
        kw["modes"] = req.modes
    if req.include_taste_map is not None:
        kw["include_taste_map"] = req.include_taste_map
    return _call(fn, compound=req.compound, conc_M=req.conc_M, n=req.n, seed=req.seed, **kw)


@app.get("/api/analysis/selectivity")
def analysis_selectivity(conc_M: float = 1e-6):
    if conc_M < 0:
        raise HTTPException(status_code=400, detail="conc_M must be >= 0")
    fn = _lazy("flylab.analysis.selectivity", "receptor_selectivity_table")
    return _call(fn, conc_M=conc_M)


@app.post("/api/analysis/selectivity-landscape")
def analysis_selectivity_landscape(req: LandscapeRequest):
    fn = _lazy("flylab.analysis.selectivity", "selectivity_landscape")
    kw: dict[str, Any] = {"assay": req.assay}
    if req.graphs:
        kw["graphs"] = req.graphs
    if req.compounds:
        kw["compounds"] = req.compounds
    return _call(fn, **kw)


@app.post("/api/analysis/circuit-si")
def analysis_circuit_si(req: CircuitSiRequest):
    fn = _lazy("flylab.analysis.selectivity", "circuit_selectivity_index")
    return _call(fn, compound=req.compound, assay=req.assay, graph=req.graph)


@app.get("/api/predictions")
def predictions(n_rep: int = 4, seed: int = 0):
    if not 1 <= n_rep <= MAX_N_REP:
        raise HTTPException(status_code=400, detail=f"n_rep must be 1..{MAX_N_REP}")
    fn = _lazy("flylab.analysis.predictions", "prediction_table")
    return _call(fn, n_rep=n_rep, seed=seed)


@app.get("/api/genotypes")
def genotypes():
    fn = _lazy("flylab.pharm.genotype", "list_genotypes")
    return _call(fn)


@app.post("/api/mixture")
def mixture(req: MixtureRequest):
    if not req.components:
        raise HTTPException(status_code=400, detail="mixture needs at least one component")
    fn = _lazy("flylab.pharm.mixtures", "mixture_assay")
    kw: dict[str, Any] = {"assay": req.assay, "model": req.model}
    if req.graph:
        kw["graph"] = req.graph
    return _call(fn, components=[c.model_dump() for c in req.components], **kw)


@app.get("/api/expression")
def expression():
    fn = _lazy("flylab.pharm.expression", "expression_table")
    return _call(fn)


@app.get("/api/validation")
def validation():
    fn = _lazy("flylab.validation.rank", "validate_all")
    return _call(fn)


# --------------------------------------------------------------------------
# dashboard / analysis (block 1-10 of the decision layer)
# --------------------------------------------------------------------------
# These routes validate their request here and then run the *same* function the
# browser build runs (``flylab.browser.bridge``), so a dashboard rendered from
# the served bench and one rendered from the static build are the same object.
# Every route that costs more than a second states an estimated runtime, and
# answers ``estimate_only=true`` with the estimate alone.
def _bridged(route: str, payload: dict[str, Any]) -> Any:
    try:
        return bridge.handle(route, payload)
    except bridge.BridgeError as exc:
        raise HTTPException(status_code=exc.status, detail=exc.detail)


#: the dashboard ladder runs on the rate engine, which covers these two assays
LadderAssay = Literal["subgraph", "taste_map"]


class DashboardRequest(BaseModel):
    compound: str = "imidacloprid"
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    graph: GraphName | None = "named"
    assay: LadderAssay = "subgraph"
    genotype: str | None = None
    #: dependence shuffles; the browser build and this default are FAST_N
    n: int = Field(bridge.FAST_DEPENDENCE_N, ge=1, le=500)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    include_dependence: bool = True
    include_ladder: bool = True
    estimate_only: bool = False


class CompareRequest(BaseModel):
    compounds: list[str] = Field(
        default_factory=lambda: ["imidacloprid", "fipronil", "deltamethrin"]
    )
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    assay: LadderAssay = "subgraph"
    graph: GraphName | None = "named"
    n: int = Field(bridge.FAST_DEPENDENCE_N, ge=1, le=500)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    include_dependence: bool = True
    estimate_only: bool = False


class ClaimsRequest(BaseModel):
    compound: str | None = "imidacloprid"
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    assay: AssayName = "subgraph"
    graph: GraphName | None = "named"
    #: a notebook the caller already has; nothing is re-run when it is supplied
    notebook: dict[str, Any] | None = None
    run_assay: bool = True


class DependenceRequest(BaseModel):
    compound: str = "imidacloprid"
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    assay: AssayName = "subgraph"
    readout: str = "auto"
    graph: GraphName | None = None
    n: int = Field(bridge.FAST_DEPENDENCE_N, ge=1, le=1000)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    modes: list[str] | None = None
    estimate_only: bool = False


class DependenceLandscapeRequest(BaseModel):
    compounds: list[str] | None = None
    concs_M: list[float] | None = None
    assay: AssayName = "subgraph"
    graph: GraphName | None = None
    n: int = Field(bridge.FAST_DEPENDENCE_N, ge=1, le=1000)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    #: a landscape is minutes of compute, so the default answer is the estimate
    estimate_only: bool = True


class AblationRequest(BaseModel):
    compound: str = "imidacloprid"
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    assay: AssayName = "subgraph"
    graph: GraphName | None = "named"
    compounds: list[str] | None = None
    concs_M: list[float] | None = None
    table: bool = False
    include_information_gain: bool = True
    estimate_only: bool = False


class StabilityRequest(BaseModel):
    fast: bool = True
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    estimate_only: bool = False


class ThresholdRequest(BaseModel):
    circuit_fracs: list[float] | None = None
    vert_limits: list[float] | None = None
    compounds: list[str] | None = None
    assay: AssayName = "subgraph"
    graph: GraphName | None = None
    estimate_only: bool = False


class GlobalUncertaintyRequest(BaseModel):
    compound: str = "imidacloprid"
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    readout: str = "mean_hz"
    graph: GraphName | None = "named"
    #: the browser build defaults to the cheap design; the paper uses 128+
    n_base: int = Field(32, ge=8, le=1024)
    n_boot: int = Field(50, ge=0, le=2000)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    estimate_only: bool = False


class VoiRequest(BaseModel):
    compound: str = "imidacloprid"
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    readout: str = "mean_hz"
    factors: list[str] | None = None
    n_base: int = Field(32, ge=8, le=1024)
    seed: int = Field(0, ge=0, le=2**31 - 1)
    estimate_only: bool = False


@app.get("/api/dashboard")
def dashboard(
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    graph: GraphName | None = "named",
    assay: LadderAssay = "subgraph",
    genotype: str | None = None,
    n: int = bridge.FAST_DEPENDENCE_N,
    seed: int = 0,
    include_dependence: bool = True,
    include_ladder: bool = True,
    estimate_only: bool = False,
):
    """One call with everything the landing page needs.

    Blocks 1-5 and 10 of the dashboard: the compound overview and headline
    tiles, the factual selectivity panel, the circuit consequence with its
    connectome-dependence verdict, the concentration ladder, the per-number
    evidence records, the trust panel and the generated explanation.
    """
    req = DashboardRequest(
        compound=compound,
        conc_M=conc_M,
        graph=graph,
        assay=assay,
        genotype=genotype,
        n=n,
        seed=seed,
        include_dependence=include_dependence,
        include_ladder=include_ladder,
        estimate_only=estimate_only,
    )
    return _bridged("/api/dashboard", req.model_dump())


@app.post("/api/compare")
def compare(req: CompareRequest):
    """The compare-compounds screen: receptor selectivity beside circuit selectivity."""
    return _bridged("/api/compare", req.model_dump())


@app.get("/api/claims")
def claims_get(
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    assay: AssayName = "subgraph",
    graph: GraphName | None = "named",
    run_assay: bool = True,
):
    """Claim provenance for a fresh run of this compound at this dose."""
    req = ClaimsRequest(
        compound=compound, conc_M=conc_M, assay=assay, graph=graph, run_assay=run_assay
    )
    return _bridged("/api/claims", req.model_dump())


@app.post("/api/claims")
def claims_post(req: ClaimsRequest):
    """Claim provenance for a notebook the caller already holds."""
    return _bridged("/api/claims", req.model_dump())


@app.post("/api/dependence")
def dependence(req: DependenceRequest):
    """Connectome-dependence profile: is this effect a wiring result?"""
    return _bridged("/api/dependence", req.model_dump())


@app.post("/api/dependence/landscape")
def dependence_landscape(req: DependenceLandscapeRequest):
    """The compound x concentration dependence landscape (estimate first)."""
    return _bridged("/api/dependence/landscape", req.model_dump())


@app.post("/api/ablation")
def ablation(req: AblationRequest):
    """The four-level model ablation ladder and what each level adds."""
    return _bridged("/api/ablation", req.model_dump())


@app.post("/api/robustness/stability")
def robustness_stability(req: StabilityRequest):
    """Every pre-registered conclusion re-derived under every admissible rule."""
    return _bridged("/api/robustness/stability", req.model_dump())


@app.post("/api/robustness/thresholds")
def robustness_thresholds(req: ThresholdRequest):
    """The amplify / buffer split recomputed across the threshold grid."""
    return _bridged("/api/robustness/thresholds", req.model_dump())


@app.post("/api/uncertainty/global")
def uncertainty_global(req: GlobalUncertaintyRequest):
    """Variance-based attribution over the model's nine uncertain factors."""
    return _bridged("/api/uncertainty/global", req.model_dump())


@app.post("/api/voi")
def voi(req: VoiRequest):
    """Which experiment would remove the most model variance."""
    return _bridged("/api/voi", req.model_dump())


class GenotypePanelRequest(BaseModel):
    compound: str = "deltamethrin"
    conc_M: float = Field(1e-6, ge=0, le=MAX_CONC_M)
    assay: AssayName = "subgraph"
    genotypes: list[str] | None = None


class IsobologramRequest(BaseModel):
    compound_a: str = "imidacloprid"
    compound_b: str = "fipronil"
    assay: str = "occupancy"
    readout: str | None = None
    effect_frac: float = Field(0.5, gt=0, le=0.999)
    n: int = Field(7, ge=3, le=25)


@app.post("/api/genotype/panel")
def genotype_panel(req: GenotypePanelRequest):
    """Wild type against each relevant resistance allele: potency and circuit."""
    return _bridged("/api/genotype/panel", req.model_dump())


@app.post("/api/mixture/isobologram")
def mixture_isobologram(req: IsobologramRequest):
    """Iso-effect concentration pairs for two compounds, plus the Loewe line."""
    return _bridged("/api/mixture/isobologram", req.model_dump())
