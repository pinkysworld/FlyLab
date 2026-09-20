"""JSON-in / JSON-out mirror of the bench HTTP API, with no web framework.

``flylab/server.py`` is FastAPI + pydantic and cannot run in WebAssembly.  This
module answers exactly the same routes with exactly the same library calls,
defaults, bounds and error contract, so the static (Pyodide) build of the bench
produces the same notebooks -- same numbers, same provenance -- as the served
build::

    from flylab.browser import bridge
    bridge.call("/api/meta")
    bridge.call("/api/assay/subgraph", {"compound": "imidacloprid", "conc_M": 1e-6})

Differences from the HTTP transport, all of them transport-level:

* errors are *returned*, never raised: ``{"error": {"status": ..., "detail": ...}}``
  where ``KeyError``/``FileNotFoundError`` -> 404, ``ValueError`` and
  out-of-bounds fields -> 400, a module that is not importable yet -> 501.
  (FastAPI answers 422 for a malformed request body; the bridge, which has no
  request objects, calls that a 400.)
* ``POST /api/experiment/csv`` returns ``{"csv": ..., "filename": ...,
  "content_type": ...}`` instead of a file download.
* query-string values arrive as strings, so every field is coerced here.

Import cost is deliberately tiny -- stdlib only at module level; numpy, the
graphs and the assay modules are imported inside the handler that needs them,
so ``import flylab.browser.bridge`` stays cheap on a cold browser tab.
"""
from __future__ import annotations

import csv
import io
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

#: must match ``flylab.server.VERSION`` (tests/test_browser_bridge.py asserts it)
VERSION = "0.5.0"

#: bounds shared with the server's pydantic models and documented in /api/meta
MAX_N_REP = 64
MAX_T_MS = 3000.0
MAX_EXPERIMENT_ROWS = 2000
MAX_CONC_M = 1.0

NOT_READY = "module not available yet"

GRAPH_NAMES = ("named", "taste_motor")
ASSAY_NAMES = ("subgraph", "spiking", "taste", "taste_map")

__all__ = ["call", "routes", "version", "BridgeError", "data_root", "VERSION"]


# --------------------------------------------------------------------------
# errors
# --------------------------------------------------------------------------
class BridgeError(Exception):
    """A status + detail the caller should see as ``{"error": {...}}``."""

    def __init__(self, status: int, detail: str):
        super().__init__(detail)
        self.status = int(status)
        self.detail = str(detail)


def _message(exc: Exception) -> str:
    if isinstance(exc, KeyError) and exc.args:
        return str(exc.args[0])
    return str(exc) or exc.__class__.__name__


# --------------------------------------------------------------------------
# data files
# --------------------------------------------------------------------------
#: set once by :func:`_prepare`
_PREPARED = False
_ROOT: Path | None = None


def _candidate_roots() -> list[Path]:
    """Directories that may contain ``data/derived`` and ``data/literature``.

    A repo checkout has them at the repo root; a wheel install has nothing --
    the data is not inside the package -- so the static build unpacks them
    beside the installed package (``<site-packages>/data/...``), which is what
    every module's ``Path(__file__).parents[2] / "data"`` already points at.
    ``FLYLAB_DATA_DIR`` overrides, and names the directory *containing* ``data``.
    """
    out: list[Path] = []
    env = os.environ.get("FLYLAB_DATA_DIR")
    if env:
        p = Path(env)
        out.append(p)
        if p.name == "data":  # tolerate being handed the data dir itself
            out.append(p.parent)
    out.append(Path.cwd())
    try:
        import flylab

        out.append(Path(flylab.__file__).resolve().parents[1])
    except Exception:  # pragma: no cover - flylab is importable by definition
        pass
    return out


def data_root() -> Path | None:
    """The directory holding ``data/derived``, or ``None`` if there is none."""
    for c in _candidate_roots():
        try:
            if (c / "data" / "derived").is_dir():
                return c
        except OSError:  # pragma: no cover - unreadable candidate
            continue
    return None


def _prepare() -> None:
    """Make every module's own data lookup succeed, without editing them.

    ``flylab.circuit.rate``, ``pharm.expression``, ``pharm.genotype``,
    ``pharm.mixtures``, ``analysis.predictions`` and ``validation.rank`` all
    search a CWD-relative ``data/...`` before giving up, so a chdir is enough
    for those.  ``flylab.maps.malecns`` resolves ``CENSUS_FALLBACK`` /
    ``GUSTATORY_SEEDS`` to absolute paths at import time with no fallback, so in
    a layout where those paths are wrong the bridge repoints the two module
    constants (it never rewrites the module).
    """
    global _PREPARED, _ROOT
    if _PREPARED:
        return
    root = data_root()
    _ROOT = root
    if root is not None:
        from flylab.circuit import rate
        from flylab.maps import malecns

        have_graph = any((d / rate.GRAPHS["named"]).exists() for d in rate.DATA_DIRS)
        if not have_graph and Path.cwd().resolve() != root.resolve():
            os.chdir(root)
        if not malecns.CENSUS_FALLBACK.exists():
            malecns.CENSUS_FALLBACK = root / "data" / "derived" / "malecns_census_v1.json"
        if not malecns.GUSTATORY_SEEDS.exists():
            malecns.GUSTATORY_SEEDS = root / "data" / "derived" / "malecns_gustatory_seeds.json"
        os.environ.setdefault("FLYLAB_DATA_DIR", str(root))
    _PREPARED = True


# --------------------------------------------------------------------------
# field coercion (the pydantic request models, without pydantic)
# --------------------------------------------------------------------------
def _missing(payload: dict[str, Any], key: str) -> bool:
    return key not in payload or payload[key] is None


def _num(payload: dict[str, Any], key: str, default: float, lo=None, hi=None, gt=None) -> float:
    if _missing(payload, key):
        value = float(default)
    else:
        try:
            value = float(payload[key])
        except (TypeError, ValueError):
            raise BridgeError(400, f"{key} must be a number")
    if lo is not None and value < lo:
        raise BridgeError(400, f"{key} must be >= {lo}")
    if gt is not None and value <= gt:
        raise BridgeError(400, f"{key} must be > {gt}")
    if hi is not None and value > hi:
        raise BridgeError(400, f"{key} must be <= {hi}")
    return value


def _int(payload: dict[str, Any], key: str, default: int, lo=None, hi=None) -> int:
    if _missing(payload, key):
        value = int(default)
    else:
        try:
            value = int(float(payload[key]))
        except (TypeError, ValueError):
            raise BridgeError(400, f"{key} must be an integer")
    if lo is not None and value < lo:
        raise BridgeError(400, f"{key} must be >= {lo}")
    if hi is not None and value > hi:
        raise BridgeError(400, f"{key} must be <= {hi}")
    return value


def _bool(payload: dict[str, Any], key: str, default: bool) -> bool:
    if _missing(payload, key):
        return bool(default)
    v = payload[key]
    if isinstance(v, str):
        low = v.strip().lower()
        if low in ("true", "1", "yes", "on"):
            return True
        if low in ("false", "0", "no", "off", ""):
            return False
        raise BridgeError(400, f"{key} must be a boolean")
    return bool(v)


def _opt_bool(payload: dict[str, Any], key: str) -> bool | None:
    return None if _missing(payload, key) else _bool(payload, key, False)


def _str(payload: dict[str, Any], key: str, default: str) -> str:
    return default if _missing(payload, key) else str(payload[key])


def _opt_str(payload: dict[str, Any], key: str) -> str | None:
    if _missing(payload, key):
        return None
    v = str(payload[key])
    return v or None


def _choice(payload: dict[str, Any], key: str, default: str | None, allowed) -> str | None:
    v = payload.get(key) if not _missing(payload, key) else default
    if v is None:
        return None
    v = str(v)
    if v not in allowed:
        raise BridgeError(400, f"{key} must be one of {tuple(allowed)}")
    return v


def _floats(payload: dict[str, Any], key: str) -> list[float] | None:
    if _missing(payload, key):
        return None
    raw = payload[key]
    if isinstance(raw, str):
        raw = [p for p in raw.replace(";", ",").split(",") if p.strip()]
    try:
        return [float(x) for x in raw]
    except (TypeError, ValueError):
        raise BridgeError(400, f"{key} must be a list of numbers")


def _strs(payload: dict[str, Any], key: str) -> list[str] | None:
    if _missing(payload, key):
        return None
    raw = payload[key]
    if isinstance(raw, str):
        raw = [p for p in raw.split(",") if p.strip()]
    return [str(x) for x in raw]


# --------------------------------------------------------------------------
# call helpers (mirrors of the server's _lazy / _call / _with_genotype / _timed)
# --------------------------------------------------------------------------
def _lazy(dotted: str, name: str) -> Callable[..., Any]:
    try:
        module = __import__(dotted, fromlist=[name])
        return getattr(module, name)
    except (ImportError, AttributeError):
        raise BridgeError(501, NOT_READY)


def _invoke(fn: Callable[..., Any], /, **kw: Any) -> Any:
    try:
        return fn(**kw)
    except TypeError as exc:
        if "unexpected keyword" in str(exc) or "required positional" in str(exc):
            raise BridgeError(501, f"module signature not compatible yet: {exc}")
        raise


def _with_genotype(fn: Callable[..., Any], genotype: str | None, /, *args: Any, **kw: Any) -> Any:
    if genotype is None:
        return fn(*args, **kw)
    try:
        return fn(*args, genotype=genotype, **kw)
    except TypeError as exc:
        if "genotype" in str(exc):
            raise BridgeError(501, "genotype not supported yet")
        raise


def _timed(fn: Callable[[], Any]) -> tuple[Any, float]:
    t0 = time.perf_counter()
    out = fn()
    return out, (time.perf_counter() - t0) * 1000.0


def _stamp(nb: Any, runtime_ms: float) -> Any:
    if isinstance(nb, dict):
        nb.setdefault("readouts", {})
        if isinstance(nb["readouts"], dict):
            nb["readouts"]["runtime_ms"] = round(runtime_ms, 2)
    return nb


# --------------------------------------------------------------------------
# shared request shapes
# --------------------------------------------------------------------------
def _assay_request(p: dict[str, Any]) -> dict[str, Any]:
    return {
        "compound": _opt_str(p, "compound"),
        "conc_M": _num(p, "conc_M", 0.0, lo=0, hi=MAX_CONC_M),
        "sugar_hz": _num(p, "sugar_hz", 150.0, lo=0, hi=1000),
        "bitter_hz": _num(p, "bitter_hz", 0.0, lo=0, hi=1000),
        "drive_hz": _num(p, "drive_hz", 40.0, lo=0, hi=1000),
        "steps": _int(p, "steps", 80, lo=1, hi=1000),
        "graph": _choice(p, "graph", None, GRAPH_NAMES),
        "genotype": _opt_str(p, "genotype"),
    }


def _graph_kw(graph: str | None, assay: str) -> dict[str, Any]:
    """The server only forwards ``graph`` to the assays that take one."""
    return {"graph": graph} if graph and assay in ("subgraph", "spiking", "taste_map") else {}


# --------------------------------------------------------------------------
# meta / library
# --------------------------------------------------------------------------
def _graph_meta() -> dict[str, Any]:
    from flylab.circuit.rate import GRAPHS, load_graph

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
    from flylab.pharm.occupancy import list_compounds

    rows = []
    for key in list_compounds(lib):
        spec = lib["compounds"][key]
        targets = [
            {
                "receptor": receptor,
                "ec50_M": float(r["ec50_M"]),
                "direction": r.get("direction", "none"),
                "evidence_tier": r.get("evidence_tier", "class_placeholder"),
                "source": r.get("source", ""),
            }
            for receptor, r in (spec.get("receptors") or {}).items()
            if receptor.startswith("insect_") and r.get("direction") not in (None, "none")
        ]
        targets.sort(key=lambda r: r["ec50_M"])
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


def _h_health(p: dict[str, Any]) -> dict[str, Any]:
    return {"ok": True, "version": VERSION}


def _h_meta(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.maps.malecns import MAP_CITATION, MAP_ID, load_gustatory_seeds
    from flylab.pharm.exposure import ROUTE_DEFAULTS
    from flylab.pharm.mechanisms import GAIN_KEYS, mechanism_table_rows
    from flylab.pharm.occupancy import (
        library_sha256,
        load_library,
        receptor_table,
        selectivity_pairs,
    )

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


def _h_census(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.maps.malecns import load_census

    return load_census()


def _h_drugs(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.taste import library_keys

    return {"compounds": library_keys()}


def _h_drug(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.pharm.occupancy import load_library

    key = p.get("key")
    if not key:
        raise BridgeError(400, "key is required")
    lib = load_library()
    k = str(key).lower().strip()
    if k not in lib["compounds"]:
        raise BridgeError(404, f"unknown compound {key!r}")
    return {"key": k, **lib["compounds"][k]}


def _h_occupancy(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.pharm.occupancy import compare_compound

    compound = _opt_str(p, "compound")
    if compound is None:
        raise BridgeError(400, "compound is required")
    conc = _num(p, "conc_M", 0.0)
    if conc < 0:
        raise BridgeError(400, "conc_M must be >= 0")
    return compare_compound(compound, conc)


def _h_occupancy_curve(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.pharm.occupancy import occupancy_curve, receptor_table

    compound = _opt_str(p, "compound")
    if compound is None:
        raise BridgeError(400, "compound is required")
    return {
        "compound": compound,
        "points": occupancy_curve(compound),
        "receptors": receptor_table(),
    }


def _h_occupancy_ci(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.pharm.uncertainty import occupancy_ci

    compound = _opt_str(p, "compound")
    if compound is None:
        raise BridgeError(400, "compound is required")
    conc = _num(p, "conc_M", 0.0)
    if conc < 0:
        raise BridgeError(400, "conc_M must be >= 0")
    n = _int(p, "n", 500)
    if not 2 <= n <= 5000:
        raise BridgeError(400, "n must be between 2 and 5000")
    return occupancy_ci(
        compound, conc, sd_log10=_num(p, "sd_log10", 0.3), n=n, seed=_int(p, "seed", 0)
    )


# --------------------------------------------------------------------------
# assays
# --------------------------------------------------------------------------
def _h_assay_taste(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.taste import run_taste_assay

    r = _assay_request(p)
    nb, ms = _timed(
        lambda: _with_genotype(
            run_taste_assay,
            r["genotype"],
            r["compound"],
            r["conc_M"],
            r["sugar_hz"],
            r["bitter_hz"],
        )
    )
    return _stamp(nb, ms)


def _h_assay_taste_dose(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.taste import dose_response

    compound = _opt_str(p, "compound")
    if compound is None:
        raise BridgeError(400, "compound is required")
    return {"compound": compound, "points": dose_response(compound, [10**e for e in range(-10, -4)])}


def _h_assay_cns(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.wholens import run_wholens_assay

    r = _assay_request(p)
    nb, ms = _timed(
        lambda: _with_genotype(run_wholens_assay, r["genotype"], r["compound"], r["conc_M"])
    )
    return _stamp(nb, ms)


def _h_assay_subgraph(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.subgraph import run_subgraph_assay

    r = _assay_request(p)
    nb, ms = _timed(
        lambda: _with_genotype(
            run_subgraph_assay,
            r["genotype"],
            r["compound"],
            r["conc_M"],
            drive_hz=r["drive_hz"],
            steps=r["steps"],
            graph=r["graph"],
        )
    )
    return _stamp(nb, ms)


def _h_assay_subgraph_dose(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.subgraph import dose_response_subgraph

    compound = _opt_str(p, "compound")
    if compound is None:
        raise BridgeError(400, "compound is required")
    return {"compound": compound, "points": dose_response_subgraph(compound)}


def _h_assay_spiking(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.spiking import run_spiking_assay

    compound = _opt_str(p, "compound")
    conc = _num(p, "conc_M", 0.0, lo=0, hi=MAX_CONC_M)
    drive_hz = _num(p, "drive_hz", 40.0, lo=0, hi=1000)
    t_ms = _num(p, "t_ms", 500.0, gt=0, hi=MAX_T_MS)
    seed = _int(p, "seed", 0, lo=0, hi=2**31 - 1)
    graph = _choice(p, "graph", None, GRAPH_NAMES)
    genotype = _opt_str(p, "genotype")
    nb, ms = _timed(
        lambda: _with_genotype(
            run_spiking_assay,
            genotype,
            compound,
            conc,
            drive_hz=drive_hz,
            t_ms=t_ms,
            seed=seed,
            graph=graph,
        )
    )
    return _stamp(nb, ms)


def _h_assay_taste_map(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.taste_map import run_taste_map_assay

    compound = _opt_str(p, "compound")
    conc = _num(p, "conc_M", 0.0, lo=0, hi=MAX_CONC_M)
    sugar_hz = _num(p, "sugar_hz", 150.0, lo=0, hi=1000)
    bitter_hz = _num(p, "bitter_hz", 0.0, lo=0, hi=1000)
    engine = _choice(p, "engine", "rate", ("rate", "lif"))
    seed = _int(p, "seed", 0, lo=0, hi=2**31 - 1)
    graph = _choice(p, "graph", "taste_motor", GRAPH_NAMES)
    genotype = _opt_str(p, "genotype")
    nb, ms = _timed(
        lambda: _with_genotype(
            run_taste_map_assay,
            genotype,
            compound,
            conc,
            sugar_hz=sugar_hz,
            bitter_hz=bitter_hz,
            engine=engine,
            seed=seed,
            graph=graph or "taste_motor",
        )
    )
    return _stamp(nb, ms)


def _h_assay_ensemble(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.ensemble import run_ensemble

    assay = _choice(p, "assay", "subgraph", ASSAY_NAMES)
    compound = _opt_str(p, "compound")
    conc = _num(p, "conc_M", 0.0, lo=0, hi=MAX_CONC_M)
    n_rep = _int(p, "n_rep", 8, lo=2, hi=MAX_N_REP)
    sd = _num(p, "ec50_sd_log10", 0.3, lo=0, hi=3)
    seed = _int(p, "seed", 0, lo=0, hi=2**31 - 1)
    graph = _choice(p, "graph", None, GRAPH_NAMES)
    genotype = _opt_str(p, "genotype")
    kw = _graph_kw(graph, assay)
    nb, ms = _timed(
        lambda: _with_genotype(
            run_ensemble,
            genotype,
            assay,
            compound,
            conc,
            n_rep=n_rep,
            ec50_sd_log10=sd,
            seed=seed,
            **kw,
        )
    )
    return _stamp(nb, ms)


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------
def _h_sensitivity(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.ensemble import sensitivity

    assay = _choice(p, "assay", "subgraph", ASSAY_NAMES)
    compound = _str(p, "compound", "imidacloprid")
    conc = _num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M)
    readout = _str(p, "readout", "mean_hz")
    factor = _num(p, "factor", 2.0, gt=1, hi=100)
    seed = _int(p, "seed", 0, lo=0, hi=2**31 - 1)
    graph = _choice(p, "graph", None, GRAPH_NAMES)
    genotype = _opt_str(p, "genotype")
    kw = _graph_kw(graph, assay)
    rows, ms = _timed(
        lambda: _with_genotype(
            sensitivity,
            genotype,
            assay,
            compound,
            conc,
            readout=readout,
            factor=factor,
            seed=seed,
            **kw,
        )
    )
    return {
        "assay": assay,
        "compound": compound,
        "concentration_M": conc,
        "readout": readout,
        "factor": factor,
        "rows": rows,
        "runtime_ms": round(ms, 2),
        "label": "model_derived",
        "warnings": [
            "One-at-a-time tornado over teaching parameters; it reports model "
            "sensitivity, not biological variability.",
        ],
    }


def _h_ic50(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.ensemble import circuit_ic50

    assay = _choice(p, "assay", "subgraph", ASSAY_NAMES)
    compound = _str(p, "compound", "imidacloprid")
    readout = _str(p, "readout", "mean_hz")
    concs = _floats(p, "concs_M")
    n_boot = _int(p, "n_boot", 200, lo=0, hi=2000)
    n_rep = _int(p, "n_rep", 4, lo=1, hi=MAX_N_REP)
    seed = _int(p, "seed", 0, lo=0, hi=2**31 - 1)
    graph = _choice(p, "graph", None, GRAPH_NAMES)
    genotype = _opt_str(p, "genotype")
    if concs is not None:
        if not concs:
            raise BridgeError(400, "concs_M must not be empty")
        if len(concs) > 64:
            raise BridgeError(400, "concs_M is limited to 64 points")
        if any(c < 0 for c in concs):
            raise BridgeError(400, "concentrations must be >= 0")
    kw = _graph_kw(graph, assay)
    out, ms = _timed(
        lambda: _with_genotype(
            circuit_ic50,
            genotype,
            assay,
            compound,
            readout=readout,
            concs=concs,
            n_boot=n_boot,
            n_rep=n_rep,
            seed=seed,
            **kw,
        )
    )
    out["runtime_ms"] = round(ms, 2)
    return out


# --------------------------------------------------------------------------
# graph viewer
# --------------------------------------------------------------------------
def _h_graph(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.analysis.layout import graph_for_viewer
    from flylab.circuit.rate import DEFAULT_GAINS, compute_gains, load_graph, rate_network

    graph = _choice(p, "graph", "named", GRAPH_NAMES) or "named"
    min_weight = _num(p, "min_weight", 5.0)
    max_edges = _int(p, "max_edges", 2000)
    compound = _opt_str(p, "compound")
    conc_M = _num(p, "conc_M", 0.0)
    drive_hz = _num(p, "drive_hz", 40.0)
    steps = _int(p, "steps", 80)
    if min_weight < 0:
        raise BridgeError(400, "min_weight must be >= 0")
    if not 1 <= max_edges <= 20000:
        raise BridgeError(400, "max_edges must be between 1 and 20000")
    if conc_M < 0:
        raise BridgeError(400, "conc_M must be >= 0")

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


def _h_graph_impact(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.analysis.impact import grn_to_mn9_paths, summarize_impact
    from flylab.circuit.rate import compute_gains, load_graph

    compound = _opt_str(p, "compound")
    conc_M = _num(p, "conc_M", 0.0, lo=0, hi=MAX_CONC_M)
    name = _choice(p, "graph", "named", GRAPH_NAMES) or "named"
    drive_hz = _num(p, "drive_hz", 40.0, lo=0, hi=1000)
    top_edges = _int(p, "top_edges", 25, lo=1, hi=200)
    top_nodes = _int(p, "top_nodes", 25, lo=1, hi=200)
    top_paths = _int(p, "top_paths", 10, lo=1, hi=100)
    max_len = _int(p, "max_len", 2, lo=1, hi=4)
    out, ms = _timed(
        lambda: summarize_impact(
            compound,
            conc_M,
            graph=name,
            drive_hz=drive_hz,
            top_edges=top_edges,
            top_nodes=top_nodes,
            top_paths=top_paths,
            max_len=max_len,
        )
    )
    if name == "taste_motor":
        g = load_graph(name)
        gains, _ = compute_gains(compound, conc_M)
        out["grn_to_mn9_paths"] = grn_to_mn9_paths(g, max_len=3, top=top_paths, gains=gains)
    out["runtime_ms"] = round(ms, 2)
    return out


# --------------------------------------------------------------------------
# exposure
# --------------------------------------------------------------------------
def _h_exposure(p: dict[str, Any]) -> dict[str, Any]:
    from flylab.pharm.exposure import exposure_profile

    compound = _str(p, "compound", "imidacloprid")
    dose = _num(p, "dose", 1.0, lo=0, hi=1e6)
    route = _choice(p, "route", "feeding", ("feeding", "topical", "bath"))
    t_h = _num(p, "t_h", 24.0, gt=0, hi=240)
    dt_h = _num(p, "dt_h", 0.05, gt=0, hi=24)
    if dt_h > t_h:
        raise BridgeError(400, "dt_h must be <= t_h")
    out, ms = _timed(lambda: exposure_profile(compound, dose, route, t_h=t_h, dt_h=dt_h))
    out["runtime_ms"] = round(ms, 2)
    return out


# --------------------------------------------------------------------------
# experiment
# --------------------------------------------------------------------------
def _experiment_request(p: dict[str, Any]) -> dict[str, Any]:
    compounds = _strs(p, "compounds")
    concs = _floats(p, "concs_M")
    return {
        "assay": _choice(p, "assay", "subgraph", ASSAY_NAMES),
        "compounds": ["imidacloprid"] if compounds is None else compounds,
        "concs_M": [1e-9, 1e-8, 1e-7, 1e-6, 1e-5] if concs is None else concs,
        "replicates": _int(p, "replicates", 1, lo=1, hi=MAX_N_REP),
        "seed": _int(p, "seed", 0, lo=0, hi=2**31 - 1),
        "readouts": _strs(p, "readouts"),
        "jitter_log10": _num(p, "jitter_log10", 0.0, lo=0, hi=3),
        "drive_hz": None if _missing(p, "drive_hz") else _num(p, "drive_hz", 0.0, lo=0, hi=1000),
        "include_vehicle": _bool(p, "include_vehicle", True),
        "graph": _choice(p, "graph", None, GRAPH_NAMES),
        "keep_notebooks": _bool(p, "keep_notebooks", False),
        "genotype": _opt_str(p, "genotype"),
    }


def _run_design(req: dict[str, Any]) -> dict[str, Any]:
    from flylab.assays.experiment import ExperimentDesign, run_experiment

    n_rows = len(req["compounds"]) * len(req["concs_M"]) * req["replicates"]
    if n_rows > MAX_EXPERIMENT_ROWS:
        raise BridgeError(
            400,
            f"design would produce {n_rows} rows; the limit is "
            f"{MAX_EXPERIMENT_ROWS}. Reduce compounds, concentrations or replicates.",
        )
    if any(c < 0 for c in req["concs_M"]):
        raise BridgeError(400, "concentrations must be >= 0")
    data = {
        k: v
        for k, v in req.items()
        if k not in ("genotype", "readouts") and v is not None
    }
    if req["readouts"]:
        data["readouts"] = req["readouts"]
    design = ExperimentDesign(**data)
    out, ms = _timed(lambda: _with_genotype(run_experiment, req["genotype"], design))
    out["runtime_ms"] = round(ms, 2)
    return out


def _h_experiment(p: dict[str, Any]) -> dict[str, Any]:
    req = _experiment_request(p)
    out = _run_design(req)
    if not req["keep_notebooks"]:
        out.pop("notebooks", None)
    return out


def _h_experiment_csv(p: dict[str, Any]) -> dict[str, Any]:
    """The HTTP route streams a file; here the same text comes back as JSON."""
    from flylab.assays.experiment import rows_to_csv

    req = _experiment_request(p)
    out = _run_design(req)
    text = out.get("csv") or rows_to_csv(out.get("rows", []) + out.get("vehicle_rows", []))
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return {
        "csv": text,
        "filename": f"flylab_experiment_{stamp}.csv",
        "content_type": "text/csv; charset=utf-8",
    }


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


def _h_live_lab(p: dict[str, Any]) -> dict[str, Any]:
    text = str(p.get("csv_text") or "")
    if not text.strip():
        raise BridgeError(400, "csv_text is empty")
    try:
        dialect = csv.Sniffer().sniff(text[:4096], delimiters=",;\t|")
    except csv.Error:
        dialect = csv.excel
    reader = csv.DictReader(io.StringIO(text), dialect=dialect)
    if not reader.fieldnames:
        raise BridgeError(400, "csv_text has no header row")
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
        raise BridgeError(400, "csv_text has a header but no data rows")

    notebook = dict(p.get("notebook") or {})
    notebook["live_lab"] = {
        "assay_name": _str(p, "assay_name", "live_lab"),
        "source": _str(p, "source", "imported CSV"),
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
# forward-compatible endpoints
# --------------------------------------------------------------------------
def _h_null(p: dict[str, Any]) -> dict[str, Any]:
    fn = _lazy("flylab.analysis.nullmodels", "null_distribution")
    graph = _choice(p, "graph", None, GRAPH_NAMES)
    kw = {"graph": graph} if graph else {}
    return _invoke(
        fn,
        assay=_choice(p, "assay", "subgraph", ASSAY_NAMES),
        compound=_str(p, "compound", "imidacloprid"),
        conc_M=_num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M),
        mode=_str(p, "mode", "sign_permute"),
        n=_int(p, "n", 20, lo=1, hi=500),
        seed=_int(p, "seed", 0, lo=0, hi=2**31 - 1),
        **kw,
    )


def _h_null_panel(p: dict[str, Any]) -> dict[str, Any]:
    fn = _lazy("flylab.analysis.nullmodels", "null_panel")
    modes = _strs(p, "modes")
    include_taste_map = _opt_bool(p, "include_taste_map")
    kw: dict[str, Any] = {}
    if modes:
        kw["modes"] = modes
    if include_taste_map is not None:
        kw["include_taste_map"] = include_taste_map
    return _invoke(
        fn,
        compound=_str(p, "compound", "imidacloprid"),
        conc_M=_num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M),
        n=_int(p, "n", 8, lo=1, hi=500),
        seed=_int(p, "seed", 0, lo=0, hi=2**31 - 1),
        **kw,
    )


def _h_selectivity(p: dict[str, Any]) -> dict[str, Any]:
    conc = _num(p, "conc_M", 1e-6)
    if conc < 0:
        raise BridgeError(400, "conc_M must be >= 0")
    fn = _lazy("flylab.analysis.selectivity", "receptor_selectivity_table")
    return _invoke(fn, conc_M=conc)


def _h_selectivity_landscape(p: dict[str, Any]) -> dict[str, Any]:
    fn = _lazy("flylab.analysis.selectivity", "selectivity_landscape")
    kw: dict[str, Any] = {"assay": _choice(p, "assay", "subgraph", ASSAY_NAMES)}
    graphs = _strs(p, "graphs")
    compounds = _strs(p, "compounds")
    if graphs:
        kw["graphs"] = graphs
    if compounds:
        kw["compounds"] = compounds
    return _invoke(fn, **kw)


def _h_circuit_si(p: dict[str, Any]) -> dict[str, Any]:
    fn = _lazy("flylab.analysis.selectivity", "circuit_selectivity_index")
    return _invoke(
        fn,
        compound=_str(p, "compound", "imidacloprid"),
        assay=_choice(p, "assay", "subgraph", ASSAY_NAMES),
        graph=_choice(p, "graph", "named", GRAPH_NAMES),
    )


def _h_predictions(p: dict[str, Any]) -> dict[str, Any]:
    n_rep = _int(p, "n_rep", 4)
    if not 1 <= n_rep <= MAX_N_REP:
        raise BridgeError(400, f"n_rep must be 1..{MAX_N_REP}")
    fn = _lazy("flylab.analysis.predictions", "prediction_table")
    return _invoke(fn, n_rep=n_rep, seed=_int(p, "seed", 0))


def _h_genotypes(p: dict[str, Any]) -> Any:
    fn = _lazy("flylab.pharm.genotype", "list_genotypes")
    return _invoke(fn)


def _h_mixture(p: dict[str, Any]) -> dict[str, Any]:
    raw = p.get("components") or []
    if not raw:
        raise BridgeError(400, "mixture needs at least one component")
    components = []
    for item in raw:
        item = dict(item or {})
        if not item.get("compound"):
            raise BridgeError(400, "each mixture component needs a compound")
        components.append(
            {
                "compound": str(item["compound"]),
                "conc_M": _num(item, "conc_M", 0.0, lo=0, hi=MAX_CONC_M),
            }
        )
    fn = _lazy("flylab.pharm.mixtures", "mixture_assay")
    kw: dict[str, Any] = {
        "assay": _choice(p, "assay", "subgraph", ASSAY_NAMES),
        "model": _choice(p, "model", "bliss", ("bliss", "loewe")),
    }
    graph = _choice(p, "graph", None, GRAPH_NAMES)
    if graph:
        kw["graph"] = graph
    return _invoke(fn, components=components, **kw)


def _h_expression(p: dict[str, Any]) -> Any:
    fn = _lazy("flylab.pharm.expression", "expression_table")
    return _invoke(fn)


def _h_validation(p: dict[str, Any]) -> Any:
    fn = _lazy("flylab.validation.rank", "validate_all")
    return _invoke(fn)


# --------------------------------------------------------------------------
# routing table
# --------------------------------------------------------------------------
#: route -> (http method, handler).  The method is documentation only: the
#: bridge has one calling convention, ``call(route, payload)``.
ROUTES: dict[str, tuple[str, Callable[[dict[str, Any]], Any]]] = {
    "/api/health": ("GET", _h_health),
    "/api/meta": ("GET", _h_meta),
    "/api/census": ("GET", _h_census),
    "/api/drugs": ("GET", _h_drugs),
    "/api/drugs/{key}": ("GET", _h_drug),
    "/api/occupancy": ("GET", _h_occupancy),
    "/api/occupancy/curve": ("GET", _h_occupancy_curve),
    "/api/occupancy/ci": ("GET", _h_occupancy_ci),
    "/api/assay/taste": ("POST", _h_assay_taste),
    "/api/assay/taste/dose-response": ("GET", _h_assay_taste_dose),
    "/api/assay/cns": ("POST", _h_assay_cns),
    "/api/assay/subgraph": ("POST", _h_assay_subgraph),
    "/api/assay/subgraph/dose-response": ("GET", _h_assay_subgraph_dose),
    "/api/assay/spiking": ("POST", _h_assay_spiking),
    "/api/assay/taste-map": ("POST", _h_assay_taste_map),
    "/api/assay/ensemble": ("POST", _h_assay_ensemble),
    "/api/analysis/sensitivity": ("POST", _h_sensitivity),
    "/api/analysis/ic50": ("POST", _h_ic50),
    "/api/graph": ("GET", _h_graph),
    "/api/graph/impact": ("POST", _h_graph_impact),
    "/api/exposure": ("POST", _h_exposure),
    "/api/experiment": ("POST", _h_experiment),
    "/api/experiment/csv": ("POST", _h_experiment_csv),
    "/api/notebook/live-lab": ("POST", _h_live_lab),
    "/api/analysis/null": ("POST", _h_null),
    "/api/analysis/null-panel": ("POST", _h_null_panel),
    "/api/analysis/selectivity": ("GET", _h_selectivity),
    "/api/analysis/selectivity-landscape": ("POST", _h_selectivity_landscape),
    "/api/analysis/circuit-si": ("POST", _h_circuit_si),
    "/api/predictions": ("GET", _h_predictions),
    "/api/genotypes": ("GET", _h_genotypes),
    "/api/mixture": ("POST", _h_mixture),
    "/api/expression": ("GET", _h_expression),
    "/api/validation": ("GET", _h_validation),
}


def routes() -> list[str]:
    """Every route the bridge answers, in the server's declaration order."""
    return list(ROUTES)


def version() -> dict[str, Any]:
    """Version / capability block, cheap enough to call before anything loads."""
    import platform
    import sys

    return {
        "version": VERSION,
        "notebook_version": "0.3",
        "backend": "browser-bridge",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "n_routes": len(ROUTES),
        "data_root": str(data_root()) if data_root() else None,
        "limits": {
            "max_n_rep": MAX_N_REP,
            "max_t_ms": MAX_T_MS,
            "max_experiment_rows": MAX_EXPERIMENT_ROWS,
            "max_conc_M": MAX_CONC_M,
        },
    }


def _split(route: str) -> tuple[str, dict[str, Any]]:
    """Split a request path into (route, query params).

    The UI builds GET URLs with a query string; the bridge accepts either that
    or an explicit payload, so ``call("/api/occupancy?compound=nicotine")`` and
    ``call("/api/occupancy", {"compound": "nicotine"})`` are the same request.
    """
    from urllib.parse import parse_qsl, urlsplit

    parts = urlsplit(str(route or ""))
    path = parts.path or "/"
    if len(path) > 1 and path.endswith("/"):
        path = path.rstrip("/")
    query = {k: v for k, v in parse_qsl(parts.query, keep_blank_values=False)}
    return path, query


def _resolve(path: str) -> tuple[Callable[[dict[str, Any]], Any], dict[str, Any]]:
    if path in ROUTES:
        return ROUTES[path][1], {}
    if path.startswith("/api/drugs/"):
        key = path[len("/api/drugs/") :]
        if key and "/" not in key:
            from urllib.parse import unquote

            return _h_drug, {"key": unquote(key)}
    raise BridgeError(404, f"no such route {path!r}")


def call(route: str, payload: dict[str, Any] | None = None) -> Any:
    """Answer one bench request.

    ``route`` is the server path (``"/api/assay/subgraph"``), optionally with a
    query string.  ``payload`` carries the body of a POST route or the query
    parameters of a GET route.  Returns the same JSON-able object the HTTP
    route returns, or ``{"error": {"status": ..., "detail": ...}}`` -- it never
    raises for a bad request, a missing file or an unknown compound.
    """
    try:
        _prepare()
        path, query = _split(route)
        handler, path_params = _resolve(path)
        data: dict[str, Any] = dict(query)
        if isinstance(payload, dict):
            data.update(payload)
        data.update(path_params)
        return handler(data)
    except BridgeError as exc:
        return {"error": {"status": exc.status, "detail": exc.detail}}
    except KeyError as exc:
        return {"error": {"status": 404, "detail": _message(exc)}}
    except FileNotFoundError as exc:
        return {"error": {"status": 404, "detail": _message(exc)}}
    except ValueError as exc:
        return {"error": {"status": 400, "detail": _message(exc)}}


def call_json(route: str, payload_json: str | None = None) -> str:
    """``call`` for hosts that can only pass strings (Pyodide's JS bridge)."""
    import json

    payload = json.loads(payload_json) if payload_json else None
    return json.dumps(call(route, payload), default=str)
