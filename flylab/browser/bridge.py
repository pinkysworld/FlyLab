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

__all__ = [
    "call",
    "handle",
    "routes",
    "version",
    "BridgeError",
    "data_root",
    "VERSION",
    "build_dashboard",
    "build_compare",
    "build_claims",
    "FAST_DEPENDENCE_N",
]


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


def _h_genotype_panel(p: dict[str, Any]) -> dict[str, Any]:
    fn = _lazy("flylab.pharm.genotype", "genotype_panel")
    genotypes = _strs(p, "genotypes")
    kw: dict[str, Any] = {}
    if genotypes:
        kw["genotypes"] = genotypes
    compound = _opt_str(p, "compound")
    if compound is None:
        raise BridgeError(400, "compound is required")
    return _invoke(
        fn,
        compound=compound,
        conc_M=_num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M),
        assay=_choice(p, "assay", "subgraph", ASSAY_NAMES),
        **kw,
    )


def _h_isobologram(p: dict[str, Any]) -> dict[str, Any]:
    fn = _lazy("flylab.pharm.mixtures", "isobologram")
    a = _opt_str(p, "compound_a")
    b = _opt_str(p, "compound_b")
    if not a or not b:
        raise BridgeError(400, "compound_a and compound_b are required")
    kw: dict[str, Any] = {}
    readout = _opt_str(p, "readout")
    if readout:
        kw["readout"] = readout
    return _invoke(
        fn,
        compound_a=a,
        compound_b=b,
        assay=_str(p, "assay", "occupancy"),
        effect_frac=_num(p, "effect_frac", 0.5, gt=0, hi=0.999),
        n=_int(p, "n", 7, lo=3, hi=25),
        **kw,
    )


def _h_expression(p: dict[str, Any]) -> Any:
    fn = _lazy("flylab.pharm.expression", "expression_table")
    return _invoke(fn)


def _h_validation(p: dict[str, Any]) -> Any:
    fn = _lazy("flylab.validation.rank", "validate_all")
    return _invoke(fn)


# --------------------------------------------------------------------------
# dashboard / decision layer (shared with flylab.server, which imports these)
# --------------------------------------------------------------------------
# The dashboard is the bench's interpretation layer: one call that answers
# "what does this compound do, at this dose, and how much of the answer is
# measured?".  The assembly lives here rather than in ``flylab/server.py``
# because the served bench and the browser build must return byte-identical
# payloads -- sharing the function is the only way to guarantee that.

#: shuffle count the browser (and every dashboard call) defaults to.  Kept in
#: step with :data:`flylab.analysis.dependence.FAST_N` by
#: ``tests/test_browser_bridge.py``.
FAST_DEPENDENCE_N = 20

#: dependence null modes used by the dashboard's "is it wiring?" verdict
DEPENDENCE_MODES = (
    "sign_permute",
    "weight_permute",
    "rewire_degree_preserving",
    "erdos_renyi",
)

#: per-compound cost of a compare row, measured on the committed cut
COMPARE_SECONDS_PER_COMPOUND = 0.6

#: assays whose contrast the rate engine can compute directly.  The dashboard's
#: concentration ladder is built on that path, so it accepts these two only --
#: the LIF assay would answer with zeros rather than refuse.
LADDER_ASSAYS = ("subgraph", "taste_map")

DASHBOARD_WARNINGS = [
    "Every number on this page is simulated. FlyLab has never dosed a fly.",
    "Engagement is not occupancy: only a Kd/Ki row reports fractional receptor "
    "occupancy, an EC50/IC50 row reports normalised functional engagement.",
    "A missing value means 'not modelled', never zero.",
    "There is no blended confidence score here on purpose: the layers fail "
    "independently and averaging them would hide which one is weak.",
]


def _chip(row: dict[str, Any]) -> str:
    """The evidence chip for one library row (see analysis/claims.CLASSIFICATIONS)."""
    value = row.get("param_value_M", row.get("ec50_M"))
    if value is None or row.get("engagement_model") == "not_modelled":
        return "NOT MODELLED"
    if row.get("evidence_tier") == "literature_order":
        return "LITERATURE"
    if row.get("evidence_tier") == "measured_fit":
        return "MODEL-DERIVED"
    return "MODEL-ASSUMPTION"


def _engagement_of(row: dict[str, Any]) -> float | None:
    value = row.get("engagement", row.get("occupancy"))
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if out == out else None


def _evidence_rows(occ: dict[str, Any]) -> list[dict[str, Any]]:
    """``compare_compound`` rows, each with its provenance record and chip."""
    from flylab.pharm.evidence import describe

    out = []
    for row in occ.get("receptors") or []:
        try:
            rec = dict(describe(row))
        except Exception:  # pragma: no cover - a malformed row must not 500
            rec = {"receptor": row.get("receptor")}
        rec["engagement"] = _engagement_of(row)
        rec["direction"] = row.get("direction")
        rec["organism"] = "insect" if str(row.get("receptor", "")).startswith("insect_") else "vertebrate"
        rec["classification"] = _chip(row)
        rec["not_modelled_reason"] = row.get("not_modelled_reason")
        out.append(rec)
    return out


def _split_rows(evidence: list[dict[str, Any]]):
    insect = [r for r in evidence if r["organism"] == "insect" and r["engagement"] is not None]
    vert = [r for r in evidence if r["organism"] == "vertebrate" and r["engagement"] is not None]
    return insect, vert


def _best_pair(occ: dict[str, Any]) -> tuple[str | None, dict[str, Any] | None]:
    """The sourced insect/vertebrate pair with the widest ratio."""
    best_name, best = None, None
    for name, pair in (occ.get("selectivity") or {}).items():
        if pair.get("placeholder") or not pair.get("comparable"):
            continue
        ratio = pair.get("ratio")
        if ratio is None:
            continue
        if best is None or float(ratio) > float(best.get("ratio") or 0.0):
            best_name, best = name, pair
    return best_name, best


def _threshold_conc(value_M: float | None, n: float | None, target: float) -> float | None:
    """Closed-form concentration at which a Hill row reaches ``target`` engagement."""
    try:
        v = float(value_M)
        hill = float(n or 1.0) or 1.0
    except (TypeError, ValueError):
        return None
    if not (v > 0) or not (0.0 < target < 1.0):
        return None
    return v * (target / (1.0 - target)) ** (1.0 / hill)


def _hill(conc: float, value_M: float | None, n: float | None) -> float | None:
    try:
        v = float(value_M)
        hill = float(n or 1.0) or 1.0
        c = float(conc)
    except (TypeError, ValueError):
        return None
    if not (v > 0) or c < 0:
        return None
    if c == 0:
        return 0.0
    cn = c**hill
    return cn / (v**hill + cn)


# -- block 1: compound overview + headline tiles ---------------------------
def _headline(occ: dict[str, Any], evidence: list[dict[str, Any]]) -> dict[str, Any]:
    insect, vert = _split_rows(evidence)
    top_insect = max(insect, key=lambda r: r["engagement"]) if insect else None
    top_vert = max(vert, key=lambda r: r["engagement"]) if vert else None
    pair_name, pair = _best_pair(occ)
    sourced = [r for r in evidence if r["classification"] == "LITERATURE"]
    tier = (
        "literature_order"
        if any(r["organism"] == "insect" for r in sourced)
        else "class_placeholder"
    )
    return {
        "insect_engagement": {
            "value": top_insect["engagement"] if top_insect else None,
            "receptor": top_insect["receptor"] if top_insect else None,
            "param_type": top_insect.get("param_type") if top_insect else None,
            "param_value_M": top_insect.get("param_value_M") if top_insect else None,
            # v0.6.1: the tile has to say what KIND of number it is showing and
            # how far its source is from the modelled target.
            "engagement_model": top_insect.get("engagement_model") if top_insect else "not_modelled",
            "evidence_distance": top_insect.get("evidence_distance") if top_insect else "E4",
            "evidence_distance_label": (
                top_insect.get("evidence_distance_label") if top_insect else "E4 (unsupported)"
            ),
            "provenance_warning": top_insect.get("provenance_warning") if top_insect else None,
            "classification": top_insect["classification"] if top_insect else "NOT MODELLED",
            "unit": "engagement 0-1",
        },
        "vertebrate_engagement": {
            "value": top_vert["engagement"] if top_vert else None,
            "receptor": top_vert["receptor"] if top_vert else None,
            "param_type": top_vert.get("param_type") if top_vert else None,
            "param_value_M": top_vert.get("param_value_M") if top_vert else None,
            "engagement_model": top_vert.get("engagement_model") if top_vert else "not_modelled",
            "evidence_distance": top_vert.get("evidence_distance") if top_vert else "E4",
            "evidence_distance_label": (
                top_vert.get("evidence_distance_label") if top_vert else "E4 (unsupported)"
            ),
            "provenance_warning": top_vert.get("provenance_warning") if top_vert else None,
            "classification": top_vert["classification"] if top_vert else "NOT MODELLED",
            "unit": "engagement 0-1",
        },
        "receptor_selectivity": {
            "pair": pair_name,
            "ratio_vert_over_insect": (pair or {}).get("ratio"),
            "log10_ratio": (pair or {}).get("log10_ratio_vert_over_insect"),
            "insect_receptor": (pair or {}).get("insect_receptor"),
            "vertebrate_receptor": (pair or {}).get("vertebrate_receptor"),
            "classification": "MODEL-DERIVED" if pair else "NOT MODELLED",
            "unit": "fold (potency ratio)",
        },
        "evidence_tier": {
            "value": tier,
            "n_sourced": len(sourced),
            "n_rows": len(evidence),
            "n_not_modelled": sum(1 for r in evidence if r["classification"] == "NOT MODELLED"),
            "classification": "LITERATURE" if sourced else "NOT MODELLED",
            "unit": "library rows",
        },
    }


# -- block 2: selectivity, as facts ----------------------------------------
def _selectivity_block(
    occ: dict[str, Any], evidence: list[dict[str, Any]], conc_M: float, occ_limit: float = 0.2
) -> dict[str, Any]:
    pair_name, pair = _best_pair(occ)
    vert_rows = []
    for row in evidence:
        if row["organism"] != "vertebrate":
            continue
        threshold = _threshold_conc(row.get("param_value_M"), row.get("n"), occ_limit)
        vert_rows.append(
            {
                "receptor": row.get("receptor"),
                "engagement_at_dose": row.get("engagement"),
                "param_type": row.get("param_type"),
                "param_value_M": row.get("param_value_M"),
                "n": row.get("n"),
                "species": row.get("species"),
                "relation": row.get("relation"),
                "source": row.get("source"),
                "doi": row.get("doi"),
                "pmid": row.get("pmid"),
                "classification": row.get("classification"),
                "conc_at_limit_M": threshold,
            }
        )
    reachable = [r for r in vert_rows if r["conc_at_limit_M"] is not None]
    reachable.sort(key=lambda r: r["conc_at_limit_M"])
    limiting = reachable[0] if reachable else None
    insect_eng = (pair or {}).get("insect_engagement")
    vert_eng = (pair or {}).get("vertebrate_engagement")
    difference = None
    if insect_eng is not None and vert_eng is not None:
        difference = float(insect_eng) - float(vert_eng)
    return {
        "occ_limit": occ_limit,
        "pair": pair_name,
        "insect_receptor": (pair or {}).get("insect_receptor"),
        "vertebrate_receptor": (pair or {}).get("vertebrate_receptor"),
        "insect_engagement": insect_eng,
        "vertebrate_engagement": vert_eng,
        "ratio_vert_over_insect": (pair or {}).get("ratio"),
        "log10_ratio": (pair or {}).get("log10_ratio_vert_over_insect"),
        "engagement_difference": difference,
        "vertebrate_limit_conc_M": (limiting or {}).get("conc_at_limit_M"),
        "limiting_vertebrate_receptor": limiting,
        "vertebrate_rows": vert_rows,
        "dose_over_vertebrate_limit": (
            float(conc_M) / float(limiting["conc_at_limit_M"])
            if limiting and limiting["conc_at_limit_M"]
            else None
        ),
        "statement": (
            "These are values, not a verdict: a ratio is not a safety margin, and the "
            "vertebrate engagement printed beside it is the number that matters at this dose."
        ),
    }


# -- block 3: circuit consequence ------------------------------------------
def _circuit_block(
    compound: str | None,
    conc_M: float,
    graph: str | None,
    genotype: str | None = None,
    drive_hz: float = 40.0,
    steps: int = 80,
) -> dict[str, Any]:
    from flylab.assays.subgraph import run_subgraph_assay

    nb, ms = _timed(
        lambda: _with_genotype(
            run_subgraph_assay,
            genotype,
            compound,
            conc_M,
            drive_hz=drive_hz,
            steps=steps,
            graph=graph,
        )
    )
    nb = _stamp(nb, ms)
    readouts = nb.get("readouts") or {}
    vehicle = readouts.get("vehicle") or {}
    out: dict[str, Any] = {
        "assay": nb.get("assay"),
        "graph": graph or "named",
        "gains": nb.get("gains") or {},
        "n_nodes": readouts.get("n_nodes"),
        "n_edges": readouts.get("n_edges"),
        "drive_hz": readouts.get("drive_hz"),
        "readouts": {},
        "notebook": nb,
        "classification": "PREDICTION",
    }
    for key, label in (("mean_hz", "mean circuit rate"), ("mn9_hz", "MN9"), ("dnp01_hz", "DNp01")):
        treated = readouts.get(key)
        base = vehicle.get(key)
        try:
            t = float(treated)
            b = float(base)
        except (TypeError, ValueError):
            continue
        pct = ((t - b) / b * 100.0) if abs(b) > 1e-12 else None
        out["readouts"][key] = {
            "label": label,
            "vehicle": b,
            "treated": t,
            "delta": t - b,
            "percent": pct,
            "direction": "suppressed" if t < b else ("enhanced" if t > b else "unchanged"),
            "unit": "Hz",
        }
    mean = out["readouts"].get("mean_hz") or {}
    out["direction"] = mean.get("direction", "unchanged")
    out["percent"] = mean.get("percent")
    out["warnings"] = list(nb.get("warnings") or [])
    return out


# -- block 3b: is the effect wiring-dependent? ------------------------------
def _dependence_verdict(cls: str | None) -> str:
    """Plain-language verdict for a dependence class, used by both routes."""
    if cls == "topology-dependent":
        return "specific wiring evidence: present (topology-dependent)"
    if cls == "composition-dominated":
        return "specific wiring evidence: weak (composition-dominated)"
    if cls == "network-insensitive":
        return "specific wiring evidence: none (network-insensitive)"
    return f"specific wiring evidence: {cls}"


def _dependence_block(
    compound: str,
    conc_M: float,
    graph: str | None,
    n: int,
    assay: str = "subgraph",
    seed: int = 0,
) -> dict[str, Any]:
    fn = _lazy("flylab.analysis.dependence", "dependence_profile")
    kw: dict[str, Any] = {}
    if graph:
        kw["graph"] = graph
    profile = _invoke(
        fn,
        compound=compound,
        conc_M=conc_M,
        assay=assay,
        n=int(n),
        seed=int(seed),
        modes=list(DEPENDENCE_MODES),
        **kw,
    )
    cls = profile.get("class")
    structural = [
        m
        for m in profile.get("modes") or []
        if m.get("mode") in ("weight_permute", "rewire_degree_preserving") and m.get("beats_null")
    ]
    verdict = _dependence_verdict(cls)
    return {
        "class": cls,
        "verdict": verdict,
        "reason": (profile.get("classification") or {}).get("reason"),
        "n": profile.get("n"),
        "p_resolution": profile.get("p_resolution"),
        "real_effect": profile.get("real_effect"),
        "readout": profile.get("readout"),
        "structural_modes_beaten": [m.get("mode") for m in structural],
        "modes": [
            {
                "mode": m.get("mode"),
                "information_kept": m.get("information_kept"),
                "p_two_sided": m.get("p_two_sided"),
                "p_resolution": m.get("p_resolution"),
                "at_resolution_floor": bool(m.get("resolution_limited")),
                "beats_null": m.get("beats_null"),
                "z": m.get("z"),
                "n": m.get("n"),
            }
            for m in profile.get("modes") or []
        ],
        "necessary_information_level": profile.get("necessary_information_level"),
        "classification": "MODEL-DERIVED",
        "warnings": list(profile.get("warnings") or []),
        "note": (
            "The headline statistic is the empirical two-sided permutation p with its "
            f"resolution 1/(n+1); z is reported second and is never a probability."
        ),
    }


# -- block 4: one concentration axis ---------------------------------------
def _ladder(
    compound: str,
    assay: str,
    graph: str | None,
    readout: str,
    concs: list[float],
    threshold_frac: float = 0.5,
) -> dict[str, Any]:
    """Circuit response across the concentration ladder, on the rate engine.

    ``analysis.nullmodels.fast_drug_effect`` is the project's own engine path
    for exactly this contrast: same arithmetic as the notebook assay, ~80x
    faster because the vehicle arm is not re-run at every rung.  The crossing
    rule below is the one in ``analysis.selectivity.circuit_threshold_conc``,
    and ``tests/test_server.py`` asserts the two agree.
    """
    import math

    fn = _lazy("flylab.analysis.nullmodels", "fast_drug_effect")
    curve: list[dict[str, Any]] = []
    rel: list[float | None] = []
    vehicle_val: float | None = None
    for c in concs:
        res = _invoke(fn, assay=assay, compound=compound, conc_M=c, readout=readout, graph=graph)
        treated, vehicle = res.get("treated"), res.get("vehicle")
        if vehicle is not None:
            vehicle_val = vehicle
        r = None
        if treated is not None and vehicle is not None and abs(vehicle) > 1e-12:
            r = abs(treated - vehicle) / abs(vehicle)
        rel.append(r)
        curve.append(
            {
                "conc_M": c,
                "treated": treated,
                "vehicle": vehicle,
                "effect": res.get("effect"),
                "rel_change": r,
            }
        )

    hit_index = next((i for i, r in enumerate(rel) if r is not None and r >= threshold_frac), None)
    finite = [(r, c) for r, c in zip(rel, concs) if r is not None]
    max_rel, max_rel_conc = max(finite, default=(None, None))
    hit_conc: float | None = None
    if hit_index == 0:
        hit_conc = concs[0]
    elif hit_index is not None:
        lo_r, hi_r = rel[hit_index - 1], rel[hit_index]
        lo_x, hi_x = math.log10(concs[hit_index - 1]), math.log10(concs[hit_index])
        if lo_r is not None and hi_r is not None and hi_r > lo_r:
            frac = (threshold_frac - lo_r) / (hi_r - lo_r)
            hit_conc = float(10.0 ** (lo_x + frac * (hi_x - lo_x)))
        else:
            hit_conc = concs[hit_index]
    return {
        "curve": curve,
        "vehicle": vehicle_val,
        "circuit_threshold_M": hit_conc,
        "crossing_index": hit_index,
        "max_rel_change": max_rel,
        "max_rel_change_conc_M": max_rel_conc,
        "threshold_frac": threshold_frac,
        "readout": readout,
    }


def _ladder_block(
    compound: str,
    conc_M: float,
    evidence: list[dict[str, Any]],
    graph: str | None,
    assay: str = "subgraph",
    occ: dict[str, Any] | None = None,
    occ_limit: float = 0.2,
) -> dict[str, Any]:
    import math

    from flylab.analysis.selectivity import DEFAULT_CONCS

    readout = "mean_hz" if assay == "subgraph" else "mn9_hz"
    concs = [float(c) for c in DEFAULT_CONCS]
    lad = _ladder(compound, assay, graph or "named", readout, concs)
    curve = lad["curve"]
    insect, vert = _split_rows(evidence)

    def _series(rows):
        out = []
        for c in concs:
            values = [
                _hill(c, r.get("param_value_M"), r.get("n"))
                for r in rows
                if r.get("param_value_M") is not None
            ]
            values = [v for v in values if v is not None]
            out.append(max(values) if values else None)
        return out

    vt = _invoke(
        _lazy("flylab.analysis.selectivity", "vertebrate_threshold_conc"),
        compound=compound,
        occ_limit=occ_limit,
    )
    c_circ, c_vert = lad["circuit_threshold_M"], vt.get("conc_M")
    circuit_si = (
        float(math.log10(c_vert) - math.log10(c_circ)) if (c_circ and c_vert) else None
    )
    _pair_name, pair = _best_pair(occ or {})
    receptor_si = (pair or {}).get("log10_ratio_vert_over_insect")
    gap = None
    if circuit_si is not None and receptor_si not in (None, float("-inf")):
        gap = float(circuit_si - float(receptor_si))

    per_receptor = [
        {
            "receptor": r.get("receptor"),
            "organism": r.get("organism"),
            "classification": r.get("classification"),
            "param_type": r.get("param_type"),
            "param_value_M": r.get("param_value_M"),
            "values": [_hill(c, r.get("param_value_M"), r.get("n")) for c in concs],
        }
        for r in insect + vert
    ]
    return {
        "concs_M": concs,
        "insect_engagement": _series(insect),
        "vertebrate_engagement": _series(vert),
        "circuit_response": [p.get("rel_change") for p in curve],
        "circuit_treated_hz": [p.get("treated") for p in curve],
        "circuit_vehicle_hz": lad["vehicle"],
        "readout": readout,
        "per_receptor": per_receptor,
        "circuit_threshold_M": c_circ,
        "threshold_frac": lad["threshold_frac"],
        "vertebrate_threshold_M": c_vert,
        "vertebrate_receptor": vt.get("receptor"),
        "vertebrate_threshold_placeholder": bool(vt.get("placeholder")),
        "circuit_si_log10": circuit_si,
        "receptor_si_log10": receptor_si,
        "si_gap_circuit_minus_receptor": gap,
        "max_rel_change": lad["max_rel_change"],
        "max_rel_change_conc_M": lad["max_rel_change_conc_M"],
        "current_conc_M": conc_M,
        "axis_note": (
            "All three series are dimensionless fractions on one axis: engagement is 0-1 "
            "and the circuit response is |treated - vehicle| / vehicle. No second y-axis."
        ),
        "classification": "MODEL-DERIVED",
        "warnings": [
            "The circuit response is a relative change of a simulated rate, not a measured "
            "dose-response.",
            "A positive circuit selectivity index means only that this simulation moves "
            "before the teaching library's vertebrate receptor fills; it is not a safety margin.",
        ],
    }


# -- block 8: what can I trust? --------------------------------------------
def _trust_block(
    evidence: list[dict[str, Any]], circuit: dict[str, Any], graph: str | None
) -> dict[str, Any]:
    sourced = [r for r in evidence if r["classification"] == "LITERATURE"]
    placeholder = [r for r in evidence if r["classification"] == "NOT MODELLED"]
    coverage, mn9_gap = _expression_coverage()
    rows = [
        {
            "area": "Receptor values",
            "status": f"{len(sourced)} of {len(evidence)} rows literature-supported",
            "classification": "LITERATURE" if sourced else "NOT MODELLED",
            "basis": (
                "Every modelled row cites a paper and carries its parameter type; "
                f"{len(placeholder)} row(s) are class placeholders and report nothing."
            ),
            "detail": {
                "n_sourced": len(sourced),
                "n_placeholder": len(placeholder),
                "receptors_not_modelled": [r["receptor"] for r in placeholder],
            },
        },
        {
            "area": "Mechanism rule (engagement to gain)",
            "status": "asserted, never fitted",
            "classification": "MODEL-ASSUMPTION",
            "basis": (
                "The occupancy-to-gain transformation is the project's central modelling "
                "choice. No experiment in the literature measures it for these receptors."
            ),
            "detail": {"gains": circuit.get("gains") or {}, "gain_floor": 0.05},
        },
        {
            "area": "Exposure prediction",
            "status": "low confidence",
            "classification": "MODEL-ASSUMPTION",
            "basis": (
                "One-compartment first-order model. Nothing connects an applied dose to the "
                "free concentration at the receptor in this species."
            ),
            "detail": {"model": "one-compartment C(t)", "fitted_to_animal_data": False},
        },
        {
            "area": "PK constants",
            "status": "placeholder",
            "classification": "MODEL-ASSUMPTION",
            "basis": "Absorption and elimination constants are teaching defaults per route.",
            "detail": {"source": "flylab/pharm/exposure.py ROUTE_DEFAULTS"},
        },
        {
            "area": "Circuit topology",
            "status": "real MaleCNS v1.0",
            "classification": "LITERATURE",
            "basis": (
                f"{circuit.get('n_nodes')} cells and {circuit.get('n_edges')} edges from the "
                "public reconstruction, hops-limited with a synapse-count floor."
            ),
            "detail": {
                "graph": graph or "named",
                "n_nodes": circuit.get("n_nodes"),
                "n_edges": circuit.get("n_edges"),
            },
        },
        {
            "area": "Transmitter identity",
            "status": "predicted",
            "classification": "PREDICTION",
            "basis": (
                "Edge signs follow MaleCNS predicted consensus transmitters, not staining. "
                "A wrong label flips a synapse's sign."
            ),
            "detail": {"source": "MaleCNS v1.0 predicted consensus neurotransmitter"},
        },
        {
            "area": "Receptor expression",
            "status": (
                f"{coverage * 100:.1f}% of cells mapped" if coverage is not None else "unmapped"
            ),
            "classification": "MODEL-ASSUMPTION",
            "basis": (
                "Gains are applied uniformly. Adult motor-neuron receptor expression -- the "
                "class MN9 belongs to -- is a confirmed gap in the literature."
            ),
            "detail": {"fraction_cells_mapped": coverage, "motor_neuron_gap": mn9_gap},
        },
        {
            "area": "Live validation",
            "status": "none",
            "classification": "NOT MODELLED",
            "basis": (
                "live_lab is null. No FlyLab code path may write a live-animal number; the "
                "slot only holds a table a human imported."
            ),
            "detail": {"live_lab": (circuit.get("notebook") or {}).get("live_lab")},
        },
    ]
    return {
        "rows": rows,
        "blended_confidence": None,
        "note": (
            "There is deliberately no single confidence percentage: these layers fail "
            "independently, and one number would hide which of them is weak."
        ),
    }


def _expression_coverage() -> tuple[float | None, dict[str, Any] | None]:
    try:
        from flylab.pharm.expression import expression_table

        table = expression_table()
    except Exception:  # pragma: no cover - optional module
        return None, None
    coverage = (table.get("coverage") or {}).get("fraction_known_overall")
    gap = table.get("motor_neuron_gap")
    slim = None
    if isinstance(gap, dict):
        slim = {"status": gap.get("status"), "description": gap.get("description")}
    try:
        coverage = float(coverage)
    except (TypeError, ValueError):
        coverage = None
    return coverage, slim


# -- block 10: why this happened -------------------------------------------
def _why_block(
    occ: dict[str, Any],
    evidence: list[dict[str, Any]],
    circuit: dict[str, Any],
    dependence: dict[str, Any] | None,
    conc_M: float,
) -> dict[str, Any]:
    """Generated from this run: saturation state, dependence class, dominant term.

    Never a per-compound string: every clause below is selected by a number the
    run produced.
    """
    insect, _vert = _split_rows(evidence)
    top = max(insect, key=lambda r: r["engagement"]) if insect else None
    theta = top["engagement"] if top else None
    receptor = (top or {}).get("receptor")
    param_type = (top or {}).get("param_type") or "potency"
    value = (top or {}).get("param_value_M")
    n = (top or {}).get("n") or 1.0
    # how much a two-fold error in the cited parameter would move engagement
    span = None
    if value:
        lo = _hill(conc_M, float(value) * 2.0, n)
        hi = _hill(conc_M, float(value) / 2.0, n)
        if lo is not None and hi is not None:
            span = abs(hi - lo)

    if theta is None:
        sat = "not modelled"
    elif theta >= 0.95:
        sat = "essentially saturated"
    elif theta >= 0.5:
        sat = "past half-maximal"
    elif theta >= 0.05:
        sat = "only partly engaged"
    else:
        sat = "barely engaged"

    potency_inert = span is not None and span < 0.05
    dominant = (
        "the occupancy-to-gain rule"
        if potency_inert
        else f"the cited {param_type} and the occupancy-to-gain rule together"
    )

    cls = (dependence or {}).get("class")
    if cls == "topology-dependent":
        structure = "the specific MaleCNS wiring pattern, which the permutation nulls do not reproduce"
        excluded = "global transmitter composition alone"
    elif cls == "composition-dominated":
        structure = "the global transmitter composition of the cut"
        excluded = ("the exact MaleCNS wiring: no structure-preserving shuffle of this cut could be told apart from it, which is weaker than saying the shuffles reproduce the effect")
    elif cls == "network-insensitive":
        structure = "the pharmacology alone"
        excluded = "the graph, which no null model distinguishes here"
    else:
        structure = "the circuit as configured"
        excluded = "any single layer"

    mean = (circuit.get("readouts") or {}).get("mean_hz") or {}
    direction = mean.get("direction", "change")
    pct = mean.get("percent")
    gains = circuit.get("gains") or {}
    moved = sorted(
        ((k, v) for k, v in gains.items() if isinstance(v, (int, float))),
        key=lambda kv: -abs(float(kv[1]) - 1.0),
    )
    driver = moved[0][0] if moved and abs(float(moved[0][1]) - 1.0) > 1e-9 else None

    conc_text = f"{conc_M:.2g} M" if conc_M else "this dose"
    parts = []
    if receptor:
        parts.append(
            f"At {conc_text} the insect {str(receptor).replace('insect_', '')} receptor is {sat} "
            + (f"(engagement {theta:.2f}" if theta is not None else "(")
            + (
                f"; a two-fold error in the cited {param_type} would move it by {span:.3f})."
                if span is not None
                else ")."
            )
        )
    else:
        parts.append(f"At {conc_text} no insect receptor in this row carries a sourced value.")
    if driver:
        parts.append(
            f"The mechanism rule turns that into {driver} = {float(gains[driver]):.3f}, "
            f"and the circuit is {direction}"
            + (f" by {abs(pct):.1f}%" if pct is not None else "")
            + "."
        )
    else:
        parts.append("No mechanism rule fires at this dose, so the circuit is unchanged.")
    parts.append(
        f"The {direction} response is therefore governed primarily by {dominant} and by "
        f"{structure}, rather than by {excluded}."
    )
    return {
        "text": " ".join(parts),
        "inputs": {
            "insect_receptor": receptor,
            "insect_engagement": theta,
            "saturation": sat,
            "param_type": param_type,
            "two_fold_engagement_span": span,
            "potency_is_inert": potency_inert,
            "dependence_class": cls,
            "dominant_uncertainty": dominant,
            "dominant_gain": driver,
            "direction": direction,
            "percent_change": pct,
        },
        "classification": "MODEL-DERIVED",
        "note": "Generated from this run, not from a stored sentence per compound.",
    }


# -- assembly ---------------------------------------------------------------
def dashboard_runtime_estimate(include_dependence: bool, n: int, graph: str | None) -> dict[str, Any]:
    """Seconds this dashboard call is expected to take, before it runs."""
    base = 0.9  # occupancy + subgraph assay + the concentration ladder
    depend = 0.0
    detail: dict[str, Any] = {"assembly_s": base}
    if include_dependence:
        try:
            fn = _lazy("flylab.analysis.dependence", "estimate_landscape_runtime")
            est = _invoke(fn, n_cells=1, n=int(n), graph=graph or "named")
            depend = float(est.get("estimate_s") or 0.0)
            detail["dependence"] = est
        except BridgeError:
            depend = 0.0
    return {
        "estimate_s": round(base + depend, 2),
        "n": int(n),
        "include_dependence": bool(include_dependence),
        "detail": detail,
        "note": "Order-of-magnitude, from the project's own per-shuffle cost model.",
    }


def build_dashboard(
    compound: str,
    conc_M: float = 1e-6,
    *,
    graph: str | None = "named",
    assay: str = "subgraph",
    genotype: str | None = None,
    n: int = FAST_DEPENDENCE_N,
    seed: int = 0,
    include_dependence: bool = True,
    include_ladder: bool = True,
    estimate_only: bool = False,
) -> dict[str, Any]:
    """Everything the landing page needs, in one call."""
    from flylab.pharm.occupancy import compare_compound, library_report, load_library

    estimate = dashboard_runtime_estimate(include_dependence, n, graph)
    if estimate_only:
        return {
            "compound": compound,
            "concentration_M": conc_M,
            "runtime_estimate": estimate,
            "estimate_only": True,
        }

    t0 = time.perf_counter()
    lib = load_library()
    key = str(compound).lower().strip()
    spec = (lib.get("compounds") or {}).get(key)
    if spec is None:
        raise BridgeError(404, f"unknown compound {compound!r}")
    occ = compare_compound(key, conc_M)
    evidence = _evidence_rows(occ)
    headline = _headline(occ, evidence)
    selectivity = _selectivity_block(occ, evidence, conc_M)
    circuit = _circuit_block(key, conc_M, graph, genotype)
    dependence = (
        _dependence_block(key, conc_M, graph, n, assay=assay, seed=seed)
        if include_dependence
        else None
    )
    ladder = (
        _ladder_block(key, conc_M, evidence, graph, assay=assay, occ=occ)
        if include_ladder
        else None
    )
    trust = _trust_block(evidence, circuit, graph)
    why = _why_block(occ, evidence, circuit, dependence, conc_M)

    from flylab.analysis.claims import claim_audit

    notebook = circuit.get("notebook") or {}
    claims = claim_audit(notebook)

    insect_targets = [
        {
            "receptor": r.get("receptor"),
            "direction": r.get("direction"),
            "param_type": r.get("param_type"),
            "param_value_M": r.get("param_value_M"),
            "classification": r.get("classification"),
        }
        for r in evidence
        if r["organism"] == "insect" and r.get("param_value_M") is not None
    ]
    insect_targets.sort(key=lambda r: r["param_value_M"])
    primary = insect_targets[0] if insect_targets else None

    warnings = list(DASHBOARD_WARNINGS)
    for source in (occ.get("notes"), circuit.get("warnings"), (ladder or {}).get("warnings")):
        for w in source or []:
            if isinstance(w, str) and w not in warnings:
                warnings.append(w)

    return {
        "compound": {
            "key": key,
            "name": spec.get("name", key),
            "class": spec.get("class"),
            "cas": spec.get("cas"),
            "target_receptor": (primary or {}).get("receptor"),
            "mode": (primary or {}).get("direction"),
            "insect_targets": insect_targets,
        },
        "concentration_M": conc_M,
        "graph": graph or "named",
        "assay": assay,
        "genotype": genotype,
        "headline": headline,
        "coverage": {
            "n_rows": len(evidence),
            "n_sourced": sum(1 for r in evidence if r["classification"] == "LITERATURE"),
            "n_not_modelled": sum(1 for r in evidence if r["classification"] == "NOT MODELLED"),
            "library": library_report(lib),
        },
        "occupancy": occ,
        "evidence": evidence,
        "selectivity": selectivity,
        "circuit": circuit,
        "dependence": dependence,
        "ladder": ladder,
        "trust": trust,
        "why": why,
        "claims": claims,
        "runtime_estimate": estimate,
        "runtime_s": round(time.perf_counter() - t0, 3),
        "label": "model_derived",
        "warnings": warnings,
        "disclaimer": occ.get("disclaimer"),
    }


def compare_runtime_estimate(
    compounds: list[str], include_dependence: bool, n: int, graph: str | None
) -> dict[str, Any]:
    base = COMPARE_SECONDS_PER_COMPOUND * max(1, len(compounds))
    depend = 0.0
    detail: dict[str, Any] = {"assembly_s": round(base, 2)}
    if include_dependence:
        try:
            fn = _lazy("flylab.analysis.dependence", "estimate_landscape_runtime")
            est = _invoke(fn, n_cells=len(compounds), n=int(n), graph=graph or "named")
            depend = float(est.get("estimate_s") or 0.0)
            detail["dependence"] = est
        except BridgeError:
            depend = 0.0
    return {
        "estimate_s": round(base + depend, 2),
        "n_compounds": len(compounds),
        "n": int(n),
        "include_dependence": bool(include_dependence),
        "detail": detail,
        "note": "Order-of-magnitude, from the project's own per-shuffle cost model.",
    }


def build_compare(
    compounds: list[str],
    conc_M: float = 1e-6,
    *,
    assay: str = "subgraph",
    graph: str | None = "named",
    n: int = FAST_DEPENDENCE_N,
    seed: int = 0,
    include_dependence: bool = True,
    estimate_only: bool = False,
) -> dict[str, Any]:
    """The decision table: receptor selectivity beside circuit selectivity."""
    from flylab.pharm.occupancy import compare_compound, load_library

    keys = [str(c).lower().strip() for c in compounds if str(c).strip()]
    if not keys:
        raise BridgeError(400, "compare needs at least one compound")
    if len(keys) > 8:
        raise BridgeError(400, "compare is limited to 8 compounds")
    estimate = compare_runtime_estimate(keys, include_dependence, n, graph)
    if estimate_only:
        return {"compounds": keys, "runtime_estimate": estimate, "estimate_only": True}

    lib = load_library()
    t0 = time.perf_counter()
    rows = []
    for key in keys:
        spec = (lib.get("compounds") or {}).get(key)
        if spec is None:
            raise BridgeError(404, f"unknown compound {key!r}")
        occ = compare_compound(key, conc_M)
        evidence = _evidence_rows(occ)
        insect, vert = _split_rows(evidence)
        top_insect = max(insect, key=lambda r: r["engagement"]) if insect else None
        top_vert = max(vert, key=lambda r: r["engagement"]) if vert else None
        csi = _ladder_block(key, conc_M, evidence, graph, assay=assay, occ=occ)
        circuit = _circuit_block(key, conc_M, graph)
        mean = (circuit.get("readouts") or {}).get("mean_hz") or {}
        dep = (
            _dependence_block(key, conc_M, graph, n, assay=assay, seed=seed)
            if include_dependence
            else None
        )
        targets = [r for r in evidence if r["organism"] == "insect" and r.get("param_value_M")]
        targets.sort(key=lambda r: r["param_value_M"])
        sourced = [r for r in evidence if r["classification"] == "LITERATURE"]
        rows.append(
            {
                "compound": key,
                "name": spec.get("name", key),
                "class": spec.get("class"),
                "target_receptor": (targets[0] if targets else {}).get("receptor"),
                "mode": (targets[0] if targets else {}).get("direction"),
                "insect_engagement": (top_insect or {}).get("engagement"),
                "insect_receptor": (top_insect or {}).get("receptor"),
                "vertebrate_engagement": (top_vert or {}).get("engagement"),
                "vertebrate_receptor": (top_vert or {}).get("receptor"),
                "receptor_si_log10": csi.get("receptor_si_log10"),
                "circuit_si_log10": csi.get("circuit_si_log10"),
                "si_gap_circuit_minus_receptor": csi.get("si_gap_circuit_minus_receptor"),
                "circuit_threshold_M": csi.get("circuit_threshold_M"),
                "vertebrate_threshold_M": csi.get("vertebrate_threshold_M"),
                "circuit_delta_hz": mean.get("delta"),
                "circuit_delta_percent": mean.get("percent"),
                "circuit_direction": mean.get("direction"),
                "topology_dependence": (dep or {}).get("class"),
                "topology_verdict": (dep or {}).get("verdict"),
                "topology_p": (dep or {}).get("modes"),
                "evidence_tier": "literature_order" if sourced else "class_placeholder",
                "n_sourced": len(sourced),
                "n_rows": len(evidence),
                "n_not_modelled": sum(1 for r in evidence if r["classification"] == "NOT MODELLED"),
                "vertebrate_threshold_placeholder": csi.get("vertebrate_threshold_placeholder"),
                "max_rel_change": csi.get("max_rel_change"),
                "ladder": csi,
            }
        )

    scored = [
        r
        for r in rows
        if r.get("receptor_si_log10") is not None and r.get("circuit_si_log10") is not None
    ]
    best_receptor = max(scored, key=lambda r: r["receptor_si_log10"]) if scored else None
    best_circuit = max(scored, key=lambda r: r["circuit_si_log10"]) if scored else None
    findings = []
    if best_receptor and best_circuit:
        if best_receptor["compound"] != best_circuit["compound"]:
            findings.append(
                f"The highest receptor selectivity ({best_receptor['name']}, "
                f"{best_receptor['receptor_si_log10']:.2f} log10) is NOT the highest circuit "
                f"selectivity ({best_circuit['name']}, {best_circuit['circuit_si_log10']:.2f} log10)."
            )
        else:
            findings.append(
                f"{best_receptor['name']} leads on both receptor and circuit selectivity in this set."
            )
    return {
        "concentration_M": conc_M,
        "assay": assay,
        "graph": graph or "named",
        "rows": rows,
        "findings": findings,
        "best_receptor_si": (best_receptor or {}).get("compound"),
        "best_circuit_si": (best_circuit or {}).get("compound"),
        "runtime_estimate": estimate,
        "runtime_s": round(time.perf_counter() - t0, 3),
        "label": "model_derived",
        "columns": [
            "compound",
            "target_receptor",
            "insect_engagement",
            "vertebrate_engagement",
            "receptor_si_log10",
            "circuit_si_log10",
            "si_gap_circuit_minus_receptor",
            "circuit_delta_percent",
            "topology_dependence",
            "evidence_tier",
        ],
        "warnings": [
            "Receptor selectivity and circuit selectivity are different quantities; a large "
            "receptor ratio does not buy a wide circuit window.",
            "Eight compounds in the library have no circuit selectivity index at all: RDL / "
            "GluCl block cannot reach the effect threshold on these cuts.",
            "Topology dependence is a documented label from permutation nulls at the stated "
            "n, not a hypothesis test with multiplicity control.",
        ],
    }


def build_claims(
    compound: str | None = None,
    conc_M: float = 1e-6,
    *,
    notebook: dict[str, Any] | None = None,
    assay: str = "subgraph",
    graph: str | None = "named",
    run_assay: bool = True,
) -> dict[str, Any]:
    """Claim provenance for a supplied notebook, or for a fresh run."""
    from flylab.analysis.claims import claim_audit

    payload = notebook
    if payload is None:
        if not compound:
            raise BridgeError(400, "claims needs a compound or a notebook")
        if run_assay:
            payload = _circuit_block(str(compound).lower().strip(), conc_M, graph).get("notebook")
        else:
            payload = {
                "compound": str(compound).lower().strip(),
                "concentration_M": conc_M,
                "assay": assay,
            }
    return claim_audit(payload)


# --------------------------------------------------------------------------
# dashboard / analysis routes
# --------------------------------------------------------------------------
def _h_dashboard(p: dict[str, Any]) -> dict[str, Any]:
    compound = _opt_str(p, "compound")
    if compound is None:
        raise BridgeError(400, "compound is required")
    return build_dashboard(
        compound,
        _num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M),
        graph=_choice(p, "graph", "named", GRAPH_NAMES),
        assay=_choice(p, "assay", "subgraph", LADDER_ASSAYS) or "subgraph",
        genotype=_opt_str(p, "genotype"),
        n=_int(p, "n", FAST_DEPENDENCE_N, lo=1, hi=500),
        seed=_int(p, "seed", 0, lo=0, hi=2**31 - 1),
        include_dependence=_bool(p, "include_dependence", True),
        include_ladder=_bool(p, "include_ladder", True),
        estimate_only=_bool(p, "estimate_only", False),
    )


def _h_compare(p: dict[str, Any]) -> dict[str, Any]:
    compounds = _strs(p, "compounds") or []
    return build_compare(
        compounds,
        _num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M),
        assay=_choice(p, "assay", "subgraph", LADDER_ASSAYS) or "subgraph",
        graph=_choice(p, "graph", "named", GRAPH_NAMES),
        n=_int(p, "n", FAST_DEPENDENCE_N, lo=1, hi=500),
        seed=_int(p, "seed", 0, lo=0, hi=2**31 - 1),
        include_dependence=_bool(p, "include_dependence", True),
        estimate_only=_bool(p, "estimate_only", False),
    )


def _h_claims(p: dict[str, Any]) -> dict[str, Any]:
    notebook = p.get("notebook")
    return build_claims(
        _opt_str(p, "compound"),
        _num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M),
        notebook=notebook if isinstance(notebook, dict) and notebook else None,
        assay=_choice(p, "assay", "subgraph", ASSAY_NAMES) or "subgraph",
        graph=_choice(p, "graph", "named", GRAPH_NAMES),
        run_assay=_bool(p, "run_assay", True),
    )


def _dependence_estimate(n_cells: int, n: int, graph: str | None) -> dict[str, Any]:
    fn = _lazy("flylab.analysis.dependence", "estimate_landscape_runtime")
    return _invoke(fn, n_cells=int(n_cells), n=int(n), graph=graph or "named")


def _h_dependence(p: dict[str, Any]) -> dict[str, Any]:
    compound = _str(p, "compound", "imidacloprid")
    conc = _num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M)
    graph = _choice(p, "graph", None, GRAPH_NAMES)
    n = _int(p, "n", FAST_DEPENDENCE_N, lo=1, hi=1000)
    estimate = _dependence_estimate(1, n, graph)
    if _bool(p, "estimate_only", False):
        return {"compound": compound, "runtime_estimate": estimate, "estimate_only": True}
    fn = _lazy("flylab.analysis.dependence", "dependence_profile")
    kw: dict[str, Any] = {}
    if graph:
        kw["graph"] = graph
    modes = _strs(p, "modes")
    out = _invoke(
        fn,
        compound=compound,
        conc_M=conc,
        assay=_choice(p, "assay", "subgraph", ASSAY_NAMES),
        readout=_str(p, "readout", "auto"),
        n=n,
        seed=_int(p, "seed", 0, lo=0, hi=2**31 - 1),
        modes=list(modes) if modes else list(DEPENDENCE_MODES),
        **kw,
    )
    out["runtime_estimate"] = estimate
    out["verdict"] = _dependence_verdict(out.get("class"))
    return out


def _h_dependence_landscape(p: dict[str, Any]) -> dict[str, Any]:
    compounds = _strs(p, "compounds")
    concs = _floats(p, "concs_M") or [1e-8, 1e-7, 1e-6, 1e-5]
    graph = _choice(p, "graph", None, GRAPH_NAMES)
    n = _int(p, "n", FAST_DEPENDENCE_N, lo=1, hi=1000)
    if compounds is None:
        from flylab.pharm.occupancy import list_compounds, load_library

        compounds = list(list_compounds(load_library()))
    if len(compounds) * len(concs) > 200:
        raise BridgeError(400, "landscape is limited to 200 cells; narrow compounds or concs_M")
    estimate = _dependence_estimate(len(compounds) * len(concs), n, graph)
    if _bool(p, "estimate_only", True):
        # a full landscape is minutes of compute: the default answer is the
        # estimate, and the caller has to ask again to actually run it.
        return {
            "compounds": compounds,
            "concs_M": concs,
            "n": n,
            "runtime_estimate": estimate,
            "estimate_only": True,
            "note": "Pass estimate_only=false to run it.",
        }
    fn = _lazy("flylab.analysis.dependence", "dependence_landscape")
    kw: dict[str, Any] = {}
    if graph:
        kw["graph"] = graph
    out = _invoke(
        fn,
        compounds=compounds,
        concs_M=concs,
        assay=_choice(p, "assay", "subgraph", ASSAY_NAMES),
        n=n,
        seed=_int(p, "seed", 0, lo=0, hi=2**31 - 1),
        **kw,
    )
    out["runtime_estimate"] = estimate
    return out


def _h_ablation(p: dict[str, Any]) -> dict[str, Any]:
    compounds = _strs(p, "compounds")
    concs = _floats(p, "concs_M") or [1e-6]
    table = _bool(p, "table", False)
    # one ablation already scores the whole library once (that is what the
    # information-gain block needs); the table repeats it per concentration.
    estimate = {
        "estimate_s": round(1.5 * len(concs), 2) if table else 1.5,
        "table": table,
        "note": (
            "The ablation ladder scores every library compound at all four levels once; "
            "the table repeats that at each concentration."
        ),
    }
    if _bool(p, "estimate_only", False):
        return {"runtime_estimate": estimate, "estimate_only": True}
    if table:
        fn = _lazy("flylab.analysis.baselines", "ablation_table")
        kw: dict[str, Any] = {
            "concs_M": concs,
            "graph": _choice(p, "graph", "named", GRAPH_NAMES),
            "assay": _choice(p, "assay", "subgraph", ASSAY_NAMES),
        }
        if compounds:
            kw["compounds"] = compounds
        out = _invoke(fn, **kw)
    else:
        fn = _lazy("flylab.analysis.baselines", "ablation")
        out = _invoke(
            fn,
            compound=_str(p, "compound", "imidacloprid"),
            conc_M=_num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M),
            graph=_choice(p, "graph", "named", GRAPH_NAMES),
            assay=_choice(p, "assay", "subgraph", ASSAY_NAMES),
            include_information_gain=_bool(p, "include_information_gain", True),
        )
    out["runtime_estimate"] = estimate
    return out


def _h_robustness_stability(p: dict[str, Any]) -> dict[str, Any]:
    fast = _bool(p, "fast", True)
    estimate = {
        "estimate_s": 20.0 if fast else 240.0,
        "fast": fast,
        "note": (
            "Re-derives every prospective conclusion under every admissible mechanism "
            "specification. The browser build defaults to the fast family."
        ),
    }
    if _bool(p, "estimate_only", False):
        return {"runtime_estimate": estimate, "estimate_only": True}
    fn = _lazy("flylab.analysis.robustness", "conclusion_stability")
    out = _invoke(
        fn,
        fast=fast,
        seed=_int(p, "seed", 0, lo=0, hi=2**31 - 1),
        conc_M=_num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M),
        n_jobs=1,
    )
    out["runtime_estimate"] = estimate
    return out


def _h_robustness_thresholds(p: dict[str, Any]) -> dict[str, Any]:
    fracs = _floats(p, "circuit_fracs") or [0.25, 0.4, 0.5]
    limits = _floats(p, "vert_limits") or [0.1, 0.2, 0.3]
    compounds = _strs(p, "compounds")
    estimate = {
        "estimate_s": round(8.0 * len(compounds or range(10)), 1),
        "note": "Recomputes the amplify / buffer split on every cell of the threshold grid.",
    }
    if _bool(p, "estimate_only", False):
        return {"runtime_estimate": estimate, "estimate_only": True}
    fn = _lazy("flylab.analysis.robustness", "threshold_sensitivity")
    kw: dict[str, Any] = {
        "circuit_fracs": fracs,
        "vert_limits": limits,
        "assay": _choice(p, "assay", "subgraph", ASSAY_NAMES),
    }
    graph = _choice(p, "graph", None, GRAPH_NAMES)
    if graph:
        kw["graph"] = graph
    if compounds:
        kw["compounds"] = compounds
    out = _invoke(fn, **kw)
    out["runtime_estimate"] = estimate
    return out


def _h_uncertainty_global(p: dict[str, Any]) -> dict[str, Any]:
    n_base = _int(p, "n_base", 32, lo=8, hi=1024)
    estimate = {
        # measured: ~0.21 s per base sample once the rate surrogate exists, plus a
        # one-time ~40 s surrogate build the first time a process asks for one.
        "estimate_s": round(n_base * 0.21 + 40.0, 1),
        "estimate_s_warm": round(n_base * 0.21, 1),
        "surrogate_build_s": 40.0,
        "n_base": n_base,
        "n_evaluations": n_base * 11,
        "note": (
            "Saltelli cross-sampling over 9 factors; cost is linear in n_base on top of a "
            "one-time rate-surrogate build. The browser build defaults to n_base=32, the "
            "paper uses 128 or more."
        ),
    }
    if _bool(p, "estimate_only", False):
        return {"runtime_estimate": estimate, "estimate_only": True}
    fn = _lazy("flylab.analysis.uncertainty_global", "sobol_analysis")
    out = _invoke(
        fn,
        compound=_str(p, "compound", "imidacloprid"),
        conc_M=_num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M),
        readout=_str(p, "readout", "mean_hz"),
        n_base=n_base,
        n_boot=_int(p, "n_boot", 50, lo=0, hi=2000),
        seed=_int(p, "seed", 0, lo=0, hi=2**31 - 1),
        graph=_choice(p, "graph", "named", GRAPH_NAMES) or "named",
    )
    budget = _lazy("flylab.analysis.uncertainty_global", "uncertainty_budget")
    out["budget"] = _invoke(budget, result=out)
    out["runtime_estimate"] = estimate
    return out


def _h_voi(p: dict[str, Any]) -> dict[str, Any]:
    n_base = _int(p, "n_base", 32, lo=8, hi=1024)
    estimate = {
        "estimate_s": round(n_base * 0.23 + 40.0, 1),
        "estimate_s_warm": round(n_base * 0.23, 1),
        "surrogate_build_s": 40.0,
        "n_base": n_base,
        "note": (
            "Value of information reuses one Sobol run, so it costs the same as the global "
            "uncertainty route at the same n_base, including the one-time surrogate build."
        ),
    }
    if _bool(p, "estimate_only", False):
        return {"runtime_estimate": estimate, "estimate_only": True}
    fn = _lazy("flylab.analysis.voi", "value_of_information")
    factors = _strs(p, "factors")
    kw: dict[str, Any] = {}
    if factors:
        kw["factors"] = factors
    out = _invoke(
        fn,
        compound=_str(p, "compound", "imidacloprid"),
        conc_M=_num(p, "conc_M", 1e-6, lo=0, hi=MAX_CONC_M),
        readout=_str(p, "readout", "mean_hz"),
        n_base=n_base,
        seed=_int(p, "seed", 0, lo=0, hi=2**31 - 1),
        **kw,
    )
    out["runtime_estimate"] = estimate
    return out


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
    "/api/genotype/panel": ("POST", _h_genotype_panel),
    "/api/mixture/isobologram": ("POST", _h_isobologram),
    "/api/expression": ("GET", _h_expression),
    "/api/validation": ("GET", _h_validation),
    "/api/dashboard": ("GET", _h_dashboard),
    "/api/compare": ("POST", _h_compare),
    "/api/dependence": ("POST", _h_dependence),
    "/api/dependence/landscape": ("POST", _h_dependence_landscape),
    "/api/ablation": ("POST", _h_ablation),
    "/api/robustness/stability": ("POST", _h_robustness_stability),
    "/api/robustness/thresholds": ("POST", _h_robustness_thresholds),
    "/api/uncertainty/global": ("POST", _h_uncertainty_global),
    "/api/voi": ("POST", _h_voi),
    "/api/claims": ("POST", _h_claims),
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


def handle(route: str, payload: dict[str, Any] | None = None) -> Any:
    """Answer one request and let errors *raise* (:class:`BridgeError`).

    ``flylab/server.py`` delegates its newer routes here after validating the
    request with pydantic, so the served bench and the browser build cannot
    drift: they run the same function on the same arguments.  The FastAPI layer
    maps :class:`BridgeError` onto its own status codes.
    """
    _prepare()
    path, query = _split(route)
    handler, path_params = _resolve(path)
    data: dict[str, Any] = dict(query)
    if isinstance(payload, dict):
        data.update(payload)
    data.update(path_params)
    return handler(data)


def call(route: str, payload: dict[str, Any] | None = None) -> Any:
    """Answer one bench request.

    ``route`` is the server path (``"/api/assay/subgraph"``), optionally with a
    query string.  ``payload`` carries the body of a POST route or the query
    parameters of a GET route.  Returns the same JSON-able object the HTTP
    route returns, or ``{"error": {"status": ..., "detail": ...}}`` -- it never
    raises for a bad request, a missing file or an unknown compound.
    """
    try:
        return handle(route, payload)
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
