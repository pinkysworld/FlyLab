"""Declarative experiment specifications: one YAML file describes a whole run.

A researcher who wants to reuse FlyLab should not have to learn its Python
API.  A spec file names the compounds, the doses, the graph, the readouts, the
engines and which of the analysis layers to run, and :func:`run_spec` turns it
into a *self-describing run directory*: a manifest with every hash and seed, the
resolved spec (defaults filled in), one notebook JSON per assay run, the tabular
readouts, the requested analysis payloads, one claim card per headline result,
and figures when matplotlib happens to be installed.

    name: fipronil_connectome_dependence
    compound: fipronil
    concentrations: [1e-8, 1e-7, 1e-6, 1e-5]
    graph: named
    readouts: [mean_hz, mn9_hz]
    engines: [rate]
    analyses:
      dependence: {permutations: 1000, correction: benjamini-hochberg}
      robustness: {specifications: default_family}
      uncertainty: {samples: 2048}

Three rules shape this module.

*Determinism.*  The same spec and the same seed must produce a byte-identical
``results.csv`` and identical artifact hashes.  Wall-clock timings and
timestamps obviously cannot be identical, so every artifact carries **two**
hashes: ``sha256`` over the bytes on disk and ``content_sha256`` over the same
payload with the declared volatile keys (:data:`VOLATILE_KEYS`) removed.  The
manifest's ``determinism.digest`` is a hash over the content hashes, and that is
what a reproduction check compares.

*No hardcoded results.*  Nothing in this module knows a numeric answer.  The
analysis entry points are discovered by name at call time
(:data:`ANALYSES`), imported lazily, and a missing function produces a message
naming the module, the names that were tried and what is actually exported.

*It must import without pydantic.*  ``ExperimentSpec`` is a pydantic model when
pydantic is importable and an equivalent dataclass otherwise, the same pattern
``flylab/assays/experiment.py`` uses, so the science core still runs in the
browser build.  Validation lives in plain functions shared by both paths, so the
error messages do not depend on which one is active.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import platform as _platform
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

__all__ = [
    "SPEC_VERSION",
    "RUN_MANIFEST_VERSION",
    "SpecError",
    "ExperimentSpec",
    "ENGINES",
    "ANALYSES",
    "VOLATILE_KEYS",
    "EXAMPLE_SPEC",
    "make_spec",
    "load_spec",
    "run_spec",
    "spec_schema",
    "estimate_runtime",
    "resolved_spec_dict",
    "spec_sha256",
    "HAVE_PYDANTIC",
]

SPEC_VERSION = "1.0"
RUN_MANIFEST_VERSION = "1.0"

#: engine name -> the assay in :mod:`flylab.assays.ensemble` that implements it.
#: ``subgraph``/``spiking`` are accepted as the historical aliases of the first
#: two so that a v0.5 design file is also a valid spec.
ENGINES: dict[str, str] = {
    "rate": "subgraph",
    "lif": "spiking",
    "taste": "taste",
    "taste_map": "taste_map",
    "subgraph": "subgraph",
    "spiking": "spiking",
}

#: the canonical engine names, in the order a schema should list them
ENGINE_NAMES: tuple[str, ...] = ("rate", "lif", "taste", "taste_map")

#: analysis name -> where to find it.  ``candidates`` is tried in order, so the
#: spec layer keeps working when the analysis layer renames its entry point.
ANALYSES: dict[str, dict[str, Any]] = {
    "dependence": {
        "module": "flylab.analysis.dependence",
        "candidates": ("dependence_profile",),
        "artifact": "dependence.json",
        "enums": {"correction": ("benjamini-hochberg", "bh", "none", "off")},
        "options": {
            "permutations": "int, permutations per null model (maps to n=)",
            "correction": "benjamini-hochberg | none, applied across the cells of this run",
            "modes": "list of null models, or null for the analysis default",
            "readout": "readout key, or 'auto'",
            "concentrations": "list of concentrations, or 'headline' (the highest in the spec) or 'all'",
            "alpha": "float, significance level",
            "n_jobs": "int, worker processes the analysis may use",
        },
    },
    "robustness": {
        "module": "flylab.analysis.robustness",
        "candidates": ("conclusion_stability",),
        "artifact": "robustness.json",
        "enums": {"specifications": ("default_family", "fast_family")},
        "options": {
            "specifications": "default_family | fast_family",
            "shuffles": "int, topology shuffles per conclusion, or null for the default",
            "n_jobs": "int, worker processes the analysis may use",
        },
    },
    "uncertainty": {
        "module": "flylab.analysis.uncertainty_global",
        "candidates": ("sobol_analysis",),
        "artifact": "uncertainty.json",
        "options": {
            "samples": "int, base samples per Saltelli matrix (maps to n_base=)",
            "readout": "readout key the variance budget is computed on",
            "bootstrap": "int, bootstrap resamples for the confidence intervals",
            "concentrations": "list of concentrations, or 'headline' or 'all'",
        },
    },
}

#: keys dropped, recursively, before an artifact's ``content_sha256`` is taken.
#: They are wall-clock or host facts: identical inputs cannot reproduce them and
#: no scientific result depends on them.
VOLATILE_KEYS: frozenset[str] = frozenset(
    {
        "created_utc",
        "generated_utc",
        "platform",
        "runtime_s",
        "runtime_ms",
        "runtime_by_spec_s",
        "elapsed_s",
        "wall_s",
        "duration_s",
        "timings_s",
    }
)

#: per-replicate wall-clock cost of each engine, seconds, measured on the
#: reference machine.  Used only by the runtime estimate; never by a result.
ENGINE_COST_S: dict[str, float] = {
    "subgraph": 0.25,
    "spiking": 0.80,
    "taste": 0.12,
    "taste_map": 0.55,
}

#: seconds per (permutation x null model) for the dependence analysis, and per
#: Sobol' model evaluation.  Order-of-magnitude planning numbers only.
DEPENDENCE_COST_S = 0.0125
SOBOL_COST_S = 0.0028
ROBUSTNESS_COST_S = 3.4  # per specification in the family

EXAMPLE_SPEC = """\
# A FlyLab experiment specification.  `flylab run this.yaml --dry-run` validates
# it and prints the resolved spec plus a runtime estimate; drop --dry-run to
# write a self-describing run directory.
name: fipronil_connectome_dependence
compound: fipronil                 # or:  compounds: [fipronil, imidacloprid]
concentrations: [1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5]
graph: named                       # named | taste_motor
readouts: [mean_hz, mn9_hz]
engines: [rate]                    # rate | lif | taste | taste_map
replicates: 1
seed: 0
analyses:
  dependence: {permutations: 1000, correction: benjamini-hochberg}
  robustness: {specifications: default_family}
  uncertainty: {samples: 2048}
"""


class SpecError(ValueError):
    """An experiment spec is not valid.  The message names the valid options."""


# --------------------------------------------------------------------------
# tiny helpers
# --------------------------------------------------------------------------
def _canonical(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def strip_volatile(obj: Any, keys: Iterable[str] = VOLATILE_KEYS) -> Any:
    """``obj`` with every :data:`VOLATILE_KEYS` entry removed, recursively."""
    drop = set(keys)
    if isinstance(obj, Mapping):
        return {k: strip_volatile(v, drop) for k, v in obj.items() if k not in drop}
    if isinstance(obj, (list, tuple)):
        return [strip_volatile(v, drop) for v in obj]
    return obj


def _suggest(value: Any, valid: Sequence[str], what: str) -> str:
    """``unknown X 'y'; valid are: ...`` with a near-miss hint when there is one."""
    text = str(value)
    near = [v for v in valid if v.lower().startswith(text.lower()[:3])] if text else []
    hint = f"  Did you mean {near[0]!r}?" if near else ""
    listed = ", ".join(valid) if valid else "(none available in this build)"
    return f"unknown {what} {value!r}. Valid options: {listed}.{hint}"


# --------------------------------------------------------------------------
# lazily discovered vocabularies -- the library, the graphs and the readouts
# all live in packages other agents own, so nothing here imports them eagerly
# --------------------------------------------------------------------------
def known_compounds() -> tuple[str, ...]:
    """Compound keys in the typed library, or ``()`` if it cannot be read."""
    try:
        from flylab.pharm.occupancy import load_library

        return tuple(sorted(load_library().get("compounds") or {}))
    except Exception:  # pragma: no cover - only when the pharm layer is broken
        return ()


def known_graphs() -> tuple[str, ...]:
    """Derived graph names FlyLab ships, or ``()`` if the map layer is absent."""
    try:
        from flylab.circuit.rate import GRAPHS

        return tuple(sorted(GRAPHS))
    except Exception:  # pragma: no cover
        return ()


def known_readouts() -> tuple[str, ...]:
    """The comparable scalar readout keys."""
    try:
        from flylab.assays.ensemble import READOUT_KEYS

        return tuple(READOUT_KEYS)
    except Exception:  # pragma: no cover - the assay layer mid-edit
        return ("mn9_hz", "dnp01_hz", "mean_hz", "g_ach", "g_gaba")


# --------------------------------------------------------------------------
# validation, shared by the pydantic and the dataclass path
# --------------------------------------------------------------------------
#: canonical field -> accepted input aliases (a v0.5 design file is a valid spec)
ALIASES: dict[str, tuple[str, ...]] = {
    "concentrations": ("concs_M", "concs", "conc_M"),
    "engines": ("engine", "assay", "assays"),
    "compounds": ("compound",),
    "seed": ("rng_seed",),
    "replicates": ("n_rep",),
}

FIELD_ORDER: tuple[str, ...] = (
    "flylab_spec_version",
    "name",
    "description",
    "compounds",
    "concentrations",
    "graph",
    "engines",
    "readouts",
    "replicates",
    "seed",
    "include_vehicle",
    "jitter_log10",
    "drive_hz",
    "analyses",
    "figures",
    "cards",
    "options",
)


def _defaults() -> dict[str, Any]:
    return {
        "flylab_spec_version": SPEC_VERSION,
        "name": "flylab_run",
        "description": None,
        "compounds": ["imidacloprid"],
        "concentrations": [1e-9, 1e-8, 1e-7, 1e-6, 1e-5],
        "graph": "named",
        "engines": ["rate"],
        "readouts": list(known_readouts()),
        "replicates": 1,
        "seed": 0,
        "include_vehicle": True,
        "jitter_log10": 0.0,
        "drive_hz": None,
        "analyses": {},
        "figures": True,
        "cards": True,
        "options": {},
    }


def _as_bool(v: Any) -> bool:
    if isinstance(v, str):
        low = v.strip().lower()
        if low in ("true", "yes", "on", "1"):
            return True
        if low in ("false", "no", "off", "0", ""):
            return False
        raise SpecError(f"expected a boolean, got {v!r}")
    return bool(v)


def _as_list(v: Any) -> list[Any]:
    if v is None:
        return []
    if isinstance(v, (str, bytes)) or isinstance(v, Mapping):
        return [v]
    if isinstance(v, (list, tuple, set)):
        return list(v)
    return [v]


def _normalise_keys(data: Mapping[str, Any]) -> dict[str, Any]:
    """Apply the aliases and reject anything that is not a spec field."""
    out: dict[str, Any] = {}
    alias_to_field = {a: f for f, aliases in ALIASES.items() for a in aliases}
    valid = set(FIELD_ORDER) | set(alias_to_field)
    for key, value in dict(data).items():
        name = str(key)
        if name in ("flylab_spec_version", "spec_version"):
            out["flylab_spec_version"] = value
            continue
        field = alias_to_field.get(name, name)
        if field not in FIELD_ORDER:
            raise SpecError(
                _suggest(name, sorted(valid), "spec key")
                + "  (a v0.5 design file is also a valid spec: assay, concs_M and "
                "compound are accepted as aliases.)"
            )
        if field in out and out[field] not in (None, [], {}):
            # `compound: fipronil` plus `compounds: [...]` -- merge rather than drop
            if field == "compounds":
                out[field] = _as_list(out[field]) + _as_list(value)
                continue
            raise SpecError(f"spec key {name!r} duplicates {field!r}; give only one")
        out[field] = value
    return out


def _validate_analyses(raw: Any) -> dict[str, dict[str, Any]]:
    if raw is None:
        return {}
    if isinstance(raw, (list, tuple, set)):  # `analyses: [dependence]`
        raw = {str(name): {} for name in raw}
    if not isinstance(raw, Mapping):
        raise SpecError(
            "analyses must be a mapping of analysis name -> options "
            f"(or a list of names). Valid analyses: {', '.join(sorted(ANALYSES))}."
        )
    out: dict[str, dict[str, Any]] = {}
    for name, opts in raw.items():
        key = str(name)
        if key not in ANALYSES:
            raise SpecError(_suggest(key, sorted(ANALYSES), "analysis"))
        if opts in (None, True):
            opts = {}
        if opts is False:
            continue  # `dependence: false` switches one off without deleting it
        if not isinstance(opts, Mapping):
            raise SpecError(
                f"options for analysis {key!r} must be a mapping; got {type(opts).__name__}. "
                f"Valid options: {', '.join(sorted(ANALYSES[key]['options']))}."
            )
        allowed = set(ANALYSES[key]["options"])
        enums = dict(ANALYSES[key].get("enums") or {})
        for opt, value in opts.items():
            if str(opt) not in allowed:
                raise SpecError(
                    _suggest(opt, sorted(allowed), f"option for analysis {key!r}")
                )
            choices = enums.get(str(opt))
            if choices and str(value).lower() not in choices:
                raise SpecError(
                    _suggest(value, list(choices), f"{opt} for analysis {key!r}")
                )
        out[key] = {str(k): v for k, v in opts.items()}
    return out


def _coerce(data: Mapping[str, Any] | None) -> dict[str, Any]:
    """Validate a raw spec mapping and return the fully resolved field dict."""
    if data is None:
        data = {}
    if not isinstance(data, Mapping):
        raise SpecError(f"a spec must be a mapping, got {type(data).__name__}")
    given = _normalise_keys(data)
    out = _defaults()

    version = given.pop("flylab_spec_version", SPEC_VERSION)
    if str(version).split(".")[0] != SPEC_VERSION.split(".")[0]:
        raise SpecError(
            f"spec version {version!r} is not supported; this build reads {SPEC_VERSION}"
        )
    out["flylab_spec_version"] = SPEC_VERSION

    if "name" in given and given["name"] is not None:
        name = str(given["name"]).strip()
        if not name:
            raise SpecError("name must not be empty")
        out["name"] = name
    if "description" in given:
        out["description"] = None if given["description"] is None else str(given["description"])

    # compounds -----------------------------------------------------------
    if "compounds" in given:
        compounds = [str(c).strip() for c in _as_list(given["compounds"]) if str(c).strip()]
        if not compounds:
            raise SpecError("a spec needs at least one compound")
        seen: list[str] = []
        for c in compounds:
            if c not in seen:
                seen.append(c)
        known = known_compounds()
        if known:
            for c in seen:
                if c not in known:
                    raise SpecError(_suggest(c, list(known), "compound"))
        out["compounds"] = seen

    # concentrations -------------------------------------------------------
    if "concentrations" in given:
        raw_concs = _as_list(given["concentrations"])
        if not raw_concs:
            raise SpecError("a spec needs at least one concentration (free molar, e.g. 1.0e-6)")
        concs: list[float] = []
        for c in raw_concs:
            try:
                value = float(c)
            except (TypeError, ValueError):
                raise SpecError(f"concentration {c!r} is not a number (free molar, e.g. 1.0e-6)")
            if value < 0:
                raise SpecError("concentrations must be >= 0 M")
            concs.append(value)
        out["concentrations"] = sorted(dict.fromkeys(concs))

    # graph ----------------------------------------------------------------
    if "graph" in given and given["graph"] is not None:
        graph = str(given["graph"])
        graphs = known_graphs()
        if graphs and graph not in graphs and not Path(graph).suffix:
            raise SpecError(_suggest(graph, list(graphs), "graph"))
        out["graph"] = graph

    # engines --------------------------------------------------------------
    if "engines" in given:
        engines = [str(e).strip() for e in _as_list(given["engines"]) if str(e).strip()]
        if not engines:
            raise SpecError(f"a spec needs at least one engine. Valid: {', '.join(ENGINE_NAMES)}.")
        canonical: list[str] = []
        for engine in engines:
            if engine not in ENGINES:
                raise SpecError(_suggest(engine, list(ENGINE_NAMES), "engine"))
            if engine not in canonical:
                canonical.append(engine)
        out["engines"] = canonical

    # readouts -------------------------------------------------------------
    if "readouts" in given:
        readouts = [str(r).strip() for r in _as_list(given["readouts"]) if str(r).strip()]
        if not readouts:
            raise SpecError(f"a spec needs at least one readout. Valid: {', '.join(known_readouts())}.")
        valid_readouts = known_readouts()
        for r in readouts:
            if r not in valid_readouts:
                raise SpecError(_suggest(r, list(valid_readouts), "readout"))
        out["readouts"] = list(dict.fromkeys(readouts))

    # scalars --------------------------------------------------------------
    if "replicates" in given:
        try:
            replicates = int(given["replicates"])
        except (TypeError, ValueError):
            raise SpecError(f"replicates must be an integer >= 1, got {given['replicates']!r}")
        if replicates < 1:
            raise SpecError("replicates must be >= 1")
        out["replicates"] = replicates
    if "seed" in given:
        try:
            out["seed"] = int(given["seed"])
        except (TypeError, ValueError):
            raise SpecError(f"seed must be an integer, got {given['seed']!r}")
    if "include_vehicle" in given:
        out["include_vehicle"] = _as_bool(given["include_vehicle"])
    if "jitter_log10" in given:
        jitter = float(given["jitter_log10"] or 0.0)
        if jitter < 0:
            raise SpecError("jitter_log10 must be >= 0 (0 disables the library Monte-Carlo)")
        out["jitter_log10"] = jitter
    if "drive_hz" in given:
        out["drive_hz"] = None if given["drive_hz"] is None else float(given["drive_hz"])
    if "figures" in given:
        out["figures"] = _as_bool(given["figures"])
    if "cards" in given:
        out["cards"] = _as_bool(given["cards"])
    if "options" in given:
        opts = given["options"] or {}
        if not isinstance(opts, Mapping):
            raise SpecError("options must be a mapping passed through to the assay")
        out["options"] = dict(opts)

    out["analyses"] = _validate_analyses(given.get("analyses"))
    return {k: out[k] for k in FIELD_ORDER}


# --------------------------------------------------------------------------
# the model itself: pydantic when available, an equivalent dataclass otherwise
# --------------------------------------------------------------------------
try:  # pragma: no cover - both paths are exercised, one per environment
    from pydantic import BaseModel, ConfigDict, Field

    HAVE_PYDANTIC = True
except Exception:  # pragma: no cover - selected on Pyodide / minimal installs
    HAVE_PYDANTIC = False


if HAVE_PYDANTIC:

    class ExperimentSpec(BaseModel):  # type: ignore[no-redef]
        """A validated, fully resolved experiment specification.

        Build one with :func:`make_spec` or :func:`load_spec`; both raise
        :class:`SpecError` with the valid options named.
        """

        model_config = ConfigDict(extra="forbid")

        flylab_spec_version: str = SPEC_VERSION
        name: str = "flylab_run"
        description: str | None = None
        compounds: list[str] = Field(default_factory=lambda: _defaults()["compounds"])
        concentrations: list[float] = Field(default_factory=lambda: _defaults()["concentrations"])
        graph: str = "named"
        engines: list[str] = Field(default_factory=lambda: ["rate"])
        readouts: list[str] = Field(default_factory=lambda: list(known_readouts()))
        replicates: int = 1
        seed: int = 0
        include_vehicle: bool = True
        jitter_log10: float = 0.0
        drive_hz: float | None = None
        analyses: dict[str, dict[str, Any]] = Field(default_factory=dict)
        figures: bool = True
        cards: bool = True
        options: dict[str, Any] = Field(default_factory=dict)

        @classmethod
        def from_mapping(cls, data: Mapping[str, Any] | None) -> "ExperimentSpec":
            return cls(**_coerce(data))

else:  # pragma: no cover - selected only when pydantic is absent

    from dataclasses import asdict, dataclass, field

    @dataclass
    class ExperimentSpec:  # type: ignore[no-redef]
        """pydantic-free twin of the model above: same fields, same defaults."""

        flylab_spec_version: str = SPEC_VERSION
        name: str = "flylab_run"
        description: str | None = None
        compounds: list[str] = field(default_factory=lambda: _defaults()["compounds"])
        concentrations: list[float] = field(default_factory=lambda: _defaults()["concentrations"])
        graph: str = "named"
        engines: list[str] = field(default_factory=lambda: ["rate"])
        readouts: list[str] = field(default_factory=lambda: list(known_readouts()))
        replicates: int = 1
        seed: int = 0
        include_vehicle: bool = True
        jitter_log10: float = 0.0
        drive_hz: float | None = None
        analyses: dict[str, dict[str, Any]] = field(default_factory=dict)
        figures: bool = True
        cards: bool = True
        options: dict[str, Any] = field(default_factory=dict)

        def __init__(self, **kw: Any) -> None:
            unknown = [k for k in kw if k not in FIELD_ORDER]
            if unknown:
                raise SpecError(_suggest(unknown[0], list(FIELD_ORDER), "spec key"))
            values = _defaults()
            values.update({k: v for k, v in kw.items()})
            for name in FIELD_ORDER:
                object.__setattr__(self, name, values[name])

        def model_dump(self, **_: Any) -> dict[str, Any]:
            return asdict(self)

        def dict(self, **kw: Any) -> dict[str, Any]:  # pragma: no cover - alias
            return self.model_dump(**kw)

        @classmethod
        def from_mapping(cls, data: Mapping[str, Any] | None) -> "ExperimentSpec":
            return cls(**_coerce(data))


def make_spec(data: Mapping[str, Any] | ExperimentSpec | None) -> ExperimentSpec:
    """Validate a raw mapping (or pass an ``ExperimentSpec`` through)."""
    if isinstance(data, ExperimentSpec):
        return data
    return ExperimentSpec.from_mapping(data)


def load_spec(path: str | Path) -> ExperimentSpec:
    """Read and validate a spec YAML (or JSON) file.

    Raises :class:`SpecError` with the valid options named for an unknown key,
    compound, graph, engine, readout or analysis.
    """
    import yaml

    p = Path(path)
    if not p.exists():
        raise SpecError(f"spec file not found: {p}")
    try:
        data = yaml.safe_load(p.read_text()) or {}
    except Exception as exc:
        raise SpecError(f"{p} is not valid YAML: {exc}") from exc
    if isinstance(data, Mapping) and isinstance(data.get("spec"), Mapping):
        data = data["spec"]  # a spec nested under `spec:` is accepted
    if isinstance(data, Mapping) and isinstance(data.get("design"), Mapping):
        data = data["design"]  # ... and so is a v0.5 design file
    if not isinstance(data, Mapping):
        raise SpecError(f"{p} must contain a mapping at the top level, got {type(data).__name__}")
    return make_spec(data)


def resolved_spec_dict(spec: ExperimentSpec) -> dict[str, Any]:
    """The spec as a plain dict, in declaration order, defaults included."""
    dumped = spec.model_dump()
    return {k: dumped[k] for k in FIELD_ORDER if k in dumped}


def spec_sha256(spec: ExperimentSpec | Mapping[str, Any]) -> str:
    """SHA-256 over the canonical JSON of the *resolved* spec."""
    payload = resolved_spec_dict(spec) if isinstance(spec, ExperimentSpec) else dict(spec)
    return _sha256_text(_canonical(payload))


# --------------------------------------------------------------------------
# schema
# --------------------------------------------------------------------------
def spec_schema() -> dict[str, Any]:
    """A machine-readable description of the spec format.

    JSON-Schema shaped (``type`` / ``properties`` / ``enum`` / ``default``) with
    two FlyLab additions: ``aliases`` records the v0.5 design keys each field
    still accepts, and ``x-analyses`` documents the per-analysis option blocks.
    Enumerations that depend on installed data (the compound library, the
    derived graphs) are filled in when they can be read and omitted otherwise.
    """
    defaults = _defaults()
    compounds = list(known_compounds())
    graphs = list(known_graphs())
    readouts = list(known_readouts())

    def prop(kind: str, description: str, *, enum: list[Any] | None = None, **extra: Any) -> dict[str, Any]:
        out: dict[str, Any] = {"type": kind, "description": description}
        if enum:
            out["enum"] = enum
        out.update(extra)
        return out

    properties: dict[str, Any] = {
        "flylab_spec_version": prop("string", "Spec format version.", default=SPEC_VERSION),
        "name": prop("string", "Run name; used for the run directory and the manifest.", default=defaults["name"]),
        "description": prop("string", "Free text recorded in the manifest.", default=None, nullable=True),
        "compounds": prop(
            "array",
            "Compound keys from the typed library. `compound: x` is accepted for a single one.",
            items={"type": "string", **({"enum": compounds} if compounds else {})},
            default=defaults["compounds"],
            minItems=1,
        ),
        "concentrations": prop(
            "array",
            "Free concentrations in molar. Sorted and de-duplicated on load.",
            items={"type": "number", "minimum": 0},
            default=defaults["concentrations"],
            minItems=1,
        ),
        "graph": prop(
            "string",
            "Derived MaleCNS cut the circuit runs on.",
            enum=graphs or None,
            default=defaults["graph"],
        ),
        "engines": prop(
            "array",
            "Circuit runtimes to run. `rate` is the deterministic leaky-rate engine; "
            "`lif` is the Shiu-style spiking engine.",
            items={"type": "string", "enum": list(ENGINE_NAMES)},
            default=defaults["engines"],
            minItems=1,
        ),
        "readouts": prop(
            "array",
            "Scalar readouts written to results.csv.",
            items={"type": "string", **({"enum": readouts} if readouts else {})},
            default=defaults["readouts"],
            minItems=1,
        ),
        "replicates": prop("integer", "Replicates per cell.", default=1, minimum=1),
        "seed": prop("integer", "Master RNG seed; replicate r uses seed + r.", default=0),
        "include_vehicle": prop("boolean", "Also run a vehicle (compound = null) cell.", default=True),
        "jitter_log10": prop(
            "number",
            "Log10 SD of the library Monte-Carlo. 0 disables it and keeps the run deterministic "
            "in the strict sense; any value > 0 is still reproducible from the seed.",
            default=0.0,
            minimum=0,
        ),
        "drive_hz": prop("number", "Sensory drive override; null uses the assay default.", default=None, nullable=True),
        "analyses": prop(
            "object",
            "Which analysis layers to run, and their options. Omit to run none.",
            default={},
            properties={
                name: {
                    "type": "object",
                    "description": f"Options for the {name} analysis.",
                    "properties": {
                        opt: {"description": text} for opt, text in meta["options"].items()
                    },
                    "additionalProperties": False,
                }
                for name, meta in ANALYSES.items()
            },
            additionalProperties=False,
        ),
        "figures": prop("boolean", "Write figures when matplotlib is installed.", default=True),
        "cards": prop("boolean", "Write one claim card per headline result.", default=True),
        "options": prop("object", "Extra keyword arguments passed through to the assay.", default={}),
    }

    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://github.com/pinkysworld/FlyLab/spec/v" + SPEC_VERSION,
        "title": "FlyLab experiment specification",
        "flylab_spec_version": SPEC_VERSION,
        "type": "object",
        "additionalProperties": False,
        "required": [],
        "properties": properties,
        "aliases": {field: list(aliases) for field, aliases in ALIASES.items()},
        "x-analyses": {
            name: {
                "module": meta["module"],
                "entry_points_tried": list(meta["candidates"]),
                "artifact": meta["artifact"],
                "options": dict(meta["options"]),
                "choices": {k: list(v) for k, v in (meta.get("enums") or {}).items()},
            }
            for name, meta in ANALYSES.items()
        },
        "x-run-directory": {
            "manifest.json": "spec hash, code version, git sha, library sha, map ids, seeds, timings, artifact hashes",
            "spec.yaml": "the spec as resolved, including the defaults that were filled in",
            "notebooks/": "one notebook JSON per assay run",
            "results.csv": "the tabular readouts",
            "dependence.json": "written only when the dependence analysis was requested",
            "uncertainty.json": "written only when the uncertainty analysis was requested",
            "robustness.json": "written only when the robustness analysis was requested",
            "cards/": "one claim card per headline result, as JSON and as Markdown",
            "figures/": "written only when matplotlib is importable",
        },
        "x-volatile-keys": sorted(VOLATILE_KEYS),
        "example": EXAMPLE_SPEC,
    }


# --------------------------------------------------------------------------
# lazily discovered analysis entry points
# --------------------------------------------------------------------------
def analysis_entry_point(name: str) -> tuple[Callable[..., Any], str]:
    """``(function, dotted_name)`` for one analysis, imported at call time.

    The analysis layer is owned by other modules and its function names move,
    so nothing here holds a reference to one.  Every candidate in
    :data:`ANALYSES` is tried in order; if none exists the error names the
    module, what was tried and what that module does export.
    """
    meta = ANALYSES.get(name)
    if meta is None:
        raise SpecError(_suggest(name, sorted(ANALYSES), "analysis"))
    module_name = str(meta["module"])
    try:
        module = importlib.import_module(module_name)
    except Exception as exc:  # the module may be mid-rewrite, or absent in a slim build
        raise SpecError(
            f"analysis {name!r} needs {module_name}, which could not be imported "
            f"({type(exc).__name__}: {exc}). Remove it from the spec's `analyses:` block "
            "to run the assays without it."
        ) from exc
    for candidate in meta["candidates"]:
        fn = getattr(module, candidate, None)
        if callable(fn):
            return fn, f"{module_name}.{candidate}"
    exported = sorted(getattr(module, "__all__", None) or [n for n in vars(module) if not n.startswith("_")])
    raise SpecError(
        f"analysis {name!r}: none of {', '.join(meta['candidates'])} is exported by "
        f"{module_name}. It currently exports: {', '.join(exported) or '(nothing public)'}. "
        "Update flylab.spec.ANALYSES if the entry point was renamed."
    )


def _call_filtered(fn: Callable[..., Any], **kw: Any) -> Any:
    """Call ``fn`` with only the keyword arguments it actually accepts.

    The analysis signatures change between versions; a spec must degrade to the
    analysis's own default rather than crash on a dropped parameter.
    """
    import inspect

    try:
        sig = inspect.signature(fn)
    except (TypeError, ValueError):  # pragma: no cover - builtins
        return fn(**kw)
    if any(p.kind is inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()):
        accepted = kw
    else:
        accepted = {k: v for k, v in kw.items() if k in sig.parameters}
    return fn(**accepted)


def _analysis_concentrations(value: Any, spec: ExperimentSpec) -> list[float]:
    """``headline`` (the highest dose), ``all``, or an explicit list."""
    concs = list(spec.concentrations)
    if value in (None, "headline", "max"):
        return [max(concs)] if concs else []
    if value == "all":
        return concs
    out: list[float] = []
    for item in _as_list(value):
        try:
            out.append(float(item))
        except (TypeError, ValueError):
            raise SpecError(
                f"concentrations for an analysis must be numbers, 'headline' or 'all'; got {item!r}"
            )
    return sorted(dict.fromkeys(out))


# --------------------------------------------------------------------------
# runtime estimate
# --------------------------------------------------------------------------
def estimate_runtime(spec: ExperimentSpec | Mapping[str, Any]) -> dict[str, Any]:
    """Order-of-magnitude wall-clock estimate for a spec, without running it.

    The per-unit costs are measured planning constants (:data:`ENGINE_COST_S`
    and friends), not results; they scale linearly and are wrong by a factor of
    a few on unfamiliar hardware.  ``--dry-run`` prints this.
    """
    s = make_spec(spec)
    n_cells = len(s.compounds) * len(s.concentrations) * s.replicates
    per_engine: dict[str, Any] = {}
    assay_s = 0.0
    runs = 0
    for engine in s.engines:
        assay = ENGINES[engine]
        cost = ENGINE_COST_S.get(assay, 0.3)
        n_runs = n_cells + (s.replicates if s.include_vehicle else 0)
        per_engine[engine] = {"assay": assay, "runs": n_runs, "seconds": round(cost * n_runs, 1)}
        assay_s += cost * n_runs
        runs += n_runs

    analyses: dict[str, Any] = {}
    for name, opts in s.analyses.items():
        if name == "dependence":
            n = int(opts.get("permutations", 1000))
            modes = opts.get("modes")
            n_modes = len(modes) if modes else 4
            cells = len(s.compounds) * len(_analysis_concentrations(opts.get("concentrations"), s))
            seconds = DEPENDENCE_COST_S * n * n_modes * max(cells, 1)
            detail: dict[str, Any] = {"cells": cells, "permutations": n, "null_models": n_modes}
            # prefer the analysis layer's own estimator when it exposes one
            try:
                fn = getattr(importlib.import_module(ANALYSES[name]["module"]), "estimate_landscape_runtime", None)
                if callable(fn):
                    own = _call_filtered(fn, n_cells=max(cells, 1), n=n, graph=s.graph)
                    if isinstance(own, Mapping):
                        for key in ("estimate_s", "seconds", "runtime_s", "total_s"):
                            if own.get(key) is not None:
                                seconds = float(own[key])
                                detail["source"] = f"{ANALYSES[name]['module']}.estimate_landscape_runtime"
                                break
            except Exception:  # pragma: no cover - the estimator is optional
                pass
            analyses[name] = {**detail, "seconds": round(seconds, 1)}
        elif name == "robustness":
            family = str(opts.get("specifications", "default_family"))
            n_specs = 9 if family == "fast_family" else 25
            seconds = ROBUSTNESS_COST_S * n_specs
            analyses[name] = {"specifications": family, "family_size_assumed": n_specs, "seconds": round(seconds, 1)}
        elif name == "uncertainty":
            samples = int(opts.get("samples", 128))
            cells = len(s.compounds) * len(_analysis_concentrations(opts.get("concentrations"), s))
            evals = samples * 11 * max(cells, 1)  # k + 2 with k = 9 factors
            analyses[name] = {
                "base_samples": samples,
                "cells": cells,
                "model_evaluations_assumed": evals,
                "seconds": round(SOBOL_COST_S * evals, 1),
            }
    total = assay_s + sum(float(a["seconds"]) for a in analyses.values())
    return {
        "n_assay_runs": runs,
        "assays": per_engine,
        "assays_seconds": round(assay_s, 1),
        "analyses": analyses,
        "estimate_s": round(total, 1),
        "estimate_human": _human_time(total),
        "note": (
            "Linear extrapolation from measured per-unit costs on the reference machine; "
            "it is a planning number, not a measurement of this run."
        ),
    }


def _human_time(seconds: float) -> str:
    if seconds < 90:
        return f"{seconds:.0f} s"
    if seconds < 5400:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.1f} h"


# --------------------------------------------------------------------------
# artifacts and hashing
# --------------------------------------------------------------------------
#: file suffixes whose payload is JSON and can therefore be hashed with the
#: volatile keys removed
_JSONISH = {".json"}


def artifact_entry(path: Path) -> dict[str, Any]:
    """``{bytes, sha256, content_sha256, volatile_fields_removed}`` for one file.

    ``sha256`` is over the bytes on disk; ``content_sha256`` is over the same
    payload with :data:`VOLATILE_KEYS` stripped, and is what two runs of the
    same spec are compared on.  They are equal for every artifact that carries
    no timestamp or wall-clock timing.
    """
    data = path.read_bytes()
    entry: dict[str, Any] = {
        "bytes": len(data),
        "sha256": _sha256_bytes(data),
        "content_sha256": _sha256_bytes(data),
        "volatile_fields_removed": [],
    }
    if path.suffix.lower() in _JSONISH:
        try:
            payload = json.loads(data.decode("utf-8"))
        except Exception:  # pragma: no cover - a non-JSON .json file
            return entry
        removed = sorted(_volatile_present(payload))
        entry["content_sha256"] = _sha256_text(_canonical(strip_volatile(payload)))
        entry["volatile_fields_removed"] = removed
    return entry


def _volatile_present(obj: Any, found: set[str] | None = None) -> set[str]:
    out = found if found is not None else set()
    if isinstance(obj, Mapping):
        for k, v in obj.items():
            if k in VOLATILE_KEYS:
                out.add(str(k))
            else:
                _volatile_present(v, out)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            _volatile_present(v, out)
    return out


def collect_artifacts(run_dir: Path, skip: Sequence[str] = ("manifest.json",)) -> dict[str, Any]:
    """Every file under ``run_dir`` (except ``skip``) with its two hashes."""
    root = Path(run_dir)
    out: dict[str, Any] = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(root).as_posix()
        if rel in set(skip):
            continue
        out[rel] = artifact_entry(path)
    return out


def determinism_digest(artifacts: Mapping[str, Any]) -> str:
    """One hash over ``path -> content_sha256``: the reproduction fingerprint."""
    lines = [f"{path}\t{meta['content_sha256']}" for path, meta in sorted(artifacts.items())]
    return _sha256_text("\n".join(lines))


# --------------------------------------------------------------------------
# the run
# --------------------------------------------------------------------------
def _conc_tag(conc: float) -> str:
    return "0" if not conc else f"{conc:.2e}".replace("-0", "-").replace("+0", "")


def _notebook_name(engine: str, key: Mapping[str, Any]) -> str:
    compound = key.get("compound") or "vehicle"
    return f"{engine}__{compound}__{_conc_tag(float(key.get('conc_M') or 0.0))}__r{int(key.get('replicate') or 0)}.json"


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def _graph_facts(name: str) -> dict[str, Any]:
    try:
        from flylab.circuit.rate import load_graph, resolve_graph

        g = load_graph(resolve_graph(name))
    except Exception as exc:  # pragma: no cover - a missing derived graph
        return {"graph": name, "available": False, "reason": f"{type(exc).__name__}: {exc}"}
    return {
        "graph": name,
        "available": True,
        "map": g.get("map"),
        "citation": g.get("citation"),
        "n_nodes": g.get("n_nodes"),
        "n_edges": g.get("n_edges"),
        "hops": g.get("hops"),
        "min_weight": g.get("min_weight"),
        "closure_min_weight": g.get("closure_min_weight"),
    }


def _code_provenance(seed: int) -> dict[str, Any]:
    out: dict[str, Any] = {
        "flylab_version": None,
        "git_sha": None,
        "library_sha256": None,
        "python": _platform.python_version(),
        "platform": _platform.platform(),
    }
    try:
        from flylab.notebook.schema import provenance

        p = provenance(seed)
        out.update(
            {
                "flylab_version": p.get("flylab_version"),
                "git_sha": p.get("git_sha"),
                "library_sha256": p.get("library_sha256"),
            }
        )
    except Exception:  # pragma: no cover
        try:
            import flylab

            out["flylab_version"] = getattr(flylab, "__version__", None)
        except Exception:
            pass
    try:
        from flylab.pharm.occupancy import library_sha256, load_library

        out["library_sha256"] = library_sha256()
        lib = load_library()
        out["library_version"] = lib.get("library_version")
        out["library_schema_version"] = lib.get("schema_version")
    except Exception:  # pragma: no cover - the pharm layer mid-edit
        pass
    return out


def _map_ids() -> dict[str, Any]:
    try:
        from flylab.maps.malecns import MAP_CITATION, MAP_ID

        return {"map_id": MAP_ID, "map_citation": MAP_CITATION}
    except Exception:  # pragma: no cover
        return {"map_id": None, "map_citation": None}


def _run_assays(spec: ExperimentSpec, run_dir: Path, say: Callable[[str], None]) -> dict[str, Any]:
    """Every engine x compound x concentration x replicate, plus the notebooks."""
    from flylab.assays.experiment import ExperimentDesign, rows_to_csv, run_experiment

    rows: list[dict[str, Any]] = []
    notebook_paths: list[str] = []
    warnings: list[str] = []
    headline: dict[tuple[str, str], dict[str, Any]] = {}
    vehicles: dict[str, dict[str, Any]] = {}
    per_engine: dict[str, Any] = {}
    nb_dir = run_dir / "notebooks"
    nb_dir.mkdir(parents=True, exist_ok=True)
    headline_conc = max(spec.concentrations) if spec.concentrations else 0.0

    for engine in spec.engines:
        assay = ENGINES[engine]
        t0 = time.perf_counter()
        design = ExperimentDesign(
            assay=assay,
            compounds=list(spec.compounds),
            concs_M=list(spec.concentrations),
            replicates=spec.replicates,
            seed=spec.seed,
            readouts=list(spec.readouts),
            jitter_log10=spec.jitter_log10,
            drive_hz=spec.drive_hz,
            include_vehicle=spec.include_vehicle,
            graph=spec.graph,
            keep_notebooks=True,
            options=dict(spec.options),
        )
        result = run_experiment(design)
        say(f"    {engine}: {result['n_rows']} rows + {len(result['vehicle_rows'])} vehicle rows")
        for row in result["rows"]:
            rows.append({"engine": engine, "condition": "drug", **row})
        for row in result["vehicle_rows"]:
            rows.append({"engine": engine, "condition": "vehicle", **row})
        for key, nb in zip(result.get("notebook_keys") or [], result.get("notebooks") or []):
            name = _notebook_name(engine, key)
            _write_json(nb_dir / name, nb)
            notebook_paths.append(f"notebooks/{name}")
            if key.get("condition") == "vehicle" and int(key.get("replicate") or 0) == 0:
                vehicles[engine] = nb
            if (
                key.get("condition") == "drug"
                and int(key.get("replicate") or 0) == 0
                and float(key.get("conc_M") or 0.0) == headline_conc
            ):
                headline[(str(key.get("compound")), engine)] = {
                    "notebook": nb,
                    "path": f"notebooks/{name}",
                    "compound": key.get("compound"),
                    "conc_M": key.get("conc_M"),
                    "engine": engine,
                    "assay": assay,
                }
        for w in result.get("warnings") or []:
            if w not in warnings:
                warnings.append(w)
        per_engine[engine] = {
            "assay": assay,
            "n_rows": result["n_rows"],
            "n_vehicle_rows": len(result["vehicle_rows"]),
            "seconds": round(time.perf_counter() - t0, 3),
        }

    for (_compound, engine), item in headline.items():
        item["vehicle"] = vehicles.get(engine)

    fields = ["engine", "condition", "compound", "conc_M", "replicate"] + list(spec.readouts)
    csv_text = rows_to_csv(rows, fields)
    (run_dir / "results.csv").write_text(csv_text)
    return {
        "rows": rows,
        "csv": csv_text,
        "notebook_paths": notebook_paths,
        "headline": headline,
        "warnings": warnings,
        "per_engine": per_engine,
        "headline_conc_M": headline_conc,
    }


def _run_dependence(spec: ExperimentSpec, opts: Mapping[str, Any], say: Callable[[str], None]) -> dict[str, Any]:
    fn, dotted = analysis_entry_point("dependence")
    n = int(opts.get("permutations", 1000))
    correction = str(opts.get("correction", "benjamini-hochberg")).lower()
    if correction not in ("benjamini-hochberg", "bh", "none", "off"):
        raise SpecError(
            _suggest(correction, ["benjamini-hochberg", "none"], "dependence correction")
        )
    concs = _analysis_concentrations(opts.get("concentrations"), spec)
    assay = ENGINES[spec.engines[0]]
    readout = opts.get("readout", "auto")
    cells: list[dict[str, Any]] = []
    for compound in spec.compounds:
        for conc in concs:
            say(f"    dependence: {compound} @ {conc:.2e} M, n={n}")
            kw: dict[str, Any] = {
                "compound": compound,
                "conc_M": conc,
                "assay": assay,
                "graph": spec.graph,
                "readout": readout,
                "n": n,
                "seed": spec.seed,
                "alpha": float(opts.get("alpha", 0.05)),
                "n_jobs": int(opts.get("n_jobs", 1)),
            }
            if opts.get("modes"):  # otherwise leave the analysis its own default set
                kw["modes"] = tuple(opts["modes"])
            cells.append(_call_filtered(fn, **kw))
    payload: dict[str, Any] = {
        "analysis": "dependence",
        "entry_point": dotted,
        "settings": {
            "permutations": n,
            "correction": correction,
            "assay": assay,
            "graph": spec.graph,
            "readout": readout,
            "concentrations": concs,
            "seed": spec.seed,
        },
        "cells": cells,
    }
    if correction in ("benjamini-hochberg", "bh"):
        payload["multiplicity_correction"] = _apply_bh(cells, float(opts.get("alpha", 0.05)))
    else:
        payload["multiplicity_correction"] = {
            "method": "none",
            "note": "no multiplicity correction was requested in the spec",
        }
    return payload


def _apply_bh(cells: Sequence[Mapping[str, Any]], alpha: float) -> dict[str, Any]:
    """Benjamini-Hochberg across every (cell, null model) permutation p."""
    labels: list[str] = []
    pvalues: list[float | None] = []
    for cell in cells:
        for mode in cell.get("modes") or []:
            labels.append(
                f"{cell.get('compound')}@{cell.get('conc_M')}:{mode.get('mode')}"
            )
            p = mode.get("p_two_sided", mode.get("p"))
            pvalues.append(float(p) if isinstance(p, (int, float)) else None)
    if not labels:
        return {"method": "benjamini-hochberg", "n_tests": 0, "note": "no p-values to correct"}
    try:
        module = importlib.import_module(ANALYSES["dependence"]["module"])
        bh = getattr(module, "benjamini_hochberg", None)
        if not callable(bh):
            raise AttributeError("benjamini_hochberg")
        out = dict(_call_filtered(bh, pvalues=pvalues, alpha=alpha))
    except Exception as exc:
        return {
            "method": "benjamini-hochberg",
            "n_tests": len(labels),
            "available": False,
            "reason": (
                f"{ANALYSES['dependence']['module']}.benjamini_hochberg is not available "
                f"({type(exc).__name__}: {exc}); the uncorrected per-mode p values stand."
            ),
            "tests": labels,
        }
    out["method"] = "benjamini-hochberg"
    out["alpha"] = alpha
    out["tests"] = labels
    out["n_tests"] = len(labels)
    return out


def _run_robustness(spec: ExperimentSpec, opts: Mapping[str, Any], say: Callable[[str], None]) -> dict[str, Any]:
    fn, dotted = analysis_entry_point("robustness")
    family = str(opts.get("specifications", "default_family"))
    if family not in ("default_family", "fast_family"):
        raise SpecError(_suggest(family, ["default_family", "fast_family"], "specification family"))
    say(f"    robustness: {family}")
    kw: dict[str, Any] = {
        "fast": family == "fast_family",
        "seed": spec.seed,
        "conc_M": max(spec.concentrations) if spec.concentrations else 1e-6,
        "n_jobs": int(opts.get("n_jobs", 1)),
    }
    if opts.get("shuffles") is not None:
        kw["n_shuffles"] = int(opts["shuffles"])
    result = _call_filtered(fn, **kw)
    return {
        "analysis": "robustness",
        "entry_point": dotted,
        "settings": {"specifications": family, "seed": spec.seed},
        "result": result,
    }


def _run_uncertainty(spec: ExperimentSpec, opts: Mapping[str, Any], say: Callable[[str], None]) -> dict[str, Any]:
    fn, dotted = analysis_entry_point("uncertainty")
    samples = int(opts.get("samples", 128))
    readout = str(opts.get("readout") or (spec.readouts[0] if spec.readouts else "mean_hz"))
    concs = _analysis_concentrations(opts.get("concentrations"), spec)
    cells: list[dict[str, Any]] = []
    for compound in spec.compounds:
        for conc in concs:
            say(f"    uncertainty: {compound} @ {conc:.2e} M, n_base={samples}")
            result = _call_filtered(
                fn,
                compound=compound,
                conc_M=conc,
                readout=readout,
                n_base=samples,
                seed=spec.seed,
                graph=spec.graph,
                n_boot=int(opts.get("bootstrap", 200)),
            )
            cell = {"compound": compound, "conc_M": conc, "sobol": result}
            try:
                module = importlib.import_module(ANALYSES["uncertainty"]["module"])
                budget = getattr(module, "uncertainty_budget", None)
                if callable(budget):
                    cell["budget"] = budget(result)
            except Exception as exc:  # pragma: no cover - optional summary
                cell["budget_unavailable"] = f"{type(exc).__name__}: {exc}"
            cells.append(cell)
    return {
        "analysis": "uncertainty",
        "entry_point": dotted,
        "settings": {"samples": samples, "readout": readout, "concentrations": concs, "seed": spec.seed},
        "cells": cells,
    }


def _write_figures(spec: ExperimentSpec, rows: Sequence[Mapping[str, Any]], run_dir: Path) -> dict[str, Any]:
    """Dose-response figures, when matplotlib happens to be installed."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        return {
            "available": False,
            "reason": f"matplotlib is not installed ({type(exc).__name__}); every other artifact was written",
            "files": [],
        }
    fig_dir = run_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    for engine in spec.engines:
        for readout in spec.readouts:
            series: dict[str, list[tuple[float, float]]] = {}
            for row in rows:
                if row.get("engine") != engine or row.get("condition") != "drug":
                    continue
                value = row.get(readout)
                if value is None:
                    continue
                try:
                    value = float(value)
                except (TypeError, ValueError):
                    continue
                series.setdefault(str(row.get("compound")), []).append((float(row["conc_M"]), value))
            if not series:
                continue
            fig, ax = plt.subplots(figsize=(4.2, 3.0), dpi=150)
            for compound in sorted(series):
                points = sorted(series[compound])
                ax.plot([p[0] for p in points], [p[1] for p in points], marker="o", label=compound)
            ax.set_xscale("log")
            ax.set_xlabel("free concentration (M)")
            ax.set_ylabel(readout)
            ax.set_title(f"{spec.name}: {readout} ({engine})", fontsize=9)
            ax.legend(fontsize=7)
            fig.tight_layout()
            name = f"{engine}_{readout}.png"
            fig.savefig(fig_dir / name, metadata={"Software": "FlyLab"})
            plt.close(fig)
            written.append(f"figures/{name}")
    return {"available": True, "files": written}


def run_spec(
    spec: ExperimentSpec | Mapping[str, Any] | str | Path,
    outdir: str | Path,
    *,
    dry_run: bool = False,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Run a spec and write a self-describing run directory.

    Returns the run manifest (the same dict that is written to
    ``manifest.json``).  With ``dry_run`` nothing is written and nothing is
    computed: the return value carries the resolved spec and the runtime
    estimate only.

    Determinism: for a given spec and seed, ``results.csv`` is byte-identical
    between runs and every artifact's ``content_sha256`` -- and therefore
    ``determinism.digest`` -- is unchanged.  ``manifest.json`` itself is not,
    because it records wall-clock timings and a timestamp.
    """
    if isinstance(spec, (str, Path)):
        spec = load_spec(spec)
    s = make_spec(spec)
    say = progress or (lambda _msg: None)
    resolved = resolved_spec_dict(s)
    estimate = estimate_runtime(s)

    if dry_run:
        return {
            "dry_run": True,
            "spec": resolved,
            "spec_sha256": spec_sha256(s),
            "runtime_estimate": estimate,
            "outdir": str(Path(outdir)),
            "would_write": sorted(
                ["manifest.json", "spec.yaml", "results.csv", "notebooks/"]
                + ([ANALYSES[a]["artifact"] for a in s.analyses])
                + (["cards/"] if s.cards else [])
                + (["figures/"] if s.figures else [])
            ),
            "warnings": [
                "Nothing was executed: --dry-run validates the spec and prints the plan.",
            ],
        }

    run_dir = Path(outdir)
    run_dir.mkdir(parents=True, exist_ok=True)
    t_start = time.perf_counter()
    timings: dict[str, float] = {}

    # the resolved spec, defaults included ---------------------------------
    import yaml

    (run_dir / "spec.yaml").write_text(
        "# FlyLab resolved experiment spec (defaults filled in).\n"
        f"# spec_sha256: {spec_sha256(s)}\n"
        + yaml.safe_dump(json.loads(json.dumps(resolved, default=str)), sort_keys=False)
    )

    say(f"  assays ({', '.join(s.engines)})")
    t0 = time.perf_counter()
    assay_out = _run_assays(s, run_dir, say)
    timings["assays_s"] = round(time.perf_counter() - t0, 3)

    analysis_results: dict[str, Any] = {}
    analysis_errors: dict[str, str] = {}
    runners = {
        "dependence": _run_dependence,
        "robustness": _run_robustness,
        "uncertainty": _run_uncertainty,
    }
    for name in ("dependence", "robustness", "uncertainty"):
        if name not in s.analyses:
            continue
        say(f"  {name}")
        t0 = time.perf_counter()
        try:
            payload = runners[name](s, s.analyses[name], say)
        except SpecError:
            raise
        except Exception as exc:  # an analysis failure must not lose the assays
            analysis_errors[name] = f"{type(exc).__name__}: {exc}"
            say(f"    {name} failed: {analysis_errors[name]}")
            continue
        analysis_results[name] = payload
        _write_json(run_dir / ANALYSES[name]["artifact"], payload)
        timings[f"{name}_s"] = round(time.perf_counter() - t0, 3)

    # claim cards ----------------------------------------------------------
    cards: list[str] = []
    if s.cards:
        t0 = time.perf_counter()
        from flylab.report.card import claim_card, to_markdown

        card_dir = run_dir / "cards"
        card_dir.mkdir(parents=True, exist_ok=True)
        for (compound, engine), item in sorted(assay_out["headline"].items()):
            card = claim_card(
                item["notebook"],
                compound=compound,
                conc_M=item["conc_M"],
                engine=engine,
                graph=s.graph,
                run_name=s.name,
                vehicle=item.get("vehicle"),
                dependence=analysis_results.get("dependence"),
                robustness=analysis_results.get("robustness"),
                requested_analyses=sorted(s.analyses),
                notebook_path=item["path"],
            )
            stem = f"{compound}__{engine}"
            _write_json(card_dir / f"{stem}.json", card)
            (card_dir / f"{stem}.md").write_text(to_markdown(card))
            cards += [f"cards/{stem}.json", f"cards/{stem}.md"]
        timings["cards_s"] = round(time.perf_counter() - t0, 3)

    # figures --------------------------------------------------------------
    figures = {"available": False, "reason": "figures were switched off in the spec", "files": []}
    if s.figures:
        t0 = time.perf_counter()
        figures = _write_figures(s, assay_out["rows"], run_dir)
        timings["figures_s"] = round(time.perf_counter() - t0, 3)

    timings["total_s"] = round(time.perf_counter() - t_start, 3)

    artifacts = collect_artifacts(run_dir)
    warnings = list(assay_out["warnings"])
    if not figures.get("available"):
        warnings.append(f"no figures were written: {figures.get('reason')}")
    for name, reason in analysis_errors.items():
        warnings.append(f"analysis {name!r} did not complete: {reason}")
    if s.jitter_log10 > 0:
        warnings.append(
            "jitter_log10 > 0: the library Monte-Carlo is active. The run is still "
            "reproducible from the seed, but the readouts are not the unjittered ones."
        )

    manifest: dict[str, Any] = {
        "flylab_run_manifest_version": RUN_MANIFEST_VERSION,
        "kind": "flylab-run",
        "name": s.name,
        "description": s.description,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "spec_sha256": spec_sha256(s),
        "spec": resolved,
        "code": _code_provenance(s.seed),
        **_map_ids(),
        "graphs": {s.graph: _graph_facts(s.graph)},
        "seeds": {
            "global": s.seed,
            "replicate_seeds": [s.seed + r for r in range(s.replicates)],
            "note": "replicate r runs with seed + r; the rate engine is deterministic anyway",
        },
        "engines": assay_out["per_engine"],
        "analyses_requested": sorted(s.analyses),
        "analyses_completed": sorted(analysis_results),
        "analysis_errors": analysis_errors,
        "figures": figures,
        "cards": cards,
        "n_rows": len(assay_out["rows"]),
        "notebooks": assay_out["notebook_paths"],
        "timings_s": timings,
        "runtime_estimate": estimate,
        "artifacts": artifacts,
        "determinism": {
            "digest": determinism_digest(artifacts),
            "volatile_keys": sorted(VOLATILE_KEYS),
            "note": (
                "content_sha256 is taken with the volatile keys above removed. Two runs of "
                "the same spec and seed agree on every content_sha256 and on this digest; "
                "manifest.json itself differs, because it records timings and a timestamp."
            ),
        },
        "warnings": warnings,
        "label": "model_derived",
        "disclaimer": (
            "Every number in this directory is simulation output on a public connectome. "
            "No value here was measured in a living fly."
        ),
    }
    _write_json(run_dir / "manifest.json", manifest)
    say(f"  -> {run_dir}  ({len(artifacts) + 1} files, digest {manifest['determinism']['digest'][:12]})")
    return manifest
