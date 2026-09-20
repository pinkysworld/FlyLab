"""Replicates, confidence intervals, sensitivity and model-derived IC50.

Three sources of replicate-to-replicate variation are supported:

1. the RNG seed (spiking assay only - the rate model is deterministic),
2. the teaching library's EC50s, jittered log-normally
   (``flylab.pharm.uncertainty.sample_library`` when available, otherwise an
   equivalent inline jitter),
3. the drive rate, jittered by +/- ``drive_jitter`` (10% by default).

Every IC50 produced here is the **model's own circuit-level IC50**: the
concentration at which the simulated readout is half way between its low-dose
and high-dose plateau.  It is *not* an animal IC50 and it is not a receptor
IC50; notebooks label it ``model_derived``.
"""
from __future__ import annotations

import contextlib
import copy
import threading
from typing import Any, Iterable, Sequence

import numpy as np

from flylab.analysis.fit import bootstrap_fit, hill4_fit
from flylab.circuit.rate import load_graph, resolve_graph
from flylab.notebook.schema import empty_notebook

ASSAYS = ("subgraph", "spiking", "taste", "taste_map")

#: the scalar readouts the ensemble layer summarises
READOUT_KEYS = ("mn9_hz", "dnp01_hz", "mean_hz", "g_ach", "g_gaba")

DEFAULT_DRIVE = {"subgraph": 40.0, "spiking": 40.0, "taste": 150.0, "taste_map": 150.0}

IC50_WARNING = (
    "IC50 is the MODEL's circuit-level IC50: the concentration at which this "
    "simulated readout is halfway between its low-dose and high-dose plateau. "
    "It is not an animal IC50, not a receptor IC50, and it is entirely a "
    "consequence of the teaching EC50 library and the gain patch rules."
)

_TASTE_LOCK = threading.Lock()


# --------------------------------------------------------------------------
# library jitter
# --------------------------------------------------------------------------
def sample_library(rng: np.random.Generator, sd_log10: float, library: dict | None = None) -> dict:
    """Log-normally jittered copy of the teaching library.

    Delegates to ``flylab.pharm.uncertainty.sample_library`` when that module
    exists (it is the pharm agent's single source of truth); the inline
    fallback multiplies every ``ec50_M`` by ``10 ** N(0, sd_log10)``.
    """
    if sd_log10 <= 0:
        return copy.deepcopy(library) if library else None
    try:  # pragma: no cover - depends on the concurrently built pharm module
        from flylab.pharm.uncertainty import sample_library as _sample

        out = _sample(rng, sd_log10)
        if out:
            return out
    except Exception:
        pass
    from flylab.pharm.occupancy import load_library

    lib = copy.deepcopy(library or load_library())
    for comp in lib.get("compounds", {}).values():
        for spec in comp.get("receptors", {}).values():
            try:
                spec["ec50_M"] = float(spec["ec50_M"]) * float(10.0 ** rng.normal(0.0, sd_log10))
            except Exception:
                continue
    return lib


@contextlib.contextmanager
def _taste_library(library: dict | None):
    """Make ``assays/taste.py`` read a jittered library without editing it."""
    if library is None:
        yield
        return
    import flylab.assays.taste as taste_mod

    with _TASTE_LOCK:
        orig = taste_mod.compare_compound
        taste_mod.compare_compound = lambda name, conc, **kw: orig(name, conc, library=library)
        try:
            yield
        finally:
            taste_mod.compare_compound = orig


# --------------------------------------------------------------------------
# one replicate
# --------------------------------------------------------------------------
def run_assay(
    assay: str,
    compound: str | None,
    conc_M: float,
    *,
    library: dict | None = None,
    drive_hz: float | None = None,
    seed: int = 0,
    rule_overrides: dict[str, float] | None = None,
    graph_obj: dict | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Run one assay replicate and return its notebook."""
    if assay not in ASSAYS:
        raise ValueError(f"unknown assay {assay!r}; expected one of {ASSAYS}")
    drive = float(drive_hz if drive_hz is not None else DEFAULT_DRIVE[assay])
    if assay == "subgraph":
        from flylab.assays.subgraph import run_subgraph_assay

        return run_subgraph_assay(
            compound, conc_M, drive_hz=drive, library=library,
            rule_overrides=rule_overrides, graph_obj=graph_obj, **kw,
        )
    if assay == "spiking":
        from flylab.assays.spiking import run_spiking_assay

        return run_spiking_assay(
            compound, conc_M, drive_hz=drive, seed=seed, library=library,
            rule_overrides=rule_overrides, graph_obj=graph_obj, **kw,
        )
    if assay == "taste_map":
        from flylab.assays.taste_map import run_taste_map_assay

        return run_taste_map_assay(
            compound, conc_M, sugar_hz=drive, seed=seed, library=library,
            rule_overrides=rule_overrides, **kw,
        )
    from flylab.assays.taste import run_taste_assay

    with _taste_library(library):
        return run_taste_assay(compound, conc_M, sugar_hz=drive, **kw)


def readouts_of(assay: str, nb: dict[str, Any]) -> dict[str, float | None]:
    """The five comparable scalar readouts, whatever the assay."""
    r = nb.get("readouts", {})
    gains = nb.get("gains", {}) or {}
    if assay == "taste":
        return {
            "mn9_hz": r.get("mn9_sugar_hz"),
            "dnp01_hz": None,
            "mean_hz": r.get("mn9_sugar_hz"),
            "g_ach": r.get("g_ach", gains.get("g_ach")),
            "g_gaba": gains.get("g_gaba", 1.0),
        }
    return {k: r.get(k, gains.get(k)) for k in READOUT_KEYS}


# --------------------------------------------------------------------------
# ensemble
# --------------------------------------------------------------------------
def _summarise(values: dict[str, list[float]], ci: tuple[float, float]) -> dict[str, Any]:
    out_ci, out_mean, out_sd = {}, {}, {}
    for k, vals in values.items():
        v = np.array([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
        if v.size == 0:
            continue
        out_mean[k] = float(v.mean())
        out_sd[k] = float(v.std(ddof=1)) if v.size > 1 else 0.0
        lo, hi = np.percentile(v, ci) if v.size > 1 else (float(v[0]), float(v[0]))
        out_ci[k] = [float(lo), float(hi)]
    return {"ci": out_ci, "mean": out_mean, "sd": out_sd}


def run_ensemble(
    assay: str,
    compound: str | None = None,
    conc_M: float = 0.0,
    n_rep: int = 8,
    ec50_sd_log10: float = 0.3,
    seed: int = 0,
    drive_jitter: float = 0.10,
    ci: tuple[float, float] = (2.5, 97.5),
    **kw: Any,
) -> dict[str, Any]:
    """Replicate an assay and attach an ``uncertainty`` block to its notebook.

    The notebook body is the *unjittered* run (so its headline numbers match a
    plain assay call); ``uncertainty`` summarises ``n_rep`` jittered replicates.
    """
    base_drive = float(kw.pop("drive_hz", DEFAULT_DRIVE[assay]))
    nb = run_assay(assay, compound, conc_M, drive_hz=base_drive, seed=seed, **kw)

    reps: list[dict[str, Any]] = []
    values: dict[str, list[float]] = {k: [] for k in READOUT_KEYS}
    for i in range(int(n_rep)):
        rng = np.random.default_rng([int(seed), i])
        lib = sample_library(rng, ec50_sd_log10) if (compound and ec50_sd_log10 > 0) else None
        drive = base_drive * float(1.0 + drive_jitter * rng.uniform(-1.0, 1.0))
        rep_nb = run_assay(
            assay, compound, conc_M, library=lib, drive_hz=drive, seed=int(seed) + i, **kw
        )
        vals = readouts_of(assay, rep_nb)
        reps.append({"replicate": i, "seed": int(seed) + i, "drive_hz": drive, **vals})
        for k in READOUT_KEYS:
            values[k].append(vals.get(k))

    summary = _summarise(values, ci)
    nb.setdefault("gains", {})
    nb["uncertainty"] = {
        "n_rep": int(n_rep),
        "method": (
            f"Monte-Carlo over the teaching library (log10 EC50 sd={ec50_sd_log10}), "
            f"+/-{int(drive_jitter * 100)}% drive jitter and the RNG seed; "
            f"percentile CI {ci[0]}-{ci[1]}%."
        ),
        "ci_percentiles": list(ci),
        "ec50_sd_log10": float(ec50_sd_log10),
        "drive_jitter": float(drive_jitter),
        "seed": int(seed),
        "ci": summary["ci"],
        "mean": summary["mean"],
        "sd": summary["sd"],
        "replicates": reps,
        "label": "model_derived",
    }
    nb["warnings"] = list(nb.get("warnings", [])) + [
        "Uncertainty is Monte-Carlo over teaching EC50s and drive, not a "
        "measurement error bar; it says how sensitive the model is, not how "
        "variable a fly is."
    ]
    return nb


# --------------------------------------------------------------------------
# sensitivity
# --------------------------------------------------------------------------
SENSITIVITY_PARAMS = ("ec50", "hill_n", "drive_hz", "gain_coef", "weight_threshold")


def _dominant_insect_receptor(compound: str, conc_M: float) -> str | None:
    """Insect receptor with the highest occupancy and a real direction."""
    from flylab.pharm.occupancy import compare_compound

    rows = compare_compound(compound, conc_M)["receptors"]
    cand = [
        r for r in rows
        if str(r["receptor"]).startswith("insect_") and r.get("direction") not in (None, "none", "unknown")
    ]
    if not cand:
        cand = [r for r in rows if str(r["receptor"]).startswith("insect_")]
    if not cand:
        return None
    return max(cand, key=lambda r: r["occupancy"])["receptor"]


def _library_with(receptor: str, field: str, factor: float) -> dict:
    from flylab.pharm.occupancy import load_library

    lib = copy.deepcopy(load_library())
    for comp in lib.get("compounds", {}).values():
        spec = comp.get("receptors", {}).get(receptor)
        if spec is None:
            continue
        spec[field] = float(spec.get(field, 1.0)) * float(factor)
    return lib


def _threshold_graph(graph: dict[str, Any], min_weight: float) -> dict[str, Any]:
    out = dict(graph)
    out["edges"] = [e for e in graph["edges"] if float(e["weight"]) >= float(min_weight)]
    out["n_edges"] = len(out["edges"])
    return out


def sensitivity(
    assay: str = "subgraph",
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    params: Sequence[str] = SENSITIVITY_PARAMS,
    factor: float = 2.0,
    readout: str = "mn9_hz",
    seed: int = 0,
    graph: str | None = None,
    **kw: Any,
) -> list[dict[str, Any]]:
    """One-at-a-time +/- ``factor`` tornado table for one readout.

    Parameters probed:

    ``ec50``              EC50 of the dominant insect receptor
    ``hill_n``            Hill coefficient of that receptor
    ``drive_hz``          input drive
    ``gain_coef``         the 0.4 / 1.6 coefficients of the nAChR agonist rule,
                          applied through ``rule_overrides`` so the default rule
                          in ``pharm/mechanisms.py`` is never edited
    ``weight_threshold``  minimum synapse count kept in the graph
    """
    base_drive = float(kw.pop("drive_hz", DEFAULT_DRIVE[assay]))
    base_nb = run_assay(assay, compound, conc_M, drive_hz=base_drive, seed=seed, **kw)
    base = readouts_of(assay, base_nb).get(readout)

    receptor = _dominant_insect_receptor(compound, conc_M)
    g = load_graph(resolve_graph(graph)) if assay in ("subgraph", "spiking") else None
    base_thr = float((g or {}).get("min_weight", 5.0) or 5.0)

    rows: list[dict[str, Any]] = []
    for p in params:
        low_kw: dict[str, Any] = {}
        high_kw: dict[str, Any] = {}
        low_set: Any
        high_set: Any
        if p == "ec50":
            if receptor is None:
                continue
            low_kw = {"library": _library_with(receptor, "ec50_M", 1.0 / factor)}
            high_kw = {"library": _library_with(receptor, "ec50_M", factor)}
            low_set, high_set = f"{receptor} EC50 / {factor}", f"{receptor} EC50 x {factor}"
        elif p == "hill_n":
            if receptor is None:
                continue
            low_kw = {"library": _library_with(receptor, "n", 1.0 / factor)}
            high_kw = {"library": _library_with(receptor, "n", factor)}
            low_set, high_set = f"{receptor} n / {factor}", f"{receptor} n x {factor}"
        elif p == "drive_hz":
            low_kw = {"drive_hz": base_drive / factor}
            high_kw = {"drive_hz": base_drive * factor}
            low_set, high_set = base_drive / factor, base_drive * factor
        elif p == "gain_coef":
            low_kw = {"rule_overrides": {"agonist_a": 0.4 / factor, "agonist_b": 1.6 / factor}}
            high_kw = {"rule_overrides": {"agonist_a": 0.4 * factor, "agonist_b": 1.6 * factor}}
            low_set, high_set = "a=0.4/f, b=1.6/f", "a=0.4*f, b=1.6*f"
        elif p == "weight_threshold":
            if g is None:
                continue
            low_kw = {"graph_obj": _threshold_graph(g, base_thr)}
            high_kw = {"graph_obj": _threshold_graph(g, base_thr * factor)}
            low_set, high_set = base_thr, base_thr * factor
        else:
            raise ValueError(f"unknown sensitivity parameter {p!r}")

        lo_nb = run_assay(assay, compound, conc_M, seed=seed, **{"drive_hz": base_drive, **kw, **low_kw})
        hi_nb = run_assay(assay, compound, conc_M, seed=seed, **{"drive_hz": base_drive, **kw, **high_kw})
        lo = readouts_of(assay, lo_nb).get(readout)
        hi = readouts_of(assay, hi_nb).get(readout)
        rows.append(
            {
                "param": p,
                "low": lo,
                "high": hi,
                "base": base,
                "readout": readout,
                "factor": float(factor),
                "low_setting": low_set,
                "high_setting": high_set,
                "span": (
                    abs((hi if hi is not None else 0.0) - (lo if lo is not None else 0.0))
                    if (lo is not None and hi is not None)
                    else None
                ),
            }
        )
    rows.sort(key=lambda r: -(r["span"] if r["span"] is not None else -1.0))
    return rows


# --------------------------------------------------------------------------
# circuit IC50
# --------------------------------------------------------------------------
DEFAULT_CONCS = [1e-10, 1e-9, 1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 1e-5, 1e-4]

#: readout that actually moves with dose, per assay.  On the seed-driven
#: ``subgraph`` rate assay MN9 is *clamped* by its own 40 Hz drive, so its rate
#: barely moves; the network mean is the informative readout there.  In the LIF
#: assay the drive is sub-threshold Poisson, so MN9 is free to move.
DEFAULT_READOUT = {"subgraph": "mean_hz", "spiking": "mn9_hz", "taste": "mn9_hz", "taste_map": "mn9_hz"}


def circuit_ic50(
    assay: str = "subgraph",
    compound: str = "imidacloprid",
    readout: str = "mn9_hz",
    concs: Iterable[float] | None = None,
    n_boot: int = 200,
    seed: int = 0,
    n_rep: int = 4,
    ec50_sd_log10: float = 0.3,
    drive_jitter: float = 0.10,
    **kw: Any,
) -> dict[str, Any]:
    """Dose ladder -> 4-parameter Hill fit -> bootstrap CI on the model IC50.

    ``readout="auto"`` picks :data:`DEFAULT_READOUT` for the assay.  Note that
    on the ``subgraph`` assay MN9 is driven at ``drive_hz`` itself, so
    ``mn9_hz`` is nearly dose-independent there and ``mean_hz`` is the readout
    that carries the dose-response.
    """
    if readout == "auto":
        readout = DEFAULT_READOUT.get(assay, "mn9_hz")
    xs = [float(c) for c in (concs if concs is not None else DEFAULT_CONCS)]
    base_drive = float(kw.pop("drive_hz", DEFAULT_DRIVE[assay]))

    points: list[dict[str, Any]] = []
    y_base: list[float] = []
    y_reps: list[list[float]] = []
    for ci_, c in enumerate(xs):
        nb = run_assay(assay, compound, c, drive_hz=base_drive, seed=seed, **kw)
        val = readouts_of(assay, nb).get(readout)
        rep_vals: list[float] = []
        for i in range(int(n_rep)):
            rng = np.random.default_rng([int(seed), ci_, i])
            lib = sample_library(rng, ec50_sd_log10) if ec50_sd_log10 > 0 else None
            drive = base_drive * float(1.0 + drive_jitter * rng.uniform(-1.0, 1.0))
            rnb = run_assay(
                assay, compound, c, library=lib, drive_hz=drive, seed=int(seed) + i, **kw
            )
            v = readouts_of(assay, rnb).get(readout)
            rep_vals.append(float(v) if v is not None else float("nan"))
        y_base.append(float(val) if val is not None else float("nan"))
        y_reps.append(rep_vals)
        arr = np.array([v for v in rep_vals if np.isfinite(v)], dtype=float)
        points.append(
            {
                "conc_M": c,
                "value": float(val) if val is not None else None,
                "replicates": rep_vals,
                "mean": float(arr.mean()) if arr.size else None,
                "sd": float(arr.std(ddof=1)) if arr.size > 1 else 0.0,
            }
        )

    warnings = [IC50_WARNING, "Hops-limited MaleCNS neighborhood; teaching EC50 library."]
    boot: dict[str, Any] = {"n_boot": int(n_boot), "n_ok": 0}
    try:
        b = bootstrap_fit(xs, y_reps, n_boot=int(n_boot), seed=int(seed))
        fit = hill4_fit(xs, y_base)
        ci_block = {k: v for k, v in (("ic50", b["ci"].get("ic50")), ("slope", b["ci"].get("slope"))) if v}
        boot = {"n_boot": int(n_boot), "n_ok": b["n_ok"], "method": b["method"]}
    except ValueError as exc:
        fit = {"error": str(exc)}
        ci_block = {}
        warnings.append(f"fit failed: {exc}")
    if fit.get("note"):
        warnings.append(fit["note"])
    if isinstance(fit.get("r2"), float) and fit["r2"] < 0.8:
        warnings.append(
            f"Hill fit is poor (r2={fit['r2']:.2f}) for readout {readout!r}; the "
            "dose-response may be non-monotonic or flat. On the seed-driven "
            "subgraph assay MN9 is clamped by its own drive - try 'mean_hz'."
        )
    if fit.get("in_range") is False:
        warnings.append(
            "Fitted IC50 lies outside the tested concentration range; it is an "
            "extrapolation, not a measurement of the ladder."
        )

    return {
        "assay": assay,
        "compound": compound,
        "readout": readout,
        "concs_M": xs,
        "points": points,
        "fit": fit,
        "ci": ci_block,
        "bootstrap": boot,
        "label": "model_derived",
        "warnings": warnings,
    }


__all__ = [
    "run_ensemble",
    "sensitivity",
    "circuit_ic50",
    "run_assay",
    "readouts_of",
    "sample_library",
    "READOUT_KEYS",
    "DEFAULT_READOUT",
    "ASSAYS",
    "SENSITIVITY_PARAMS",
    "IC50_WARNING",
]
