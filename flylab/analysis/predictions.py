"""Pre-registered predictions and the sample size they would need.

Every row this module produces is a **software prediction**: what the FlyLab
model says should happen if the hypothesis is tested in a real fly, written
down *before* any live table exists.  ``live_result`` is always ``None`` here;
it is filled in by hand from a real experiment, never by this code
(``HANDOFF.md``: "do not invent live PER/climbing numbers").

The effect sizes are ``treated - vehicle`` on a model readout, replicated with
the ensemble layer's Monte-Carlo (teaching EC50 jitter + drive jitter + RNG
seed).  The resulting ``d = mean / sd`` is therefore **model-internal
variability**: it says how stable this simulation is under its own parameter
uncertainty.  It is *not* biological variance and it must never be reported as
one, so the suggested sample sizes cap ``d`` at a large-but-plausible
biological effect before doing the power calculation.

Hypotheses H1-H5 are the ones in ``docs/RESEARCH_MAP.md``.  H6 and H7 are new
here: H6 is the map-path prediction the ``taste_motor`` extract made possible,
H7 is the non-selective negative control that the paper needs in order to
claim anything about selectivity at all.
"""
from __future__ import annotations

import math
import time
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

__all__ = [
    "HYPOTHESES",
    "hypothesis",
    "predict",
    "power_for_per",
    "power_for_continuous",
    "prediction_table",
    "to_markdown",
    "behavioral_baseline",
    "DEFAULT_BASELINE_PER",
    "DEFAULT_BASELINE_SD",
]

#: documented fallbacks when no literature YAML is available
DEFAULT_BASELINE_PER = 0.8
DEFAULT_BASELINE_SD = 0.15

#: cap on the model-internal effect size used for power (Cohen's "large")
D_CAP = 1.0

LIT_PATHS = [
    Path("data/literature/behavioral_assays.yaml"),
    Path(__file__).resolve().parents[2] / "data/literature/behavioral_assays.yaml",
]

BASE_WARNINGS = [
    "status software_prediction: these are model predictions written before any "
    "live experiment, not results. live_result stays null until a real table is "
    "imported by hand.",
    "d is model-internal variability (teaching EC50 jitter, drive jitter, RNG "
    "seed). It is NOT biological variance and the model has no access to any.",
    "model_derived: the effect sizes come from the teaching EC50 library, the "
    "gain patch rules and one hops-limited MaleCNS neighborhood.",
]


# --------------------------------------------------------------------------
# the pre-registered hypotheses
# --------------------------------------------------------------------------
def _h(**kw: Any) -> dict[str, Any]:
    row = {
        "assay": "subgraph",
        "compound": None,
        "conc_M": 1e-6,
        "readout": "mean_hz",
        "contrast": "treated_minus_vehicle",
        "direction": "decrease",
        "endpoint": "continuous",
        "graph": None,
        "source": "docs/RESEARCH_MAP.md",
        "live_protocol": None,
    }
    row.update(kw)
    return row


HYPOTHESES: dict[str, dict[str, Any]] = {
    "H1": _h(
        statement=(
            "At 1 uM, imidacloprid occupancy on insect nAChR is >= 0.9 while "
            "vertebrate a4b2 occupancy stays <= 0.2 (dual scorecard gap > 0)."
        ),
        assay="occupancy",
        compound="imidacloprid",
        conc_M=1e-6,
        readout="occupancy_gap_insect_minus_vertebrate",
        contrast="occupancy_gap",
        direction="increase",
        endpoint="receptor_occupancy",
        live_protocol=None,
    ),
    "H2": _h(
        statement=(
            "The same dose lowers model drive through ACh edges in the "
            "MN9/DNp01 neighborhood: g_ach falls and the network mean rate "
            "falls with it."
        ),
        compound="imidacloprid",
        conc_M=1e-6,
        readout="mean_hz",
        direction="decrease",
        endpoint="per",
        live_protocol="protocols/per.md",
    ),
    "H3": _h(
        statement=(
            "Fipronil at 1 uM lowers g_gaba (RDL block) and disinhibits the "
            "neighborhood: the network mean rate rises above vehicle."
        ),
        compound="fipronil",
        conc_M=1e-6,
        readout="mean_hz",
        direction="increase",
        endpoint="climbing",
        live_protocol="protocols/climbing.md",
    ),
    "H4": _h(
        statement=(
            "Diazepam occupies vertebrate GABA-A and does not occupy insect "
            "RDL in this library, so it moves no fly readout: the predicted "
            "circuit effect is exactly zero."
        ),
        compound="diazepam",
        conc_M=1e-6,
        readout="mean_hz",
        direction="none",
        endpoint="per",
        live_protocol="protocols/per.md",
    ),
    "H5": _h(
        statement=(
            "Bitter input vetoes sugar-driven MN9 on the map-extracted "
            "labellar-GRN path (Shiu 2024 direction): adding bitter drive "
            "lowers MN9 below the sugar-only rate in vehicle."
        ),
        assay="taste_map",
        compound=None,
        conc_M=0.0,
        readout="mn9_sugar_bitter_hz - mn9_sugar_hz",
        contrast="within",
        direction="decrease",
        endpoint="per",
        graph="taste_motor",
        live_protocol="protocols/per.md",
    ),
    "H6": _h(
        statement=(
            "Fipronil at 1 uM relieves the bitter veto on the map path: the "
            "bitter veto ratio (MN9 sugar+bitter / MN9 sugar) rises above its "
            "vehicle value because RDL block removes the inhibition that "
            "carries the veto."
        ),
        assay="taste_map",
        compound="fipronil",
        conc_M=1e-6,
        readout="bitter_veto_ratio",
        direction="increase",
        endpoint="per",
        graph="taste_motor",
        source="FlyLab v0.5, new prediction on the taste_motor extract",
        live_protocol="protocols/per.md",
    ),
    "H7": _h(
        statement=(
            "Picrotoxin is the non-selective control: it perturbs the fly "
            "neighborhood (RDL block raises the mean rate) while its receptor "
            "selectivity index is ~0, so any insect-selectivity claim must "
            "fail for it."
        ),
        compound="picrotoxin",
        conc_M=1e-6,
        readout="mean_hz",
        direction="increase",
        endpoint="climbing",
        source="FlyLab v0.5, negative control for the selectivity claim",
        live_protocol="protocols/climbing.md",
    ),
}
for _k, _v in HYPOTHESES.items():
    _v["h_id"] = _k


def hypothesis(h_id: str) -> dict[str, Any]:
    """One pre-registered hypothesis by id (``"H1"`` ... ``"H7"``)."""
    key = str(h_id).upper().strip()
    if key not in HYPOTHESES:
        raise KeyError(f"unknown hypothesis {h_id!r}; known: {', '.join(HYPOTHESES)}")
    return dict(HYPOTHESES[key])


# --------------------------------------------------------------------------
# literature baseline (written by another agent; schema unknown)
# --------------------------------------------------------------------------
def _walk(obj: Any, path: tuple[str, ...] = ()):  # pragma: no cover - trivial
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield path + (str(k),), v
            yield from _walk(v, path + (str(k),))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            yield path + (str(i),), v
            yield from _walk(v, path + (str(i),))


def _as_proportion(v: Any) -> float | None:
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(x) or x < 0:
        return None
    if 0.0 < x <= 1.0:
        return x
    if 1.0 < x <= 100.0:
        return x / 100.0
    return None


def behavioral_baseline(path: str | Path | None = None) -> dict[str, Any]:
    """Baseline PER statistics, read *defensively* from the literature YAML.

    ``data/literature/behavioral_assays.yaml`` is written by another agent and
    its schema is not fixed, so this walks whatever it finds and accepts the
    first plausible PER baseline proportion (and SD) it can recognise, treating
    values in ``(1, 100]`` as percentages.  When nothing usable is present it
    falls back to the documented defaults (PER 0.8, SD 0.15) and says so in
    ``warnings``.
    """
    warnings: list[str] = []
    candidates = [Path(path)] if path else list(LIT_PATHS)
    found_path = next((p for p in candidates if p.exists()), None)
    out = {
        "baseline_per": float(DEFAULT_BASELINE_PER),
        "sd": float(DEFAULT_BASELINE_SD),
        "source": "documented default (HANDOFF/protocols): baseline PER 0.8, SD 0.15",
        "path": None,
        "loaded": False,
        "label": "model_derived",
        "warnings": warnings,
    }
    if found_path is None:
        warnings.append(
            "data/literature/behavioral_assays.yaml not found; using the "
            "documented defaults baseline PER 0.8, SD 0.15."
        )
        return out
    out["path"] = str(found_path)
    try:
        import yaml

        data = yaml.safe_load(found_path.read_text())
    except Exception as exc:  # pragma: no cover - defensive
        warnings.append(f"could not parse {found_path}: {exc}; using documented defaults.")
        return out
    if not isinstance(data, (dict, list)):
        warnings.append(f"{found_path} has no mapping at its root; using documented defaults.")
        return out

    per_val: float | None = None
    per_key: str | None = None
    sd_val: float | None = None
    sd_key: str | None = None
    candidates: list[dict[str, Any]] = []
    for keypath, value in _walk(data):
        leaf = keypath[-1].lower()
        ctx = " ".join(k.lower() for k in keypath)
        if "per" not in ctx and "proboscis" not in ctx:
            continue
        v = _as_proportion(value)
        if v is not None and len(candidates) < 12 and any(
            t in leaf for t in ("rate", "fraction", "proportion", "probability", "responder", "baseline", "p0")
        ):
            candidates.append({"key": ".".join(keypath), "value": v})
        if per_val is None and any(
            t in leaf for t in ("baseline", "p0", "control", "vehicle", "proportion", "response_rate", "rate", "mean")
        ) and "non_respond" not in leaf and "non-respond" not in leaf:
            if v is not None:
                per_val, per_key = v, ".".join(keypath)
        if sd_val is None and any(t in leaf for t in ("sd", "std", "stdev", "sigma")):
            if v is not None:
                sd_val, sd_key = v, ".".join(keypath)
    out["candidates"] = candidates
    if per_val is not None:
        out["baseline_per"] = per_val
        out["loaded"] = True
        out["source"] = f"{found_path}:{per_key}"
    else:
        msg = (
            f"{found_path} exists but no PER baseline proportion was "
            f"recognised; using the documented default {DEFAULT_BASELINE_PER}."
        )
        if candidates:
            msg += (
                " Proportions that were found but not used: "
                + ", ".join(f"{c['key']}={c['value']:.3g}" for c in candidates[:6])
                + " -- pass baseline_p explicitly to use one of them."
            )
        warnings.append(msg)
    if sd_val is not None:
        out["sd"] = sd_val
        out["sd_source"] = f"{found_path}:{sd_key}"
    else:
        warnings.append("no PER SD recognised; using the documented default 0.15.")
    return out


# --------------------------------------------------------------------------
# power
# --------------------------------------------------------------------------
def _z(p: float) -> float:
    """Inverse standard-normal CDF by bisection on ``math.erf`` (no scipy)."""
    p = min(max(float(p), 1e-12), 1 - 1e-12)
    lo, hi = -12.0, 12.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if 0.5 * (1.0 + math.erf(mid / math.sqrt(2.0))) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def power_for_per(
    effect_frac: float,
    baseline_p: float | None = None,
    alpha: float = 0.05,
    power: float = 0.8,
    baseline: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Two-proportion z-test sample size for a PER experiment.

    PER is binary per fly, so a predicted *relative* drop ``effect_frac`` in the
    response proportion becomes ``p2 = baseline_p * (1 - effect_frac)`` and

    ``n = (z_(1-a/2) * sqrt(2 p_bar q_bar) + z_power * sqrt(p1 q1 + p2 q2))^2 / (p1 - p2)^2``

    per group (normal approximation, no continuity correction, no scipy).
    ``baseline_p`` defaults to whatever :func:`behavioral_baseline` finds.
    """
    base = baseline or behavioral_baseline()
    p1 = float(baseline_p if baseline_p is not None else base["baseline_per"])
    warnings = list(base.get("warnings", []))
    frac = float(effect_frac)
    p1 = min(max(p1, 1e-6), 1 - 1e-6)
    p2 = min(max(p1 * (1.0 - frac), 1e-6), 1 - 1e-6)
    delta = p1 - p2
    if abs(delta) < 1e-9:
        warnings.append("no predicted difference in PER: the sample size is undefined.")
        return {
            "n_per_group": None, "p1": p1, "p2": p2, "effect_frac": frac,
            "alpha": float(alpha), "power": float(power),
            "baseline_source": base.get("source"),
            "method": "two-proportion z-test (normal approximation)",
            "label": "model_derived", "warnings": warnings,
        }
    za = _z(1.0 - float(alpha) / 2.0)
    zb = _z(float(power))
    pbar = 0.5 * (p1 + p2)
    n = ((za * math.sqrt(2.0 * pbar * (1.0 - pbar)) + zb * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2))) ** 2) / (delta ** 2)
    warnings.append(
        "Sample size uses the MODEL's predicted effect as if it were the "
        "biological effect; a real power calculation needs a pilot PER table."
    )
    return {
        "n_per_group": int(math.ceil(n)),
        "n_exact": float(n),
        "p1": p1,
        "p2": p2,
        "effect_frac": frac,
        "alpha": float(alpha),
        "power": float(power),
        "baseline_source": base.get("source"),
        "baseline_loaded": bool(base.get("loaded")),
        "method": "two-proportion z-test (normal approximation, no continuity correction)",
        "label": "model_derived",
        "warnings": warnings,
    }


def power_for_continuous(
    d: float | None,
    alpha: float = 0.05,
    power: float = 0.8,
    d_cap: float = D_CAP,
) -> dict[str, Any]:
    """Two-sample t-test sample size, normal approximation: ``n = 2 (za+zb)^2 / d^2``.

    ``d`` from :func:`predict` is model-internal, so it is capped at ``d_cap``
    (default 1.0, a "large" biological effect) before the calculation; without
    the cap a near-deterministic simulation would suggest n = 2.
    """
    warnings: list[str] = []
    if d is None or not math.isfinite(float(d)) or abs(float(d)) < 1e-9:
        warnings.append("effect size is zero or undefined: the sample size is undefined.")
        return {
            "n_per_group": None, "d": d, "d_used": None, "alpha": float(alpha),
            "power": float(power), "method": "two-sample t (normal approximation)",
            "label": "model_derived", "warnings": warnings,
        }
    d_raw = abs(float(d))
    d_used = min(d_raw, float(d_cap))
    if d_used < d_raw:
        warnings.append(
            f"model-internal d={d_raw:.2f} was capped at {d_cap:.2f} for the "
            "power calculation: model stability is not biological effect size."
        )
    za = _z(1.0 - float(alpha) / 2.0)
    zb = _z(float(power))
    n = 2.0 * (za + zb) ** 2 / (d_used ** 2)
    return {
        "n_per_group": int(math.ceil(n)),
        "n_exact": float(n),
        "d": float(d),
        "d_used": float(d_used),
        "alpha": float(alpha),
        "power": float(power),
        "method": "two-sample t-test, normal approximation n = 2 (z_a + z_b)^2 / d^2",
        "label": "model_derived",
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# prediction
# --------------------------------------------------------------------------
def _plan(n_rep: int, seed: int, compound: str | None, ec50_sd_log10: float,
          drive_jitter: float, base_drive: float) -> list[dict[str, Any]]:
    """The replicate plan ``run_ensemble`` uses, reproduced so the vehicle arm
    can be run at exactly the same jittered drive."""
    from flylab.assays.ensemble import sample_library

    out = []
    for i in range(int(n_rep)):
        rng = np.random.default_rng([int(seed), i])
        lib = sample_library(rng, ec50_sd_log10) if (compound and ec50_sd_log10 > 0) else None
        drive = float(base_drive) * float(1.0 + drive_jitter * rng.uniform(-1.0, 1.0))
        out.append({"replicate": i, "seed": int(seed) + i, "library": lib, "drive_hz": drive})
    return out


def _value(nb: dict[str, Any], readout: str) -> float | None:
    from flylab.analysis.nullmodels import readout_value

    return readout_value(nb, readout)


def _occupancy_gap(compound: str, conc_M: float, library: dict[str, Any] | None) -> float | None:
    """insect - vertebrate occupancy for the compound's best receptor pair."""
    from flylab.analysis.selectivity import _best_pair, _compound_rows

    occ = _compound_rows(compound, conc_M, library)
    _, pair = _best_pair(occ)
    if not pair:
        return None
    return float(pair["insect_occupancy"] - pair["vertebrate_occupancy"])


def _ci(values: np.ndarray) -> list[float] | None:
    if values.size == 0:
        return None
    if values.size == 1:
        return [float(values[0]), float(values[0])]
    sd = float(values.std(ddof=1))
    mean = float(values.mean())
    half = 1.959963985 * sd / math.sqrt(values.size)
    return [mean - half, mean + half]


def predict(
    h_id: str,
    n_rep: int = 8,
    seed: int = 0,
    ec50_sd_log10: float = 0.3,
    drive_jitter: float = 0.10,
    **kw: Any,
) -> dict[str, Any]:
    """Predicted effect, 95% CI and standardised effect size for one hypothesis.

    The treated arm is :func:`flylab.assays.ensemble.run_ensemble` (teaching
    EC50 jitter + drive jitter + RNG seed); the vehicle arm is re-run at each
    replicate's *own* jittered drive so the difference is paired.  ``d`` is
    ``mean / sd`` across replicates and is labelled ``model_internal``: it
    measures the stability of the simulation, not the variability of a fly.
    """
    from flylab.assays.ensemble import DEFAULT_DRIVE, run_assay, run_ensemble

    h = hypothesis(h_id)
    t0 = time.perf_counter()
    assay, compound, conc = h["assay"], h["compound"], float(h["conc_M"])
    readout, contrast = h["readout"], h["contrast"]
    warnings = list(BASE_WARNINGS)
    kw = dict(kw)
    if h.get("graph") and assay != "occupancy":
        kw.setdefault("graph", h["graph"])

    effects: list[float] = []
    treated_vals: list[float] = []
    vehicle_vals: list[float] = []
    ensemble_block: dict[str, Any] | None = None

    if contrast == "occupancy_gap":
        from flylab.assays.ensemble import sample_library

        for i in range(int(n_rep)):
            rng = np.random.default_rng([int(seed), i])
            lib = sample_library(rng, ec50_sd_log10) if ec50_sd_log10 > 0 else None
            gap = _occupancy_gap(compound, conc, lib)
            if gap is None:
                continue
            effects.append(gap)
            treated_vals.append(gap)
        warnings.append(
            "H1 is a receptor-level prediction: its replicates are Monte-Carlo "
            "over the teaching log10 EC50s, with no circuit involved."
        )
    else:
        base_drive = float(kw.pop("drive_hz", DEFAULT_DRIVE[assay]))
        plan = _plan(n_rep, seed, compound, ec50_sd_log10, drive_jitter, base_drive)
        ens = run_ensemble(
            assay, compound, conc, n_rep=int(n_rep), seed=int(seed),
            ec50_sd_log10=ec50_sd_log10, drive_jitter=drive_jitter,
            drive_hz=base_drive, **kw,
        )
        ensemble_block = ens.get("uncertainty")
        for w in ens.get("warnings", []):
            if w not in warnings:
                warnings.append(w)
        reps = (ensemble_block or {}).get("replicates") or []
        for i, item in enumerate(plan):
            rep = reps[i] if i < len(reps) else {}
            if contrast == "within":
                nb = run_assay(
                    assay, compound, conc, library=item["library"],
                    drive_hz=item["drive_hz"], seed=item["seed"], **kw,
                )
                a, b = [s.strip() for s in readout.split("-")]
                va, vb = _value(nb, a), _value(nb, b)
                if va is None or vb is None:
                    continue
                effects.append(float(va - vb))
                treated_vals.append(float(va))
                vehicle_vals.append(float(vb))
                continue
            t_val = rep.get(readout) if readout in rep else None
            if t_val is None:
                nb = run_assay(
                    assay, compound, conc, library=item["library"],
                    drive_hz=item["drive_hz"], seed=item["seed"], **kw,
                )
                t_val = _value(nb, readout)
            veh_nb = run_assay(
                assay, None, 0.0, drive_hz=item["drive_hz"], seed=item["seed"], **kw
            )
            v_val = _value(veh_nb, readout)
            if t_val is None or v_val is None:
                continue
            treated_vals.append(float(t_val))
            vehicle_vals.append(float(v_val))
            effects.append(float(t_val) - float(v_val))

    arr = np.array([e for e in effects if np.isfinite(e)], dtype=float)
    mean = float(arr.mean()) if arr.size else None
    sd = float(arr.std(ddof=1)) if arr.size > 1 else 0.0
    d: float | None = None
    if mean is not None and sd > 0:
        d = float(mean / sd)
    elif mean is not None:
        warnings.append(
            "the model gives an identical effect in every replicate (sd = 0), "
            "so a standardised effect size is undefined."
        )
    veh_mean = float(np.mean(vehicle_vals)) if vehicle_vals else None
    rel = None
    if mean is not None and veh_mean not in (None, 0.0) and abs(veh_mean) > 1e-12:
        rel = float(abs(mean) / abs(veh_mean))

    observed = "none"
    if mean is not None:
        observed = "increase" if mean > 1e-9 else ("decrease" if mean < -1e-9 else "none")
    return {
        "h_id": h["h_id"],
        "statement": h["statement"],
        "assay": assay,
        "compound": compound,
        "conc_M": conc,
        "readout": readout,
        "contrast": contrast,
        "graph": h.get("graph"),
        "endpoint": h["endpoint"],
        "live_protocol": h.get("live_protocol"),
        "predicted_direction": h["direction"],
        "observed_direction": observed,
        "direction_matches": bool(observed == h["direction"]),
        "n_rep": int(n_rep),
        "seed": int(seed),
        "predicted_effect": mean,
        "sd": sd,
        "ci": _ci(arr),
        "ci_method": "mean +/- 1.96 * sd / sqrt(n_rep) over model replicates",
        "d": d,
        "d_meaning": "model_internal_variability",
        "treated_mean": float(np.mean(treated_vals)) if treated_vals else None,
        "vehicle_mean": veh_mean,
        "relative_effect": rel,
        "replicate_effects": [float(x) for x in arr],
        "uncertainty": ensemble_block,
        "status": "software_prediction",
        "live_result": None,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }


#: the model happily predicts a near-total ablation of a readout; a behavioural
#: endpoint that moves by more than this is not a plausible design assumption.
REL_EFFECT_CAP = 0.5


def _suggested_n(pred: dict[str, Any], baseline: dict[str, Any]) -> dict[str, Any]:
    """Endpoint-appropriate sample size for one prediction.

    Binary (PER) endpoints use the two-proportion test on the model's
    *relative* effect, capped at :data:`REL_EFFECT_CAP`.  Continuous endpoints
    (climbing, receptor occupancy) use the t approximation on the capped
    model-internal ``d``.  When a PER-relative effect cannot be formed (the
    vehicle readout is zero) the continuous approximation is used instead.
    """
    if pred["endpoint"] == "per":
        rel = pred.get("relative_effect")
        if rel is None or not math.isfinite(float(rel)) or rel <= 0:
            out = power_for_continuous(pred.get("d"))
            out["warnings"] = list(out.get("warnings", [])) + [
                "the vehicle readout is zero, so no relative PER effect could "
                "be formed; the continuous approximation was used instead."
            ]
            return out
        capped = min(max(float(rel), 0.05), REL_EFFECT_CAP)
        out = power_for_per(capped, baseline=baseline)
        if capped < float(rel):
            out["warnings"] = list(out.get("warnings", [])) + [
                f"the model's relative effect ({rel:.0%}) was capped at "
                f"{REL_EFFECT_CAP:.0%} for the power calculation: near-total "
                "ablation of a model readout is not a plausible behavioural "
                "effect size."
            ]
        return out
    return power_for_continuous(pred.get("d"))


def prediction_table(n_rep: int = 6, seed: int = 0, h_ids: Iterable[str] | None = None,
                     **kw: Any) -> dict[str, Any]:
    """Run every pre-registered hypothesis and return the pre-registration table."""
    t0 = time.perf_counter()
    ids = list(h_ids) if h_ids is not None else list(HYPOTHESES)
    baseline = behavioral_baseline()
    rows: list[dict[str, Any]] = []
    warnings = list(BASE_WARNINGS) + list(baseline.get("warnings", []))
    warnings.append(f"PER baseline source: {baseline.get('source')}")
    for h_id in ids:
        pred = predict(h_id, n_rep=n_rep, seed=seed, **kw)
        power = _suggested_n(pred, baseline)
        rows.append(
            {
                "h_id": pred["h_id"],
                "statement": pred["statement"],
                "assay": pred["assay"],
                "compound": pred["compound"],
                "conc_M": pred["conc_M"],
                "readout": pred["readout"],
                "predicted_direction": pred["predicted_direction"],
                "observed_direction": pred["observed_direction"],
                "predicted_effect": pred["predicted_effect"],
                "ci": pred["ci"],
                "d": pred["d"],
                "d_meaning": pred["d_meaning"],
                "endpoint": pred["endpoint"],
                "live_protocol": pred["live_protocol"],
                "suggested_n_per_group": power.get("n_per_group"),
                "power_method": power.get("method"),
                "status": "software_prediction",
                "live_result": None,
            }
        )
        for w in pred["warnings"] + power.get("warnings", []):
            if w not in warnings:
                warnings.append(w)
    return {
        "n_rep": int(n_rep),
        "seed": int(seed),
        "baseline": baseline,
        "n_rows": len(rows),
        "rows": rows,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }


def _fmt(x: Any, nd: int = 3) -> str:
    if x is None:
        return "-"
    if isinstance(x, (list, tuple)):
        return "[" + ", ".join(_fmt(v, nd) for v in x) + "]"
    if isinstance(x, float):
        if x != 0 and (abs(x) < 1e-3 or abs(x) >= 1e5):
            return f"{x:.2e}"
        return f"{x:.{nd}f}"
    return str(x)


def to_markdown(rows: dict[str, Any] | Sequence[dict[str, Any]]) -> str:
    """Markdown table of a :func:`prediction_table` result (or its rows)."""
    data = rows["rows"] if isinstance(rows, dict) else list(rows)
    head = (
        "| H | assay | compound | conc (M) | readout | predicted | effect | 95% CI | d (model) | n/group | status | live |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    )
    lines = list(head)
    for r in data:
        lines.append(
            "| {h} | {assay} | {cmp} | {conc} | {ro} | {dir} | {eff} | {ci} | {d} | {n} | {st} | {live} |".format(
                h=r.get("h_id"),
                assay=r.get("assay"),
                cmp=r.get("compound") or "-",
                conc=_fmt(r.get("conc_M")),
                ro=r.get("readout"),
                dir=r.get("predicted_direction"),
                eff=_fmt(r.get("predicted_effect")),
                ci=_fmt(r.get("ci")),
                d=_fmt(r.get("d"), 2),
                n=_fmt(r.get("suggested_n_per_group")),
                st=r.get("status", "software_prediction"),
                live="None" if r.get("live_result") is None else _fmt(r.get("live_result")),
            )
        )
    lines.append("")
    lines.append(
        "d is model-internal variability (teaching EC50 jitter, drive jitter, "
        "RNG seed), not biological variance; suggested n treats the model's "
        "effect as if it were the biological one and caps d at "
        f"{D_CAP:.1f}. status is software_prediction: no live fly was tested."
    )
    return "\n".join(lines)
