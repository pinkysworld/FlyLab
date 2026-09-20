"""One-compartment exposure model: dose -> C(t) -> occupancy(t).

This is a **teaching-order** model. Absorption and elimination rates are not
fitted to any fly ADME dataset; they are round numbers chosen so that the
three routes differ in the direction the literature describes (bath fastest,
topical slowest and least complete). Every profile carries that warning.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from flylab.pharm.occupancy import compare_compound, hill_occupancy, load_library

__all__ = ["ROUTE_DEFAULTS", "FLY_HAEMOLYMPH_VOLUME_L", "exposure_profile"]

#: Nominal haemolymph volume of one adult *Drosophila*, in litres (~80 nL).
#: Teaching constant: adult flies are commonly quoted at 60-100 nL haemolymph
#: for a ~1 mg animal. Used only to turn a dose into a concentration.
FLY_HAEMOLYMPH_VOLUME_L = 8.0e-8

#: Per-route first-order constants. ``f`` is the fraction of the dose that
#: reaches the haemolymph, ``ka``/``ke`` are absorption and elimination rate
#: constants in 1/h. All are ``class_placeholder`` teaching defaults.
ROUTE_DEFAULTS: dict[str, dict[str, float]] = {
    "feeding": {"f": 0.10, "ka": 0.7, "ke": 0.35},
    "topical": {"f": 0.02, "ka": 0.3, "ke": 0.25},
    "bath": {"f": 0.50, "ka": 3.0, "ke": 0.50},
}

_TEACHING_WARNING = (
    "Teaching-order exposure model: one compartment, first-order absorption "
    "and elimination, parameters are class placeholders and are NOT fitted to "
    "Drosophila ADME data."
)


def _bateman(t_h: np.ndarray, c0_M: float, ka: float, ke: float) -> np.ndarray:
    """Bateman function for first-order in / first-order out."""
    if abs(ka - ke) < 1e-9:  # limiting case, flip-flop kinetics
        return c0_M * ke * t_h * np.exp(-ke * t_h)
    return c0_M * (ka / (ka - ke)) * (np.exp(-ke * t_h) - np.exp(-ka * t_h))


def exposure_profile(
    compound: str,
    dose: float,
    route: str,
    t_h: float = 24.0,
    dt_h: float = 0.05,
    params: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Concentration and receptor occupancy over time for one dose.

    Args:
        compound: library key.
        dose: amount delivered to one fly, in **nanomoles**.
        route: ``feeding``, ``topical`` or ``bath``.
        t_h: length of the profile in hours.
        dt_h: time step in hours.
        params: overrides for ``f``, ``ka``, ``ke`` (1/h) and ``volume_L``.

    Returns:
        ``{t_h, conc_M, occupancy: {receptor: [...]}, auc, cmax, tmax,
        warnings}``. ``auc`` is the trapezoidal integral of ``conc_M`` in M*h.
    """
    if route not in ROUTE_DEFAULTS:
        raise ValueError(f"unknown route {route!r}; use one of {sorted(ROUTE_DEFAULTS)}")
    if dose < 0:
        raise ValueError("dose must be >= 0")
    if t_h <= 0 or dt_h <= 0 or dt_h > t_h:
        raise ValueError("need 0 < dt_h <= t_h")

    p = dict(ROUTE_DEFAULTS[route])
    p.update(params or {})
    volume_L = float(p.get("volume_L", FLY_HAEMOLYMPH_VOLUME_L))
    ka, ke, f = float(p["ka"]), float(p["ke"]), float(p["f"])
    if ka <= 0 or ke <= 0 or volume_L <= 0:
        raise ValueError("ka, ke and volume_L must be > 0")

    t = np.arange(0.0, t_h + 0.5 * dt_h, dt_h)
    c0 = (dose * 1e-9 * f) / volume_L  # nmol -> mol -> M in the haemolymph
    conc = np.clip(_bateman(t, c0, ka, ke), 0.0, None)

    lib = load_library()
    entry = lib["compounds"][compound.lower().strip()] if compound.lower().strip() in lib["compounds"] else None
    if entry is None:
        raise KeyError(f"unknown compound {compound!r}")
    occupancy = {
        receptor: [hill_occupancy(float(c), float(spec["ec50_M"]), float(spec.get("n", 1.0))) for c in conc]
        for receptor, spec in entry["receptors"].items()
    }

    i_max = int(np.argmax(conc))
    warnings = [
        _TEACHING_WARNING,
        f"Route defaults for {route!r}: f={f}, ka={ka}/h, ke={ke}/h (evidence_tier class_placeholder).",
        f"Dose is interpreted as nanomoles per fly, diluted into a nominal {volume_L * 1e9:.0f} nL haemolymph volume.",
        "No metabolism, no tissue binding, no protein binding, no excretion lag.",
        "Occupancy(t) assumes instantaneous equilibrium at each time point.",
    ]
    if compound.lower().strip() in {"thiamethoxam"}:
        warnings.append("This compound is a pro-insecticide in vivo; the model does not form its active metabolite.")

    return {
        "compound": entry["name"],
        "key": compound.lower().strip(),
        "route": route,
        "dose_nmol": float(dose),
        "params": {"f": f, "ka_per_h": ka, "ke_per_h": ke, "volume_L": volume_L},
        "t_h": [float(x) for x in t],
        "conc_M": [float(x) for x in conc],
        "occupancy": occupancy,
        "auc": float(np.trapezoid(conc, t)) if hasattr(np, "trapezoid") else float(np.trapz(conc, t)),
        "cmax": float(conc[i_max]),
        "tmax": float(t[i_max]),
        "warnings": warnings,
    }
