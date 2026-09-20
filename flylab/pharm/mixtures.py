"""Combination pharmacology: Bliss independence, Loewe additivity, isobolograms.

Null models (both written out in ``data/literature/mixtures.yaml``
``reference_models``):

* **Bliss independence** -- Bliss CI (1939) *Ann Appl Biol* 26:585.
  ``E_AB = E_A + E_B - E_A*E_B`` on fractional effects, equivalently
  ``S_AB = S_A * S_B`` on survivals.  The right null when the two agents act on
  **disjoint** targets.
* **Loewe additivity** -- Loewe S & Muischnek H (1926) *Naunyn-Schmiedebergs
  Arch* 114:313.  ``CI = d_A/D_A + d_B/D_B``; ``CI = 1`` is the additivity
  isobole.  The right null when the two agents act at the **same** target.
* Greco WR, Bravo G, Parsons JC (1995) *Pharmacol Rev* 47:331 -- the reason a
  synergy claim must always name its null model.  FlyLab therefore reports the
  model it used with every verdict.

Why same-target agonists come out additive here
-----------------------------------------------
With competitive occupancy at one site (Gaddum 1937, via
:func:`flylab.pharm.binding.competitive_occupancy`), writing ``a_i =
([C_i]/EC50_i)^n``, each component holds ``theta_i = a_i / (1 + sum_j a_j)``
and the site's total occupancy is ``sum(a_i) / (1 + sum(a_i))``.  That is
exactly the occupancy a single agent would produce at the summed
equi-effective dose, i.e. ``CI = 1`` identically.  The model therefore predicts
**no interaction** for imidacloprid + clothianidin, which is what Zhu et al.
(2017) observed in honeybee -- the retrospective test
:func:`validate_against_published` runs.
"""

from __future__ import annotations

import math
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from flylab.pharm.binding import competitive_occupancy
from flylab.pharm.mechanisms import gains_from_occupancy
from flylab.pharm.occupancy import hill_occupancy, load_library

__all__ = [
    "BLISS_TOL",
    "LOEWE_TOL",
    "ASSAYS",
    "bliss_expected",
    "bliss_expected_many",
    "inverse_hill",
    "loewe_index",
    "mixture_occupancy",
    "mixture_assay",
    "isobologram",
    "validate_against_published",
]

#: Additivity band on Bliss: |observed - expected| <= 0.05 fractional effect.
BLISS_TOL = 0.05
#: Additivity band on the Loewe combination index: 1/1.25 <= CI <= 1.25.
LOEWE_TOL = 0.25
#: Below this fractional effect a single agent is treated as inactive.
NULL_EFFECT_FLOOR = 1e-6

ASSAYS = ("occupancy", "subgraph", "spiking")

_ACTIVATING = {"agonist", "partial_agonist", "positive_modulator"}
_BLOCKING = {"antagonist", "negative_modulator", "inhibitor"}

MIXTURE_WARNINGS = [
    "Combination result. The null model is named per receptor and per verdict: "
    "Loewe additivity where components share a FlyLab receptor key, Bliss "
    "independence where the targets are disjoint (mixtures.yaml "
    "flylab_recommendation; Greco et al. 1995).",
    "FlyLab combines compounds at the RECEPTOR level at a stated free "
    "concentration. Published synergism ratios are whole-animal numbers "
    "dominated by metabolism and are not occupancy multipliers.",
]


# ---------------------------------------------------------------------------
# null models
# ---------------------------------------------------------------------------
def bliss_expected(e1: float, e2: float) -> float:
    """Bliss (1939) independence: ``E = e1 + e2 - e1*e2`` on fractions 0..1."""
    a, b = float(e1), float(e2)
    for value in (a, b):
        if not 0.0 <= value <= 1.0:
            raise ValueError("fractional effects must lie in [0, 1]")
    return a + b - a * b


def bliss_expected_many(effects: Iterable[float]) -> float:
    """Bliss independence for n agents: ``1 - prod(1 - e_i)``."""
    survival = 1.0
    for e in effects:
        e = float(e)
        if not 0.0 <= e <= 1.0:
            raise ValueError("fractional effects must lie in [0, 1]")
        survival *= 1.0 - e
    return 1.0 - survival


def inverse_hill(effect: float, ec50_M: float, n: float = 1.0) -> float:
    """Concentration of one agent alone producing ``effect`` (0..1).

    Inverse of the Hill equation: ``C = EC50 * (E/(1-E))**(1/n)``.  This is the
    ``D_A`` of the Loewe combination index.
    """
    e = float(effect)
    if not 0.0 < e < 1.0:
        raise ValueError("effect must lie strictly in (0, 1) to invert the Hill curve")
    if ec50_M <= 0 or n <= 0:
        raise ValueError("ec50_M and n must be > 0")
    return float(ec50_M) * (e / (1.0 - e)) ** (1.0 / float(n))


def loewe_index(
    doses: Sequence[float],
    ec50s: Sequence[float],
    effect: float,
    hills: Sequence[float] | None = None,
) -> float:
    """Loewe (1926) combination index ``CI = sum_i d_i / D_i``.

    ``D_i`` is read off component ``i``'s own dose-response curve at the
    observed ``effect`` (:func:`inverse_hill`).  ``CI < 1`` synergy, ``CI = 1``
    additivity, ``CI > 1`` antagonism.
    """
    if len(doses) != len(ec50s):
        raise ValueError("doses and ec50s must have the same length")
    ns = list(hills) if hills is not None else [1.0] * len(doses)
    if len(ns) != len(doses):
        raise ValueError("hills must have the same length as doses")
    total = 0.0
    for d, ec50, n in zip(doses, ec50s, ns):
        if d <= 0:
            continue
        total += float(d) / inverse_hill(effect, ec50, n)
    return float(total)


def _loewe_verdict(ci: float | None, tol: float = LOEWE_TOL) -> str:
    if ci is None or not math.isfinite(ci) or ci <= 0:
        return "inconclusive"
    if ci < 1.0 / (1.0 + tol):
        return "synergistic"
    if ci > 1.0 + tol:
        return "antagonistic"
    return "additive"


# ---------------------------------------------------------------------------
# occupancy of a mixture
# ---------------------------------------------------------------------------
def _components(components: Iterable[Mapping[str, Any]], conc_scale: float = 1.0) -> list[dict[str, Any]]:
    out = []
    for c in components:
        key = str(c["compound"]).lower().strip()
        conc = float(c.get("conc_M", c.get("conc", 0.0))) * float(conc_scale)
        if conc < 0:
            raise ValueError("concentrations must be >= 0")
        out.append({"compound": key, "conc_M": conc})
    if len(out) < 1:
        raise ValueError("a mixture needs at least one component")
    return out


def _specs(components: list[dict[str, Any]], library: dict[str, Any] | None) -> dict[str, dict[str, Any]]:
    lib = library or load_library()
    known = lib.get("compounds") or {}
    missing = [c["compound"] for c in components if c["compound"] not in known]
    if missing:
        raise KeyError(
            "compound(s) not in the library: "
            + ", ".join(missing)
            + ". FlyLab does not substitute an analogue silently."
        )
    return {c["compound"]: known[c["compound"]] for c in components}


def _family(direction: str) -> str | None:
    if direction in _ACTIVATING:
        return "activating"
    if direction in _BLOCKING:
        return "blocking"
    return None


def mixture_occupancy(
    components: Iterable[Mapping[str, Any]],
    library: dict[str, Any] | None = None,
    conc_scale: float = 1.0,
) -> dict[str, Any]:
    """Per-receptor combined occupancy for a mixture.

    ``components`` is ``[{"compound": key, "conc_M": float}, ...]``.

    Same-target components are combined by competitive occupancy
    (:func:`flylab.pharm.binding.competitive_occupancy`), which is Loewe-additive
    by construction; components at disjoint receptors are independent (Bliss).
    Each receptor row records ``model`` and ``why``.

    Returns ``{"components", "receptors": [...], "rows": [...], "models_used",
    "shared_receptors", "warnings"}`` where ``rows`` is directly consumable by
    :func:`flylab.pharm.mechanisms.gains_from_occupancy`.
    """
    comps = _components(components, conc_scale)
    specs = _specs(comps, library)

    by_receptor: dict[str, list[dict[str, Any]]] = {}
    for comp in comps:
        entry = specs[comp["compound"]]
        for receptor, spec in (entry.get("receptors") or {}).items():
            direction = str(spec.get("direction", "none"))
            if direction == "none":
                continue
            ec50 = float(spec["ec50_M"])
            n = float(spec.get("n", 1.0))
            by_receptor.setdefault(receptor, []).append(
                {
                    "compound": comp["compound"],
                    "conc_M": comp["conc_M"],
                    "ec50_M": ec50,
                    "n": n,
                    "direction": direction,
                    "family": _family(direction),
                    "efficacy": spec.get("efficacy"),
                    "a": (comp["conc_M"] / ec50) ** n if comp["conc_M"] > 0 else 0.0,
                    "occupancy_alone": hill_occupancy(comp["conc_M"], ec50, n),
                    "source": spec.get("source"),
                    "evidence_tier": spec.get("evidence_tier"),
                }
            )

    warnings = list(MIXTURE_WARNINGS)
    receptor_rows: list[dict[str, Any]] = []
    gain_rows: list[dict[str, Any]] = []
    shared: list[str] = []
    models_used: dict[str, str] = {}

    for receptor in sorted(by_receptor):
        parts = by_receptor[receptor]
        # Per-component occupancy comes from the shared binding module. One
        # competitor slot carries every other component: with
        # competitor_kd_M = 1, ``competitive_occupancy`` applies Gaddum's
        # Kd_app = Kd_i * (1 + competitor_M). Passing
        #   competitor_M = (1 + sum_j!=i a_j)**(1/n) - 1,   a_j = ([C_j]/EC50_j)^n
        # makes the Hill exponent act on each component separately, which is the
        # Langmuir competition generalised to a Hill exponent; for n = 1 it
        # reduces to the textbook Kd_app = Kd_i * (1 + [B]/Kb).
        for idx, p in enumerate(parts):
            others = sum(q["a"] for j, q in enumerate(parts) if j != idx)
            competitor = (1.0 + others) ** (1.0 / p["n"]) - 1.0
            p["occupancy_in_mixture"] = float(
                competitive_occupancy(p["conc_M"], p["ec50_M"], competitor, 1.0, n=p["n"])
            )
        theta_total = min(1.0, sum(p["occupancy_in_mixture"] for p in parts))

        if len(parts) == 1:
            model = "independent_single_agent"
            why = (
                "only one component acts at this receptor, so its occupancy is "
                "unchanged by the mixture; across disjoint receptors the "
                "combination follows Bliss independence."
            )
        else:
            shared.append(receptor)
            model = "loewe_competitive"
            why = (
                "two or more components share this receptor key, so they compete "
                "for one site: each component's occupancy is "
                "binding.competitive_occupancy with every other component in the "
                "competitor term (Gaddum 1937), and the site total is their sum. "
                "For n = 1 that is exactly the occupancy of a single agent at the "
                "summed equi-effective dose, i.e. Loewe additivity with CI = 1."
            )
            ns = {round(p["n"], 6) for p in parts}
            if ns != {1.0}:
                k = len(parts)
                deviation = k ** (1.0 - 1.0 / max(ns))
                warnings.append(
                    f"{receptor}: Hill n = {sorted(ns)} is not 1, so even PURE competition "
                    f"gives a dose-additivity index of about {deviation:.2f} for {k} equipotent "
                    "components instead of exactly 1. That offset is an artefact of the Hill "
                    "exponent, not an interaction, and is inside the default Loewe tolerance."
                )
            if len(ns) > 1:
                warnings.append(
                    f"{receptor}: components have different Hill coefficients {sorted(ns)}; "
                    "the competitive expression is exact only for equal n, so the combined "
                    "occupancy here is an approximation."
                )

        models_used[receptor] = model
        families = {p["family"] for p in parts if p["family"]}
        row = {
            "receptor": receptor,
            "occupancy": theta_total,
            "model": model,
            "why": why,
            "n_components": len(parts),
            "components": [
                {
                    "compound": p["compound"],
                    "conc_M": p["conc_M"],
                    "ec50_M": p["ec50_M"],
                    "n": p["n"],
                    "direction": p["direction"],
                    "occupancy_alone": p["occupancy_alone"],
                    "occupancy_in_mixture": p["occupancy_in_mixture"],
                }
                for p in parts
            ],
            "directions": sorted({p["direction"] for p in parts}),
        }

        if len(families) > 1:
            by_family = {
                fam: sum(p["occupancy_in_mixture"] for p in parts if p["family"] == fam)
                for fam in families
            }
            dominant = max(by_family, key=lambda f: by_family[f])
            suppressed = [f for f in by_family if f != dominant]
            row["model"] = models_used[receptor] = "competitive_mixed_direction"
            row["occupancy_by_family"] = by_family
            row["dominant_family"] = dominant
            row["why"] = (
                why + " Components pull this receptor in opposite directions "
                "(activating vs blocking); the circuit patch uses the family "
                "holding the larger share of the site and the other is reported "
                "but not applied."
            )
            warnings.append(
                f"{receptor}: mixed-direction combination ({sorted(row['directions'])}). "
                f"The gain patch uses the {dominant} family (occupancy "
                f"{by_family[dominant]:.3f}); {', '.join(suppressed)} "
                f"({', '.join(f'{by_family[f]:.3f}' for f in suppressed)}) is reported only. "
                "FlyLab has no rule for a receptor held partly open and partly blocked."
            )
            direction = next(p["direction"] for p in parts if p["family"] == dominant)
            occupancy_for_gain = by_family[dominant]
        else:
            direction = parts[0]["direction"]
            occupancy_for_gain = theta_total

        receptor_rows.append(row)
        gain_rows.append(
            {
                "receptor": receptor,
                "occupancy": occupancy_for_gain,
                "direction": direction,
                "source": "mixture of " + " + ".join(p["compound"] for p in parts),
                "evidence_tier": "literature_order",
                "n": parts[0]["n"],
            }
        )

    if not shared:
        warnings.append(
            "No receptor is shared by two components: the targets are disjoint, "
            "so Bliss independence is the applicable null model."
        )
    return {
        "components": comps,
        "receptors": receptor_rows,
        "rows": gain_rows,
        "shared_receptors": shared,
        "models_used": models_used,
        "warnings": warnings,
    }


def _loewe_receptor_index(mix: dict[str, Any]) -> list[dict[str, Any]]:
    """Analytic combination index at each shared receptor."""
    out = []
    for row in mix["receptors"]:
        if row["n_components"] < 2:
            continue
        theta = row["occupancy"]
        if not 0.0 < theta < 1.0:
            out.append({"receptor": row["receptor"], "combination_index": None,
                        "verdict": "inconclusive",
                        "why": "combined occupancy is 0 or 1; the curve cannot be inverted"})
            continue
        doses = [c["conc_M"] for c in row["components"]]
        ec50s = [c["ec50_M"] for c in row["components"]]
        ns = [c["n"] for c in row["components"]]
        ci = loewe_index(doses, ec50s, theta, ns)
        out.append(
            {
                "receptor": row["receptor"],
                "combination_index": ci,
                "verdict": _loewe_verdict(ci),
                "effect_used": theta,
                "why": "CI from the components' own Hill curves at the mixture's occupancy.",
            }
        )
    return out


# ---------------------------------------------------------------------------
# circuit runs
# ---------------------------------------------------------------------------
_VEHICLE_CACHE: dict[tuple, dict[str, float | None]] = {}


def _rate_readouts(gains: Mapping[str, float], graph: str | None, drive_hz: float, steps: int):
    from flylab.assays.subgraph import named_mean, named_readout
    from flylab.circuit.rate import rate_network

    net = rate_network(graph)
    drive = net.drive_vector(net.seed_drive(drive_hz))
    r = net.run(drive, dict(gains), steps=steps)
    named = named_readout(net, r)
    return {
        "mn9_hz": named_mean(named, "MN9"),
        "dnp01_hz": named_mean(named, "DNp01"),
        "mean_hz": float(r.mean()),
        "max_hz": float(r.max()),
    }, np.asarray(r, dtype=float)


def _spiking_readouts(
    gains: Mapping[str, float], graph: str | None, drive_hz: float, t_ms: float, seed: int
):
    from flylab.assays.spiking import _default_drive, _named, _named_mean
    from flylab.circuit.lif import BACKGROUND_FRACTION_DEFAULT, lif_network
    from flylab.circuit.rate import load_graph, resolve_graph

    g = load_graph(resolve_graph(graph))
    net = lif_network(g, seed=seed)
    drive = _default_drive(g, drive_hz, BACKGROUND_FRACTION_DEFAULT * drive_hz)
    result = net.run(drive, dict(gains), t_ms=t_ms)
    named = _named(net, result.rates_hz)
    return {
        "mn9_hz": _named_mean(named, "MN9"),
        "dnp01_hz": _named_mean(named, "DNp01"),
        "mean_hz": float(result.rates_hz.mean()),
        "max_hz": float(result.rates_hz.max()),
    }, np.asarray(result.rates_hz, dtype=float)


def _circuit_readouts(
    gains: Mapping[str, float],
    assay: str,
    graph: str | None,
    drive_hz: float,
    steps: int,
    t_ms: float,
    seed: int,
):
    """``(readouts, per-node rate vector)``; the vector never enters a notebook."""
    if assay == "subgraph":
        return _rate_readouts(gains, graph, drive_hz, steps)
    if assay == "spiking":
        return _spiking_readouts(gains, graph, drive_hz, t_ms, seed)
    raise ValueError(f"unknown circuit assay {assay!r}")


def _vehicle_readouts(assay: str, graph: str | None, drive_hz: float, steps: int, t_ms: float, seed: int):
    from flylab.circuit.rate import DEFAULT_GAINS

    key = (assay, str(graph), drive_hz, steps, t_ms, seed)
    if key not in _VEHICLE_CACHE:
        _VEHICLE_CACHE[key] = _circuit_readouts(DEFAULT_GAINS, assay, graph, drive_hz, steps, t_ms, seed)
    readouts, rates = _VEHICLE_CACHE[key]
    return dict(readouts), rates


def circuit_disruption(rates: np.ndarray, vehicle_rates: np.ndarray) -> float:
    """Fraction of the vehicle network activity displaced by the drug.

    ``E = mean(|r_treated - r_vehicle|) / mean(r_vehicle)``.  This is the
    default effect metric for a combination verdict because it is
    **non-negative and direction-free**: FlyLab's nAChR patch is deliberately
    non-monotone (it desensitises at high occupancy), so a signed readout such
    as ``mn9_hz`` can move up for one agent and down for another at the same
    concentration, and Bliss independence is undefined for opposing effects.
    Values above 1 mean the drug displaced more activity than the vehicle had.
    """
    veh = float(np.mean(vehicle_rates))
    if veh <= 1e-12:
        return float("nan")
    return float(np.mean(np.abs(np.asarray(rates) - np.asarray(vehicle_rates))) / veh)


def _relative(treated: float | None, vehicle: float | None) -> float | None:
    if treated is None or vehicle is None or abs(vehicle) < 1e-12:
        return None
    return (float(treated) - float(vehicle)) / float(vehicle)


def _fraction(rel: float | None) -> float | None:
    if rel is None:
        return None
    return min(1.0, abs(float(rel)))


# ---------------------------------------------------------------------------
# mixture assay
# ---------------------------------------------------------------------------
def mixture_assay(
    components: Iterable[Mapping[str, Any]],
    assay: str = "subgraph",
    conc_scale: float = 1.0,
    graph: str | None = None,
    model: str = "bliss",
    library: dict[str, Any] | None = None,
    readout: str = "circuit_disruption",
    drive_hz: float = 40.0,
    steps: int = 80,
    t_ms: float = 500.0,
    seed: int = 0,
    bliss_tol: float = BLISS_TOL,
    loewe_tol: float = LOEWE_TOL,
    **kw: Any,
) -> dict[str, Any]:
    """Run the circuit with the combined gains and judge the interaction.

    The combined gains come from :func:`mixture_occupancy` ->
    :func:`flylab.pharm.mechanisms.gains_from_occupancy`, so the same patch
    rules apply as for a single compound.  Each single agent is run on its own
    for comparison, and the observed combined effect is scored against the
    named null model.

    Tolerances (explicit, not adaptive): Bliss additivity is
    ``|observed - expected| <= bliss_tol`` (default 0.05 in fractional-effect
    units); Loewe additivity is ``1/(1+loewe_tol) <= CI <= 1 + loewe_tol``
    (default CI in [0.80, 1.25]).

    Returns a notebook with ``assay == "mixture"``.
    """
    from flylab.notebook.schema import empty_notebook, set_map
    from flylab.circuit.rate import load_graph, resolve_graph

    if assay not in ("subgraph", "spiking"):
        raise ValueError(f"mixture_assay needs a circuit assay ('subgraph'|'spiking'), got {assay!r}")
    if model not in ("bliss", "loewe"):
        raise ValueError("model must be 'bliss' or 'loewe'")

    comps = _components(components, conc_scale)
    mix = mixture_occupancy(comps, library=library)
    gains = gains_from_occupancy(mix["rows"])

    vehicle, veh_rates = _vehicle_readouts(assay, graph, drive_hz, steps, t_ms, seed)
    combined, comb_rates = _circuit_readouts(gains, assay, graph, drive_hz, steps, t_ms, seed)

    def _effect(readouts: dict[str, Any], rates: np.ndarray) -> tuple[float | None, float | None]:
        """``(signed relative change, fractional effect)`` for the chosen readout."""
        if readout == "circuit_disruption":
            e = circuit_disruption(rates, veh_rates)
            return (None, None) if e != e else (e, e)
        rel = _relative(readouts.get(readout), vehicle.get(readout))
        return rel, _fraction(rel)

    singles = []
    for comp in comps:
        single_mix = mixture_occupancy([comp], library=library)
        single_gains = gains_from_occupancy(single_mix["rows"])
        r, r_rates = _circuit_readouts(single_gains, assay, graph, drive_hz, steps, t_ms, seed)
        rel_s, e_s = _effect(r, r_rates)
        singles.append(
            {
                "compound": comp["compound"],
                "conc_M": comp["conc_M"],
                "gains": single_gains,
                "readouts": r,
                "relative_change": rel_s,
                "effect_fraction": e_s,
            }
        )

    rel_obs, e_obs = _effect(combined, comb_rates)
    warnings = list(mix["warnings"])
    bliss_block, verdict, reason = _judge_bliss(singles, rel_obs, e_obs, bliss_tol, warnings)

    loewe_receptors = _loewe_receptor_index(mix)
    loewe_block = {
        "receptor_level": loewe_receptors,
        "shared_receptors": mix["shared_receptors"],
        "note": (
            "Receptor-level CI from the components' Hill curves. For a shared "
            "site this is 1 by construction (competitive occupancy is "
            "Loewe-additive), which is the model's prediction of NO "
            "interaction, not a fitted result."
        ),
    }
    if loewe_receptors:
        cis = [r["combination_index"] for r in loewe_receptors if r["combination_index"]]
        loewe_block["combination_index"] = float(np.mean(cis)) if cis else None
        loewe_block["verdict"] = _loewe_verdict(loewe_block["combination_index"], loewe_tol)
    else:
        loewe_block["combination_index"] = None
        loewe_block["verdict"] = "not_applicable"
        loewe_block["why"] = "no receptor is shared by two components; use Bliss."

    if model == "loewe":
        if loewe_block["verdict"] in ("not_applicable", "inconclusive"):
            final, final_reason = "inconclusive", (
                "Loewe requested but no shared receptor gives an invertible curve; "
                "Bliss verdict was " + verdict
            )
        else:
            final, final_reason = loewe_block["verdict"], "Loewe combination index at the shared receptor(s)."
    else:
        final, final_reason = verdict, reason

    path = resolve_graph(graph)
    g = load_graph(path)
    nb = empty_notebook("mixture", seed=seed)
    set_map(nb, g["map"], path.stem, g["citation"])
    nb["compound"] = " + ".join(c["compound"] for c in comps)
    nb["concentration_M"] = None
    nb["occupancy"] = mix["rows"]
    nb["gains"] = dict(gains)
    nb["readouts"] = {
        "engine": assay,
        "readout": readout,
        "combined": combined,
        "vehicle": vehicle,
        "singles": singles,
        "relative_change": rel_obs,
        "effect_fraction": e_obs,
        "n_nodes": g["n_nodes"],
        "n_edges": g["n_edges"],
        "drive_hz": float(drive_hz),
    }
    nb["mixture"] = {
        "components": comps,
        "conc_scale": float(conc_scale),
        "receptors": mix["receptors"],
        "models_used": mix["models_used"],
        "shared_receptors": mix["shared_receptors"],
        "bliss": bliss_block,
        "loewe": loewe_block,
        "null_model": model,
        "tolerances": {"bliss_abs_effect": bliss_tol, "loewe_ci_relative": loewe_tol},
    }
    nb["synergy"] = {"verdict": final, "model": model, "why": final_reason}
    nb["warnings"] = warnings
    return nb


#: a single agent at or above this fractional effect is at the readout ceiling
SATURATION = 0.9


def _judge_bliss(
    singles: list[dict[str, Any]],
    rel_obs: float | None,
    e_obs: float | None,
    tol: float,
    warnings: list[str],
) -> tuple[dict[str, Any], str, str]:
    raw = [s["effect_fraction"] for s in singles]
    effects = [None if e is None else min(1.0, float(e)) for e in raw]
    rels = [s["relative_change"] for s in singles]
    e_obs_raw = e_obs
    e_obs = None if e_obs is None else min(1.0, float(e_obs))
    block: dict[str, Any] = {
        "single_effect_fractions": effects,
        "single_effect_fractions_unclipped": raw,
        "observed_effect_fraction": e_obs,
        "observed_effect_fraction_unclipped": e_obs_raw,
        "tolerance": tol,
    }
    if effects and e_obs is not None and any(e is not None and e >= SATURATION for e in effects):
        block.update({"expected_effect_fraction": None, "delta": None, "verdict": "inconclusive"})
        warnings.append(
            f"At least one single agent already displaces >= {SATURATION:.0%} of the vehicle "
            "readout, so the Bliss expectation is pinned near 1 and cannot discriminate "
            "synergy from additivity. Lower the concentration."
        )
        return block, "inconclusive", "a single agent is at the readout ceiling."
    if any(e is None for e in effects) or e_obs is None:
        block.update({"expected_effect_fraction": None, "delta": None, "verdict": "inconclusive"})
        return block, "inconclusive", "a readout was undefined (vehicle rate is zero or the cell is missing)."
    if all(e < NULL_EFFECT_FLOOR for e in effects):
        block.update({"expected_effect_fraction": 0.0, "delta": e_obs, "verdict": "inconclusive"})
        warnings.append(
            "Neither single agent moves the readout at these concentrations, so "
            "no interaction can be scored: raise the concentrations or change the readout."
        )
        return block, "inconclusive", "no single agent produced a measurable effect."

    signs = {int(np.sign(r)) for r in rels if r is not None and abs(r) > NULL_EFFECT_FLOOR}
    if len(signs) > 1:
        block.update({"expected_effect_fraction": None, "delta": None, "verdict": "inconclusive"})
        warnings.append(
            "The single agents move the readout in OPPOSITE directions; Bliss and "
            "Loewe are both undefined for opposing effects, so no synergy verdict is given."
        )
        return block, "inconclusive", "single agents act in opposite directions on this readout."

    expected = bliss_expected_many(effects)
    delta = e_obs - expected
    common_sign = signs.pop() if signs else 0
    obs_sign = int(np.sign(rel_obs)) if rel_obs is not None else 0
    block.update({"expected_effect_fraction": expected, "delta": delta})
    if common_sign and obs_sign and obs_sign != common_sign:
        block["verdict"] = "antagonistic"
        return block, "antagonistic", (
            "the combination moves the readout in the opposite direction to both single agents."
        )
    if delta > tol:
        block["verdict"] = "synergistic"
        return block, "synergistic", f"observed {e_obs:.3f} exceeds the Bliss expectation {expected:.3f} by more than {tol}."
    if delta < -tol:
        block["verdict"] = "antagonistic"
        return block, "antagonistic", f"observed {e_obs:.3f} falls short of the Bliss expectation {expected:.3f} by more than {tol}."
    block["verdict"] = "additive"
    return block, "additive", f"observed {e_obs:.3f} is within {tol} of the Bliss expectation {expected:.3f}."


# ---------------------------------------------------------------------------
# isobologram
# ---------------------------------------------------------------------------
def _target_receptor(compound: str, other: str | None, library: dict[str, Any] | None) -> str:
    lib = library or load_library()
    rows = (lib["compounds"][compound].get("receptors") or {})
    active = {r: s for r, s in rows.items() if str(s.get("direction", "none")) != "none"}
    if other is not None:
        other_rows = (lib["compounds"][other].get("receptors") or {})
        shared = [
            r
            for r in active
            if str((other_rows.get(r) or {}).get("direction", "none")) != "none"
        ]
        if shared:
            active = {r: active[r] for r in shared}
    if not active:
        raise ValueError(f"{compound} has no active receptor row")
    insect = {r: s for r, s in active.items() if r.startswith("insect_")} or active
    return min(insect, key=lambda r: float(insect[r]["ec50_M"]))


def _effect_of(
    components: list[dict[str, Any]],
    assay: str,
    readout: str | None,
    library: dict[str, Any] | None,
    **kw: Any,
) -> float:
    """Fractional effect (0..1) of a component list, for the iso-effect search."""
    if assay == "occupancy":
        mix = mixture_occupancy(components, library=library)
        target = readout
        rows = {r["receptor"]: r for r in mix["receptors"]}
        if target is None or target not in rows:
            target = max(rows, key=lambda r: rows[r]["occupancy"]) if rows else None
        return float(rows[target]["occupancy"]) if target else 0.0
    nb = mixture_assay(components, assay=assay, library=library, readout=readout or "mn9_hz", **kw)
    return float(nb["readouts"]["effect_fraction"] or 0.0)


def _iso_conc(
    compound: str,
    effect_frac: float,
    assay: str,
    readout: str | None,
    library: dict[str, Any] | None,
    lo: float = 1e-14,
    hi: float = 1e-1,
    iters: int = 45,
    **kw: Any,
) -> float | None:
    """Concentration of one compound alone producing ``effect_frac`` (bisection)."""
    f_hi = _effect_of([{"compound": compound, "conc_M": hi}], assay, readout, library, **kw)
    if f_hi < effect_frac:
        return None
    for _ in range(iters):
        mid = math.sqrt(lo * hi)
        if _effect_of([{"compound": compound, "conc_M": mid}], assay, readout, library, **kw) < effect_frac:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def isobologram(
    compound_a: str,
    compound_b: str,
    assay: str = "occupancy",
    readout: str | None = None,
    effect_frac: float = 0.5,
    n: int = 7,
    library: dict[str, Any] | None = None,
    iters: int = 40,
    **kw: Any,
) -> dict[str, Any]:
    """Iso-effect concentration pairs, plus the Loewe additivity line.

    For ``n`` mixing ratios the total dose is scaled until the mixture produces
    ``effect_frac``; the resulting ``(conc_a, conc_b)`` pairs are the isobole.
    The additivity line joins ``(D_A, 0)`` and ``(0, D_B)``: points on it mean
    ``CI = 1``, below it synergy, above it antagonism.

    ``assay="occupancy"`` (default) is analytic and fast; ``"subgraph"`` runs
    the rate network at every bisection step and is much slower.
    """
    a = compound_a.lower().strip()
    b = compound_b.lower().strip()
    if readout is None and assay == "occupancy":
        readout = _target_receptor(a, b, library)

    d_a = _iso_conc(a, effect_frac, assay, readout, library, iters=iters, **kw)
    d_b = _iso_conc(b, effect_frac, assay, readout, library, iters=iters, **kw)
    warnings = list(MIXTURE_WARNINGS)
    if d_a is None or d_b is None:
        warnings.append(
            "One component cannot reach the requested effect alone at any "
            "concentration up to 0.1 M, so no isobologram exists for it."
        )
        return {
            "compound_a": a, "compound_b": b, "assay": assay, "readout": readout,
            "effect_frac": effect_frac, "d_a": d_a, "d_b": d_b,
            "points": [], "additivity_line": [], "warnings": warnings,
        }

    points: list[dict[str, Any]] = []
    for i in range(int(n)):
        f = i / (int(n) - 1) if n > 1 else 0.0  # fraction of the dose contributed by B
        # ratio fixed at (1-f)*d_a : f*d_b, scaled until the effect is reached
        base_a, base_b = (1.0 - f) * d_a, f * d_b
        if base_a <= 0 and base_b <= 0:  # pragma: no cover - guarded by f in [0,1]
            continue
        lo, hi = 1e-4, 1e4
        def eff(scale: float) -> float:
            comps = [c for c in (
                {"compound": a, "conc_M": base_a * scale},
                {"compound": b, "conc_M": base_b * scale},
            ) if c["conc_M"] > 0]
            return _effect_of(comps, assay, readout, library, **kw)

        if eff(hi) < effect_frac:
            continue
        for _ in range(int(iters)):
            mid = math.sqrt(lo * hi)
            if eff(mid) < effect_frac:
                lo = mid
            else:
                hi = mid
        scale = math.sqrt(lo * hi)
        ca, cb = base_a * scale, base_b * scale
        points.append(
            {
                "fraction_b": f,
                "conc_a_M": ca,
                "conc_b_M": cb,
                "effect": eff(scale),
                "combination_index": (ca / d_a) + (cb / d_b),
            }
        )

    return {
        "compound_a": a,
        "compound_b": b,
        "assay": assay,
        "readout": readout,
        "effect_frac": float(effect_frac),
        "d_a": d_a,
        "d_b": d_b,
        "points": points,
        "additivity_line": [
            {"conc_a_M": d_a, "conc_b_M": 0.0},
            {"conc_a_M": 0.0, "conc_b_M": d_b},
        ],
        "model": "Loewe additivity isobole (Loewe & Muischnek 1926); CI = 1 on the line.",
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# retrospective validation
# ---------------------------------------------------------------------------
#: honeybee binary mixtures of Zhu et al. 2017 (mixtures.yaml
#: ``imidacloprid_binary_mixtures_honeybee``).  The two partners FlyLab has no
#: library entry for are mapped to the nearest library member ONLY when the
#: caller asks for it, and the result then carries ``substituted``.
SUBSTITUTES = {
    "lambda-cyhalothrin": ("deltamethrin", "same class (type II pyrethroid, insect_Nav); "
                                           "lambda-cyhalothrin is not in the FlyLab library"),
    "acephate": ("chlorpyrifos_oxon", "acephate is an organophosphate pro-insecticide acting "
                                      "through AChE inhibition; the library's only AChE "
                                      "inhibitor is chlorpyrifos-oxon"),
}

#: published classification -> FlyLab verdict vocabulary
_PUBLISHED_TO_VERDICT = {
    "no_interaction": "additive",
    "additive": "additive",
    "synergistic": "synergistic",
    "antagonistic": "antagonistic",
}


def _mixtures_dataset() -> dict[str, Any]:
    import yaml
    from pathlib import Path

    for d in (Path(__file__).resolve().parents[2] / "data" / "literature", Path("data/literature")):
        p = d / "mixtures.yaml"
        if p.exists():
            return yaml.safe_load(p.read_text()) or {}
    return {}


def validate_against_published(
    conc_M: float = 1e-9,
    assay: str = "subgraph",
    substitute: bool = True,
    library: dict[str, Any] | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Run the three honeybee binary mixtures and compare with the publication.

    Targets (Zhu YC, Yao J, Adamczyk J, Luttrell R (2017) *PLoS ONE*
    12(5):e0176837, recorded in ``mixtures.yaml``):

    ======================================  ==================
    imidacloprid + clothianidin             no interaction
    imidacloprid + lambda-cyhalothrin       additive
    imidacloprid + acephate                 additive
    ======================================  ==================

    Only clothianidin is in the FlyLab library.  With ``substitute=False`` the
    other two are reported as ``skipped`` with the missing compound named.  With
    ``substitute=True`` (default) the nearest library member is used and the row
    carries ``substituted`` plus a warning: a substituted row tests the *class*
    prediction, not the published pair.

    The concentration is a FlyLab modelling choice (default 1e-9 M each, equal
    parts).  The publication dosed imidacloprid at its honeybee LC20 in a
    whole-animal mortality assay; no free CNS concentration is available, so the
    test is of the *interaction class*, never of the mortality number.  1e-9 M
    is chosen because it keeps every agent on the rising limb of the model's own
    dose-response: above roughly 1e-8 M the nAChR patch passes its
    desensitisation peak, the rate readout stops being monotone in occupancy and
    a Bliss verdict is not interpretable.  The concentration dependence is a
    property of the patch rule and is reported, not hidden - pass another
    ``conc_M`` to see it.
    """
    data = _mixtures_dataset()
    entry = next(
        (e for e in data.get("empirical", []) if e.get("id") == "imidacloprid_binary_mixtures_honeybee"),
        None,
    )
    wanted = {"clothianidin (neonicotinoid)", "lambda-cyhalothrin (pyrethroid)", "acephate (organophosphate)"}
    results: list[dict[str, Any]] = []
    known = set((library or load_library())["compounds"])
    warnings = [
        "Retrospective validation against a HONEYBEE whole-animal mortality "
        "study (Zhu et al. 2017). FlyLab runs a Drosophila MaleCNS circuit at a "
        "stated free concentration: agreement on the interaction class is "
        "meaningful, agreement on magnitude would not be.",
        "The published imidacloprid dose was its LC20 in bees; the concentration "
        "used here is a FlyLab modelling choice, not the published dose.",
    ]

    rows = (entry or {}).get("results", []) if entry else []
    for row in rows:
        partner = str(row.get("partner", ""))
        if partner not in wanted:
            continue
        published = str(row.get("classification", ""))
        key = row.get("library_key") or partner.split(" (")[0].strip().lower()
        substituted = None
        if key not in known:
            sub = SUBSTITUTES.get(partner.split(" (")[0].strip().lower())
            if not (substitute and sub):
                results.append(
                    {
                        "partner": partner,
                        "published": published,
                        "skipped": True,
                        "reason": f"{key} is not in the FlyLab library and no substitution was requested",
                        "missing_compounds": [key],
                    }
                )
                continue
            substituted = {"requested": key, "used": sub[0], "reason": sub[1]}
            warnings.append(
                f"SUBSTITUTION: {key} is not in the library; {sub[0]} was used instead "
                f"({sub[1]}). This row tests the class prediction, not the published pair."
            )
            key = sub[0]

        components = [
            {"compound": "imidacloprid", "conc_M": conc_M},
            {"compound": key, "conc_M": conc_M},
        ]
        mix = mixture_occupancy(components, library=library)
        model = "loewe" if mix["shared_receptors"] else "bliss"
        nb = mixture_assay(components, assay=assay, model=model, library=library, **kw)
        verdict = nb["synergy"]["verdict"]
        expected = _PUBLISHED_TO_VERDICT.get(published, published)
        results.append(
            {
                "partner": partner,
                "components": components,
                "published": published,
                "published_as_flylab_verdict": expected,
                "model_verdict": verdict,
                "null_model": model,
                "matches": verdict == expected,
                "shared_receptors": mix["shared_receptors"],
                "bliss": nb["mixture"]["bliss"],
                "loewe": {k: nb["mixture"]["loewe"][k] for k in ("combination_index", "verdict")},
                "readout": nb["readouts"]["readout"],
                "effect_fraction": nb["readouts"]["effect_fraction"],
                "substituted": substituted,
                "skipped": False,
                "why": nb["synergy"]["why"],
            }
        )

    evaluated = [r for r in results if not r["skipped"]]
    n_match = sum(1 for r in evaluated if r["matches"])
    return {
        "source": (entry or {}).get("source"),
        "species": (entry or {}).get("species"),
        "conc_M": conc_M,
        "assay": assay,
        "results": results,
        "summary": {
            "n_targets": len(results),
            "n_evaluated": len(evaluated),
            "n_skipped": len(results) - len(evaluated),
            "n_matching": n_match,
            "n_substituted": sum(1 for r in evaluated if r.get("substituted")),
            "all_match": bool(evaluated) and n_match == len(evaluated),
        },
        "warnings": warnings,
        "honesty_note": (
            "Mismatches are reported, not tuned away. The library is not adjusted "
            "to make a published interaction class come out right."
        ),
    }
