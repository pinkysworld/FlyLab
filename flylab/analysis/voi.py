"""Value of information: which experiment would remove the most model variance?

``VOI_j = Var(Y) - E[Var(Y | X_j)]`` is the variance of the readout that would
disappear if factor ``j`` were resolved -- learned exactly, at its true value --
while every other assumption stayed as uncertain as it is now.  By the law of
total variance,

.. math::

    E[Var(Y | X_j)] = Var(Y)\\,(1 - S_j)
    \\quad\\Longrightarrow\\quad
    VOI_j = S_j \\cdot Var(Y)

so the quantity is the **first-order Sobol' index scaled back into the units of
the readout's variance**, and it is computed from exactly the same sample as
:func:`flylab.analysis.uncertainty_global.sobol_analysis` -- no extra model
evaluations.  ``S_j`` is the right index here and ``T_j`` is not: resolving one
factor removes its main effect but leaves its interactions with everything that
stays unknown.  ``T_j * Var(Y)`` is reported alongside as
``voi_upper_bound_var``: the variance that would go if resolving ``j`` also
resolved every interaction it takes part in, which is what you would get by
resolving ``j`` *last*.

Each factor is then mapped to a **concrete experiment**, so the ranking answers
the question a reviewer actually asks: *what should be measured next?*

The honest caveat, repeated in every result's ``warnings``: this is variance of
a **model output** under **assumed** input ranges.  A factor whose range was
declared generously will look informative, and a real measurement will not
collapse a factor to a point.  Nothing here is a statement about biological
variability between flies.
"""
from __future__ import annotations

import time
from typing import Any, Mapping, Sequence

from flylab.analysis.uncertainty_global import FACTOR_NAMES, sobol_analysis

__all__ = [
    "EXPERIMENTS",
    "VOI_WARNING",
    "voi",
    "value_of_information",
    "to_markdown",
]

VOI_WARNING = (
    "VOI here is over MODEL variance, not biological variance: it says how much "
    "of the spread of a simulated readout would disappear if an assumption were "
    "resolved, under the input ranges declared in "
    "flylab.analysis.uncertainty_global.FACTORS. It is not an expected gain in "
    "predictive accuracy about a living fly, and a real experiment reduces a "
    "factor's range rather than collapsing it to a point, so every VOI below is "
    "an upper bound on what that experiment would actually buy."
)

#: factor -> the experiment that would resolve it.  ``blocking`` marks the
#: experiments FlyLab's own roadmap already names as gates.
EXPERIMENTS: dict[str, dict[str, Any]] = {
    "potency": {
        "experiment": "electrophysiological dose-response on the stated receptor",
        "design": (
            "Two-electrode voltage clamp or patch on heterologously expressed "
            "Drosophila receptor of the dominant target (for the neonicotinoids, "
            "a subunit-resolved nAChR; for fipronil, RDL), 8-10 concentrations "
            "spanning 4 log units, n >= 5 cells, fitted with a 4-parameter Hill."
        ),
        "resolves": "the sourced potency value and its parameter type (EC50 vs IC50 vs Kd)",
        "why_now": (
            "at a saturating dose the model is exactly insensitive to this "
            "parameter, so a single high dose constrains nothing; only a ladder does"
        ),
        "cost": "medium",
        "feasibility": "standard method, published protocols exist",
        "blocking_gate": None,
    },
    "hill_n": {
        "experiment": "the same dose-response, read for slope rather than midpoint",
        "design": (
            "The Hill coefficient comes free with the dose-response above, but "
            "needs denser sampling around the midpoint (>= 4 points within one "
            "log unit) and per-cell rather than pooled fits."
        ),
        "resolves": "the cooperativity of the engagement curve",
        "why_now": "shares its cost with the potency experiment; never worth its own study",
        "cost": "low (rides along with the potency ladder)",
        "feasibility": "standard",
        "blocking_gate": None,
    },
    "gain_transform": {
        "experiment": "calibration of receptor engagement against synaptic gain",
        "design": (
            "Paired recording or optogenetically evoked postsynaptic current at a "
            "known cholinergic (or GABAergic) synapse while the compound is "
            "applied at 5-6 concentrations; plot normalised PSC amplitude against "
            "the receptor engagement predicted for that concentration. The shape "
            "of that curve -- monotone, saturating, or rise-then-fall -- IS the "
            "gain rule, which FlyLab currently asserts rather than measures."
        ),
        "resolves": "the functional form of the engagement -> gain transformation",
        "why_now": (
            "the reviewer's central objection: the nicotinic 'buffering' result "
            "follows in large part from the deliberately biphasic nAChR rule, so "
            "this is the one experiment that could convert a modelling choice "
            "into a measurement"
        ),
        "cost": "high",
        "feasibility": "hard but published precedent exists for evoked-PSC pharmacology in Drosophila",
        "blocking_gate": "the single fitted parameter the project allows",
    },
    "gain_coef": {
        "experiment": "the same calibration, read for amplitude rather than shape",
        "design": (
            "The magnitude of the gain change at a fixed, known engagement; one "
            "concentration, many synapses, is enough to set the coefficient once "
            "the shape is fixed."
        ),
        "resolves": "the size of the gain-rule coefficients",
        "why_now": (
            "the one-at-a-time tornado already names this the dominant parameter; "
            "it is the project's one allowed fitted parameter"
        ),
        "cost": "medium",
        "feasibility": "easier than resolving the shape",
        "blocking_gate": "the single fitted parameter the project allows",
    },
    "drive": {
        "experiment": "in vivo baseline firing rates of the driven cell classes",
        "design": (
            "Extracellular or calcium recordings of labellar GRNs and MN9 in the "
            "quiescent and stimulated fly, to set the operating point the model "
            "currently assumes (40 Hz seeds; 150 Hz sugar)."
        ),
        "resolves": "the external drive rate the network is simulated at",
        "why_now": "sets the operating point; no result's direction currently depends on it",
        "cost": "medium",
        "feasibility": "published rates exist for some classes; MN9 is the gap",
        "blocking_gate": None,
    },
    "expression": {
        "experiment": "FISH / scRNA-seq for receptor expression in the named cell class",
        "design": (
            "HCR-FISH for nAChR subunits, Rdl and GluCl on identified adult "
            "labellar motor neurons (MN9) and the labellar GRN types, or a "
            "MaleCNS-registered single-cell atlas that carries superclass labels."
        ),
        "resolves": (
            "which cells actually carry the receptor the drug acts on; FlyLab's "
            "expression coverage is 0.21 and adult motor-neuron receptor "
            "expression -- MN9, the main readout -- has no source at all"
        ),
        "why_now": "a confirmed literature gap, and a data problem rather than a modelling one",
        "cost": "medium",
        "feasibility": "standard method; the MN9 result would be new",
        "blocking_gate": "v0.6 gate 1 (expression coverage)",
    },
    "weight_threshold": {
        "experiment": "no experiment: a reconstruction-confidence analysis",
        "design": (
            "Recompute the readout against MaleCNS synapse-confidence strata, or "
            "against an independent netlist (FlyWire) for the overlapping cells, "
            "instead of against a hand-chosen minimum synapse count."
        ),
        "resolves": "how many weak edges belong in the cut",
        "why_now": "answerable from existing public data, at no bench cost",
        "cost": "none (compute only)",
        "feasibility": "immediate",
        "blocking_gate": None,
    },
    "transmitter": {
        "experiment": "immunostaining or ground-truth transmitter labels",
        "design": (
            "Immunohistochemistry (ChAT / VGlut / GABA) or in-situ labels on the "
            "cell types that carry the most signed weight in the cut, compared "
            "against the MaleCNS predicted consensus transmitter to estimate the "
            "per-class error rate directly."
        ),
        "resolves": "the probability that a predicted transmitter -- hence an edge's sign -- is wrong",
        "why_now": (
            "every edge sign in the model comes from a prediction; the error rate "
            "is currently assumed, not measured"
        ),
        "cost": "medium",
        "feasibility": "standard; published validation sets already exist for some classes",
        "blocking_gate": None,
    },
    "lif_seed": {
        "experiment": "none: increase the number of simulated replicates",
        "design": (
            "Pure simulation noise. Averaging more LIF seeds removes it at no "
            "experimental cost; on the deterministic rate engine it is exactly zero."
        ),
        "resolves": "nothing about biology",
        "why_now": "a compute decision, not an experiment",
        "cost": "none (compute only)",
        "feasibility": "immediate",
        "blocking_gate": None,
    },
}


def voi(
    factors: Sequence[str] | None = None,
    result: Mapping[str, Any] | None = None,
    compound: str = "imidacloprid",
    conc_M: float = 1e-6,
    readout: str = "mean_hz",
    n_base: int = 128,
    seed: int = 0,
    **kw: Any,
) -> dict[str, Any]:
    """Rank the named factors by the model variance a decisive experiment would remove.

    ``result`` reuses an existing
    :func:`flylab.analysis.uncertainty_global.sobol_analysis` output instead of
    re-running it (the VOI needs no extra model evaluations); otherwise one is
    run with the arguments given.  ``factors`` restricts the ranking.

    Returns the ranked table, the recommended next experiment in plain
    language, and the standing :data:`VOI_WARNING`.
    """
    t0 = time.perf_counter()
    res = dict(result) if result is not None else sobol_analysis(
        compound=compound, conc_M=conc_M, readout=readout, n_base=n_base, seed=seed, **kw
    )
    var = float(res["output_variance"])
    sd = float(res.get("output_sd") or (var ** 0.5 if var > 0 else 0.0))
    wanted = set(factors) if factors is not None else set(FACTOR_NAMES)

    rows: list[dict[str, Any]] = []
    for r in res["rows"]:
        name = r["factor"]
        if name not in wanted:
            continue
        exp = EXPERIMENTS.get(name, {})
        v = float(r["first_order"]) * var
        v_hi = float(r["total_order"]) * var
        ci = r.get("first_order_ci")
        rows.append(
            {
                "factor": name,
                "voi_var": v,
                "voi_fraction": float(r["first_order"]),
                "voi_ci_var": [ci[0] * var, ci[1] * var] if ci else None,
                "voi_upper_bound_var": v_hi,
                "voi_upper_bound_fraction": float(r["total_order"]),
                # how much of the readout's SD the experiment would remove
                "residual_sd": float(max(0.0, var * (1.0 - float(r["first_order"])))) ** 0.5,
                "sd_reduction": sd - float(max(0.0, var * (1.0 - float(r["first_order"])))) ** 0.5,
                "experiment": exp.get("experiment"),
                "design": exp.get("design"),
                "resolves": exp.get("resolves"),
                "why_now": exp.get("why_now"),
                "cost": exp.get("cost"),
                "feasibility": exp.get("feasibility"),
                "blocking_gate": exp.get("blocking_gate"),
            }
        )
    rows.sort(key=lambda r: -r["voi_var"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i

    bench = [r for r in rows if r.get("cost") not in ("none (compute only)",)]
    top = bench[0] if bench else (rows[0] if rows else None)
    free = [r for r in rows if r.get("cost") == "none (compute only)" and r["voi_fraction"] > 0]

    recommendation = None
    if top is not None:
        recommendation = (
            f"Next experiment: {top['experiment']}. Resolving {top['factor']} would "
            f"remove about {100 * top['voi_fraction']:.0f}% of the model variance of "
            f"{res['readout']} for {res['compound']} at {res['conc_M']:g} M "
            f"({top['voi_var']:.3g} of {var:.3g} Hz^2; the readout's SD would fall "
            f"from {sd:.3g} to {top['residual_sd']:.3g} Hz). {top['why_now']}."
        )
        if free:
            recommendation += (
                " Before spending bench time, note that "
                + ", ".join(f"{r['factor']} ({100 * r['voi_fraction']:.0f}%)" for r in free)
                + " can be reduced from existing data or compute alone."
            )

    warnings = [VOI_WARNING] + list(res.get("warnings", []))
    warnings.append(
        "VOI_j uses the FIRST-ORDER index: resolving one factor leaves its "
        "interactions with the factors that stay unknown. voi_upper_bound_var "
        "uses the total-order index and is what you would get by resolving that "
        "factor last, after everything else."
    )
    if res.get("interaction_share", 0.0) > 0.5:
        warnings.append(
            f"interactions carry {100 * res['interaction_share']:.0f}% of the "
            "variance, so no single experiment removes much on its own; the "
            "ranking is about ordering, not about absolute gains."
        )

    return {
        "compound": res["compound"],
        "conc_M": res["conc_M"],
        "readout": res["readout"],
        "engine": res.get("engine"),
        "n_base": res.get("n_base"),
        "n_evaluations": res.get("n_evaluations"),
        "output_variance": var,
        "output_sd": sd,
        "interaction_share": res.get("interaction_share"),
        "rows": rows,
        "recommendation": recommendation,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }


def to_markdown(result: Mapping[str, Any]) -> str:
    """The VOI ranking as a Markdown table plus the recommendation."""
    lines = [
        f"### Value of information - {result['compound']} {result['conc_M']:g} M, "
        f"readout {result['readout']} (Var(Y) = {result['output_variance']:.4g})",
        "",
        "| # | factor | VOI (share of Var Y) | VOI (Hz^2) | upper bound | experiment that resolves it | cost |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in result["rows"]:
        lines.append(
            f"| {r['rank']} | {r['factor']} | {r['voi_fraction']:+.3f} | "
            f"{r['voi_var']:+.3g} | {r['voi_upper_bound_fraction']:+.3f} | "
            f"{r['experiment']} | {r['cost']} |"
        )
    if result.get("recommendation"):
        lines += ["", f"**{result['recommendation']}**"]
    lines += ["", f"_{VOI_WARNING}_"]
    return "\n".join(lines)


#: Alias exported by :mod:`flylab.analysis` under a name that does not shadow
#: this module.  ``from flylab.analysis import voi`` would otherwise bind the
#: function over the submodule and break ``flylab.analysis.voi.to_markdown``.
value_of_information = voi
