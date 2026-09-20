"""Mechanism-rule robustness: is a conclusion a property of the pharmacology
and the connectome, or of the gain function we happened to write down?

Why this module exists
----------------------
The one-at-a-time tornado in :func:`flylab.assays.ensemble.sensitivity` reports
that ``gain_coef`` dominates every other parameter at 1 uM (span 2.24 Hz
against 0.69 Hz for the drive and exactly 0 Hz for the EC50 and the Hill
coefficient).  The gain rules in :mod:`flylab.pharm.mechanisms` -- ``1 + 0.4t -
1.6t^2``, ``1 + 1.5t``, ``1 + 2t`` -- are modelling choices, not measured
transformations.  A reviewer is therefore right to ask whether the headline
phenomena are properties of *those curves* rather than of the
pharmacology/connectome integration.

This module answers that question the only way it can be answered: by
re-running the principal qualitative conclusions across a **prespecified**
family of alternative gain rules and reporting, per conclusion, the fraction of
admissible specifications that retain it -- the **Conclusion Stability
Matrix**.  The paper can then write "X persists across 88% of 25 prespecified
model specifications, failing only under <named specs>" instead of "our model
predicts X".

Design constraints honoured here
--------------------------------
* :mod:`flylab.pharm.mechanisms` is **never edited**.  The default rules stay
  numerically frozen; ``tests/test_pharm_mechanisms.py`` still passes.  An
  alternative specification is installed for the duration of a ``with`` block
  by rebinding the ``compute_gains`` symbol that the assay modules imported
  (the same trick ``flylab.assays.ensemble._taste_library`` already uses for
  the library).  The default member of the family *delegates* to
  :func:`flylab.pharm.mechanisms.gains_from_occupancy`, so "the default
  specification reproduces the current numbers" is true by construction, not by
  luck.
* Engagement (occupancy) of ``None`` means **not modelled**, never 0.  Such a
  row is skipped and recorded in ``not_modelled``; a conclusion whose dominant
  receptor is unmodelled evaluates to ``None`` (undecidable), not ``False``.
* Every family member is bounded (``[0.05, 10.0]``), respects the 0.05 floor
  and reduces to 1.0 as engagement -> 0, so "vehicle" always means "no change".

Runtime
-------
``conclusion_stability()`` with the full 25-member family takes roughly
**2-4 minutes** on the committed graphs (about 5-6 s per specification:
2 rate assays, 1 taste-map assay, 2 rewiring nulls at ``n_shuffles`` shuffles
and 10 circuit-selectivity ladders).  ``conclusion_stability(fast=True)`` uses
the 9-member subsample and reduced compound sets and takes roughly 30-50 s.
"""
from __future__ import annotations

import contextlib
import math
import threading
import time
from typing import Any, Callable, Iterable, Mapping, Sequence

__all__ = [
    "GAIN_FLOOR",
    "GAIN_CEIL",
    "SHAPES",
    "COEFFICIENT_SCALES",
    "RULE_SPECS",
    "MECHANISM_FAMILY",
    "DEFAULT_SPEC_NAME",
    "FAST_SHAPES",
    "FAST_COEFFICIENT_SCALES",
    "NICOTINIC_SET",
    "AMPLIFY_SET",
    "CONCLUSIONS",
    "receptor_si",
    "vertebrate_threshold",
    "subsample_family",
    "shape_function",
    "make_spec",
    "spec_gains",
    "spec_by_name",
    "default_spec",
    "family_table",
    "mechanism_spec",
    "SpecRun",
    "conclusion_stability",
    "threshold_sensitivity",
    "to_markdown",
]

# --------------------------------------------------------------------------
# bounds
# --------------------------------------------------------------------------
#: the 0.05 floor of :mod:`flylab.pharm.mechanisms` (a dead network is not
#: informative).  Every family member respects it.
GAIN_FLOOR = 0.05

#: hard ceiling.  The default rules top out at ``ach_tone = 3.0``; the most
#: aggressive family member (``1 + 3.0*t`` for AChE tone at coefficient 1.5x)
#: reaches 4.0, so 10.0 never binds on the default rules and only guards
#: against a runaway specification.
GAIN_CEIL = 10.0


def _clip(x: float) -> float:
    return min(max(float(x), GAIN_FLOOR), GAIN_CEIL)


# --------------------------------------------------------------------------
# the prespecified transformation family
# --------------------------------------------------------------------------
#: Half-saturation engagement of the ``saturating`` shape.  0.25 means the
#: transformation is half-way to its endpoint at 25% receptor engagement, which
#: is the shape a strongly amplifying transduction step would have.
SATURATING_K = 0.25

#: Alternative monotone / biphasic transformations, prespecified before any
#: result was looked at.  Each maps engagement ``t`` in [0, 1] and a rule's
#: ``(sign s, linear coefficient a, quadratic coefficient b)`` onto a gain.
#:
#: All five reduce to ``g = 1`` at ``t = 0`` (vehicle is vehicle), are bounded
#: in ``[GAIN_FLOOR, GAIN_CEIL]``, and are continuous in ``t``.
SHAPES: dict[str, dict[str, Any]] = {
    "linear": {
        "name": "linear",
        "formula": "g = 1 + s*a*t",
        "rationale": (
            "The null transformation: gain moves in proportion to receptor "
            "engagement. For a blocker this is g = 1 - a*t, the textbook "
            "'a competitive blocker removes that fraction of drive'."
        ),
        "biphasic": False,
    },
    "saturating": {
        "name": "saturating",
        "formula": f"g = 1 + s*a*(1+k)*t/(k+t), k={SATURATING_K}",
        "rationale": (
            "Receptor -> synaptic gain passes through a saturating "
            "transduction step, so most of the effect is spent by the time "
            "engagement is moderate. Normalised to the same endpoint as "
            "'linear' at t = 1 so shape, not magnitude, is what varies."
        ),
        "biphasic": False,
    },
    "weak_biphasic": {
        "name": "weak_biphasic",
        "formula": "g = 1 + s*a*t - 0.5*b_ref*t^2",
        "rationale": (
            "Half of FlyLab's desensitisation/depolarisation-block term. For "
            "the nAChR agonist rule the gain still turns over, but later and "
            "less deeply; for monotone rules the quadratic term acts as a mild "
            "supra-linear steepening rather than a true rise-then-fall."
        ),
        "biphasic": True,
    },
    "flylab_biphasic": {
        "name": "flylab_biphasic",
        "formula": "g = 1 + s*a*t - b*t^2   (FlyLab default; b = 0 for monotone rules)",
        "rationale": (
            "The shipped FlyLab rule set, verbatim: nAChR agonist "
            "1 + 0.4t - 1.6t^2, every other rule monotone (b = 0). Included as "
            "a member of the family so the default is testable, not privileged."
        ),
        "biphasic": True,
    },
    "strong_biphasic": {
        "name": "strong_biphasic",
        "formula": "g = 1 + s*a*t - 2*b_ref*t^2",
        "rationale": (
            "Twice FlyLab's desensitisation term: the nAChR agonist rule turns "
            "over early and hits the 0.05 floor well before full engagement, "
            "which is what a strongly depolarisation-blocking neonicotinoid "
            "would look like."
        ),
        "biphasic": True,
    },
}

#: Crossed with every shape.  1.0 is the shipped coefficient set.
COEFFICIENT_SCALES: tuple[float, ...] = (0.5, 0.75, 1.0, 1.25, 1.5)

#: The subsample used by ``fast=True``: the two monotone extremes, the default,
#: and the coefficient extremes plus the default.
FAST_SHAPES: tuple[str, ...] = ("linear", "saturating", "flylab_biphasic")
FAST_COEFFICIENT_SCALES: tuple[float, ...] = (0.5, 1.0, 1.5)

#: ``(receptor, direction-class) -> (gain key, sign, a, b)`` decomposition of
#: every rule in :data:`flylab.pharm.mechanisms.MECHANISM_TABLE`.
#:
#: ``sign`` is +1 for a rule that raises the gain and -1 for one that lowers
#: it; ``a`` is the linear coefficient and ``b`` the quadratic
#: (desensitisation) coefficient of the shipped rule.  ``b_ref`` is the
#: magnitude the ``weak_biphasic`` / ``strong_biphasic`` shapes scale: FlyLab's
#: own ``b`` where there is one, otherwise the linear coefficient, so that the
#: biphasic axis is non-degenerate for every rule.
RULE_SPECS: dict[tuple[str, str], dict[str, Any]] = {
    ("insect_nAChR", "activating"): {"gain": "g_ach", "s": +1.0, "a": 0.4, "b": 1.6},
    ("insect_nAChR", "blocking"): {"gain": "g_ach", "s": -1.0, "a": 1.0, "b": 0.0},
    ("insect_RDL", "blocking"): {"gain": "g_gaba", "s": -1.0, "a": 1.0, "b": 0.0},
    ("insect_RDL", "activating"): {"gain": "g_gaba", "s": +1.0, "a": 0.4, "b": 0.0},
    ("insect_GluCl", "activating"): {"gain": "g_glu", "s": +1.0, "a": 0.8, "b": 0.0},
    ("insect_GluCl", "blocking"): {"gain": "g_glu", "s": -1.0, "a": 1.0, "b": 0.0},
    ("insect_Nav", "activating"): {"gain": "g_nav", "s": +1.0, "a": 1.5, "b": 0.0},
    ("insect_Nav", "blocking"): {"gain": "g_nav", "s": -1.0, "a": 1.0, "b": 0.0},
    ("insect_OctR", "activating"): {"gain": "g_oct", "s": +1.0, "a": 0.5, "b": 0.0},
    ("insect_OctR", "blocking"): {"gain": "g_oct", "s": -1.0, "a": 1.0, "b": 0.0},
    ("insect_AChE", "blocking"): {"gain": "ach_tone", "s": +1.0, "a": 2.0, "b": 0.0},
}

_ACTIVATING = {"agonist", "partial_agonist", "positive_modulator"}
_BLOCKING = {"antagonist", "negative_modulator", "inhibitor"}


def _direction_class(direction: str | None) -> str | None:
    d = str(direction or "none")
    if d in _ACTIVATING:
        return "activating"
    if d in _BLOCKING:
        return "blocking"
    return None


def shape_function(
    shape: str,
    s: float,
    a: float,
    b: float,
    coef: float = 1.0,
    b_ref: float | None = None,
) -> Callable[[float], float]:
    """Unclipped transformation ``t -> g`` for one rule under one shape.

    ``coef`` scales both the linear and the quadratic coefficient (the
    coefficient ensemble axis).  ``b_ref`` defaults to ``b`` when the rule has a
    quadratic term and to ``a`` when it does not, which keeps the biphasic axis
    meaningful for monotone rules.

    The returned function is *not* clipped; :func:`spec_gains` applies the
    floor and ceiling exactly where :mod:`flylab.pharm.mechanisms` does.
    """
    if shape not in SHAPES:
        raise ValueError(f"unknown shape {shape!r}; expected one of {tuple(SHAPES)}")
    aa = float(a) * float(coef)
    bb = float(b) * float(coef)
    ref = (float(b) if float(b) > 0.0 else float(a)) if b_ref is None else float(b_ref)
    ref *= float(coef)
    sgn = float(s)

    if shape == "linear":
        return lambda t: 1.0 + sgn * aa * t
    if shape == "saturating":
        k = SATURATING_K
        return lambda t: 1.0 + sgn * aa * (1.0 + k) * t / (k + t)
    if shape == "weak_biphasic":
        return lambda t: 1.0 + sgn * aa * t - 0.5 * ref * t * t
    if shape == "flylab_biphasic":
        return lambda t: 1.0 + sgn * aa * t - bb * t * t
    # strong_biphasic
    return lambda t: 1.0 + sgn * aa * t - 2.0 * ref * t * t


def _rule_formula(receptor: str, cls: str, shape: str, coef: float) -> str:
    spec = RULE_SPECS[(receptor, cls)]
    a = spec["a"] * coef
    b = spec["b"] * coef
    ref = (spec["b"] if spec["b"] > 0 else spec["a"]) * coef
    sgn = "+" if spec["s"] > 0 else "-"
    key = spec["gain"]
    if shape == "linear":
        body = f"1 {sgn} {a:g}t"
    elif shape == "saturating":
        body = f"1 {sgn} {a * (1 + SATURATING_K):g}t/({SATURATING_K:g}+t)"
    elif shape == "weak_biphasic":
        body = f"1 {sgn} {a:g}t - {0.5 * ref:g}t^2"
    elif shape == "flylab_biphasic":
        body = f"1 {sgn} {a:g}t" + (f" - {b:g}t^2" if b else "")
    else:
        body = f"1 {sgn} {a:g}t - {2 * ref:g}t^2"
    return f"{key} = max({GAIN_FLOOR}, {body})"


def _make_spec(shape: str, coef: float) -> dict[str, Any]:
    is_default = shape == "flylab_biphasic" and coef == 1.0
    rules = {
        f"{receptor}:{cls}": _rule_formula(receptor, cls, shape, coef)
        for (receptor, cls) in RULE_SPECS
    }
    return {
        "name": f"{shape}@{coef:.2f}x",
        "shape": shape,
        "coef_scale": float(coef),
        "formula": SHAPES[shape]["formula"],
        "rationale": SHAPES[shape]["rationale"],
        "biphasic": bool(SHAPES[shape]["biphasic"]),
        "params": {
            "coef_scale": float(coef),
            "saturating_k": SATURATING_K,
            "floor": GAIN_FLOOR,
            "ceiling": GAIN_CEIL,
        },
        "rules": rules,
        "is_default": is_default,
        "delegates_to_default": is_default,
    }


def make_spec(shape: str, coef: float = 1.0) -> dict[str, Any]:
    """One family member built on demand (``shape`` x ``coef``)."""
    return _make_spec(shape, float(coef))


#: The prespecified specification family: 5 transformations x 5 coefficient
#: scales = 25 admissible specifications, exposed as data so the paper can
#: table it directly (``family_table()``).
MECHANISM_FAMILY: list[dict[str, Any]] = [
    _make_spec(shape, coef) for shape in SHAPES for coef in COEFFICIENT_SCALES
]

DEFAULT_SPEC_NAME = "flylab_biphasic@1.00x"


def default_spec() -> dict[str, Any]:
    """The shipped FlyLab rule set as a member of the family."""
    return spec_by_name(DEFAULT_SPEC_NAME)


def spec_by_name(name: str, family: Sequence[Mapping[str, Any]] | None = None) -> dict[str, Any]:
    for spec in family if family is not None else MECHANISM_FAMILY:
        if spec["name"] == name:
            return dict(spec)
    raise KeyError(f"unknown specification {name!r}")


def subsample_family(
    shapes: Sequence[str] = FAST_SHAPES,
    coefs: Sequence[float] = FAST_COEFFICIENT_SCALES,
) -> list[dict[str, Any]]:
    """A documented subsample of :data:`MECHANISM_FAMILY` (the ``--fast`` set).

    The default subsample keeps the two monotone extremes and the shipped
    biphasic rule, crossed with the coefficient extremes and 1.0x: 9 of the 25
    specifications, chosen to span the corners of the family rather than to
    sample it uniformly.
    """
    return [
        s for s in MECHANISM_FAMILY if s["shape"] in tuple(shapes) and s["coef_scale"] in tuple(coefs)
    ]


def family_table(family: Sequence[Mapping[str, Any]] | None = None) -> list[dict[str, Any]]:
    """:data:`MECHANISM_FAMILY` flattened for a paper table."""
    fam = family if family is not None else MECHANISM_FAMILY
    return [
        {
            "spec": s["name"],
            "shape": s["shape"],
            "coefficient_scale": s["coef_scale"],
            "formula": s["formula"],
            "nachr_agonist_rule": s["rules"]["insect_nAChR:activating"],
            "rdl_antagonist_rule": s["rules"]["insect_RDL:blocking"],
            "nav_rule": s["rules"]["insect_Nav:activating"],
            "is_default": s["is_default"],
            "rationale": s["rationale"],
        }
        for s in fam
    ]


# --------------------------------------------------------------------------
# gains under a specification
# --------------------------------------------------------------------------
_GAIN_KEYS = ("g_ach", "g_gaba", "g_glu", "g_oct", "g_nav", "ach_tone")


def _engagement(row: Mapping[str, Any]) -> float | None:
    """Receptor engagement of a row, or ``None`` when it is *not modelled*.

    ``None`` is never silently turned into 0: an unsupported row means the
    model has nothing to say about that receptor, which is different from
    saying the drug does not engage it.
    """
    value = row.get("occupancy", row.get("engagement"))
    if value is None:
        return None
    try:
        t = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(t):
        return None
    return min(max(t, 0.0), 1.0)


def spec_gains(
    rows: Iterable[Mapping[str, Any]],
    spec: Mapping[str, Any] | None = None,
    *,
    force_generic: bool = False,
) -> dict[str, float]:
    """Gain dict for occupancy ``rows`` under one specification.

    The default specification delegates to
    :func:`flylab.pharm.mechanisms.gains_from_occupancy`, so it is exact by
    construction.  ``force_generic=True`` runs the generic family code on the
    default parameters instead, which is what the freeze test compares against.
    """
    from flylab.pharm.mechanisms import default_gains, gains_from_occupancy

    spec = dict(spec or default_spec())
    rows = list(rows or [])
    if spec.get("delegates_to_default") and not force_generic:
        return gains_from_occupancy(rows)

    shape = spec["shape"]
    coef = float(spec["coef_scale"])
    g = default_gains()
    ache_t: float | None = None

    for row in rows:
        receptor = str(row.get("receptor", ""))
        cls = _direction_class(row.get("direction"))
        if cls is None or (receptor, cls) not in RULE_SPECS:
            continue
        t = _engagement(row)
        if t is None:  # not modelled -> no patch, never t = 0 by default
            continue
        rule = RULE_SPECS[(receptor, cls)]
        if receptor == "insect_AChE":
            ache_t = t if ache_t is None else max(ache_t, t)
            continue
        f = shape_function(shape, rule["s"], rule["a"], rule["b"], coef)
        g[rule["gain"]] = _clip(f(t))

    if ache_t is not None:
        tone_rule = RULE_SPECS[("insect_AChE", "blocking")]
        tone = _clip(shape_function(shape, tone_rule["s"], tone_rule["a"], tone_rule["b"], coef)(ache_t))
        agon = RULE_SPECS[("insect_nAChR", "activating")]
        curve = shape_function(shape, agon["s"], agon["a"], agon["b"], coef)(ache_t)
        g["ach_tone"] = tone
        g["g_ach"] = _clip(g["g_ach"] * tone * curve)

    return {k: float(g[k]) for k in _GAIN_KEYS}


# --------------------------------------------------------------------------
# receptor-side quantities, None-safe under library schema v3
# --------------------------------------------------------------------------
# Schema v3 reports ``engagement: None`` and ``param_value_M: None`` for rows
# with no sourced evidence, and refuses a selectivity ratio whenever either
# side of a pair is unsupported.  The helpers below therefore recompute the
# receptor half of the selectivity index from the rows themselves, skipping
# not-modelled rows instead of letting a missing value become a number.  The
# circuit half still comes from
# :func:`flylab.analysis.selectivity.circuit_threshold_conc`, which is the part
# that actually depends on the gain rules.


def _param_value(row: Mapping[str, Any]) -> float | None:
    for key in ("param_value_M", "value_M", "ec50_M"):
        v = row.get(key)
        if v is None:
            continue
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if math.isfinite(f) and f > 0:
            return f
    return None


def _rows_by_receptor(compound: str, conc_M: float, library: Mapping[str, Any] | None = None):
    from flylab.pharm.occupancy import compare_compound

    occ = (
        compare_compound(compound, conc_M, library=dict(library))
        if library
        else compare_compound(compound, conc_M)
    )
    return {str(r.get("receptor")): r for r in occ["receptors"]}, occ


def receptor_si(
    compound: str,
    conc_M: float = 1e-6,
    library: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """``log10(potency_vertebrate / potency_insect)`` of the best sourced pair.

    A pair is used only when *both* sides carry a sourced value; pairs whose
    two sides are different kinds of parameter (an insect IC50 against a
    vertebrate Kd, say) are still reported but flagged ``comparable: False``
    and lose to a type-matched pair.  Returns ``si=None`` when no pair is
    modellable, which is different from an SI of 0.
    """
    from flylab.pharm.occupancy import load_library, selectivity_pairs

    rows, _ = _rows_by_receptor(compound, conc_M, library)
    pairs = selectivity_pairs(dict(library) if library else load_library())
    best: dict[str, Any] | None = None
    for pair_name, pair in pairs.items():
        ins = rows.get(pair["insect"])
        vert = rows.get(pair["vertebrate"])
        if ins is None or vert is None:
            continue
        iv, vv = _param_value(ins), _param_value(vert)
        if iv is None or vv is None:
            continue
        cand = {
            "pair": pair_name,
            "si": float(math.log10(vv / iv)),
            "insect_receptor": pair["insect"],
            "vertebrate_receptor": pair["vertebrate"],
            "insect_value_M": iv,
            "vertebrate_value_M": vv,
            "comparable": ins.get("param_type") == vert.get("param_type"),
        }
        if best is None:
            best = cand
        elif (cand["comparable"], cand["si"]) > (best["comparable"], best["si"]):
            best = cand
    if best is None:
        return {"pair": None, "si": None, "comparable": None, "reason": "no modellable insect/vertebrate pair"}
    return best


def vertebrate_threshold(
    compound: str,
    occ_limit: float = 0.2,
    library: Mapping[str, Any] | None = None,
) -> float | None:
    """Lowest concentration at which any sourced vertebrate row reaches
    ``occ_limit`` engagement (closed form ``C = v * (t/(1-t))**(1/n)``)."""
    rows, _ = _rows_by_receptor(compound, 1e-6, library)
    t = min(max(float(occ_limit), 1e-9), 1.0 - 1e-9)
    best: float | None = None
    for name, row in rows.items():
        if not name.startswith("vertebrate_"):
            continue
        v = _param_value(row)
        if v is None:
            continue
        try:
            n = float(row.get("n", 1.0)) or 1.0
        except (TypeError, ValueError):
            n = 1.0
        c = v * (t / (1.0 - t)) ** (1.0 / n)
        best = c if best is None else min(best, c)
    return best


_SPEC_LOCK = threading.Lock()

#: assay modules that imported ``compute_gains`` at module scope
_PATCH_TARGETS = (
    "flylab.assays.subgraph",
    "flylab.assays.taste_map",
    "flylab.assays.spiking",
)


@contextlib.contextmanager
def mechanism_spec(spec: Mapping[str, Any] | None):
    """Install an alternative gain specification for the duration of a block.

    :mod:`flylab.pharm.mechanisms` is not touched.  The ``compute_gains``
    symbol that the assay modules imported is rebound to a wrapper that maps
    occupancy rows through ``spec`` instead; everything downstream (rate
    assay, taste map, spiking, null models, selectivity ladders) then runs
    under that specification.  Passing ``None`` or the default specification is
    a no-op, so the shipped numbers are reproduced byte-for-byte.

    Not re-entrant across threads: a module-level lock is held for the block.
    """
    import importlib

    if spec is None or dict(spec).get("delegates_to_default"):
        yield
        return

    from flylab.circuit.rate import DEFAULT_GAINS

    def patched(compound, conc_M=0.0, library=None, rule_overrides=None):
        if not compound:
            return dict(DEFAULT_GAINS), None
        from flylab.pharm.occupancy import compare_compound

        occ = (
            compare_compound(compound, conc_M, library=library)
            if library
            else compare_compound(compound, conc_M)
        )
        return spec_gains(occ["receptors"], spec), occ

    modules = [importlib.import_module(m) for m in _PATCH_TARGETS]
    with _SPEC_LOCK:
        originals = [getattr(m, "compute_gains") for m in modules]
        for m in modules:
            m.compute_gains = patched
        try:
            yield
        finally:
            for m, orig in zip(modules, originals):
                m.compute_gains = orig


# --------------------------------------------------------------------------
# one specification's runs
# --------------------------------------------------------------------------
#: nicotinic compounds used for the "selectivity buffering" conclusion
NICOTINIC_SET: tuple[str, ...] = (
    "imidacloprid",
    "nicotine",
    "acetamiprid",
    "clothianidin",
    "nitenpyram",
    "spinosad",
)

#: Nav / AChE compounds used for the "selectivity amplification" conclusion
AMPLIFY_SET: tuple[str, ...] = ("ddt", "deltamethrin", "permethrin", "chlorpyrifos_oxon")

FAST_NICOTINIC_SET: tuple[str, ...] = ("imidacloprid", "nicotine", "clothianidin")
FAST_AMPLIFY_SET: tuple[str, ...] = ("deltamethrin", "chlorpyrifos_oxon")

PAPER_CONC_M = 1e-6


class SpecRun:
    """Every assay one specification needs, computed once and cached.

    All work happens inside :func:`mechanism_spec`, so a ``SpecRun`` is a
    complete re-derivation of the manuscript's readouts under one alternative
    gain rule.
    """

    def __init__(
        self,
        spec: Mapping[str, Any] | None = None,
        conc_M: float = PAPER_CONC_M,
        seed: int = 0,
        n_shuffles: int = 6,
        nicotinic: Sequence[str] = NICOTINIC_SET,
        amplify: Sequence[str] = AMPLIFY_SET,
        concs: Sequence[float] | None = None,
        threshold_frac: float = 0.5,
        vert_occ_limit: float = 0.2,
    ):
        self.spec = dict(spec or default_spec())
        self.conc_M = float(conc_M)
        self.seed = int(seed)
        self.n_shuffles = int(n_shuffles)
        self.nicotinic = tuple(nicotinic)
        self.amplify = tuple(amplify)
        self.concs = tuple(concs) if concs is not None else None
        self.threshold_frac = float(threshold_frac)
        self.vert_occ_limit = float(vert_occ_limit)
        self._cache: dict[tuple, Any] = {}
        self.notes: list[str] = []

    # -- primitives ------------------------------------------------------
    def _memo(self, key: tuple, fn: Callable[[], Any]) -> Any:
        if key not in self._cache:
            self._cache[key] = fn()
        return self._cache[key]

    def gains(self, compound: str) -> dict[str, float]:
        def go():
            from flylab.pharm.occupancy import compare_compound

            occ = compare_compound(compound, self.conc_M)
            return spec_gains(occ["receptors"], self.spec)

        return self._memo(("gains", compound), go)

    def subgraph(self, compound: str | None) -> dict[str, Any]:
        def go():
            from flylab.assays.subgraph import run_subgraph_assay

            with mechanism_spec(self.spec):
                return run_subgraph_assay(compound, self.conc_M if compound else 0.0)

        return self._memo(("subgraph", compound), go)

    def taste_map(self, compound: str | None) -> dict[str, Any]:
        def go():
            from flylab.assays.taste_map import run_taste_map_assay

            with mechanism_spec(self.spec):
                return run_taste_map_assay(compound, self.conc_M if compound else 0.0, seed=self.seed)

        return self._memo(("taste_map", compound), go)

    def topology_z(self, compound: str, mode: str = "rewire_degree_preserving") -> float | None:
        """|z| of the drug effect against a degree-preserving rewiring null."""

        def go():
            from flylab.analysis.nullmodels import null_distribution

            with mechanism_spec(self.spec):
                r = null_distribution(
                    "subgraph", compound, self.conc_M, mode=mode,
                    n=self.n_shuffles, seed=self.seed, readout="auto",
                )
            z = r.get("z")
            return None if z is None else abs(float(z))

        return self._memo(("topo", compound, mode), go)

    def circuit_threshold(self, compound: str) -> float | None:
        """Lowest concentration whose relative ``mean_hz`` change reaches
        ``threshold_frac`` -- the only half of the selectivity index that the
        gain rules can move."""

        def go():
            from flylab.analysis.selectivity import circuit_threshold_conc

            with mechanism_spec(self.spec):
                r = circuit_threshold_conc(
                    "subgraph", compound, readout="auto",
                    threshold_frac=self.threshold_frac, concs=self.concs,
                )
            return r.get("conc_M")

        return self._memo(("circ", compound), go)

    def si_gap(self, compound: str) -> float | None:
        """``circuit SI - receptor SI`` for one compound (``None`` if unscored).

        ``circuit SI = log10(C_vert) - log10(C_circuit)``; a positive gap means
        the neighborhood *amplifies* the receptor-level selectivity margin, a
        negative one that it *buffers* it.
        """

        def go():
            c_circ = self.circuit_threshold(compound)
            c_vert = vertebrate_threshold(compound, self.vert_occ_limit)
            rec = receptor_si(compound)
            if not c_circ or not c_vert or rec.get("si") is None:
                return None
            return float(math.log10(c_vert) - math.log10(c_circ) - float(rec["si"]))

        return self._memo(("si", compound), go)

    def mean_gap(self, compounds: Sequence[str]) -> float | None:
        vals = [v for v in (self.si_gap(c) for c in compounds) if v is not None]
        if not vals:
            return None
        return float(sum(vals) / len(vals))


# --------------------------------------------------------------------------
# the manuscript's principal qualitative claims, as predicates
# --------------------------------------------------------------------------
def _c_nicotinic_suppression(run: SpecRun) -> bool | None:
    nb = run.subgraph("imidacloprid")
    r = nb["readouts"]
    treated, vehicle = r.get("mean_hz"), (r.get("vehicle") or {}).get("mean_hz")
    if treated is None or vehicle is None:
        return None
    return bool(treated < vehicle)


def _c_rdl_disinhibition(run: SpecRun) -> bool | None:
    nb = run.subgraph("fipronil")
    r = nb["readouts"]
    treated, vehicle = r.get("mean_hz"), (r.get("vehicle") or {}).get("mean_hz")
    if treated is None or vehicle is None:
        return None
    return bool(treated > vehicle)


def _c_imidacloprid_topology_weak(run: SpecRun) -> bool | None:
    z = run.topology_z("imidacloprid")
    if z is None:
        return None
    return bool(z < 2.0)


def _c_fipronil_topology_exceeds(run: SpecRun) -> bool | None:
    z_fip = run.topology_z("fipronil")
    z_imi = run.topology_z("imidacloprid")
    if z_fip is None or z_imi is None:
        return None
    return bool(z_fip > z_imi)


def _c_nicotinic_buffering(run: SpecRun) -> bool | None:
    gap = run.mean_gap(run.nicotinic)
    if gap is None:
        return None
    return bool(gap < 0.0)


def _c_nav_ache_amplification(run: SpecRun) -> bool | None:
    gap = run.mean_gap(run.amplify)
    if gap is None:
        return None
    return bool(gap > 0.0)


def _c_map_bitter_veto(run: SpecRun) -> bool | None:
    nb = run.taste_map("fipronil")
    ratio = nb["readouts"].get("bitter_veto_ratio")
    if ratio is None:
        return None
    return bool(float(ratio) < 1.0)


#: The manuscript's principal qualitative conclusions, each a predicate that
#: returns ``True`` (retained), ``False`` (reversed / lost) or ``None``
#: (undecidable under this specification -- e.g. the readout is silenced at the
#: 0.05 floor, so a ratio is undefined).
CONCLUSIONS: dict[str, dict[str, Any]] = {
    "C1_nicotinic_suppression": {
        "claim": "a nicotinic agonist (imidacloprid, 1 uM) suppresses network activity",
        "readout": "subgraph mean_hz, treated < vehicle",
        "evaluate": _c_nicotinic_suppression,
    },
    "C2_rdl_disinhibition": {
        "claim": "an RDL antagonist (fipronil, 1 uM) disinhibits the network",
        "readout": "subgraph mean_hz, treated > vehicle",
        "evaluate": _c_rdl_disinhibition,
    },
    "C3_imidacloprid_topology_weak": {
        "claim": "imidacloprid's mean-rate effect is not a wiring result",
        "readout": "|z| < 2 against a degree-preserving rewiring null",
        "evaluate": _c_imidacloprid_topology_weak,
    },
    "C4_fipronil_topology_exceeds": {
        "claim": "fipronil's topology dependence exceeds imidacloprid's",
        "readout": "|z_fipronil| > |z_imidacloprid| on the same null",
        "evaluate": _c_fipronil_topology_exceeds,
    },
    "C5_nicotinic_buffering": {
        "claim": "the neighborhood buffers nicotinic receptor selectivity",
        "readout": "mean (circuit SI - receptor SI) < 0 over the nicotinic set",
        "evaluate": _c_nicotinic_buffering,
    },
    "C6_nav_ache_amplification": {
        "claim": "the neighborhood amplifies Nav / AChE receptor selectivity",
        "readout": "mean (circuit SI - receptor SI) > 0 over the Nav/AChE set",
        "evaluate": _c_nav_ache_amplification,
    },
    "C7_map_bitter_veto": {
        "claim": "on the map path, bitter drive lowers MN9 under fipronil",
        "readout": "taste_map bitter_veto_ratio < 1",
        "evaluate": _c_map_bitter_veto,
    },
}

BASE_WARNINGS = [
    "model_derived: a stability fraction says how robust a conclusion is to the "
    "chosen gain transformation, not how likely it is to be true in a fly.",
    "The specification family is prespecified and finite; the fraction is a "
    "fraction of THAT family, not a posterior probability.",
    "Every member of the family shares the same connectome cut, the same "
    "teaching EC50 library and the same rate/LIF runtimes; those assumptions "
    "are held fixed here and are attacked elsewhere (null models, "
    "uncertainty_global).",
]


# --------------------------------------------------------------------------
# the Conclusion Stability Matrix
# --------------------------------------------------------------------------
def conclusion_stability(
    conclusions: Mapping[str, Mapping[str, Any]] | None = None,
    family: Sequence[Mapping[str, Any]] | None = None,
    n_jobs: int = 1,
    seed: int = 0,
    fast: bool = False,
    conc_M: float = PAPER_CONC_M,
    n_shuffles: int | None = None,
    nicotinic: Sequence[str] | None = None,
    amplify: Sequence[str] | None = None,
    concs: Sequence[float] | None = None,
    progress: Callable[[str], None] | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Re-derive every conclusion under every admissible specification.

    Returns the **Conclusion Stability Matrix**: for each conclusion, the
    fraction of specifications that retain it and the *names* of the ones that
    do not.

    ``fast=True`` uses :func:`subsample_family` (9 of 25 specifications),
    reduced compound sets and 4 rewiring shuffles: about 30-50 s instead of
    2-4 minutes.

    ``n_jobs`` is accepted for API symmetry but the analysis runs serially: a
    specification is installed by rebinding a process-global symbol, so two
    specifications cannot be in flight at once.  A value > 1 is recorded in
    ``warnings``.
    """
    t0 = time.perf_counter()
    fam = list(family) if family is not None else (subsample_family() if fast else MECHANISM_FAMILY)
    concl = dict(conclusions) if conclusions is not None else CONCLUSIONS
    shuffles = int(n_shuffles if n_shuffles is not None else (4 if fast else 6))
    nic = tuple(nicotinic) if nicotinic is not None else (FAST_NICOTINIC_SET if fast else NICOTINIC_SET)
    amp = tuple(amplify) if amplify is not None else (FAST_AMPLIFY_SET if fast else AMPLIFY_SET)

    warnings = list(BASE_WARNINGS)
    if int(n_jobs) > 1:
        warnings.append(
            f"n_jobs={n_jobs} requested but ignored: an alternative mechanism "
            "specification is installed process-globally, so specifications "
            "are evaluated serially."
        )

    matrix: dict[str, dict[str, bool | None]] = {name: {} for name in concl}
    diagnostics: dict[str, dict[str, Any]] = {}
    errors: list[dict[str, str]] = []
    per_spec_runtime: dict[str, float] = {}

    for spec in fam:
        ts = time.perf_counter()
        if progress:
            progress(spec["name"])
        run = SpecRun(
            spec, conc_M=conc_M, seed=seed, n_shuffles=shuffles,
            nicotinic=nic, amplify=amp, concs=concs, **kw,
        )
        for name, c in concl.items():
            try:
                verdict = c["evaluate"](run)
            except Exception as exc:  # pragma: no cover - defensive
                verdict = None
                errors.append({"spec": spec["name"], "conclusion": name, "error": repr(exc)})
            matrix[name][spec["name"]] = verdict
        diagnostics[spec["name"]] = {
            "gains_imidacloprid": run.gains("imidacloprid"),
            "gains_fipronil": run.gains("fipronil"),
            "mean_gap_nicotinic": run.mean_gap(nic),
            "mean_gap_nav_ache": run.mean_gap(amp),
        }
        per_spec_runtime[spec["name"]] = float(time.perf_counter() - ts)

    rows: list[dict[str, Any]] = []
    for name, c in concl.items():
        verdicts = matrix[name]
        n_true = sum(1 for v in verdicts.values() if v is True)
        n_false = sum(1 for v in verdicts.values() if v is False)
        n_undec = sum(1 for v in verdicts.values() if v is None)
        n = len(verdicts) or 1
        decidable = n_true + n_false
        rows.append(
            {
                "conclusion": name,
                "claim": c["claim"],
                "readout": c.get("readout"),
                "n_specs": len(verdicts),
                "n_retained": n_true,
                "n_lost": n_false,
                "n_undecidable": n_undec,
                # headline: undecidable counts against the conclusion
                "fraction_retained": float(n_true) / float(n),
                # secondary: restricted to specifications where it is decidable
                "fraction_retained_decidable": (float(n_true) / decidable) if decidable else None,
                "failing_specs": sorted(s for s, v in verdicts.items() if v is False),
                "undecidable_specs": sorted(s for s, v in verdicts.items() if v is None),
                "fragile": bool(n_true < n),
            }
        )
    rows.sort(key=lambda r: (r["fraction_retained"], r["conclusion"]))

    if any(r["fragile"] for r in rows):
        warnings.append(
            "fragile conclusions (not retained by every specification): "
            + ", ".join(r["conclusion"] for r in rows if r["fragile"])
        )
    if errors:
        warnings.append(f"{len(errors)} conclusion evaluations raised and were recorded as undecidable")

    return {
        "family_size": len(fam),
        "fast": bool(fast),
        "n_shuffles": shuffles,
        "conc_M": float(conc_M),
        "seed": int(seed),
        "nicotinic_set": list(nic),
        "amplify_set": list(amp),
        "specs": [s["name"] for s in fam],
        "default_spec": DEFAULT_SPEC_NAME,
        "rows": rows,
        "matrix": matrix,
        "diagnostics": diagnostics,
        "errors": errors,
        "runtime_s": float(time.perf_counter() - t0),
        "runtime_by_spec_s": per_spec_runtime,
        "label": "model_derived",
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# threshold sensitivity for the amplify / buffer split
# --------------------------------------------------------------------------
MECHANISM_CLASSES: dict[str, tuple[str, ...]] = {
    "nicotinic": NICOTINIC_SET,
    "nav_ache": AMPLIFY_SET,
}


def threshold_sensitivity(
    circuit_fracs: Sequence[float] = (0.25, 0.40, 0.50),
    vert_limits: Sequence[float] = (0.10, 0.20, 0.30),
    compounds: Sequence[str] | None = None,
    assay: str = "subgraph",
    graph: str | None = None,
    concs: Sequence[float] | None = None,
    spec: Mapping[str, Any] | None = None,
    classes: Mapping[str, Sequence[str]] | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Recompute the amplify / buffer split on the threshold grid.

    The manuscript classifies a compound as *amplified* when its circuit
    selectivity index exceeds its receptor selectivity index and *buffered*
    when it does not.  Both thresholds that enter that comparison are
    conventions: the 50% relative circuit change that defines ``C_circuit`` and
    the 20% vertebrate engagement that defines ``C_vert``.  This recomputes the
    split on the full ``circuit_fracs x vert_limits`` grid and reports whether
    the *mechanism-level* conclusion ("nicotinic buffered, Nav/AChE amplified")
    survives every cell.

    The circuit ladder is computed once per (compound, circuit_frac) and reused
    across ``vert_limits``, which only shift ``C_vert`` by a closed-form factor.
    """
    from flylab.analysis.selectivity import circuit_threshold_conc

    t0 = time.perf_counter()
    cls = {k: tuple(v) for k, v in (classes or MECHANISM_CLASSES).items()}
    names = list(compounds) if compounds is not None else sorted({c for v in cls.values() for c in v})
    rec_si = {name: receptor_si(name) for name in names}

    # closed-form vertebrate thresholds, one per (compound, limit)
    vert: dict[tuple[str, float], float | None] = {}
    for name in names:
        for lim in vert_limits:
            vert[(name, float(lim))] = vertebrate_threshold(name, float(lim))

    # circuit thresholds, one ladder per (compound, frac)
    circ: dict[tuple[str, float], float | None] = {}
    with mechanism_spec(spec):
        for name in names:
            for frac in circuit_fracs:
                r = circuit_threshold_conc(
                    assay, name, readout="auto", threshold_frac=float(frac),
                    concs=concs, graph=graph, **kw,
                )
                circ[(name, float(frac))] = r["conc_M"]

    grid: list[dict[str, Any]] = []
    per_compound: list[dict[str, Any]] = []
    for frac in circuit_fracs:
        for lim in vert_limits:
            gaps: dict[str, list[float]] = {k: [] for k in cls}
            n_amp = n_buf = n_unscored = 0
            for name in names:
                c_circ = circ[(name, float(frac))]
                c_vert = vert[(name, float(lim))]
                si_rec = rec_si[name].get("si")
                gap = None
                if c_circ and c_vert and si_rec is not None and math.isfinite(si_rec):
                    gap = float(math.log10(c_vert) - math.log10(c_circ) - si_rec)
                if gap is None:
                    n_unscored += 1
                else:
                    n_amp += int(gap > 0)
                    n_buf += int(gap <= 0)
                    for k, members in cls.items():
                        if name in members:
                            gaps[k].append(gap)
                per_compound.append(
                    {
                        "compound": name,
                        "circuit_frac": float(frac),
                        "vert_limit": float(lim),
                        "circuit_threshold_M": c_circ,
                        "vertebrate_threshold_M": c_vert,
                        "receptor_si_log10": si_rec,
                        "si_gap_circuit_minus_receptor": gap,
                        "verdict": None if gap is None else ("amplify" if gap > 0 else "buffer"),
                    }
                )
            cell = {
                "circuit_frac": float(frac),
                "vert_limit": float(lim),
                "n_amplify": n_amp,
                "n_buffer": n_buf,
                "n_unscored": n_unscored,
            }
            for k in cls:
                vals = gaps[k]
                cell[f"mean_gap_{k}"] = float(sum(vals) / len(vals)) if vals else None
                cell[f"n_scored_{k}"] = len(vals)
                cell[f"verdict_{k}"] = (
                    None if not vals else ("amplify" if sum(vals) / len(vals) > 0 else "buffer")
                )
            grid.append(cell)

    expected = {"nicotinic": "buffer", "nav_ache": "amplify"}
    checks = []
    for k, want in expected.items():
        if k not in cls:
            continue
        got = [c[f"verdict_{k}"] for c in grid]
        checks.append(
            {
                "class": k,
                "expected": want,
                "n_cells": len(got),
                "n_matching": sum(1 for v in got if v == want),
                "n_undecided": sum(1 for v in got if v is None),
                "stable": all(v == want for v in got),
                "failing_cells": [
                    f"circuit_frac={c['circuit_frac']:.2f},vert_limit={c['vert_limit']:.2f}"
                    for c in grid
                    if c[f"verdict_{k}"] != want
                ],
            }
        )

    stable = all(c["stable"] for c in checks) if checks else None
    warnings = list(BASE_WARNINGS) + [
        "Changing vert_limit shifts every compound's C_vert by the same "
        "closed-form factor when the vertebrate Hill coefficients agree, so "
        "that axis moves the whole selectivity scale rather than re-ordering "
        "compounds; the circuit_frac axis is the one that can re-order them.",
        "A compound with no circuit threshold inside the tested ladder is "
        "reported as unscored, never as buffered.",
    ]
    if stable is False:
        warnings.append(
            "the mechanism-level amplify/buffer verdict is NOT stable across "
            "the threshold grid: "
            + "; ".join(f"{c['class']} fails at {', '.join(c['failing_cells'])}" for c in checks if not c["stable"])
        )

    return {
        "assay": assay,
        "graph": graph,
        "spec": (dict(spec or default_spec()))["name"],
        "circuit_fracs": [float(f) for f in circuit_fracs],
        "vert_limits": [float(v) for v in vert_limits],
        "compounds": names,
        "classes": {k: list(v) for k, v in cls.items()},
        "grid": grid,
        "per_compound": per_compound,
        "class_checks": checks,
        "stable": stable,
        "runtime_s": float(time.perf_counter() - t0),
        "label": "model_derived",
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------
def to_markdown(result: Mapping[str, Any]) -> str:
    """Markdown for a :func:`conclusion_stability` or
    :func:`threshold_sensitivity` result."""
    if "rows" in result and "matrix" in result:
        out = [
            f"### Conclusion Stability Matrix ({result['family_size']} prespecified specifications)",
            "",
            "| conclusion | claim | retained | % | failing specifications |",
            "|---|---|---|---|---|",
        ]
        for r in result["rows"]:
            fail = ", ".join(r["failing_specs"]) or "-"
            if r["undecidable_specs"]:
                fail += f" (undecidable: {', '.join(r['undecidable_specs'])})"
            out.append(
                f"| {r['conclusion']} | {r['claim']} | {r['n_retained']}/{r['n_specs']} "
                f"| {100 * r['fraction_retained']:.0f}% | {fail} |"
            )
        return "\n".join(out)
    out = [
        "### Threshold sensitivity of the amplify / buffer split",
        "",
        "| circuit frac | vertebrate limit | amplify | buffer | unscored | nicotinic | Nav/AChE |",
        "|---|---|---|---|---|---|---|",
    ]
    for c in result["grid"]:
        out.append(
            f"| {c['circuit_frac']:.0%} | {c['vert_limit']:.0%} | {c['n_amplify']} | "
            f"{c['n_buffer']} | {c['n_unscored']} | {c.get('verdict_nicotinic')} "
            f"| {c.get('verdict_nav_ache')} |"
        )
    return "\n".join(out)
