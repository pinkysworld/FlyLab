#!/usr/bin/env python3
"""Regenerate every number, figure and table in the FlyLab IJRC paper.

One command, committed artifacts only::

    python scripts/reproduce_paper.py            # full run  (~8-10 min)
    python scripts/reproduce_paper.py --fast     # reduced   (~2 min)
    python scripts/reproduce_paper.py --only dependence
    flylab reproduce-paper                       # same thing through the CLI

Nothing here downloads anything and nothing needs the 1.1 GB MaleCNS weight
matrix: every input is either ``data/derived/`` (two committed neighbourhood
cuts, the census fallback and the gustatory seed list), ``data/literature/``
(the six sourced YAML datasets) or ``flylab/pharm/library.yaml``.

Outputs (all under ``--outdir``, default ``papers/``):

===========================  ==================================================
``figures/F*.png`` / ``.svg``  300 dpi figures, colour-blind-safe, UI palette
``figures/captions.md``        one caption per figure
``tables/T*.csv`` / ``.md``    the paper and supplement tables
``results.json``               every scalar quoted in either document, keyed,
                               with a provenance block (flylab version, git
                               sha, library sha256, map ids, seeds)
``IJRC_FlyLab_draft.md``       rendered from ``IJRC_FlyLab_draft.md.in`` by
``SUPPLEMENT.md``              substituting ``{{key}}`` from ``results.json``
===========================  ==================================================

``--only`` is **incremental**: when ``<outdir>/results.json`` already exists the
run keeps every value, figure, table and note produced by the steps it is not
re-running, so a partial run cannot silently truncate the record. The document
then carries ``"partial": true`` -- except for ``--only paper``, which only
re-renders the two documents from whatever the last full run stored and so
leaves the record's status alone.

``--fast`` changes only the *statistical* effort, never the model:

=========================  ===============  ===============
knob                       default          ``--fast``
=========================  ===============  ===============
dependence permutations    1000 (named)     20
                           300 (taste map)  20
                           100 (landscape)  20
specification family       25 specs         9 (subsample)
Sobol' base sample         1024 (11264      64
                           evaluations)
Ishigami validation        16384            1024
selectivity ladder         13 half-decades  7 decades
dose-response bootstrap    200 resamples    40
dose-response replicates   4                2
prediction replicates      6                2
=========================  ===============  ===============

A ``--fast`` permutation p has a resolution of 1/21 and a ``--fast`` Sobol'
index sits close to its own noise floor; the point estimates that do not depend
on replicate or permutation count are identical.

Parallelism (``--jobs``, default: the core count capped at 4) changes the wall
clock only.  The two analyses that dominate a default run -- the Sobol' sample
and the conclusion-stability matrix -- depend on nothing the other steps
produce, so they are started in worker processes when the run begins and joined
at their own step; the permutation analyses use the remaining workers.  With
``--jobs 1`` everything runs in process and the same numbers come out, slower.

``OPENBLAS_NUM_THREADS=1`` (and its siblings, set at the top of this module
before numpy is imported) matters more than the worker count: the matrices are
small enough that a BLAS thread pool costs an order of magnitude more than it
saves.  Run this as a subprocess to get the pinning; importing it into a
process that has already imported numpy will not.

Determinism: every RNG in the pipeline is seeded from ``--seed`` (default 0)
and the rate model is deterministic, so two runs on the same commit produce
byte-identical ``results.json`` values apart from the ``runtime_s`` fields.
"""
from __future__ import annotations

import os

# BLAS threads: the matrices here are small (1126 and 1841 nodes) and the
# thread pool costs far more than it saves -- a 1126-node step goes from ~0.4 s
# to ~0.05 s single-threaded.  Must happen before numpy is imported, which is
# why it sits at the top of the module rather than inside main().
_THREAD_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)
for _var in _THREAD_VARS:
    os.environ.setdefault(_var, "1")
# If numpy was already imported by the host process these have no effect and
# the run is simply slower; run the script as a subprocess to get the pinning.

import argparse  # noqa: E402
import csv  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from dataclasses import dataclass, field  # noqa: E402
from datetime import datetime, timezone  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Callable, Iterable  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

RESULTS_SCHEMA = "flylab-paper-results/1"

# --------------------------------------------------------------------------
# palette -- mirrors flylab/static/styles.css so figures and UI agree
# --------------------------------------------------------------------------
INSECT = "#2a78d6"
VERTEBRATE = "#eb6834"
NT_COLOURS = {
    "acetylcholine": "#2a78d6",
    "gaba": "#eb6834",
    "glutamate": "#1baf7a",
    "octopamine": "#eda100",
    "dopamine": "#e87ba4",
    "serotonin": "#008300",
    "histamine": "#4a3aa7",
    "unclear": "#898781",
}
SEQ = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
GOOD = "#0ca30c"
WARN = "#fab219"
CRIT = "#d03b3b"
TEXT1 = "#0b0b0b"
TEXT2 = "#52514e"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"

#: tier -> hatch, so an unsourced (placeholder) bar is visible in greyscale
TIER_HATCH = {"literature_order": "", "measured_fit": "..", "class_placeholder": "///"}

#: the four-compound panel of F3 / T1 highlights (one per insect target family)
PANEL4 = ("imidacloprid", "fipronil", "deltamethrin", "ivermectin")
PAPER_CONC = 1e-6

#: The statistical effort of a **default** (non-``--fast``) run, in one place.
#: The referee's B4 was that the committed record had been produced at a
#: shuffle count the shipped code no longer used and that nothing detected it,
#: so ``tests/test_reproduce.py`` asserts that every one of these knobs, as
#: recorded in ``papers/results.json``, equals the value here -- and that the
#: ones that mirror a library default still equal that default.
DEFAULT_EFFORT: dict[str, int] = {
    "dep_n": 1000,
    "dep_n_taste": 300,
    "dep_n_landscape": 1000,
    "dep_n_landscape_taste": 300,
    "bal_n": 1000,
    "val_recovery_n": 200,
    "val_power_replicates": 6,
    "stab_n_shuffles": 100,
    "stab_family_size": 25,
    "abl_ref_n_draws": 50,
    "unc_n_base": 1024,
    "ish_n_base": 16384,
    "ic50_n_boot": 200,
    "pred_n_rep": 6,
}

#: knobs above that must also equal a constant shipped inside ``flylab``
LIBRARY_DEFAULTS: dict[str, tuple[str, str]] = {
    "dep_n": ("flylab.analysis.dependence", "PAPER_N"),
    "dep_n_landscape": ("flylab.analysis.dependence", "PAPER_N"),
    "stab_n_shuffles": ("flylab.analysis.robustness", "DEFAULT_SHUFFLES"),
}


#: null mode -> the key a dependence result files its probability under.
#: ``sign_permute_weight_matched`` is the ladder's rank-3 rung; the plain
#: ``sign_permute`` beside it is a joint target-set-and-sign null (T23).
_P_KEY: dict[str, str] = {
    "sign_permute": "p_sign",
    "sign_permute_weight_matched": "p_sign_wm",
    "weight_permute": "p_weight",
    "rewire_degree_preserving": "p_degree",
    "erdos_renyi": "p_ER",
}


# --------------------------------------------------------------------------
# run context
# --------------------------------------------------------------------------
@dataclass
class Ctx:
    outdir: Path
    fast: bool = False
    seed: int = 0
    quiet: bool = False
    partial: bool = False
    jobs: int = 1
    background: dict[str, Any] = field(default_factory=dict)
    values: dict[str, dict[str, Any]] = field(default_factory=dict)
    steps: dict[str, dict[str, Any]] = field(default_factory=dict)
    figures: list[dict[str, Any]] = field(default_factory=list)
    tables: list[dict[str, Any]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    _step: str = "?"

    # -- directories ------------------------------------------------------
    @property
    def figdir(self) -> Path:
        return self.outdir / "figures"

    @property
    def tabdir(self) -> Path:
        return self.outdir / "tables"

    def ensure_dirs(self) -> None:
        for d in (self.outdir, self.figdir, self.tabdir):
            d.mkdir(parents=True, exist_ok=True)

    # -- logging ----------------------------------------------------------
    def say(self, msg: str) -> None:
        if not self.quiet:
            print(msg, flush=True)

    def note(self, msg: str) -> None:
        """A finding worth printing and recording (e.g. a contradicted claim)."""
        self.notes.append(f"[{self._step}] {msg}")
        self.say(f"    ! {msg}")

    # -- recording --------------------------------------------------------
    def put(
        self,
        key: str,
        value: Any,
        text: str | None = None,
        unit: str | None = None,
        note: str | None = None,
    ) -> Any:
        """Record one paper-quotable scalar.

        ``text`` is the string the draft template substitutes for ``{{key}}``;
        when omitted it is derived from ``value`` with a sane default format.

        The key must be an identifier: the templates substitute
        ``{{name}}`` by regex, and a key with a space in it would be
        unquotable and would silently never render.
        """
        if not key.isidentifier():
            raise ValueError(f"results key {key!r} is not an identifier")
        self.values[key] = {
            "value": value,
            "text": text if text is not None else _auto_text(value),
            "unit": unit,
            "note": note,
            "step": self._step,
        }
        return value

    def get(self, key: str, default: Any = None) -> Any:
        row = self.values.get(key)
        return default if row is None else row["value"]

    def text(self, key: str, default: str = "n/a") -> str:
        row = self.values.get(key)
        return default if row is None else str(row["text"])


def _auto_text(value: Any) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return str(value) if abs(value) < 10_000 else f"{value:,}".replace(",", "\u202f")
    if isinstance(value, float):
        a = abs(value)
        if a != 0 and (a < 1e-3 or a >= 1e5):
            return f"{value:.2e}"
        if a >= 100:
            return f"{value:.0f}"
        if a >= 10:
            return f"{value:.1f}"
        return f"{value:.3f}"
    return str(value)


def _fmt_M(x: float | None) -> str:
    """Concentration as a readable molar string (3.5e-08 M -> 35 nM)."""
    if x is None:
        return "n/a"
    for scale, unit in ((1e-3, "mM"), (1e-6, "\u00b5M"), (1e-9, "nM"), (1e-12, "pM")):
        if abs(x) >= scale:
            v = x / scale
            return f"{v:.3g} {unit}"
    return f"{x:.2e} M"


# --------------------------------------------------------------------------
# matplotlib helpers
# --------------------------------------------------------------------------
def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 100,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.facecolor": SURFACE,
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "axes.edgecolor": TEXT2,
            "axes.labelcolor": TEXT1,
            "text.color": TEXT1,
            "xtick.color": TEXT2,
            "ytick.color": TEXT2,
            "grid.color": GRID,
            "legend.frameon": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    return plt


def save_fig(ctx: Ctx, fig, name: str, caption: str, svg: bool = False) -> None:
    """Write ``name.png`` (300 dpi) and optionally ``name.svg``; record both."""
    import matplotlib.pyplot as plt

    ctx.ensure_dirs()
    png = ctx.figdir / f"{name}.png"
    fig.savefig(png)
    out = [png]
    if svg:
        s = ctx.figdir / f"{name}.svg"
        fig.savefig(s)
        out.append(s)
    plt.close(fig)
    for p in out:
        ctx.figures.append(
            {
                "name": p.name,
                "path": str(p.relative_to(ctx.outdir)),
                "bytes": p.stat().st_size,
                "step": ctx._step,
                "figure": name,
                "caption": caption if p.suffix == ".png" else None,
            }
        )
    ctx.say(f"    fig {png.name} ({png.stat().st_size / 1024:.0f} kB)")


def save_table(
    ctx: Ctx,
    name: str,
    rows: list[dict[str, Any]],
    fields: Iterable[str],
    caption: str,
    md_fields: Iterable[str] | None = None,
) -> None:
    """Write ``name.csv`` (all fields) and ``name.md`` (a readable subset)."""
    ctx.ensure_dirs()
    fields = list(fields)
    csv_path = ctx.tabdir / f"{name}.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: _csv_cell(r.get(k)) for k in fields})

    mf = list(md_fields) if md_fields else fields
    lines = ["| " + " | ".join(mf) + " |", "|" + "|".join(["---"] * len(mf)) + "|"]
    for r in rows:
        lines.append("| " + " | ".join(_md_cell(r.get(k)) for k in mf) + " |")
    md_path = ctx.tabdir / f"{name}.md"
    md_path.write_text(f"**{name}.** {caption}\n\n" + "\n".join(lines) + "\n")

    for p in (csv_path, md_path):
        ctx.tables.append(
            {
                "name": p.name,
                "path": str(p.relative_to(ctx.outdir)),
                "bytes": p.stat().st_size,
                "rows": len(rows),
                "step": ctx._step,
            }
        )
    ctx.say(f"    tab {csv_path.name} ({len(rows)} rows)")


def _csv_cell(v: Any) -> Any:
    if isinstance(v, (list, tuple)):
        return "; ".join(str(x) for x in v)
    if isinstance(v, dict):
        return json.dumps(v, default=str)
    return "" if v is None else v


def _md_cell(v: Any) -> str:
    if v is None:
        return "--"
    if isinstance(v, float):
        return _auto_text(v)
    if isinstance(v, (list, tuple)):
        return "; ".join(str(x) for x in v)
    return str(v).replace("|", "/").replace("\n", " ")


# --------------------------------------------------------------------------
# background compute
# --------------------------------------------------------------------------
# Two analyses dominate the default run: the Sobol' sample (n_base x (k+2)
# model evaluations) and the conclusion-stability matrix (every conclusion
# re-derived under every specification).  Neither needs anything the other
# steps produce, so both are started in worker processes when the run begins
# and joined at their own step.  The functions below are module-level and take
# only primitives, because a fork-safe payload has to pickle.  If the pool is
# unavailable the step computes the same thing in process: the numbers are
# identical, only the wall clock differs.
def _bg_sobol(fast: bool, seed: int) -> dict[str, Any]:
    from flylab.analysis.uncertainty_global import sobol_analysis

    return sobol_analysis(
        compound="imidacloprid",
        conc_M=PAPER_CONC,
        readout="mean_hz",
        n_base=64 if fast else DEFAULT_EFFORT["unc_n_base"],
        n_boot=50 if fast else 200,
        seed=seed,
    )


def _bg_stability(fast: bool, seed: int) -> dict[str, Any]:
    from flylab.analysis.robustness import conclusion_stability, threshold_sensitivity

    return {
        "stability": conclusion_stability(fast=fast, seed=seed),
        "threshold": threshold_sensitivity(),
    }


def _background(ctx: Ctx, name: str, fn: Callable[..., dict[str, Any]], *args: Any) -> dict[str, Any]:
    """Result of a backgrounded job, or the same computation in process."""
    fut = ctx.background.pop(name, None)
    if fut is not None:
        try:
            t0 = time.perf_counter()
            out = fut.result()
            ctx.say(f"    joined background job {name} after {time.perf_counter() - t0:.1f}s")
            return out
        except Exception as exc:  # a dead worker must not lose the step
            ctx.note(f"background job {name} failed ({type(exc).__name__}: {exc}); recomputing in process")
    return fn(*args)


# ==========================================================================
# STEPS
# ==========================================================================
def step_provenance(ctx: Ctx) -> None:
    """Versions, hashes, map ids and seeds; nothing scientific."""
    from flylab.circuit.rate import load_graph, resolve_graph
    from flylab.maps.malecns import MAP_CITATION, MAP_ID
    from flylab.notebook.schema import NOTEBOOK_VERSION, provenance
    from flylab.pharm.occupancy import library_sha256, load_library

    p = provenance(ctx.seed)
    lib = load_library()
    named = load_graph(resolve_graph("named"))
    taste = load_graph(resolve_graph("taste_motor"))

    ctx.values["_provenance"] = {
        "value": {
            "flylab_version": p["flylab_version"],
            "git_sha": p["git_sha"],
            "library_sha256": library_sha256(),
            "library_version": lib.get("library_version"),
            "library_schema_version": lib.get("schema_version"),
            "notebook_schema_version": NOTEBOOK_VERSION,
            "platform": p["platform"],
            "map_id": MAP_ID,
            "map_citation": MAP_CITATION,
            "graphs": {
                "named": {
                    "n_nodes": named["n_nodes"],
                    "n_edges": named["n_edges"],
                    "hops": named.get("hops"),
                    "min_weight": named.get("min_weight"),
                    "closure_min_weight": named.get("closure_min_weight"),
                },
                "taste_motor": {
                    "n_nodes": taste["n_nodes"],
                    "n_edges": taste["n_edges"],
                    "hops": taste.get("hops"),
                    "min_weight": taste.get("min_weight"),
                    "closure_min_weight": taste.get("closure_min_weight"),
                },
            },
            "seeds": {"global": ctx.seed},
            "fast": ctx.fast,
        },
        "text": "",
        "unit": None,
        "note": "provenance block, not a paper scalar",
        "step": ctx._step,
    }
    ctx.put("flylab_version", p["flylab_version"])
    ctx.put("git_sha", p["git_sha"], text=str(p["git_sha"])[:12])
    ctx.put("library_sha256", library_sha256(), text=library_sha256()[:12])
    ctx.put("map_id", MAP_ID)
    ctx.put("named_nodes", int(named["n_nodes"]))
    ctx.put("named_edges", int(named["n_edges"]))
    ctx.put("taste_nodes", int(taste["n_nodes"]))
    ctx.put("taste_edges", int(taste["n_edges"]))
    ctx.put("named_min_weight", int(named.get("min_weight") or 0))
    ctx.put("taste_closure_min_weight", int(taste.get("closure_min_weight") or 0))
    ctx.say(
        f"    flylab {p['flylab_version']} @ {str(p['git_sha'])[:8]} "
        f"lib {library_sha256()[:8]} map {MAP_ID}"
    )


def step_census(ctx: Ctx) -> None:
    """MaleCNS census + the Gate-1 proof that labellar GRN types are present."""
    from flylab.maps.malecns import load_census_fallback

    c = load_census_fallback()
    nt = c["neurotransmitter_counts"]
    ctx.put("census_traced", int(c["n_traced"]))
    ctx.put("census_ach", int(nt.get("acetylcholine", 0)))
    ctx.put("census_glu", int(nt.get("glutamate", 0)))
    ctx.put("census_gaba", int(nt.get("gaba", 0)))
    ctx.put("census_unclear", int(nt.get("unclear", 0)))
    ctx.put("census_gustatory", int(c["gustatory"]["n_traced"]))
    ctx.put("census_labellar_bristle", int(c["gustatory"]["subclass_counts"].get("labellar bristle", 0)))

    tc = c["gustatory"]["type_counts"]
    lb_types = {k: v for k, v in tc.items() if re.fullmatch(r"LB\d[a-e]?", k)}
    sweet = {t: int(tc.get(t, 0)) for t in ("LB3b", "LB3c")}
    bitter = {t: int(tc.get(t, 0)) for t in ("LB1a", "LB1b", "LB1c", "LB1d")}
    ctx.put("lb_type_count", len(lb_types))
    ctx.put("lb_cell_count", int(sum(lb_types.values())))
    ctx.put("sweet_grn_cells", int(sum(sweet.values())))
    ctx.put("bitter_grn_cells", int(sum(bitter.values())))
    ctx.put(
        "sweet_grn_breakdown",
        sweet,
        text=", ".join(f"{k} {v}" for k, v in sweet.items()),
    )
    ctx.put(
        "bitter_grn_breakdown",
        bitter,
        text=", ".join(f"{k} {v}" for k, v in bitter.items()),
    )
    ctx.put("gate1_proven", all(v > 0 for v in {**sweet, **bitter}.values()))

    rows = [
        {"type": t, "cells": int(n), "modality_hypothesis": _modality(t)}
        for t, n in sorted(lb_types.items(), key=lambda kv: (-kv[1], kv[0]))
    ]
    save_table(
        ctx,
        "T0_labellar_types",
        rows,
        ["type", "cells", "modality_hypothesis"],
        "Labellar gustatory receptor-neuron types present in MaleCNS v1.0 `type` "
        "strings, with the sweet/bitter assignment used by FlyLab (a working "
        "hypothesis from the Cell 2026 nomenclature, not a MaleCNS annotation).",
    )
    ctx.say(
        f"    census {c['n_traced']} traced, {len(lb_types)} LB* types, "
        f"sweet {ctx.get('sweet_grn_cells')} / bitter {ctx.get('bitter_grn_cells')} cells"
    )


def _modality(t: str) -> str:
    if t in ("LB3b", "LB3c"):
        return "sweet (used)"
    if t in ("LB1a", "LB1b", "LB1c", "LB1d"):
        return "bitter (used)"
    return "unassigned"


def step_scorecard(ctx: Ctx) -> None:
    """F2 + T1: the dual scorecard, and the sourced library with tiers."""
    import numpy as np

    from flylab.pharm.occupancy import compare_compound, list_compounds, load_library

    plt = _plt()
    lib = load_library()
    pairs = list((lib.get("selectivity_pairs") or {}).keys())

    # --- numbers ------------------------------------------------------
    imi = compare_compound("imidacloprid", PAPER_CONC)
    fip = compare_compound("fipronil", PAPER_CONC)
    dia = compare_compound("diazepam", PAPER_CONC)
    s_imi = imi["selectivity"]["nAChR"]
    s_fip = fip["selectivity"]["GABA_A"]

    ctx.put("imi_insect_occ", float(s_imi["insect_occupancy"]))
    ctx.put("imi_vert_occ", float(s_imi["vertebrate_occupancy"]))
    ctx.put("imi_occ_gap", float(s_imi["occupancy_difference"]))
    ctx.put(
        "imi_ec50_ratio",
        float(s_imi["ec50_ratio_vert_over_insect"]),
        text=f"{s_imi['ec50_ratio_vert_over_insect']:.0f}",
    )
    ctx.put("imi_insect_ec50", float(s_imi["insect_ec50_M"]), text=_fmt_M(s_imi["insect_ec50_M"]))
    ctx.put("imi_vert_ec50", float(s_imi["vertebrate_ec50_M"]), text=_fmt_M(s_imi["vertebrate_ec50_M"]))
    ctx.put("fip_insect_occ", float(s_fip["insect_occupancy"]))
    ctx.put("fip_vert_occ", float(s_fip["vertebrate_occupancy"]))
    ctx.put("fip_ec50_ratio", float(s_fip["ec50_ratio_vert_over_insect"]), text=f"{s_fip['ec50_ratio_vert_over_insect']:.0f}")
    ctx.put("fip_insect_ec50", float(s_fip["insect_ec50_M"]), text=_fmt_M(s_fip["insect_ec50_M"]))
    ctx.put("fip_vert_ec50", float(s_fip["vertebrate_ec50_M"]), text=_fmt_M(s_fip["vertebrate_ec50_M"]))
    dia_rdl = next(r for r in dia["receptors"] if r["receptor"] == "insect_RDL")
    dia_gaba = next(r for r in dia["receptors"] if r["receptor"] == "vertebrate_GABA_A")
    # Schema v3: diazepam's insect RDL row carries no value, so the pipeline
    # reports N/A instead of the 1e-4 the peer review objected to.
    dia_rdl_eng = dia_rdl.get("engagement", dia_rdl.get("occupancy"))
    ctx.put(
        "dia_insect_rdl_occ",
        None if dia_rdl_eng is None else float(dia_rdl_eng),
        text="not modelled (N/A)" if dia_rdl_eng is None else None,
    )
    ctx.put("dia_vert_gabaa_occ", float(dia_gaba.get("engagement", dia_gaba.get("occupancy"))))
    ctx.put("dia_insect_rdl_tier", dia_rdl["evidence_tier"])

    # --- T1 -----------------------------------------------------------
    rows: list[dict[str, Any]] = []
    tier_counts: dict[str, int] = {}
    for key in list_compounds():
        entry = lib["compounds"][key]
        for receptor, spec in entry["receptors"].items():
            tier = spec.get("evidence_tier", "class_placeholder")
            tier_counts[tier] = tier_counts.get(tier, 0) + 1
            rows.append(
                {
                    "compound": key,
                    "name": entry["name"],
                    "class": entry.get("class"),
                    "cas": entry.get("cas"),
                    "receptor": receptor,
                    "organism": (lib["receptors"].get(receptor) or {}).get("organism"),
                    "param_type": spec.get("param_type"),
                    "value_M": spec.get("value_M"),
                    "ec50_M": spec.get("ec50_M"),  # deprecated alias; None = not modelled
                    "relation": spec.get("relation"),
                    "species": spec.get("species"),
                    "hill_n": spec.get("n", 1.0),
                    "direction": spec.get("direction"),
                    "efficacy": spec.get("efficacy"),
                    "evidence_tier": tier,
                    "source": spec.get("source"),
                }
            )
    ctx.put("library_compounds", len(list_compounds()))
    ctx.put("library_receptor_classes", len(lib["receptors"]))
    ctx.put("library_rows", len(rows))
    ctx.put("library_rows_sourced", int(tier_counts.get("literature_order", 0)))
    ctx.put("library_rows_placeholder", int(tier_counts.get("class_placeholder", 0)))
    ctx.put("library_selectivity_pairs", len(pairs))
    save_table(
        ctx,
        "T1_library",
        rows,
        [
            "compound",
            "name",
            "class",
            "cas",
            "receptor",
            "organism",
            "param_type",
            "value_M",
            "ec50_M",
            "relation",
            "species",
            "hill_n",
            "direction",
            "efficacy",
            "evidence_tier",
            "source",
        ],
        "The sourced compound library (schema v3). Every row carries an evidence "
        "tier, a named source, the KIND of parameter the source reported "
        "(`param_type`: Kd/Ki are binding constants, EC50/IC50/Kb are functional "
        "potencies) and how far that source sits from this compound/receptor/species "
        "(`relation`); `class_placeholder` rows carry no number at all and report "
        "N/A rather than a small response.",
        md_fields=["compound", "class", "receptor", "param_type", "value_M", "direction",
                   "evidence_tier"],
    )

    # --- F2 -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.5), sharey=True)
    for ax, occ, title in (
        (axes[0], imi, "Imidacloprid"),
        (axes[1], fip, "Fipronil"),
    ):
        x = np.arange(len(pairs))
        ins = [occ["selectivity"][p]["insect_occupancy"] for p in pairs]
        ver = [occ["selectivity"][p]["vertebrate_occupancy"] for p in pairs]
        ins_t = [occ["selectivity"][p]["evidence_tier"] for p in pairs]
        holder = [occ["selectivity"][p]["placeholder"] for p in pairs]
        for i in x:
            hatch = TIER_HATCH.get(ins_t[i], "")
            if holder[i]:
                ax.bar(
                    i, 1.0, 0.86, color="white", hatch="///", edgecolor=GRID,
                    linewidth=0.8, zorder=0,
                )
                ax.annotate(
                    "not modelled (no sourced value)",
                    xy=(i, 0.5),
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=TEXT2,
                    zorder=3,
                )
            # schema v3: a not-modelled side has no bar at all (None), never a
            # bar of height 0 that reads as "no effect measured here".
            if ins[i] is not None:
                ax.bar(i - 0.2, ins[i], 0.38, color=INSECT, hatch=hatch, edgecolor="white", linewidth=0.6, zorder=2)
            if ver[i] is not None:
                ax.bar(i + 0.2, ver[i], 0.38, color=VERTEBRATE, hatch=hatch, edgecolor="white", linewidth=0.6, zorder=2)
            if not holder[i]:
                for xo, val in ((i - 0.2, ins[i]), (i + 0.2, ver[i])):
                    ax.annotate(
                        f"{val:.2f}", xy=(xo, val), xytext=(0, 3), textcoords="offset points",
                        ha="center", fontsize=7, color=TEXT2,
                    )
        ax.set_xticks(x)
        ax.set_xticklabels(pairs, rotation=20, ha="right")
        ax.set_ylim(0, 1.08)
        ax.axhline(0.9, color=GRID, lw=0.8, zorder=0)
        ax.set_title(f"{title} at {_fmt_M(PAPER_CONC)}")
        ax.grid(axis="y", lw=0.5, zorder=0)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("receptor engagement (0-1); occupancy only for Kd/Ki rows)")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=INSECT, label="insect target"),
        plt.Rectangle((0, 0), 1, 1, color=VERTEBRATE, label="vertebrate counterpart"),
        plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=TEXT2, hatch="///", label="class placeholder (no sourced EC50)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12))
    fig.suptitle("Dual scorecard: same compound, same concentration, two organisms", y=1.02)
    save_fig(
        ctx,
        fig,
        "F2_dual_scorecard",
        f"Dual scorecard for imidacloprid and fipronil at {_fmt_M(PAPER_CONC)}. Bars are "
        "Hill occupancy of the insect target (blue) and its vertebrate counterpart "
        "(orange) for each of the five receptor pairs declared in "
        "`library.yaml:selectivity_pairs`. Hatched bars rest on a class placeholder: "
        "the compound has no sourced EC50 at that pair and the row is inert. "
        f"Imidacloprid reaches {ctx.get('imi_insect_occ'):.3f} at insect nAChR against "
        f"{ctx.get('imi_vert_occ'):.3f} at vertebrate a4b2; fipronil reaches "
        f"{ctx.get('fip_insect_occ'):.3f} at insect RDL against {ctx.get('fip_vert_occ'):.3f} "
        "at native heteromeric vertebrate GABA-A. Teaching-tier EC50s: these are "
        "literature-order values, not measured constants.",
    )


def step_curves(ctx: Ctx) -> None:
    """F3: occupancy curves and the insect/vertebrate selectivity window."""
    import numpy as np

    from flylab.pharm.occupancy import compare_compound, hill_occupancy, load_library

    plt = _plt()
    lib = load_library()
    concs = np.logspace(-11, -3, 161)

    fig, axes = plt.subplots(2, 2, figsize=(8.8, 6.0), sharex=True, sharey=True)
    panel_rows = []
    for ax, key in zip(axes.ravel(), PANEL4):
        occ = compare_compound(key, PAPER_CONC)
        best = max(
            (v for v in occ["selectivity"].values() if not v.get("placeholder")),
            key=lambda v: v["log10_ec50_ratio_vert_over_insect"],
            default=None,
        )
        if best is None:
            ax.set_title(f"{occ['compound']} (no sourced pair)")
            continue
        for name, colour, style in (
            (best["insect_receptor"], INSECT, "-"),
            (best["vertebrate_receptor"], VERTEBRATE, "--"),
        ):
            spec = lib["compounds"][key]["receptors"][name]
            y = [hill_occupancy(float(c), float(spec["ec50_M"]), float(spec.get("n", 1.0))) for c in concs]
            ax.plot(concs, y, style, color=colour, lw=2, label=f"{name} ({_fmt_M(spec['ec50_M'])})")
        lo, hi = best["insect_ec50_M"], best["vertebrate_ec50_M"]
        if hi > lo:
            ax.axvspan(lo, hi, color=SEQ[0], alpha=0.65, zorder=0)
            ax.annotate(
                f"selectivity window\n{best['log10_ec50_ratio_vert_over_insect']:.2f} log10 units",
                xy=(float(np.sqrt(lo * hi)), 0.55),
                ha="center",
                fontsize=8,
                color=TEXT2,
            )
        ax.axvline(PAPER_CONC, color=TEXT2, lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_ylim(0, 1.02)
        ax.set_title(f"{occ['compound']} ({occ['class']})")
        ax.legend(fontsize=7.5, loc="upper left")
        ax.grid(lw=0.5)
        ax.set_axisbelow(True)
        panel_rows.append(
            {
                "compound": key,
                "pair": f"{best['insect_receptor']} / {best['vertebrate_receptor']}",
                "insect_ec50_M": lo,
                "vertebrate_ec50_M": hi,
                "log10_window": best["log10_ec50_ratio_vert_over_insect"],
            }
        )
        ctx.put(
            f"window_log10_{key}",
            float(best["log10_ec50_ratio_vert_over_insect"]),
            text=f"{best['log10_ec50_ratio_vert_over_insect']:.2f}",
        )
    for ax in axes[-1]:
        ax.set_xlabel("free concentration (M)")
    for ax in axes[:, 0]:
        ax.set_ylabel("fractional occupancy")
    fig.suptitle("Occupancy curves and insect-over-vertebrate selectivity windows", y=0.98)
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F3_occupancy_curves",
        "Hill occupancy against free concentration for the four-compound panel, one "
        "insect target (solid, blue) and its vertebrate counterpart (dashed, orange) "
        "per compound. The shaded band is the selectivity window between the two "
        "EC50s, annotated in log10 units; the dotted vertical line is the "
        f"{_fmt_M(PAPER_CONC)} working concentration used elsewhere in the paper. "
        "Ivermectin's window is the GluCl/GlyR pair, deltamethrin's the Nav pair. "
        "Curves are the teaching library, not measured concentration-response data.",
    )
    ctx.put("panel4", list(PANEL4), text=", ".join(PANEL4))


def step_mechanisms(ctx: Ctx) -> None:
    """T2, T21: the mechanism -> gain rules and the transmitter sign table.

    T21 exists because five of the seven signed weights the rate engine
    iterates are asserted magnitudes that appeared in no manuscript table, are
    varied by no member of the specification family and are not bounded by the
    label-permuting null. Publishing them is the minimum; the table says of
    each row whether it is a convention or an assertion.
    """
    from flylab.circuit.rate import NORMALISATION_NOTES, sign_table_rows
    from flylab.pharm.mechanisms import GAIN_KEYS, mechanism_table_rows

    rows = mechanism_table_rows()
    ctx.put("mechanism_rules", len(rows))
    ctx.put("gain_keys", list(GAIN_KEYS), text=", ".join(GAIN_KEYS))

    signs = sign_table_rows()
    asserted = [r for r in signs if r["evidence"] == "asserted"]
    ctx.put("sign_rows", len(signs))
    ctx.put("sign_asserted", len(asserted))
    ctx.put(
        "sign_asserted_transmitters",
        [r["transmitter"] for r in asserted],
        text=", ".join(f"{r['transmitter']} {r['sign']:+g}" for r in asserted),
    )
    for r in signs:
        ctx.put(f"sign_{r['transmitter']}", float(r["sign"]), text=f"{r['sign']:+g}")
    ctx.put(
        "sign_in_specification_family",
        any(r["in_specification_family"] for r in signs),
    )
    ctx.put("norm_note_row_abs", NORMALISATION_NOTES["row_abs"])
    ctx.put("norm_note_none", NORMALISATION_NOTES["none"])
    ctx.put("norm_note_degree", NORMALISATION_NOTES["degree"])
    if asserted:
        ctx.note(
            "%d of the %d transmitter sign magnitudes the rate engine iterates (%s) "
            "are asserted rather than measured, are exempt from the specification "
            "family and are not bounded by the label-permuting null. They are now "
            "published in T21 rather than left in the source."
            % (
                len(asserted),
                len(signs),
                ", ".join(f"{r['transmitter']} {r['sign']:+g}" for r in asserted),
            )
        )
    save_table(
        ctx,
        "T21_transmitter_signs",
        signs,
        ["transmitter", "sign", "role", "evidence", "in_specification_family", "rationale"],
        "The signed weight every transmitter contributes in the rate engine "
        "(`flylab.circuit.rate.SIGN_TABLE`). Acetylcholine at +1 is the "
        "normalisation and GABA at -1 the convention it is scaled against; the other "
        "five magnitudes are **asserted**, carry no citation, are not varied by any "
        "member of the specification family, and are not bounded by the "
        "label-permuting null (which permutes assignments and leaves the table "
        "itself untouched). Glutamate in particular is signed inhibitory at 40 % of "
        "the stated strength, which is a hedge against the cells whose glutamatergic "
        "output is excitatory rather than a measurement.",
        md_fields=["transmitter", "sign", "role", "evidence", "in_specification_family"],
    )
    save_table(
        ctx,
        "T2_mechanisms",
        rows,
        ["receptor", "direction", "gain", "formula", "applies_to", "rationale"],
        "Mechanism to gain-patch rules. `th` is fractional occupancy at the named "
        "receptor; every gain is a dimensionless multiplier with 1.0 = vehicle and "
        "a floor of 0.05. This table is generated from "
        "`flylab.pharm.mechanisms.MECHANISM_TABLE`, the single source of truth every "
        "assay calls.",
        md_fields=["receptor", "direction", "gain", "formula"],
    )


def step_graph(ctx: Ctx) -> None:
    """F4 + T19: the MaleCNS cuts, their structure, GRN seeds and GRN -> MN9 paths.

    The structural census (T19) is here because the referee's B3 is a fact
    about the substrate, not about the analysis: a degree-preserving rewire of
    an in-star is close to the identity, so a negative topology verdict on such
    a cut is weak evidence. The paper has to be able to print these numbers.
    """
    import numpy as np

    from flylab.analysis.baselines import graph_census
    from flylab.analysis.dependence import cut_census
    from flylab.analysis.impact import grn_to_mn9_paths
    from flylab.analysis.layout import graph_for_viewer
    from flylab.analysis.selectivity import graph_composition
    from flylab.circuit.rate import load_graph, resolve_graph

    plt = _plt()
    g = load_graph(resolve_graph("taste_motor"))
    view = graph_for_viewer(g, min_weight=25, max_edges=1200)
    paths = grn_to_mn9_paths(g, max_len=3, top=80)

    ctx.put("paths_sweet", int(paths["n_sweet"]))
    ctx.put("paths_bitter", int(paths["n_bitter"]))
    ctx.put("paths_max_len", int(paths["max_len"]))
    for gname in ("named", "taste_motor"):
        comp = graph_composition(gname)
        ctx.put(f"{gname}_frac_ach_syn", float(comp["frac_ach_synapses"]))
        ctx.put(f"{gname}_frac_gaba_syn", float(comp["frac_gaba_synapses"]))
        ctx.put(f"{gname}_frac_glu_syn", float(comp["frac_glu_synapses"]))
        ctx.put(f"{gname}_ei_ratio_syn", float(comp["ei_ratio_synapses"]))
    ctx.put("taste_view_edges_drawn", len(view["edges"]))
    ctx.put("taste_view_min_weight", float(view["min_weight"]))

    # ---- T19: how much wiring each cut actually has ------------------
    census_rows: list[dict[str, Any]] = []
    for gname in ("named", "taste_motor"):
        cc = cut_census(gname)
        nts = graph_census(gname)
        total_nodes = sum(nts.values()) or 1
        ctx.put(f"{gname}_mean_degree", float(cc["mean_degree"]))
        ctx.put(f"{gname}_density", float(cc["density"]), text=f"{cc['density']:.2e}")
        ctx.put(f"{gname}_n_seed_nodes", int(cc["n_seed_nodes"]))
        ctx.put(f"{gname}_edges_onto_seeds", int(cc["edges_onto_seeds"]))
        ctx.put(
            f"{gname}_share_onto_seeds",
            float(cc["share_onto_seeds"]),
            text=f"{cc['share_onto_seeds']:.0%}",
        )
        ctx.put(f"{gname}_nodes_with_input", int(cc["n_nodes_with_in_degree"]))
        ctx.put(
            f"{gname}_share_recurrent_edges",
            float(cc["share_edges_from_nodes_with_input"]),
            text=f"{cc['share_edges_from_nodes_with_input']:.0%}",
        )
        ctx.put(f"{gname}_max_in_degree", int(cc["max_in_degree"]))
        ctx.put(f"{gname}_in_star", bool(cc["in_star"]))
        for nt in ("acetylcholine", "gaba", "glutamate", "unclear", "octopamine", "histamine"):
            ctx.put(f"{gname}_cells_{nt}", int(nts.get(nt, 0)))
        ctx.put(
            f"{gname}_share_unclear_cells",
            float(nts.get("unclear", 0) / total_nodes),
            text=f"{nts.get('unclear', 0) / total_nodes:.0%}",
        )
        census_rows.append(
            {
                "cut": gname,
                "n_nodes": cc["n_nodes"],
                "n_edges": cc["n_edges"],
                "mean_degree": round(float(cc["mean_degree"]), 3),
                "density": cc["density"],
                "n_seed_nodes": cc["n_seed_nodes"],
                "edges_onto_seeds": cc["edges_onto_seeds"],
                "share_onto_seeds": round(float(cc["share_onto_seeds"]), 4),
                "nodes_with_in_degree": cc["n_nodes_with_in_degree"],
                "share_edges_from_nodes_with_input": round(
                    float(cc["share_edges_from_nodes_with_input"]), 4
                ),
                "max_in_degree": cc["max_in_degree"],
                "in_star": cc["in_star"],
                **{f"cells_{k}": v for k, v in sorted(nts.items())},
            }
        )
    if ctx.get("named_in_star"):
        ctx.note(
            "the `named` cut is an in-star: %s of its %d edges terminate on its %d seed "
            "cells, only %d of %d nodes receive any input and only %s of edges leave a "
            "node that itself receives input (mean degree %.2f). A degree-preserving "
            "rewire of such a graph is close to the identity, so a negative topology "
            "verdict on this cut is weak evidence. The `taste_motor` cut (mean degree "
            "%.2f) is the substrate to generalise from."
            % (
                ctx.text("named_share_onto_seeds"),
                ctx.get("named_edges") or 0,
                ctx.get("named_n_seed_nodes") or 0,
                ctx.get("named_nodes_with_input") or 0,
                ctx.get("named_nodes") or 0,
                ctx.text("named_share_recurrent_edges"),
                ctx.get("named_mean_degree") or 0.0,
                ctx.get("taste_motor_mean_degree") or 0.0,
            )
        )
    if not ctx.get("named_cells_octopamine"):
        ctx.note(
            "the `named` cut contains no octopaminergic cell at all, so a compound "
            "whose only insect mechanism acts through g_oct (chlordimeform) cannot "
            "move this readout for substrate reasons rather than pharmacological ones."
        )
    save_table(
        ctx,
        "T19_cut_census",
        census_rows,
        [
            "cut",
            "n_nodes",
            "n_edges",
            "mean_degree",
            "density",
            "n_seed_nodes",
            "edges_onto_seeds",
            "share_onto_seeds",
            "nodes_with_in_degree",
            "share_edges_from_nodes_with_input",
            "max_in_degree",
            "in_star",
            "cells_acetylcholine",
            "cells_gaba",
            "cells_glutamate",
            "cells_histamine",
            "cells_unclear",
        ],
        "Structure of the two committed cuts. `share_onto_seeds` is the fraction of "
        "edges terminating on a seed cell and `share_edges_from_nodes_with_input` the "
        "fraction whose source itself receives input -- the recurrence budget, and the "
        "only edges that can carry a path longer than one hop. A cut where most edges "
        "point at a few hubs and most nodes have no input is an in-star, and a "
        "degree-preserving rewire of it destroys very little; the transmitter census "
        "is by cell, not by synapse, and a transmitter absent from a cut cannot be "
        "perturbed on it.",
        md_fields=[
            "cut",
            "n_nodes",
            "n_edges",
            "mean_degree",
            "share_onto_seeds",
            "nodes_with_in_degree",
            "share_edges_from_nodes_with_input",
            "in_star",
        ],
    )

    seeds = view["seeds"]
    sweet_ids = set(seeds.get("LB3b", []) + seeds.get("LB3c", []))
    bitter_ids = set(sum((seeds.get(t, []) for t in ("LB1a", "LB1b", "LB1c", "LB1d")), []))
    mn9_ids = set(seeds.get("MN9", []))
    dnp_ids = set(seeds.get("DNp01", []))

    fig = plt.figure(figsize=(11.6, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.22, 1.0], wspace=0.30)
    ax = fig.add_subplot(gs[0, 0])
    # The layout places each MaleCNS superclass in its own column; a small
    # deterministic in-column offset (from the bodyId, no RNG) spreads the
    # somata so the columns read as layers instead of single lines.
    pos = {
        int(n["id"]): (float(n["x"]) + ((int(n["id"]) % 89) / 89.0 - 0.5) * 210.0, float(n["y"]))
        for n in view["nodes"]
    }
    col_label: dict[float, dict[str, int]] = {}
    for n in view["nodes"]:
        col_label.setdefault(float(n["x"]), {}).setdefault(n.get("superclass") or "unknown", 0)
        col_label[float(n["x"])][n.get("superclass") or "unknown"] += 1
    for e in view["edges"]:
        a, b = pos.get(int(e["source"])), pos.get(int(e["target"]))
        if not a or not b:
            continue
        ax.plot(
            [a[0], b[0]],
            [a[1], b[1]],
            color=NT_COLOURS.get(e.get("nt") or "unclear", NT_COLOURS["unclear"]),
            lw=0.18 + 0.7 * min(1.0, float(e["weight"]) / 400.0),
            alpha=0.22,
            zorder=1,
        )
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    cols = [NT_COLOURS.get(n.get("nt") or "unclear", NT_COLOURS["unclear"]) for n in view["nodes"]]
    ax.scatter(xs, ys, s=3.0, c=cols, linewidths=0, alpha=0.6, zorder=2)
    for ids, colour, marker, label, size in (
        (sweet_ids, GOOD, "^", "sweet GRN seeds (LB3b/c)", 42),
        (bitter_ids, CRIT, "v", "bitter GRN seeds (LB1a-d)", 42),
        (dnp_ids, "#4a3aa7", "s", "DNp01", 70),
        (mn9_ids, TEXT1, "*", "MN9", 190),
    ):
        pts = [pos[i] for i in ids if i in pos]
        if pts:
            ax.scatter(
                [p[0] for p in pts],
                [p[1] for p in pts],
                s=size,
                marker=marker,
                c=colour,
                edgecolors="white",
                linewidths=0.6,
                zorder=4,
                label=label,
            )
    ytop = max(ys) if ys else 0.0
    ybot = min(ys) if ys else 0.0
    for cx, counts in col_label.items():
        name = max(counts, key=counts.get)
        ax.text(
            cx,
            ybot - 90,
            f"{name} ({sum(counts.values())})",
            ha="right",
            va="top",
            fontsize=6.4,
            color=TEXT2,
            rotation=24,
        )
    ax.set_ylim(ybot - 900, ytop + 620)
    ax.set_xlabel("")
    ax.set_ylabel("cells, transmitter-ordered within column")
    ax.set_title(
        f"MaleCNS taste-motor cut: {g['n_nodes']} cells, {g['n_edges']} edges\n"
        f"(drawn: edges >= {view['min_weight']:.0f} synapses, strongest {len(view['edges'])})",
        pad=10,
    )
    seed_legend = ax.legend(fontsize=7.0, loc="upper left", markerscale=0.9)
    ax.add_artist(seed_legend)
    ax.legend(
        handles=[
            plt.Line2D([], [], marker="o", ls="", ms=5, color=NT_COLOURS[k], label=k)
            for k in ("acetylcholine", "gaba", "glutamate", "unclear")
        ],
        fontsize=7.0,
        loc="upper left",
        bbox_to_anchor=(0.0, 0.84),
        title="predicted transmitter",
        title_fontsize=7.0,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    ax2 = fig.add_subplot(gs[0, 1])
    bars, labels, colours = [], [], []
    for modality, colour in (("sweet", GOOD), ("bitter", CRIT)):
        seen_types: set[tuple[str, ...]] = set()
        for p in _distinct_type_paths(paths[modality], 5, seen_types):
            bars.append(abs(float(p["product_before"])))
            labels.append(" -> ".join(str(t) for t in p["types"]))
            colours.append(colour)
    y = np.arange(len(bars))[::-1]
    ax2.barh(y, bars, color=colours, height=0.72)
    ax2.set_yticks(y)
    ax2.set_yticklabels(labels, fontsize=7)
    ax2.set_xscale("log")
    ax2.set_xlabel("|product of signed synapse weights| along the path")
    ax2.set_title(
        "Strongest 5 sweet and 5 bitter GRN -> MN9 routes (distinct type sequences)\n"
        f"{paths['n_sweet']} sweet and {paths['n_bitter']} bitter paths exist at <= "
        f"{paths['max_len']} synapses"
    )
    ax2.grid(axis="x", lw=0.5)
    ax2.set_axisbelow(True)
    save_fig(
        ctx,
        fig,
        "F4_taste_motor_graph",
        "Left: the committed MaleCNS v1.0 taste-motor cut "
        f"({g['n_nodes']} cells, {g['n_edges']} edges; 1 hop around MN9, DNp01 and the "
        "labellar GRN types, plus induced partner-partner edges of at least "
        f"{g.get('closure_min_weight')} synapses). Node and edge colour is the predicted "
        "consensus transmitter; layout is the deterministic superclass-column layout "
        "from `flylab.analysis.layout`. Seeds are marked: sweet GRNs (LB3b/c, up "
        "triangles), bitter GRNs (LB1a-d, down triangles), DNp01 (squares) and MN9 "
        "(stars). Right: the five strongest sweet and five strongest bitter paths of at "
        f"most {paths['max_len']} synapses from a labellar GRN to MN9, ranked by the "
        f"absolute product of signed synapse weights; there are {paths['n_sweet']} such "
        f"sweet and {paths['n_bitter']} such bitter paths in total. The sweet/bitter "
        "assignment of LB types is a working hypothesis, not a MaleCNS annotation.",
    )

    save_table(
        ctx,
        "T7_top_paths",
        [
            {
                "modality": p["modality"],
                "rank": i + 1,
                "types": " -> ".join(str(t) for t in p["types"]),
                "bodyIds": " -> ".join(str(b) for b in p["path"]),
                "transmitters": " -> ".join(str(t) for t in p["nts"]),
                "length": p["length"],
                "min_weight": p["min_weight"],
                "signed_product": p["product_before"],
            }
            for modality in ("sweet", "bitter")
            for i, p in enumerate(paths[modality][:6])
        ],
        ["modality", "rank", "types", "bodyIds", "transmitters", "length", "min_weight", "signed_product"],
        "Strongest labellar-GRN to MN9 paths in the taste-motor cut, by modality "
        "hypothesis. Signed product uses the predicted transmitter sign of each "
        "presynaptic cell.",
        md_fields=["modality", "rank", "types", "length", "min_weight", "signed_product"],
    )


def _distinct_type_paths(rows: list[dict[str, Any]], k: int, seen: set) -> list[dict[str, Any]]:
    """Strongest instance of each distinct *type* sequence (paths repeat by cell)."""
    out = []
    for p in rows:
        key = tuple(str(t) for t in p["types"])
        if key in seen:
            continue
        seen.add(key)
        out.append(p)
        if len(out) >= k:
            break
    return out


def step_veto(ctx: Ctx) -> None:
    """F5: the map-based bitter veto on the rate and LIF engines."""
    import numpy as np

    from flylab.assays.taste_map import run_taste_map_assay

    plt = _plt()
    arms = [("vehicle", None, 0.0), ("fipronil", "fipronil", PAPER_CONC), ("imidacloprid", "imidacloprid", PAPER_CONC)]
    data: dict[str, dict[str, dict[str, Any]]] = {}
    for engine in ("rate", "lif"):
        data[engine] = {}
        for label, compound, conc in arms:
            nb = run_taste_map_assay(compound, conc, engine=engine, seed=ctx.seed)
            r = nb["readouts"]
            data[engine][label] = r
            tag = f"{engine}_{label}"
            ctx.put(f"veto_{tag}_sugar", _f(r["mn9_sugar_hz"]), unit="Hz")
            ctx.put(f"veto_{tag}_sugar_bitter", _f(r["mn9_sugar_bitter_hz"]), unit="Hz")
            ctx.put(f"veto_{tag}_ratio", _f(r["bitter_veto_ratio"]))
    ctx.put("veto_n_sweet_grn", int(data["rate"]["vehicle"]["n_sweet_grn"]))
    ctx.put("veto_n_bitter_grn", int(data["rate"]["vehicle"]["n_bitter_grn"]))
    ctx.put("veto_sugar_hz", float(data["rate"]["vehicle"]["sugar_hz"]), unit="Hz")

    if data["rate"]["imidacloprid"]["bitter_veto_ratio"] is None:
        ctx.note(
            "imidacloprid at 1 uM silences MN9 on the map path on both engines "
            "(sugar-driven MN9 = 0 Hz), so its bitter-veto ratio is undefined, not 0. "
            "g_ach hits the 0.05 floor, which removes the cholinergic GRN->MN9 drive "
            "itself. This is a property of the agonist patch rule, and it is why the "
            "taste-map arm of the null panel is undefined for imidacloprid."
        )
    ctx.put("veto_imi_mn9_silenced", data["rate"]["imidacloprid"]["mn9_sugar_hz"] == 0.0)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
    for ax, engine in zip(axes, ("rate", "lif")):
        labels = ["vehicle", "fipronil", "imidacloprid"]
        x = np.arange(len(labels))
        sugar = [_f(data[engine][k]["mn9_sugar_hz"]) or 0.0 for k in labels]
        both = [_f(data[engine][k]["mn9_sugar_bitter_hz"]) or 0.0 for k in labels]
        ax.bar(x - 0.2, sugar, 0.38, color=SEQ[2], label="sugar only")
        ax.bar(x + 0.2, both, 0.38, color=SEQ[4], label="sugar + bitter")
        for i, k in enumerate(labels):
            ratio = data[engine][k]["bitter_veto_ratio"]
            ax.annotate(
                "veto ratio\nundefined" if ratio is None else f"veto ratio\n{ratio:.2f}",
                xy=(i, max(sugar[i], both[i])),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                fontsize=7.5,
                color=TEXT2,
            )
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("MN9 rate (Hz)")
        ax.set_title(f"{engine} engine")
        ax.set_ylim(0, max(sugar + both) * 1.3 + 0.1)
        ax.grid(axis="y", lw=0.5)
        ax.set_axisbelow(True)
        ax.legend(fontsize=8)
    fig.suptitle(
        "Map-extracted bitter veto of sugar-driven MN9 (MaleCNS taste-motor cut)", y=1.03
    )
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F5_bitter_veto",
        "MN9 firing on the map-extracted labellar-GRN path, driven through the real "
        f"MaleCNS sweet ({ctx.get('veto_n_sweet_grn')}) and bitter "
        f"({ctx.get('veto_n_bitter_grn')}) GRN body IDs at {ctx.get('veto_sugar_hz'):.0f} Hz, "
        "never by driving MN9 directly. Left: leaky rate engine. Right: Shiu-style LIF "
        "engine with a background Poisson drive. In vehicle, adding bitter drive takes "
        "MN9 to zero on both engines (veto ratio 0.00). Fipronil at "
        f"{_fmt_M(PAPER_CONC)} blocks RDL, raises the sugar-only rate and partially "
        f"relieves the veto (ratio {ctx.get('veto_rate_fipronil_ratio'):.2f} rate, "
        f"{ctx.get('veto_lif_fipronil_ratio'):.2f} LIF). Imidacloprid silences MN9 "
        "altogether, so its ratio is undefined rather than zero. Absolute rates differ "
        "between engines by construction and are not calibrated to any recording.",
    )


def _plural(n: int, singular: str, plural: str | None = None) -> str:
    """``"1 model inference"`` / ``"4 model inferences"``.

    The templates substitute a rendered string for each ``{{key}}``, so a count
    whose noun has to agree with it is recorded as the whole phrase rather than
    as a bare number with an ``s`` hard-coded in the prose.
    """
    return f"{n} {singular if n == 1 else (plural or singular + 's')}"


def _binomial_two_sided(k: int, n: int, p: float = 0.5) -> float | None:
    """Exact two-sided binomial tail for ``k`` successes in ``n`` trials.

    Used for the literature-concordance count, where a mean over comparisons
    half of which order only two compounds says nothing: counting concordant
    against discordant orderings and giving the tail does.
    """
    import math

    if n <= 0:
        return None
    probs = [math.comb(n, i) * p**i * (1.0 - p) ** (n - i) for i in range(n + 1)]
    target = probs[k] * (1.0 + 1e-9)
    return float(min(1.0, sum(q for q in probs if q <= target)))


def _f(x: Any) -> float | None:
    try:
        return None if x is None else float(x)
    except (TypeError, ValueError):
        return None


def step_evidence(ctx: Ctx) -> None:
    """F11 + T10, T20: the typed evidence census of the compound library.

    RQ1. Every library row records *what its source measured* (``param_type``)
    and *how far that source is from this compound/receptor/species*, graded
    E0-E4 (``evidence_distance``); the permitted transformation is a function
    of that pair, materialised in ``TRANSFORMATION_TABLE``. The census is the
    auditable summary of that typing, and T20 is the rule itself.
    """
    from flylab.pharm.evidence import (
        DISTANCE_LABELS,
        EngagementModel,
        EvidenceDistance,
        as_param_type,
        as_relation,
        model_for,
        transformation_table_rows,
    )
    from flylab.pharm.occupancy import (
        ENGAGEMENT_IS_NOT_OCCUPANCY,
        evidence_distance_table,
        library_report,
        load_library,
        spec_value_M,
        validate_library,
    )

    plt = _plt()
    rep = library_report()
    lib = load_library()
    dist = evidence_distance_table(lib)

    ctx.put("ev_schema_version", int(rep["schema_version"]))
    ctx.put("ev_rows", int(rep["n_rows"]))
    ctx.put("ev_compounds", int(rep["n_compounds"]))
    ctx.put("ev_receptor_keys", int(rep["n_receptor_keys"]))
    ctx.put("ev_rows_modelled", int(rep["n_rows_modelled"]))
    ctx.put("ev_rows_not_modelled", int(rep["n_rows_not_modelled"]))
    by_param = rep["by_param_type"]
    by_rel = rep["by_relation"]
    by_model = rep["by_engagement_model"]
    for name in ("EC50", "IC50", "Kd", "Ki", "Kb", "unknown", "class_order", "relative_potency"):
        ctx.put(f"ev_param_{name}", int(by_param.get(name, 0)))
    for name in (
        "exact_compound_exact_receptor_exact_species",
        "exact_compound_exact_receptor_other_species",
        "exact_compound_related_receptor",
        "class_extrapolation",
        "unsupported",
    ):
        ctx.put(f"ev_rel_{name}", int(by_rel.get(name, 0)))
    for model in EngagementModel:
        ctx.put(f"ev_model_{model.value}", int(by_model.get(model.value, 0)))
    ctx.put("ev_binding_rows", int(by_model.get(EngagementModel.binding_occupancy.value, 0)))
    ctx.put("ev_functional_rows", int(by_model.get(EngagementModel.functional_engagement.value, 0)))
    ctx.put("ev_not_modelled_rows", int(by_model.get(EngagementModel.not_modelled.value, 0)))
    ctx.put(
        "ev_proxy_rows",
        int(
            by_model.get(EngagementModel.functional_engagement_proxy.value, 0)
            + by_model.get(EngagementModel.binding_engagement_proxy.value, 0)
        ),
    )

    # the evidence-distance grading: the ordered form of `relation`, and what
    # the transformation rule is actually written against
    by_dist = rep["by_evidence_distance"]
    for grade in EvidenceDistance:
        ctx.put(f"ev_dist_{grade.value}", int(by_dist.get(grade.value, 0)))
    ctx.put(
        "ev_distance_labels",
        [DISTANCE_LABELS[g] for g in EvidenceDistance],
        text="; ".join(DISTANCE_LABELS[g] for g in EvidenceDistance),
    )
    ctx.put("ev_n_binding_occupancy", int(dist["n_binding_occupancy"]))

    # every row's permitted model, from the type system rather than re-derived
    binding: list[str] = []
    proxies: list[str] = []
    rows: list[dict[str, Any]] = []
    for compound, entry in (lib.get("compounds") or {}).items():
        for receptor, spec in (entry.get("receptors") or {}).items():
            param_type = as_param_type(spec.get("param_type"))
            relation = as_relation(spec.get("relation"))
            value = spec_value_M(spec)
            model = (
                model_for(param_type, relation)
                if value is not None
                else EngagementModel.not_modelled
            )
            named = f"{compound} at {receptor} ({param_type.value} {value:.2g} M)" if value else ""
            if model is EngagementModel.binding_occupancy:
                binding.append(named)
            elif model is EngagementModel.binding_engagement_proxy:
                proxies.append(named)
            rows.append(
                {
                    "compound": compound,
                    "receptor": receptor,
                    "param_type": param_type.value,
                    "value_M": value,
                    "hill_n": spec.get("n"),
                    "direction": spec.get("direction"),
                    "relation": relation.value,
                    "evidence_distance": spec.get("evidence_distance")
                    or {
                        "exact_compound_exact_receptor_exact_species": "E0",
                        "exact_compound_exact_receptor_other_species": "E1",
                        "exact_compound_related_receptor": "E2",
                        "class_extrapolation": "E3",
                        "unsupported": "E4",
                    }[relation.value],
                    "species": spec.get("species"),
                    "evidence_tier": spec.get("evidence_tier"),
                    "engagement_model": model.value,
                    "source": spec.get("source"),
                }
            )
    ctx.put("ev_binding_row_names", binding, text="; ".join(binding) or "none")
    ctx.put("ev_n_binding_named", len(binding))
    ctx.put("ev_binding_proxy_row_names", proxies, text="; ".join(proxies) or "none")
    ctx.put("ev_n_binding_proxy_named", len(proxies))

    # T20: the rule itself, as (parameter type x evidence distance) -> model
    save_table(
        ctx,
        "T20_transformation_rule",
        transformation_table_rows(),
        [
            "param_type",
            "evidence_distance",
            "evidence_distance_label",
            "relation",
            "engagement_model",
            "model_strength",
            "engagement_model_note",
            "param_type_note",
        ],
        "The transformation rule of the evidence type system, materialised. The "
        "permitted transformation is a function of two facts: the parameter the "
        "source measured and the ordered distance E0-E4 between that source and "
        "this compound at this receptor in this species. A binding constant yields "
        "physical fractional occupancy only at E0; transferred, it becomes a labelled "
        "proxy carrying a warning that names the gap. Anything at E4 is `not_modelled` "
        "and returns N/A. `tests/test_pharm_evidence.py` checks the soundness "
        "invariant exhaustively over every library row and every entry point: no "
        "number reaching a circuit readout originates from a row whose "
        "(param_type, distance) pair does not admit one.",
        md_fields=[
            "param_type",
            "evidence_distance_label",
            "engagement_model",
            "model_strength",
        ],
    )

    problems = validate_library()
    ctx.put("ev_library_problems", problems, text=", ".join(problems) or "none")
    if problems:
        ctx.note(f"library fails its own v3 validation in {len(problems)} place(s)")
    ctx.put("ev_note", ENGAGEMENT_IS_NOT_OCCUPANCY, text=ENGAGEMENT_IS_NOT_OCCUPANCY)

    save_table(
        ctx,
        "T10_evidence_types",
        rows,
        [
            "compound",
            "receptor",
            "param_type",
            "value_M",
            "hill_n",
            "direction",
            "relation",
            "evidence_distance",
            "species",
            "evidence_tier",
            "engagement_model",
            "source",
        ],
        "Typed pharmacological evidence (library schema v3). `param_type` is what the "
        "cited source measured; `evidence_distance` grades how far that source sits "
        "from this compound at this receptor in this species (E0 on-target to E4 "
        "unsupported); `engagement_model` is the only transformation that pair "
        "permits (T20). `binding_occupancy` requires a Kd or Ki measured at E0; a "
        "transferred binding constant becomes `binding_engagement_proxy`; "
        "`functional_engagement` is a normalised Hill response and not an occupancy; "
        "`not_modelled` rows return N/A and are excluded from every numeric result.",
        md_fields=[
            "compound",
            "receptor",
            "param_type",
            "value_M",
            "evidence_distance",
            "engagement_model",
        ],
    )

    # F11: what the library actually rests on
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))
    panels = (
        ("parameter type", by_param),
        (
            "evidence distance",
            {DISTANCE_LABELS[g]: by_dist.get(g.value, 0) for g in EvidenceDistance},
        ),
        ("engagement model", by_model),
    )
    for ax, (title, data) in zip(axes, panels):
        items = sorted(data.items(), key=lambda kv: -kv[1])
        names = [k for k, _ in items]
        vals = [v for _, v in items]
        colours = [
            CRIT if n in ("unknown", "unsupported", "not_modelled", "E4 (unsupported)") else
            (
                GOOD
                if n in ("Kd", "Ki", "binding_occupancy", "E0 (on-target)")
                else (WARN if str(n).endswith("_proxy") else SEQ[2])
            )
            for n in names
        ]
        y = list(range(len(names)))[::-1]
        ax.barh(y, vals, color=colours)
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=7.5)
        ax.set_xlabel("library rows")
        ax.set_title(title, fontsize=9.5)
        ax.grid(axis="x", lw=0.5)
        ax.set_axisbelow(True)
        for yy, v in zip(y, vals):
            ax.text(v, yy, f" {v}", va="center", fontsize=7, color=TEXT2)
    fig.suptitle(
        f"Typed evidence census: {rep['n_rows']} rows, {rep['n_compounds']} compounds, "
        f"{rep['n_receptor_keys']} receptor keys (schema v{rep['schema_version']})",
        y=1.04,
    )
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F11_evidence_types",
        "What the compound library rests on, row by row. Left: the parameter each "
        "cited source actually measured. Middle: the ordered evidence distance between "
        "that source and this compound at this receptor in this species. Right: the "
        "transformation the type system then permits, which is a function of the two. "
        f"Green marks the {ctx.get('ev_dist_E0')} on-target (E0) rows and the "
        f"{ctx.get('ev_n_binding_occupancy')} of them that may use a physical "
        "binding-occupancy model; amber marks the proxy models, where the same Hill "
        "expression is computed but reported as a transfer with a warning naming the "
        "gap; red marks rows that carry no modellable evidence and return N/A rather "
        "than a small number. The type system, not a convention, enforces this: asking "
        "for an engagement from an untyped row raises rather than returning a value, "
        "and a soundness invariant is checked exhaustively over every row and every "
        "entry point.",
    )


def step_dependence(ctx: Ctx) -> None:
    """F6, F12, F16 + T6, T11, T12, T12b, T22, T23: dependence analysis (RQ2).

    Four things changed after the second referee round and all of them move
    numbers rather than wording.

    *The ladder's rank-3 rung is now a weight-matched transmitter null.* Plain
    label permutation does not preserve the weighted excitation/inhibition
    balance, so it moves the drug's target set and the sign matrix at once;
    it is still run, and reported, as a joint null (T23 measures the
    difference rather than asserting it).

    *Non-rejection is no longer equivalence.* Every mode carries a three-way
    verdict against a prespecified margin, and the landscape is
    Benjamini-Hochberg corrected across its structural tests with both counts
    reported.

    *The landscape is repeated on a cut that has topology.* The ``named`` cut
    is an in-star (T19); ``taste_motor`` has mean degree 10.4, and the two do
    not agree.

    *The instrument is tested.* A known topology-dependent effect is planted
    in a synthetic cut of the same size and density and the ladder is asked to
    find it, with a power surface over planted strength and permutation count
    (T22, F16).
    """
    import numpy as np

    from flylab.analysis.dependence import (
        DEFAULT_MODES,
        INFORMATION_LADDER,
        MODE_INFORMATION,
        balance_report,
        dependence_landscape,
        dependence_profile,
        ladder_power,
        ladder_recovery,
    )
    from flylab.analysis.nullmodels import MODES, null_distribution

    plt = _plt()
    # leave a core for each job still running in the background
    jobs = max(1, ctx.jobs - len(ctx.background))
    n = 20 if ctx.fast else DEFAULT_EFFORT["dep_n"]
    n_taste = 20 if ctx.fast else DEFAULT_EFFORT["dep_n_taste"]
    n_land = 20 if ctx.fast else DEFAULT_EFFORT["dep_n_landscape"]
    n_land_taste = 20 if ctx.fast else DEFAULT_EFFORT["dep_n_landscape_taste"]
    #: fractions of n at which the permutation p is recomputed on the same
    #: draws -- a free permutation-count sweep (25, 50, 100, 200, 400, n).
    sweep = (0.025, 0.05, 0.1, 0.2, 0.4, 1.0)
    modes = tuple(DEFAULT_MODES)
    ctx.put("dep_n", n)
    ctx.put("dep_n_taste", n_taste)
    ctx.put("dep_n_landscape", n_land)
    ctx.put("dep_n_landscape_taste", n_land_taste)
    ctx.put("dep_p_floor", 1.0 / (n + 1), text=f"{1.0 / (n + 1):.4f}")
    ctx.put("dep_p_floor_landscape", 1.0 / (n_land + 1), text=f"{1.0 / (n_land + 1):.4f}")
    ctx.put("dep_modes", list(modes), text=", ".join(modes))
    ctx.put("dep_ladder", list(INFORMATION_LADDER), text=" < ".join(INFORMATION_LADDER))

    draws: dict[tuple[str, str], list[float]] = {}
    rows: list[dict[str, Any]] = []
    sweep_rows: list[dict[str, Any]] = []

    def arm(compound: str, assay: str, readout: str, n_i: int, tag: str) -> dict[str, Any]:
        """One profile, plus the raw draws the histogram and the sweep need.

        The profile is authoritative -- it carries the equivalence margin, the
        three-way verdicts, the class and the ladder. ``null_distribution`` is
        called again only for the four ``nullmodels`` degradations, to keep
        their draws for F6 and their checkpoint prefixes for T11; shuffle *i*
        depends only on ``(seed, i)``, so the two agree exactly.
        """
        prof = dependence_profile(
            compound,
            PAPER_CONC,
            assay=assay,
            readout=readout,
            modes=modes,
            n=n_i,
            seed=ctx.seed,
            n_jobs=jobs,
            checkpoints=sweep,
        )
        by_mode = {r["mode"]: r for r in prof["modes"]}
        for mode in modes:
            row = by_mode.get(mode)
            if row is None:
                continue
            info = MODE_INFORMATION[mode]
            rows.append(
                {
                    "compound": compound,
                    "assay": assay,
                    "readout": readout,
                    "mode": mode,
                    "information_kept": info["keeps"],
                    "rank": info["rank"],
                    "joint_null": bool(info.get("joint")),
                    "on_ladder": mode in INFORMATION_LADDER,
                    "n": row["n"],
                    "n_ok": row["n_ok"],
                    "real_effect": row["real_effect"],
                    "p_two_sided": row["p_two_sided"],
                    "p_resolution": row["p_resolution"],
                    "verdict": row.get("verdict"),
                    "abs_gap_from_null_median": row.get("abs_gap_from_null_median"),
                    "delta": row.get("delta"),
                    "within_tolerance": row.get("within_tolerance"),
                    "n_stabilised": row["n_stabilised"],
                    "z": row["z"],
                    "null_mean": row["null_mean"],
                    "null_median": row.get("null_median"),
                    "null_sd": row["null_sd"],
                }
            )
            ctx.put(f"dep_{tag}_{compound}_p_{mode}", _f(row["p_two_sided"]))
            ctx.put(f"dep_{tag}_{compound}_z_{mode}", _f(row["z"]))
            ctx.put(f"dep_{tag}_{compound}_verdict_{mode}", row.get("verdict"))
            ctx.put(f"dep_{tag}_{compound}_gap_{mode}", _f(row.get("abs_gap_from_null_median")))
        cls = prof["classification"]
        lvl = prof["necessary_information_level"]
        ctx.put(
            f"dep_{tag}_{compound}_effect",
            _f(prof["real_effect"]),
            unit="Hz" if assay == "subgraph" else None,
        )
        ctx.put(f"dep_{tag}_{compound}_vehicle", _f(prof.get("real_vehicle")))
        ctx.put(f"dep_{tag}_{compound}_class", cls["class"])
        ctx.put(f"dep_{tag}_{compound}_level", lvl.get("level"))
        ctx.put(f"dep_{tag}_{compound}_level_verdict", lvl.get("verdict"))
        ctx.put(f"dep_{tag}_{compound}_delta", _f(prof.get("delta")), unit="Hz")
        ctx.put(
            f"dep_{tag}_{compound}_delta_frac",
            _f(prof.get("delta_frac")),
            text=None if prof.get("delta_frac") is None else f"{prof['delta_frac']:.0%}",
        )
        ctx.put(f"dep_{tag}_{compound}_confirmatory", bool(prof.get("confirmatory")))
        ctx.put(f"dep_{tag}_{compound}_design", str(prof.get("design")))
        ctx.put(
            f"dep_{tag}_{compound}_equivalent_modes",
            cls["equivalent_modes"],
            text=", ".join(cls["equivalent_modes"]) or "none",
        )
        ctx.put(
            f"dep_{tag}_{compound}_indeterminate_modes",
            cls["indeterminate_modes"],
            text=", ".join(cls["indeterminate_modes"]) or "none",
        )
        zs = [abs(r["z"]) for r in prof["modes"] if r["z"] is not None]
        ctx.put(f"dep_{tag}_{compound}_max_abs_z", float(max(zs)) if zs else None)
        struct = [
            _f(by_mode[m]["p_two_sided"])
            for m in ("weight_permute", "rewire_degree_preserving")
            if m in by_mode
        ]
        ctx.put(
            f"dep_{tag}_{compound}_min_structural_p",
            min([p for p in struct if p is not None], default=None),
        )

        # draws and prefix sweep, for F6 and T11 only
        for mode in MODES:
            d = null_distribution(
                assay,
                compound,
                PAPER_CONC,
                mode,
                n=n_i,
                seed=ctx.seed,
                readout=readout,
                n_jobs=jobs,
                checkpoints=sweep,
            )
            for cp in (d["convergence"] or {}).get("checkpoints") or []:
                sweep_rows.append(
                    {
                        "compound": compound,
                        "assay": assay,
                        "mode": mode,
                        "n": cp["n"],
                        "p_two_sided": cp["p"],
                        "p_resolution": cp["p_resolution"],
                        "z": cp["z"],
                        "significant_at_0.05": cp["beats_null"],
                    }
                )
                ctx.put(f"sweep_{tag}_{compound}_{mode}_n{cp['n']}", _f(cp["p"]))
            if assay == "subgraph":
                draws[(compound, mode)] = [v for v in d["null_effects"] if v is not None]
        return prof

    named = {c: arm(c, "subgraph", "mean_hz", n, "named") for c in ("imidacloprid", "fipronil")}
    taste = arm("fipronil", "taste_map", "bitter_veto_ratio", n_taste, "taste")

    imi = named["imidacloprid"]
    imi_modes = {r["mode"]: r for r in imi["modes"]}
    gap_plain = (imi_modes.get("sign_permute") or {}).get("abs_gap_from_null_median")
    gap_matched = (imi_modes.get("sign_permute_weight_matched") or {}).get(
        "abs_gap_from_null_median"
    )
    if gap_plain is not None and gap_matched is not None:
        ctx.note(
            "the transmitter null matters to the equivalence claim and not to the "
            "probability: imidacloprid's gap from the null median is %.2f Hz under the "
            "plain label permutation and %.2f Hz under the weight-matched null, against "
            "a prespecified margin of %s Hz. The plain permutation is a joint "
            "target-set-and-sign null and is reported as one."
            % (gap_plain, gap_matched, ctx.text("dep_named_imidacloprid_delta"))
        )
    if taste["class"] != "topology-dependent":
        ctx.note(
            "the fipronil bitter-veto ratio on the taste-motor path is a negative at "
            "n=%d shuffles (max |z| = %s, p resolution %.4f): the veto ratio is model "
            "behaviour, not evidence about the wiring, and the verdicts say which "
            "modes are equivalent within tolerance (%s) and which are merely "
            "indeterminate (%s)."
            % (
                n_taste,
                ctx.text("dep_taste_fipronil_max_abs_z"),
                1.0 / (n_taste + 1),
                ctx.text("dep_taste_fipronil_equivalent_modes"),
                ctx.text("dep_taste_fipronil_indeterminate_modes"),
            )
        )

    save_table(
        ctx,
        "T6_dependence_profile",
        rows,
        [
            "compound",
            "assay",
            "readout",
            "mode",
            "information_kept",
            "rank",
            "on_ladder",
            "joint_null",
            "n",
            "n_ok",
            "real_effect",
            "p_two_sided",
            "p_resolution",
            "verdict",
            "abs_gap_from_null_median",
            "delta",
            "within_tolerance",
            "n_stabilised",
            "z",
            "null_mean",
            "null_median",
            "null_sd",
        ],
        "Connectome-dependence profile. `real_effect` is treated minus vehicle on the "
        "real MaleCNS cut; the null is the same contrast on `n` degraded copies of that "
        "cut, with the seed block and node order held fixed. `p_two_sided` is the "
        "empirical permutation probability (k+1)/(n+1) and is the statistic to read "
        "first; `p_resolution` is its floor, and `z` is a standardised distance from a "
        "usually non-normal null. `verdict` is the only claim the design supports: "
        "`distinguishable` when the test rejects, `equivalent_within_tolerance` when "
        "the gap from the null median is below the prespecified margin `delta`, and "
        "`indeterminate` otherwise -- a failure to reject is never evidence of "
        "equivalence. `sign_permute` is a **joint** target-set-and-sign null and is "
        "off the ladder (`on_ladder = False`); the rank-3 rung is the weight-matched "
        "transmitter null, which holds each transmitter's share of total outgoing "
        "weight fixed (T23).",
        md_fields=[
            "compound",
            "assay",
            "mode",
            "on_ladder",
            "real_effect",
            "p_two_sided",
            "verdict",
            "abs_gap_from_null_median",
            "delta",
        ],
    )

    save_table(
        ctx,
        "T11_permutation_sweep",
        sweep_rows,
        [
            "compound",
            "assay",
            "mode",
            "n",
            "p_two_sided",
            "p_resolution",
            "z",
            "significant_at_0.05",
        ],
        "Permutation-count sweep. The same draws are re-read as prefixes, so every row "
        "is a valid smaller permutation sample and no extra circuit runs were needed. "
        "It separates the count a *verdict* needs from the count a quotable *p* needs. "
        "`significant_at_0.05` is a rejection of the null at that prefix, not a "
        "statement that the shuffled graph reproduces the effect.",
        md_fields=["compound", "assay", "mode", "n", "p_two_sided", "significant_at_0.05"],
    )

    # where does the verdict settle, and where does the value settle?
    stab = [r["n_stabilised"] for r in rows if r["n_stabilised"]]
    ctx.put("dep_verdict_stable_from_n", int(min(stab)) if stab else None)
    ctx.put("dep_verdict_stable_to_n", int(max(stab)) if stab else None)
    mids = [
        r
        for r in sweep_rows
        if r["p_two_sided"] is not None and 0.05 < r["p_two_sided"] < 0.95 and r["n"] >= max(4, n // 40)
    ]
    settle: list[int] = []
    for key in {(r["compound"], r["assay"], r["mode"]) for r in mids}:
        series = sorted(
            [r for r in mids if (r["compound"], r["assay"], r["mode"]) == key], key=lambda r: r["n"]
        )
        if len(series) < 2:
            continue
        final = series[-1]["p_two_sided"]
        for r in series:
            if abs(r["p_two_sided"] - final) <= 0.02:
                settle.append(int(r["n"]))
                break
    ctx.put("dep_p_within_0p02_from_n", int(min(settle)) if settle else None)
    ctx.put("dep_p_within_0p02_worst_n", int(max(settle)) if settle else None)

    # ---- T23: what each transmitter null does to the weighted E/I balance --
    bal = balance_report(
        graph="named",
        n=20 if ctx.fast else DEFAULT_EFFORT["bal_n"],
        seed=ctx.seed,
    )
    bal_rows: list[dict[str, Any]] = []
    for mode, per_nt in bal["modes"].items():
        for nt, r in per_nt.items():
            bal_rows.append({"mode": mode, "transmitter": nt, **r})
    plain_ach = bal["modes"]["sign_permute"]["acetylcholine"]
    matched_ach = bal["modes"]["sign_permute_weight_matched"]["acetylcholine"]
    ctx.put("bal_n", int(bal["n"]))
    ctx.put("bal_real_ach_share", float(bal["real"]["acetylcholine"]))
    ctx.put("bal_plain_ach_mean", float(plain_ach["null_mean"]))
    ctx.put("bal_plain_ach_sd", float(plain_ach["null_sd"]))
    ctx.put("bal_plain_ach_min", float(plain_ach["null_min"]))
    ctx.put("bal_plain_ach_max", float(plain_ach["null_max"]))
    ctx.put("bal_plain_ach_percentile", float(plain_ach["percentile_of_real"]))
    ctx.put("bal_matched_max_deviation", float(matched_ach["max_abs_deviation"]))
    ctx.put("bal_tol", float(bal["tol"]))
    ctx.put("bal_plain_preserves", bool(bal["preserves_weighted_balance"]["sign_permute"]))
    ctx.put(
        "bal_matched_preserves",
        bool(bal["preserves_weighted_balance"]["sign_permute_weight_matched"]),
    )
    if not bal["preserves_weighted_balance"]["sign_permute"]:
        ctx.note(
            "plain transmitter-label permutation does NOT preserve the weighted "
            "excitation/inhibition balance: the cholinergic share of total outgoing "
            "synaptic weight is %.3f on the real `named` cut against %.3f +- %.3f over "
            "%d permutations, putting the real graph at the %.0fth percentile of its "
            "own null. The weight-matched null holds every tracked share within %.3f."
            % (
                bal["real"]["acetylcholine"],
                plain_ach["null_mean"],
                plain_ach["null_sd"],
                bal["n"],
                plain_ach["percentile_of_real"],
                float(matched_ach["max_abs_deviation"]),
            )
        )
    save_table(
        ctx,
        "T23_transmitter_null_balance",
        bal_rows,
        [
            "mode",
            "transmitter",
            "real",
            "null_mean",
            "null_sd",
            "null_min",
            "null_max",
            "max_abs_deviation",
            "percentile_of_real",
            "within_tol",
        ],
        "Share of total outgoing synaptic weight carried by each transmitter on the "
        "`named` cut, on the real graph and across the draws of each transmitter null. "
        "A gain patch acts *through* transmitter identity, so this share is what the "
        "drug sees. Plain label permutation moves it -- the real graph sits outside "
        "its own null on the cholinergic share -- which makes that mode a joint "
        "target-set-and-sign null rather than a null about transmitter identity. The "
        "weight-matched null is a constrained shuffle that holds every tracked share "
        "inside the tolerance band and is the ladder's rank-3 rung.",
        md_fields=[
            "mode",
            "transmitter",
            "real",
            "null_mean",
            "null_sd",
            "percentile_of_real",
            "within_tol",
        ],
    )

    # ---- T22 + F16: does the instrument work? ------------------------
    rec = ladder_recovery(
        strengths=(0.0, 0.5, 1.0, 2.0),
        n=20 if ctx.fast else DEFAULT_EFFORT["val_recovery_n"],
        seed=ctx.seed,
    )
    pw = ladder_power(
        strengths=(0.0, 0.25, 0.5, 1.0),
        ns=(20,) if ctx.fast else (50, 200, 1000),
        replicates=2 if ctx.fast else DEFAULT_EFFORT["val_power_replicates"],
        seed=ctx.seed,
    )
    ctx.put("val_recovery_n", int(rec["n"]))
    ctx.put("val_recovery_nodes", int(rec["n_nodes"]))
    ctx.put("val_recovery_edges", int(rec["n_edges"]))
    ctx.put("val_recovery_detected", int(rec["n_detected"]))
    ctx.put("val_recovery_positive", int(rec["n_positive"]))
    ctx.put("val_recovery_false_positive", bool(rec["false_positive_on_control"]))
    ctx.put("val_recovery_smallest", _f(rec["smallest_detected_strength"]))
    ctx.put("val_recovery_statement", str(rec["statement"]))
    ctx.put("val_power_replicates", int(pw["replicates"]))
    ctx.put("val_power_ns", list(pw["ns"]), text=", ".join(str(x) for x in pw["ns"]))
    ctx.put(
        "val_power_strengths",
        list(pw["strengths"]),
        text=", ".join(f"{x:g}" for x in pw["strengths"]),
    )
    ctx.put("val_power_fpr", _f(pw["false_positive_rate"]))
    controls = [g for g in pw["grid"] if g["is_control"]]
    ctx.put("val_power_control_hits", int(sum(g["n_detected"] for g in controls)))
    ctx.put("val_power_control_runs", int(sum(g["replicates"] for g in controls)))
    for g in pw["grid"]:
        ctx.put(
            f"val_power_s{str(g['loop_strength']).replace('.', 'p')}_n{g['n']}",
            _f(g["detection_rate"]),
        )
    by_strength: dict[float, list[float]] = {}
    for g in pw["grid"]:
        if g["is_control"]:
            continue
        by_strength.setdefault(float(g["loop_strength"]), []).append(
            float(g["detection_rate"] or 0.0)
        )
    full = sorted(s for s, rates in by_strength.items() if min(rates) >= 1.0)
    weak = sorted(s for s, rates in by_strength.items() if max(rates) <= 0.25)
    ctx.put("val_power_full_strengths", full, text=", ".join(f"{x:g}" for x in full) or "none")
    ctx.put("val_power_weak_strengths", weak, text=", ".join(f"{x:g}" for x in weak) or "none")
    ctx.put(
        "val_power_by_strength",
        {f"{s:g}": [round(min(r), 3), round(max(r), 3)] for s, r in sorted(by_strength.items())},
        text="; ".join(
            f"strength {s:g}: {min(r):.2f}-{max(r):.2f}" for s, r in sorted(by_strength.items())
        ),
    )
    best_n = {}
    for g in pw["grid"]:
        if g["is_control"]:
            continue
        best_n.setdefault(g["n"], []).append(g["detection_rate"] or 0.0)
    ctx.put(
        "val_power_by_n",
        {str(k): round(sum(v) / len(v), 3) for k, v in sorted(best_n.items())},
        text=", ".join(
            f"n={k}: {sum(v) / len(v):.2f}" for k, v in sorted(best_n.items())
        ),
    )
    ctx.note(
        "instrument validation: %s Detection over the planted grid by permutation "
        "count is %s and by planted strength %s; the empirical false-positive rate on "
        "the unplanted control is %d of %d. Effect size dominates the permutation "
        "count: full power at planted strength %s, detection at or below one in four "
        "at %s."
        % (
            rec["statement"],
            ctx.text("val_power_by_n"),
            ctx.text("val_power_by_strength"),
            int(sum(g["n_detected"] for g in controls)),
            int(sum(g["replicates"] for g in controls)),
            ctx.text("val_power_full_strengths"),
            ctx.text("val_power_weak_strengths"),
        )
    )
    save_table(
        ctx,
        "T22_instrument_validation",
        [
            {
                "experiment": "recovery",
                "loop_strength": r["loop_strength"],
                "n": r["n"],
                "replicates": 1,
                "class": r["class"],
                "necessary_level": r["necessary_information_level"],
                "real_effect": r["real_effect"],
                "detection_rate": 1.0 if r["topology_dependent"] else 0.0,
                "p_weight": (r["p"] or {}).get("p_weight"),
                "p_degree": (r["p"] or {}).get("p_degree"),
            }
            for r in rec["rows"]
        ]
        + [
            {
                "experiment": "power",
                "loop_strength": g["loop_strength"],
                "n": g["n"],
                "replicates": g["replicates"],
                "class": None,
                "necessary_level": None,
                "real_effect": None,
                "detection_rate": g["detection_rate"],
                "p_weight": None,
                "p_degree": None,
            }
            for g in pw["grid"]
        ],
        [
            "experiment",
            "loop_strength",
            "n",
            "replicates",
            "class",
            "necessary_level",
            "real_effect",
            "detection_rate",
            "p_weight",
            "p_degree",
        ],
        "Ground-truth recovery and power for the dependence ladder. A recurrent "
        "cholinergic cycle of known strength is planted in a synthetic cut of the same "
        "size, density and transmitter composition as the `named` cut, and a gain patch "
        "that collapses `g_ach` removes an amplification that exists only while the "
        "cycle is intact -- so the drug *contrast*, not merely the rate, depends on the "
        "wiring. `loop_strength = 0` plants nothing and is the negative control, whose "
        "detection rate is the empirical false-positive rate. Nothing here touches the "
        "connectome or the compound library: the experiment is about the instrument.",
        md_fields=["experiment", "loop_strength", "n", "replicates", "detection_rate", "class"],
    )

    # ---- the landscapes ----------------------------------------------
    def landscape(graph: str, n_i: int, tag: str) -> dict[str, Any]:
        land = dependence_landscape(graph=graph, n=n_i, seed=ctx.seed, n_jobs=jobs, modes=modes)
        s = land["summary"]
        counts = s["class_counts"]
        counts_raw = s["class_counts_raw"]
        counts_abs = s["class_counts_absolute_floor"]
        ctx.put(f"dep_{tag}_cells", int(land["n_cells"]))
        ctx.put(f"dep_{tag}_compounds", len(land["compounds"]))
        ctx.put(f"dep_{tag}_concs", len(land["concs_M"]))
        ctx.put(f"dep_{tag}_graph", str(land["graph"]))
        for label in (
            "topology-dependent",
            "mixed",
            "composition-dominated",
            "no-effect",
            "undefined",
        ):
            key = label.replace("-", "_")
            ctx.put(f"dep_{tag}_{key}", int(counts.get(label, 0)))
            ctx.put(f"dep_{tag}_{key}_raw", int(counts_raw.get(label, 0)))
            ctx.put(f"dep_{tag}_{key}_absfloor", int(counts_abs.get(label, 0)))
        ctx.put(f"dep_{tag}_design", str(s["design"]))
        ctx.put(f"dep_{tag}_statement", str(s["statement"]))
        ctx.put(f"dep_{tag}_fdr_alpha", _f(s["fdr_alpha"]))
        ctx.put(f"dep_{tag}_n_structural_tests", int(s["n_structural_tests"]))
        ctx.put(f"dep_{tag}_n_structural_rejected", int(s["n_structural_rejected_fdr"]))
        ctx.put(f"dep_{tag}_fdr_can_reject", bool(s["fdr_can_reject"]))
        ctx.put(f"dep_{tag}_non_monotone", int(s["n_non_monotone_ladders"]))
        ctx.put(f"dep_{tag}_below_relative_floor", int(s["n_below_relative_effect_floor"]))
        ctx.put(
            f"dep_{tag}_below_relative_floor_cells",
            [f"{r['compound']}@{r['conc_M']:.0e}" for r in s["below_relative_effect_floor"]],
            text=", ".join(
                f"{r['compound']} {_fmt_M(r['conc_M'])} ({r['real_effect']:.3g} Hz)"
                for r in s["below_relative_effect_floor"]
            )
            or "none",
        )
        ctx.put(
            f"dep_{tag}_non_monotone_cells",
            [f"{r['compound']}@{r['conc_M']:.0e}" for r in s["non_monotone_cells"]],
            text=", ".join(
                f"{r['compound']} {_fmt_M(r['conc_M'])}" for r in s["non_monotone_cells"]
            )
            or "none",
        )
        topo = list(s["topology_dependent_compounds"])
        ctx.put(f"dep_{tag}_topology_compounds", topo, text=", ".join(topo) or "none")
        ctx.put(
            f"dep_{tag}_topology_compounds_raw",
            list(s["topology_dependent_compounds_raw"]),
            text=", ".join(s["topology_dependent_compounds_raw"]) or "none",
        )
        ctx.put(
            f"dep_{tag}_topology_cells",
            [
                f"{c['compound']}@{c['conc_M']:.0e}"
                for c in land["cells"]
                if c["class"] == "topology-dependent"
            ],
            text=", ".join(
                f"{c['compound']} {_fmt_M(c['conc_M'])}"
                for c in land["cells"]
                if c["class"] == "topology-dependent"
            )
            or "none",
        )
        verdicts = s["verdict_counts"]
        for v, k in verdicts.items():
            ctx.put(f"dep_{tag}_verdicts_{v}", int(k))
        deltas = land.get("delta_range") or [None, None]
        ctx.put(f"dep_{tag}_delta_min", _f(deltas[0]), unit="Hz")
        ctx.put(f"dep_{tag}_delta_max", _f(deltas[1]), unit="Hz")
        # the cells the paper names one by one
        for cell in land["cells"]:
            if abs(cell["conc_M"] - PAPER_CONC) > 1e-18:
                continue
            if cell["compound"] not in ("imidacloprid", "fipronil"):
                continue
            c = cell["compound"]
            ctx.put(f"dep_{tag}_{c}_class", cell["class"])
            ctx.put(f"dep_{tag}_{c}_class_raw", cell["class_raw"])
            ctx.put(f"dep_{tag}_{c}_q", _f(cell.get("q_value")))
            ctx.put(f"dep_{tag}_{c}_effect", _f(cell.get("real_effect")))
            for mode in modes:
                ctx.put(
                    f"dep_{tag}_{c}_p_{mode}",
                    _f((cell.get("p") or {}).get(_P_KEY.get(mode, f"p_{mode}"))),
                )
        return land

    land = landscape("named", n_land, "land")
    land_taste = landscape("taste_motor", n_land_taste, "land_taste")
    cells = land["cells"]
    counts = land["summary"]["class_counts"]

    if land["summary"]["n_topology_dependent"] != land["summary"]["n_topology_dependent_raw"]:
        ctx.note(
            "multiplicity changes the headline count on the `named` cut: %d cells are "
            "topology-dependent on the raw permutation p and %d survive "
            "Benjamini-Hochberg across the %d structural tests of the landscape. The "
            "corrected count is the one the paper quotes."
            % (
                land["summary"]["n_topology_dependent_raw"],
                land["summary"]["n_topology_dependent"],
                land["summary"]["n_structural_tests"],
            )
        )
    if land_taste["summary"]["class_counts"].get("composition-dominated", 0) == 0:
        ctx.note(
            "the central RQ2 result REVERSES on a cut that has topology: repeated on "
            "`taste_motor` (%d nodes, %d edges, mean degree %.1f) the landscape returns "
            "%d composition-dominated cells of %d, against %d of %d on the in-star "
            "`named` cut. The composition-dominated majority is a property of the "
            "substrate, not of receptor perturbation on connectomes."
            % (
                ctx.get("taste_nodes") or 0,
                ctx.get("taste_edges") or 0,
                ctx.get("taste_motor_mean_degree") or 0.0,
                land_taste["summary"]["class_counts"].get("composition-dominated", 0),
                land_taste["n_cells"],
                counts.get("composition-dominated", 0),
                land["n_cells"],
            )
        )

    def landscape_rows(res: dict[str, Any]) -> list[dict[str, Any]]:
        out = []
        for c in res["cells"]:
            by_mode = {r["mode"]: r for r in c["modes"]}
            out.append(
                {
                    "compound": c["compound"],
                    "conc_M": c["conc_M"],
                    "real_effect": c.get("real_effect"),
                    "vehicle": c.get("real_vehicle"),
                    "effect_floor": c.get("effect_floor"),
                    "delta": c.get("delta"),
                    "class": c["class"],
                    "class_raw": c.get("class_raw"),
                    "class_absolute_floor": c.get("class_absolute_floor"),
                    "q_value": c.get("q_value"),
                    "necessary_level": (c.get("necessary_information_level") or {}).get("level"),
                    "necessary_level_verdict": (
                        c.get("necessary_information_level") or {}
                    ).get("verdict"),
                    "non_monotone": bool(
                        (c.get("necessary_information_level") or {}).get("non_monotone")
                    ),
                    "p_sign": (c.get("p") or {}).get("p_sign"),
                    "p_sign_wm": (c.get("p") or {}).get("p_sign_wm"),
                    "p_weight": (c.get("p") or {}).get("p_weight"),
                    "p_degree": (c.get("p") or {}).get("p_degree"),
                    "p_ER": (c.get("p") or {}).get("p_ER"),
                    "q_weight": (by_mode.get("weight_permute") or {}).get("q_value"),
                    "q_degree": (by_mode.get("rewire_degree_preserving") or {}).get("q_value"),
                    "verdict_weight": c["verdicts"].get("weight_permute"),
                    "verdict_degree": c["verdicts"].get("rewire_degree_preserving"),
                }
            )
        return out

    landscape_fields = [
        "compound",
        "conc_M",
        "real_effect",
        "vehicle",
        "effect_floor",
        "delta",
        "class",
        "class_raw",
        "class_absolute_floor",
        "q_value",
        "necessary_level",
        "necessary_level_verdict",
        "non_monotone",
        "p_sign",
        "p_sign_wm",
        "p_weight",
        "p_degree",
        "p_ER",
        "q_weight",
        "q_degree",
        "verdict_weight",
        "verdict_degree",
    ]
    landscape_md = [
        "compound",
        "conc_M",
        "real_effect",
        "class",
        "class_raw",
        "q_value",
        "necessary_level",
        "non_monotone",
    ]
    save_table(
        ctx,
        "T12_dependence_landscape",
        landscape_rows(land),
        landscape_fields,
        f"Connectome-dependence landscape on the `named` cut: every compound of the "
        f"library at {len(land['concs_M'])} concentrations, {n_land} shuffles per mode, "
        "the same shuffled graphs reused by every cell (a paired design). `class` is "
        "decided on the Benjamini-Hochberg adjusted probabilities across the "
        "structural tests of the whole landscape and `class_raw` on the uncorrected "
        "ones; `class_absolute_floor` drops the prespecified relative effect floor (1 % "
        "of the vehicle readout) so the cost of that convention is visible. "
        "`necessary_level` is the weakest graph model on the information ladder that "
        "this test could **not** distinguish from the real cut -- never a claim that "
        "the model reproduces the effect -- and `non_monotone` marks the cells where a "
        "richer model on the ladder *was* distinguishable, so the level must be read as "
        "'the cheapest graph model this test cannot tell apart from the real one' and "
        "not as 'everything above it is indistinguishable too'.",
        md_fields=landscape_md,
    )
    save_table(
        ctx,
        "T12b_dependence_landscape_taste_motor",
        landscape_rows(land_taste),
        landscape_fields,
        f"The same landscape repeated on the denser `taste_motor` cut "
        f"({ctx.get('taste_nodes')} nodes, {ctx.get('taste_edges')} edges, mean degree "
        f"{ctx.get('taste_motor_mean_degree') or 0.0:.1f}) at {n_land_taste} shuffles per "
        "mode. The `named` cut is an in-star (T19) and a degree-preserving rewire of it "
        "is close to the identity; this table is what the same instrument says about a "
        "cut that has wiring to destroy. The two do not agree, and that disagreement "
        "is the paper's main RQ2 result.",
        md_fields=landscape_md,
    )

    # ---- F6: the profile itself --------------------------------------
    fig, axes = plt.subplots(2, len(MODES), figsize=(12.0, 5.4))
    for row, compound in enumerate(("imidacloprid", "fipronil")):
        sub = {r["mode"]: r for r in rows if r["assay"] == "subgraph" and r["compound"] == compound}
        for col, mode in enumerate(MODES):
            ax = axes[row, col]
            vals = np.asarray(draws.get((compound, mode), []), dtype=float)
            real = sub[mode]["real_effect"]
            if vals.size:
                ax.hist(vals, bins=min(30, max(5, vals.size // 12)), color=SEQ[1], edgecolor="white", linewidth=0.4)
            if real is not None:
                ax.axvline(float(real), color=CRIT, lw=2, label="real MaleCNS cut")
            p = sub[mode]["p_two_sided"]
            verdict = str(sub[mode].get("verdict") or "")
            ax.set_title(
                f"{mode}\np = {'n/a' if p is None else f'{p:.3f}'}\n{verdict.replace('_', ' ')}",
                fontsize=7.6,
                color=CRIT if verdict == "distinguishable" else (
                    GOOD if verdict == "equivalent_within_tolerance" else TEXT2
                ),
            )
            ax.tick_params(labelsize=7)
            ax.grid(axis="y", lw=0.4)
            ax.set_axisbelow(True)
            if col == 0:
                ax.set_ylabel(f"{compound}\nshuffles")
            if row == 1:
                ax.set_xlabel("effect on mean rate (Hz)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle(
        f"Connectome-dependence profile, {n} permutations per mode "
        f"(network mean rate, {_fmt_M(PAPER_CONC)}); p floor 1/(n+1) = {1.0 / (n + 1):.4f}",
        y=1.02,
    )
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F6_dependence_profile",
        "Histograms are the drug effect (treated minus vehicle on the network mean "
        f"rate) measured on {n} degraded copies of the 1-hop MN9/DNp01 `named` cut; the "
        "red line is the same contrast on the real cut. Columns run in increasing order "
        "of destruction: transmitter labels permuted, synapse weights permuted, "
        "degree-preserving double-edge rewiring, Erdos-Renyi. Titles carry the "
        "empirical two-sided permutation p and the three-way verdict against the "
        f"prespecified equivalence margin (delta = "
        f"{ctx.text('dep_named_imidacloprid_delta')} Hz, 5 % of the vehicle readout): "
        "`distinguishable` means the test rejected, `equivalent within tolerance` means "
        "the gap from the null median is below the margin, and `indeterminate` means "
        "neither -- a failure to reject is not evidence that the degraded graph gives "
        "the same effect. Plain label permutation is shown because it is informative, "
        "but it moves the weighted excitation/inhibition balance as well as transmitter "
        "identity (T23), so it is a joint null and not a rung of the ladder. "
        f"Imidacloprid is {imi['class']} (necessary level: "
        f"{imi['necessary_information_level'].get('level')}, "
        f"{imi['necessary_information_level'].get('verdict')}); fipronil is "
        f"{named['fipronil']['class']} (necessary level: "
        f"{named['fipronil']['necessary_information_level'].get('level')}).",
    )

    # ---- F12: the two landscapes side by side -------------------------
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    class_idx = {
        "no-effect": 0,
        "composition-dominated": 1,
        "mixed": 2,
        "topology-dependent": 3,
        "undefined": 4,
    }
    cmap = ListedColormap([GRID, SEQ[1], WARN, CRIT, SURFACE])
    order = sorted({c["compound"] for c in cells})
    concs = sorted({c["conc_M"] for c in cells})
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 7.4), sharey=True)
    for ax, (res, title, n_i) in zip(
        axes,
        (
            (land, f"`named` cut (mean degree {ctx.get('named_mean_degree') or 0.0:.2f})", n_land),
            (
                land_taste,
                f"`taste_motor` cut (mean degree {ctx.get('taste_motor_mean_degree') or 0.0:.1f})",
                n_land_taste,
            ),
        ),
    ):
        grid = np.full((len(order), len(concs)), np.nan)
        for c in res["cells"]:
            if c["compound"] in order and c["conc_M"] in concs:
                grid[order.index(c["compound"]), concs.index(c["conc_M"])] = class_idx.get(
                    c["class"], np.nan
                )
        ax.imshow(grid, cmap=cmap, norm=BoundaryNorm(list(range(6)), cmap.N), aspect="auto")
        ax.set_xticks(range(len(concs)))
        ax.set_xticklabels([_fmt_M(c) for c in concs], fontsize=8)
        ax.set_xlabel("free concentration")
        ax.set_xticks([x - 0.5 for x in range(1, len(concs))], minor=True)
        ax.set_yticks([y - 0.5 for y in range(1, len(order))], minor=True)
        ax.grid(which="minor", color=SURFACE, linewidth=1.2)
        ax.tick_params(which="minor", length=0)
        ax.set_title(f"{title}\n{n_i} permutations per mode", fontsize=9)
    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels(order, fontsize=7.5)
    axes[1].legend(
        handles=[
            Patch(facecolor=GRID, label="no effect (below the relative floor)"),
            Patch(facecolor=SEQ[1], label="composition-dominated"),
            Patch(facecolor=WARN, label="mixed"),
            Patch(facecolor=CRIT, label="topology-dependent"),
        ],
        fontsize=7,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
    )
    fig.suptitle(
        "Connectome-dependence landscape on two cuts of the same connectome", y=0.99
    )
    fig.tight_layout()
    ct = land_taste["summary"]["class_counts"]
    save_fig(
        ctx,
        fig,
        "F12_dependence_landscape",
        "Every compound of the library at four concentrations, on both committed cuts, "
        "classified by the weakest graph model the permutation test could **not** "
        "distinguish from the real cut -- never by a model that reproduces the effect. "
        "Classes are decided on Benjamini-Hochberg adjusted probabilities across each "
        "landscape's structural tests, with a prespecified relative effect floor of 1 % "
        f"of the vehicle readout. Left, the 1-hop `named` cut: "
        f"{counts.get('topology-dependent', 0)} of {land['n_cells']} cells need the "
        f"wiring pattern, {counts.get('composition-dominated', 0)} are not "
        "distinguishable from any structure-preserving degradation, and "
        f"{counts.get('no-effect', 0)} do not move the readout at all. Right, the "
        f"denser `taste_motor` cut: {ct.get('topology-dependent', 0)} topology-dependent, "
        f"{ct.get('composition-dominated', 0)} composition-dominated, "
        f"{ct.get('no-effect', 0)} no-effect. The composition-dominated majority does "
        "not survive the change of substrate, and the `named` cut is an in-star whose "
        "degree-preserving rewire is close to the identity (T19).",
    )

    # ---- F16: the instrument's own power ------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), gridspec_kw={"width_ratios": [1.0, 1.1]})
    ax = axes[0]
    for i, n_i in enumerate(pw["ns"]):
        xs = [g["loop_strength"] for g in pw["grid"] if g["n"] == n_i]
        ys = [g["detection_rate"] for g in pw["grid"] if g["n"] == n_i]
        ax.plot(xs, ys, "o-", color=SEQ[min(i + 1, 4)], lw=1.8, ms=7, label=f"n = {n_i}")
    ax.axhline(0.05, color=CRIT, ls="--", lw=1)
    ax.text(0.02, 0.06, "alpha = 0.05", fontsize=7, color=CRIT, transform=ax.get_yaxis_transform())
    ax.set_xlabel("planted loop strength (0 = nothing planted)")
    ax.set_ylabel("share classified topology-dependent")
    ax.set_ylim(-0.05, 1.08)
    ax.set_title(
        f"Detection rate, {pw['replicates']} synthetic cuts per cell", fontsize=9.5
    )
    ax.legend(fontsize=7.5)
    ax.grid(lw=0.5)
    ax.set_axisbelow(True)

    ax2 = axes[1]
    labels = [f"{r['loop_strength']:g}" for r in rec["rows"]]
    effects = [abs(float(r["real_effect"] or 0.0)) for r in rec["rows"]]
    cols = [CRIT if r["topology_dependent"] else SEQ[1] for r in rec["rows"]]
    ax2.bar(range(len(labels)), effects, color=cols)
    ax2.set_xticks(range(len(labels)))
    ax2.set_xticklabels(labels)
    ax2.set_xlabel("planted loop strength")
    ax2.set_ylabel("|drug effect| on the synthetic cut (Hz)")
    ax2.set_title(
        f"Ground-truth recovery at n = {rec['n']}\n"
        f"{rec['n_detected']} of {rec['n_positive']} planted effects recovered, "
        + ("a false positive" if rec["false_positive_on_control"] else "no false positive")
        + " on the control",
        fontsize=9.5,
    )
    ax2.grid(axis="y", lw=0.5)
    ax2.set_axisbelow(True)
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F16_instrument_validation",
        "The dependence ladder pointed at a graph whose answer is known in advance. A "
        "recurrent cholinergic cycle of known strength is planted in a synthetic cut of "
        "the same size, density and transmitter composition as the `named` cut; a gain "
        "patch that collapses `g_ach` then removes an amplification that exists only "
        "while the cycle is intact, so the drug contrast depends on the wiring by "
        "construction. Left: the share of independently generated graphs the ladder "
        "classifies topology-dependent, against planted strength and permutation count. "
        "Right: the recovery experiment itself, with red marking the strengths the "
        "ladder recovered. The unplanted control is the empirical false-positive rate "
        f"({ctx.get('val_power_control_hits')} of "
        f"{ctx.get('val_power_control_runs')}). Without this panel, a negative from the "
        "ladder would be indistinguishable from an under-powered test.",
    )

def step_landscape(ctx: Ctx) -> None:
    """F7: receptor selectivity index versus circuit selectivity index."""
    import numpy as np

    from flylab.analysis.selectivity import DEFAULT_CONCS, selectivity_landscape

    plt = _plt()
    concs = tuple(10.0 ** (-10 + i) for i in range(7)) if ctx.fast else DEFAULT_CONCS
    ctx.put("landscape_n_concs", len(concs))
    land = selectivity_landscape(concs=concs)
    rows = land["rows"]
    scored = [r for r in rows if r.get("circuit_si_log10") is not None]
    ctx.put("landscape_rows", len(rows))
    ctx.put("landscape_scored", len(scored))
    ctx.put("landscape_compounds_scored", len({r["compound"] for r in scored}))

    scored = [r for r in scored if r.get("si_gap_circuit_minus_receptor") is not None]
    amplify = [r for r in scored if r["si_gap_circuit_minus_receptor"] > 0]
    buffer = [r for r in scored if r["si_gap_circuit_minus_receptor"] <= 0]
    a_gap = float(np.mean([r["si_gap_circuit_minus_receptor"] for r in amplify])) if amplify else float("nan")
    b_gap = float(np.mean([r["si_gap_circuit_minus_receptor"] for r in buffer])) if buffer else float("nan")
    ctx.put("landscape_amplify_rows", len(amplify))
    ctx.put("landscape_buffer_rows", len(buffer))
    ctx.put("landscape_amplify_compounds", sorted({r["compound"] for r in amplify}), text=", ".join(sorted({r["compound"] for r in amplify})))
    ctx.put("landscape_buffer_compounds", sorted({r["compound"] for r in buffer}), text=", ".join(sorted({r["compound"] for r in buffer})))
    ctx.put("landscape_amplify_mean_gap", a_gap)
    ctx.put("landscape_buffer_mean_gap", b_gap)
    ctx.put("landscape_gap_split", float(a_gap - b_gap))

    unscored = sorted({r["compound"] for r in rows if r.get("circuit_si_log10") is None and not r.get("skipped_reason")})
    ctx.put("landscape_unscored_compounds", unscored, text=", ".join(unscored))
    if unscored:
        ctx.note(
            "no circuit selectivity index exists for "
            + ", ".join(unscored)
            + ": on this neighbourhood their largest relative change in the mean rate "
            "never reaches the 50% threshold (RDL/GluCl blockers are bounded by the "
            "GABA/glutamate share of the cut and by the 0.05 gain floor)."
        )

    save_table(
        ctx,
        "T8_landscape",
        rows,
        [
            "compound",
            "class",
            "graph",
            "receptor_pair",
            "receptor_si_log10",
            "circuit_threshold_M",
            "vertebrate_threshold_M",
            "circuit_si_log10",
            "si_gap_circuit_minus_receptor",
            "max_rel_change",
            "frac_ach_synapses",
            "frac_gaba_synapses",
            "skipped_reason",
        ],
        "Selectivity landscape: the receptor index log10(EC50_vert / EC50_insect) "
        "against the circuit index log10(C_vert) - log10(C_circuit) on each committed "
        "graph. A negative gap means the neighbourhood buffers what the receptor "
        "numbers promise.",
        md_fields=["compound", "class", "graph", "receptor_si_log10", "circuit_si_log10", "si_gap_circuit_minus_receptor"],
    )

    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    lim = [-1.5, 4.0]
    ax.plot(lim, lim, color=TEXT2, lw=1, ls="--", zorder=1, label="circuit = receptor")
    ax.fill_between(lim, lim, [lim[1], lim[1]], color=GOOD, alpha=0.06, zorder=0)
    ax.fill_between(lim, [lim[0], lim[0]], lim, color=CRIT, alpha=0.06, zorder=0)
    seen = set()
    for r in scored:
        amp = r["si_gap_circuit_minus_receptor"] > 0
        ax.scatter(
            r["receptor_si_log10"],
            r["circuit_si_log10"],
            s=64 if r["graph"] == "named" else 40,
            marker="o" if r["graph"] == "named" else "D",
            c=GOOD if amp else CRIT,
            edgecolors="white",
            linewidths=0.7,
            zorder=3,
        )
        if r["graph"] == "named" and r["compound"] not in seen:
            seen.add(r["compound"])
            ax.annotate(
                r["compound"],
                (r["receptor_si_log10"], r["circuit_si_log10"]),
                textcoords="offset points",
                xytext=(6, -3),
                fontsize=7.5,
                color=TEXT2,
            )
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("receptor selectivity index  log10(EC50 vertebrate / EC50 insect)")
    ax.set_ylabel("circuit selectivity index  log10(C vertebrate 20% occ / C circuit 50% change)")
    ax.set_title(
        "Does the neighbourhood amplify or buffer receptor selectivity?\n"
        f"amplify (green, n={len(amplify)}): mean gap {a_gap:+.2f}   |   "
        f"buffer (red, n={len(buffer)}): mean gap {b_gap:+.2f} log10 units"
    )
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=TEXT2, label="1-hop MN9/DNp01 cut"),
        plt.Line2D([], [], marker="D", ls="", color=TEXT2, label="taste-motor cut"),
        plt.Line2D([], [], marker="o", ls="", color=GOOD, label="circuit amplifies"),
        plt.Line2D([], [], marker="o", ls="", color=CRIT, label="circuit buffers"),
    ]
    ax.legend(handles=handles, fontsize=8, loc="upper left")
    ax.grid(lw=0.5)
    ax.set_axisbelow(True)
    save_fig(
        ctx,
        fig,
        "F7_selectivity_landscape",
        "Each point is one compound on one committed graph. The x axis is the receptor "
        "selectivity index from the library alone; the y axis is the circuit index, the "
        "log10 width of the window between the concentration that moves the simulated "
        "network mean rate by 50% and the concentration at which the most potent "
        "vertebrate target reaches 20% occupancy. Points below the dashed identity line "
        "are buffered by the circuit, points above are amplified. The split is by "
        f"mechanism, not by potency: {ctx.text('landscape_amplify_compounds')} (Nav "
        f"modulators and the AChE inhibitor) amplify, mean gap {a_gap:+.2f} log10 units, "
        f"while the nAChR-acting set ({ctx.text('landscape_buffer_compounds')}) buffers, "
        f"mean gap {b_gap:+.2f}. Compounds with no circuit index "
        f"({ctx.text('landscape_unscored_compounds')}) never reach the 50% "
        "threshold on this cut and are omitted.",
    )


def step_dose(ctx: Ctx) -> None:
    """F8: model dose-response with bootstrap CI, plus the sensitivity tornado."""
    import numpy as np

    from flylab.assays.ensemble import circuit_ic50, sensitivity

    plt = _plt()
    n_boot = 40 if ctx.fast else DEFAULT_EFFORT["ic50_n_boot"]
    n_rep = 2 if ctx.fast else 4
    ctx.put("ic50_n_boot", n_boot)
    ctx.put("ic50_n_rep", n_rep)

    ic = circuit_ic50(
        "subgraph", "imidacloprid", readout="auto", n_boot=n_boot, n_rep=n_rep, seed=ctx.seed
    )
    fit, ci = ic["fit"], ic["ci"]
    ctx.put("ic50_M", float(fit["ic50"]), text=_fmt_M(fit["ic50"]))
    ctx.put("ic50_lo_M", float(ci["ic50"][0]), text=_fmt_M(ci["ic50"][0]))
    ctx.put("ic50_hi_M", float(ci["ic50"][1]), text=_fmt_M(ci["ic50"][1]))
    ctx.put("ic50_slope", float(fit["slope"]))
    ctx.put("ic50_r2", float(fit["r2"]), text=f"{fit['r2']:.4f}")
    ctx.put("ic50_top_hz", float(fit["top"]), unit="Hz")
    ctx.put("ic50_bottom_hz", float(fit["bottom"]), unit="Hz")
    ctx.put("ic50_in_range", bool(fit["in_range"]))
    ctx.put("ic50_readout", ic["readout"])
    ctx.put(
        "ic50_vs_receptor_ec50",
        float(fit["ic50"] / ctx.get("imi_insect_ec50", 2e-8)),
        text=f"{fit['ic50'] / ctx.get('imi_insect_ec50', 2e-8):.1f}",
    )

    tor = sensitivity("subgraph", "imidacloprid", PAPER_CONC, readout="mean_hz", seed=ctx.seed)
    tor = sorted(tor, key=lambda r: abs(float(r["span"])), reverse=True)
    for r in tor:
        ctx.put(f"tornado_span_{r['param']}", float(r["span"]), unit="Hz")
    ctx.put("tornado_top_param", tor[0]["param"])
    zero_span = [r["param"] for r in tor if abs(float(r["span"])) < 1e-9]
    ctx.put("tornado_zero_span_params", zero_span, text=", ".join(zero_span) or "none")
    if "ec50" in zero_span:
        ctx.note(
            "at 1 uM imidacloprid the model is exactly insensitive to a 2x change in "
            "the insect_nAChR EC50 and in the Hill coefficient (span 0 Hz): occupancy "
            "is saturated, so the gain-rule coefficient is the only parameter that "
            "still moves the readout. A dose-response experiment, not a single high "
            "dose, is what constrains the EC50."
        )

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), gridspec_kw={"width_ratios": [1.25, 1.0]})
    ax = axes[0]
    xs = np.asarray([p["conc_M"] for p in ic["points"]], dtype=float)
    mean = np.asarray([p["mean"] for p in ic["points"]], dtype=float)
    sd = np.asarray([p["sd"] for p in ic["points"]], dtype=float)
    ax.errorbar(xs, mean, yerr=sd, fmt="o", color=SEQ[3], capsize=3, ms=5, label=f"model replicates (n={n_rep})")
    grid = np.logspace(np.log10(xs.min()), np.log10(xs.max()), 200)
    yy = fit["bottom"] + (fit["top"] - fit["bottom"]) / (1.0 + (grid / fit["ic50"]) ** fit["slope"])
    ax.plot(grid, yy, color=CRIT, lw=2, label="4-parameter Hill fit")
    ax.axvspan(ci["ic50"][0], ci["ic50"][1], color=CRIT, alpha=0.14, label="bootstrap 95% CI on IC50")
    ax.axvline(fit["ic50"], color=CRIT, lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("imidacloprid concentration (M)")
    ax.set_ylabel(f"network {ic['readout']} (Hz)")
    ax.set_title(
        f"Model IC50 = {_fmt_M(fit['ic50'])}\n"
        f"[{_fmt_M(ci['ic50'][0])}, {_fmt_M(ci['ic50'][1])}], slope {fit['slope']:.2f}, r2 {fit['r2']:.4f}"
    )
    ax.legend(fontsize=7.5)
    ax.grid(lw=0.5)
    ax.set_axisbelow(True)

    ax2 = axes[1]
    names = [r["param"] for r in tor][::-1]
    base = float(tor[0]["base"])
    y = np.arange(len(names))
    for i, r in enumerate(tor[::-1]):
        lo, hi = float(r["low"]), float(r["high"])
        ax2.barh(i, lo - base, left=base, height=0.6, color=SEQ[1])
        ax2.barh(i, hi - base, left=base, height=0.6, color=SEQ[3])
    ax2.axvline(base, color=TEXT1, lw=1)
    ax2.set_yticks(y)
    ax2.set_yticklabels(names)
    ax2.set_xlabel(f"network mean rate (Hz); treated baseline at 1 uM = {base:.2f} Hz")
    ax2.set_title("One-at-a-time sensitivity\n(each parameter halved and doubled)")
    ax2.grid(axis="x", lw=0.5)
    ax2.set_axisbelow(True)
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F8_dose_response_sensitivity",
        "Left: the model's own dose-response for imidacloprid on the 1-hop MN9/DNp01 "
        f"cut. Points are {n_rep} replicates per concentration with log-normal EC50 "
        "jitter (sd 0.3 log10 units) and drive jitter; the curve is a 4-parameter Hill "
        f"fit, and the band is a {n_boot}-resample bootstrap CI on the fitted IC50 "
        f"({_fmt_M(fit['ic50'])}). This is the *model's* circuit IC50, not a receptor "
        "IC50 and not an animal IC50. Right: one-at-a-time sensitivity of the same "
        f"readout at {_fmt_M(PAPER_CONC)}, each parameter halved and doubled. The "
        "gain-rule coefficient dominates; the EC50 and Hill coefficient have zero span "
        "because occupancy is already saturated at this dose.",
    )


def step_genotype(ctx: Ctx) -> None:
    """F9: resistance-allele panels and their per-compound specificity."""
    import numpy as np

    from flylab.pharm.genotype import genotype_panel, list_genotypes

    plt = _plt()
    ctx.put("n_genotypes", len(list_genotypes()) - 1)  # minus wild type

    panels = {
        c: genotype_panel(c, PAPER_CONC, assay=None)
        for c in ("fipronil", "imidacloprid", "deltamethrin", "ddt", "gaba")
    }
    rows: list[dict[str, Any]] = []
    for compound, p in panels.items():
        for r in p["rows"]:
            rows.append(
                {
                    "compound": compound,
                    "genotype": r["genotype"],
                    "gene": r.get("gene"),
                    "allele": r.get("allele"),
                    "species": r.get("species"),
                    "receptor": r.get("receptor"),
                    "fold_shift": r.get("fold_shift"),
                    "shift_kind": r.get("shift_kind"),
                    "ec50_M": r.get("ec50_M"),
                    "occupancy": r.get("occupancy"),
                    "delta_occupancy_vs_wt": r.get("delta_occupancy_vs_wt"),
                    "shifted": r.get("shifted"),
                    "note": r.get("note"),
                }
            )
            key = f"geno_{compound}_{r['genotype']}"
            occ, fold = _f(r.get("occupancy")), _f(r.get("fold_shift"))
            ctx.put(f"{key}_occ", occ)
            ctx.put(f"{key}_fold", fold, text=("n/a" if fold is None else f"{fold:g}"))
    save_table(
        ctx,
        "T9_genotype",
        rows,
        [
            "compound",
            "genotype",
            "gene",
            "allele",
            "species",
            "receptor",
            "fold_shift",
            "shift_kind",
            "ec50_M",
            "occupancy",
            "delta_occupancy_vs_wt",
            "shifted",
            "note",
        ],
        "Resistance-allele panel. A shift is applied per compound, never per receptor: "
        "an allele with no sourced number for a compound leaves that compound "
        "unchanged and says so.",
        md_fields=["compound", "genotype", "receptor", "fold_shift", "shift_kind", "occupancy"],
    )

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
    for ax, (compound, title) in zip(
        axes,
        (
            ("fipronil", "Fipronil (insect RDL)"),
            ("imidacloprid", "Imidacloprid (insect nAChR)"),
            (None, "Specificity check: para alleles"),
        ),
    ):
        if compound is None:
            gts = ["wt", "para_L1014F_kdr", "para_M918T_superkdr"]
            labels = ["wild type", "L1014F kdr", "M918T super-kdr"]
            x = np.arange(len(gts))
            for off, (c, colour) in zip((-0.19, 0.19), (("deltamethrin", SEQ[3]), ("ddt", SEQ[1]))):
                by = {r["genotype"]: r for r in panels[c]["rows"]}
                vals = [_f((by.get(gt) or {}).get("occupancy")) or 0.0 for gt in gts]
                ax.bar(x + off, vals, 0.36, color=colour, label=c)
                for i, gt in enumerate(gts):
                    fold = _f((by.get(gt) or {}).get("fold_shift"))
                    if fold:
                        ax.annotate(
                            f"{fold:g}x", xy=(i + off, vals[i]), xytext=(0, 3),
                            textcoords="offset points", ha="center", fontsize=6.5, color=TEXT2,
                        )
            ax.set_xticks(x)
            ax.set_xticklabels(labels, fontsize=7, rotation=15, ha="right")
            ax.legend(fontsize=7.5)
        else:
            rws = panels[compound]["rows"]
            labels = [r["genotype"].replace("_dmel", "").replace("_dsim", "").replace("_nlug", "").replace("_mper", "") for r in rws]
            vals = [_f(r["occupancy"]) or 0.0 for r in rws]
            cols = [TEXT2 if r["genotype"] == "wt" else (INSECT if r.get("shifted") else GRID) for r in rws]
            x = np.arange(len(labels))
            ax.bar(x, vals, color=cols)
            for i, r in enumerate(rws):
                f = r.get("fold_shift")
                ax.annotate(
                    "no sourced\nshift" if not r.get("shifted") and r["genotype"] != "wt" else (f"{f:.0f}x" if f and f != 1.0 else ""),
                    xy=(i, vals[i]),
                    xytext=(0, 5),
                    textcoords="offset points",
                    ha="center",
                    fontsize=7,
                    color=TEXT2,
                )
            ax.set_xticks(x)
            ax.set_xticklabels(labels, fontsize=7, rotation=20, ha="right")
        ax.set_ylim(0, 1.15)
        ax.set_ylabel("target occupancy at 1 uM")
        ax.set_title(title, fontsize=9.5)
        ax.grid(axis="y", lw=0.5)
        ax.set_axisbelow(True)
    fig.suptitle("Resistance alleles shift compounds, not receptors", y=1.04)
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F9_genotype_panel",
        "Target-receptor occupancy at 1 uM for wild type and each published resistance "
        "allele that has a sourced fold-shift for that compound. Grey bars are alleles "
        "with no sourced number for the compound: they are left identical to wild type "
        "and reported, never extrapolated. Left and centre: fipronil and imidacloprid. "
        "Right: the specificity check that motivates the per-compound rule. *para* "
        f"M918T (super-kdr) shifts deltamethrin {ctx.text('geno_deltamethrin_para_M918T_superkdr_fold')}-fold "
        f"but leaves DDT at {ctx.text('geno_ddt_para_M918T_superkdr_fold')}-fold, and *Rdl* "
        f"A301S shifts GABA {ctx.text('geno_gaba_rdl_A301S_nlug_fold')}-fold while leaving "
        f"fipronil at {ctx.text('geno_fipronil_rdl_A301S_nlug_fold')}-fold. A toggle that "
        "shifted a whole receptor uniformly would be wrong in both cases.",
    )


_CITE_RE = re.compile(r"10\.\d{4,9}/[^\s,;)\]]+|PMID[:\s]*(\d+)", re.I)


def _cites(text: str | None) -> set[str]:
    """DOIs and PMIDs mentioned in a source string, normalised."""
    out: set[str] = set()
    for m in _CITE_RE.finditer(str(text or "")):
        out.add((m.group(1) or m.group(0)).lower().rstrip(". "))
    return out


def _source_overlap(validation: dict[str, Any]) -> set[str]:
    """Rank-order entries whose citation the library already uses.

    An entry counts as shared-source when one of the publications it comes from
    is also cited by a library row for one of the compounds that entry orders.
    Those comparisons cannot be out-of-sample validation, whatever their rho.
    """
    from flylab.pharm.occupancy import load_library

    lib = load_library()
    by_compound: dict[str, set[str]] = {}
    for compound, entry in (lib.get("compounds") or {}).items():
        ids: set[str] = set()
        for spec in (entry.get("receptors") or {}).values():
            ids |= _cites(spec.get("source"))
        by_compound[compound] = ids
    shared: set[str] = set()
    for res in validation["assays"].values():
        for r in res["rows"]:
            ids = _cites(r.get("source"))
            if not ids:
                continue
            names = [c.get("compound") for c in (r.get("compounds") or [])]
            for name in names:
                if ids & by_compound.get(str(name), set()):
                    shared.add(r["id"])
                    break
    return shared


def step_validation(ctx: Ctx) -> None:
    """F10 + T4: model ordering against published orderings."""
    import numpy as np

    from flylab.validation.rank import KNOWN_DISCREPANCIES, validate_all

    plt = _plt()
    v = validate_all(PAPER_CONC)
    s = v["summary"]
    ctx.put("rank_entries", int(s["n_entries"]))
    ctx.put("rank_evaluated", int(s["n_evaluated"]))
    ctx.put("rank_skipped", int(s["n_skipped"]))
    ctx.put("rank_mean_rho", float(s["mean_rho"]))

    # I13: half of the evaluable comparisons order two compounds, where rho can
    # only be +-1 and the two-sided p is 1.0. A mean over them is a mean over
    # coin flips, so the degenerate entries are counted, excluded from the mean
    # and replaced by a concordance count with its binomial tail.
    scored = [
        r
        for res in v["assays"].values()
        for r in res["rows"]
        if not r.get("skipped") and r.get("spearman_rho") is not None
    ]
    degenerate = [r for r in scored if int(r.get("n_compounds") or 0) <= 2]
    informative = [r for r in scored if int(r.get("n_compounds") or 0) > 2]
    ctx.put("rank_n_degenerate", len(degenerate))
    ctx.put("rank_n_informative", len(informative))
    ctx.put(
        "rank_mean_rho_informative",
        (sum(float(r["spearman_rho"]) for r in informative) / len(informative))
        if informative
        else None,
    )
    ctx.put(
        "rank_min_n_compounds",
        min((int(r.get("n_compounds") or 0) for r in scored), default=None),
    )
    ctx.put(
        "rank_max_n_compounds",
        max((int(r.get("n_compounds") or 0) for r in scored), default=None),
    )
    concordant = sum(1 for r in scored if float(r["spearman_rho"]) > 0)
    discordant = sum(1 for r in scored if float(r["spearman_rho"]) < 0)
    ctx.put("rank_concordant", concordant)
    ctx.put("rank_discordant", discordant)
    ctx.put("rank_binomial_p", _binomial_two_sided(concordant, concordant + discordant))
    if degenerate:
        ctx.note(
            "the literature-concordance mean is not quotable bare: %d of the %d "
            "evaluable comparisons order exactly two compounds, where Spearman's rho "
            "can only be +-1 and the two-sided p is 1.0 whatever the model does. Over "
            "the %d comparisons with more than two compounds the mean is %s; as a "
            "concordance count it is %d concordant against %d discordant orderings "
            "(two-sided binomial p = %s)."
            % (
                len(degenerate),
                len(scored),
                len(informative),
                ctx.text("rank_mean_rho_informative"),
                concordant,
                discordant,
                ctx.text("rank_binomial_p"),
            )
        )
    ctx.put("rank_mean_rho_occupancy", float(v["assays"]["occupancy"]["mean_rho"]))
    ctx.put("rank_mean_rho_subgraph", float(v["assays"]["subgraph"]["mean_rho"]))
    ctx.put("rank_exact_occupancy", int(v["assays"]["occupancy"]["n_exact"]))
    ctx.put("rank_known_discrepancies", len(KNOWN_DISCREPANCIES))
    ctx.put("rank_species_drosophila", sum(
        1
        for res in v["assays"].values()
        for r in res["rows"]
        if not r.get("skipped") and "drosophila" in str(r.get("species", "")).lower()
    ))
    ctx.put(
        "rank_discrepancy_ids",
        [d["id"] for d in KNOWN_DISCREPANCIES],
        text=", ".join(d["id"] for d in KNOWN_DISCREPANCIES),
    )

    # --- is a comparison independent of the evidence the library was built on?
    # A concordance test is only out-of-sample when the published ordering does
    # not come from a publication the library already cites for one of the
    # compounds being ordered. This is computed from the DOIs/PMIDs in both
    # files, not asserted.
    overlap_ids = _source_overlap(v)
    ctx.put("rank_shared_source_entries", sorted(overlap_ids), text=", ".join(sorted(overlap_ids)) or "none")
    ctx.put("rank_shared_source", len(overlap_ids))
    evaluated_ids = {
        r["id"] for res in v["assays"].values() for r in res["rows"] if not r.get("skipped")
    }
    ctx.put("rank_source_disjoint", len(evaluated_ids - overlap_ids))
    ctx.put("rank_entries_evaluable", len(evaluated_ids))
    n_shared = len(overlap_ids)
    ctx.put("rank_shared_source_verb", n_shared, text="comes" if n_shared == 1 else "come")
    ctx.put("rank_shared_source_be", n_shared, text="is" if n_shared == 1 else "are")
    n_disjoint = len(evaluated_ids - overlap_ids)
    ctx.put("rank_source_disjoint_be", n_disjoint, text="is" if n_disjoint == 1 else "are")

    rows = []
    for assay, res in v["assays"].items():
        for r in res["rows"]:
            rows.append(
                {
                    "assay": assay,
                    "id": r["id"],
                    "species": r.get("species"),
                    "published_assay": r.get("assay_published"),
                    "use_in_flylab": r.get("use_in_flylab"),
                    "independence": (
                        "shared-source (not out-of-sample)"
                        if r["id"] in overlap_ids
                        else "source-disjoint"
                    ),
                    "n_compounds": r.get("n_compounds"),
                    "skipped": r.get("skipped"),
                    "spearman_rho": r.get("spearman_rho"),
                    "kendall_tau": r.get("kendall_tau"),
                    "exact_match": r.get("exact_match"),
                    "receptor": r.get("receptor"),
                    "source": r.get("source"),
                }
            )
    save_table(
        ctx,
        "T4_rank_validation",
        rows,
        [
            "assay",
            "id",
            "species",
            "published_assay",
            "use_in_flylab",
            "independence",
            "n_compounds",
            "skipped",
            "spearman_rho",
            "kendall_tau",
            "exact_match",
            "receptor",
            "source",
        ],
        "Literature concordance of the model's ordering with published orderings in "
        "`data/literature/published_rank_orders.yaml`. Read `spearman_rho` beside its "
        "`n_compounds`: an entry ordering two compounds can only return +-1 and its "
        "two-sided p is 1.0 whatever the model does, so those entries are excluded "
        "from any mean and the aggregate the paper quotes is a count of concordant "
        "against discordant orderings with its binomial tail. `independence` is "
        "computed by "
        "comparing DOIs/PMIDs: `shared-source` means the published ordering comes from "
        "a publication the library already cites for one of the compounds it orders, so "
        "agreement there is internal consistency and not out-of-sample validation. "
        "`skipped` entries are those the library cannot cover; they are reported, not "
        "dropped.",
        md_fields=["assay", "id", "species", "independence", "n_compounds", "skipped", "spearman_rho", "exact_match"],
    )

    occ_rows = [r for r in v["assays"]["occupancy"]["rows"] if not r.get("skipped") and r.get("spearman_rho") is not None]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), gridspec_kw={"width_ratios": [1.1, 1.0]})

    ax = axes[0]
    markers = ["o", "s", "^", "D", "v", "P", "X"]
    n_e = max(1, len(occ_rows))
    for i, r in enumerate(occ_rows):
        off = (i - (n_e - 1) / 2.0) * 0.10
        pub = [float(c["published_rank"]) + off for c in r["compounds"]]
        mod = [m + off for m in _model_ranks(r["compounds"])]
        ax.plot(
            pub, mod, markers[i % len(markers)], ms=8, mfc="none", mew=1.6,
            color=SEQ[min(i, 4)], ls="",
            label=f"{r['id'][:30]} (rho {r['spearman_rho']:.2f})",
        )
    lim = [0.5, max([c["published_rank"] for r in occ_rows for c in r["compounds"]] + [5]) + 0.5]
    ax.plot(lim, lim, ls="--", color=TEXT2, lw=1, zorder=0)
    ax.set_xlabel("published rank (1 = most potent)")
    ax.set_ylabel("model rank (1 = most potent)")
    ax.set_title(f"Model vs published ordering\nmean Spearman rho = {s['mean_rho']:.2f} over {s['n_evaluated']} entries")
    ax.legend(fontsize=6.5, loc="lower right")
    ax.grid(lw=0.5)
    ax.set_axisbelow(True)

    ax2 = axes[1]
    disc = (v.get("known_discrepancies") or {}).get("discrepancies") or []
    seen_c: set[str] = set()
    names, ys, cols = [], [], []
    for d in disc:
        for row in d.get("model_rows") or []:
            if row["compound"] in seen_c:
                continue
            seen_c.add(row["compound"])
            names.append(str(row["compound"]))
            ys.append(float(row["ec50_M"]))
            cols.append(CRIT if row["compound"] in ("nitenpyram", "clothianidin") else SEQ[2])
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    names = [names[i] for i in order]
    ys = [ys[i] for i in order]
    cols = [cols[i] for i in order]
    x = np.arange(len(names))
    ax2.bar(x, ys, color=cols)
    ax2.set_yscale("log")
    ax2.set_xticks(x)
    ax2.set_xticklabels(names, fontsize=7.5, rotation=25, ha="right")
    ax2.set_ylabel("library insect_nAChR EC50 (M), lower = more potent")
    ax2.set_title(
        "Two recorded inversions (red)\n"
        "nitenpyram: least potent in vivo; clothianidin: more potent in vivo",
        fontsize=9.5,
    )
    ax2.grid(axis="y", lw=0.5)
    ax2.set_axisbelow(True)
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F10_rank_validation",
        "Left: the model's compound ordering against each published ordering that the "
        "teaching library can cover, with Spearman's rho per entry and the number of "
        f"compounds it orders. {s['n_evaluated']} entry-assay combinations are "
        f"evaluable and {s['n_skipped']} are skipped because the library cannot cover "
        f"them; {ctx.get('rank_n_degenerate')} of the evaluable ones order exactly two "
        "compounds, where rho can only be +-1, so no bare mean is quoted: over the "
        f"{ctx.get('rank_n_informative')} entries with more than two compounds the mean "
        f"is {ctx.text('rank_mean_rho_informative')}, and as a concordance count the "
        f"result is {ctx.get('rank_concordant')} concordant against "
        f"{ctx.get('rank_discordant')} discordant orderings (two-sided binomial p = "
        f"{ctx.text('rank_binomial_p')}). Right: the two inversions the validation module records "
        "rather than fixes. Nitenpyram sits near the top of the library's insect nAChR "
        "potency ordering and is the least potent neonicotinoid in *Drosophila* "
        "whole-animal bioassays, and clothianidin/imidacloprid are ordered the opposite "
        "way from the bioassays. Both are receptor-EC50 versus whole-animal-potency "
        "differences: the whole-animal number contains uptake, penetration and "
        "metabolism, which FlyLab does not model.",
    )


def _model_ranks(compounds: list[dict[str, Any]]) -> list[float]:
    scores = [c.get("model_score") for c in compounds]
    order = sorted(range(len(scores)), key=lambda i: -(scores[i] if scores[i] is not None else -1e18))
    ranks = [0.0] * len(scores)
    for pos, i in enumerate(order):
        ranks[i] = pos + 1.0
    return ranks


def step_predictions(ctx: Ctx) -> None:
    """T3: the prospective predictions with effect sizes and a planning minimum."""
    from flylab.analysis.predictions import prediction_table

    n_rep = 2 if ctx.fast else DEFAULT_EFFORT["pred_n_rep"]
    ctx.put("pred_n_rep", n_rep)
    tab = prediction_table(n_rep=n_rep, seed=ctx.seed)
    rows = []
    for r in tab["rows"]:
        rows.append(
            {
                "id": r["h_id"],
                "statement": r["statement"],
                "assay": r["assay"],
                "compound": r.get("compound"),
                "conc_M": r.get("conc_M"),
                "readout": r["readout"],
                "predicted_direction": r["predicted_direction"],
                "observed_direction": r.get("observed_direction"),
                "predicted_effect": r.get("predicted_effect"),
                "ci_low": (r.get("ci") or [None, None])[0],
                "ci_high": (r.get("ci") or [None, None])[1],
                "d_model_internal": r.get("d"),
                "endpoint": r.get("endpoint"),
                "suggested_n_per_group": r.get("suggested_n_per_group"),
                "live_protocol": r.get("live_protocol"),
                "status": r.get("status"),
                "live_result": r.get("live_result"),
            }
        )
        ctx.put(f"pred_{r['h_id']}_effect", _f(r.get("predicted_effect")))
        ctx.put(f"pred_{r['h_id']}_n", r.get("suggested_n_per_group"))
        ctx.put(f"pred_{r['h_id']}_dir_ok", r.get("predicted_direction") == r.get("observed_direction"))
    ctx.put("pred_n_hypotheses", len(rows))
    ctx.put("pred_baseline_per", _f((tab.get("baseline") or {}).get("p")))
    ns = [r["suggested_n_per_group"] for r in rows if r["suggested_n_per_group"]]
    ctx.put("pred_max_n", int(max(ns)) if ns else None)
    ctx.put("pred_live_results", sum(1 for r in rows if r["live_result"] is not None))
    mismatched = [r["id"] for r in rows if r["observed_direction"] and r["predicted_direction"] != r["observed_direction"] and r["predicted_direction"] != "none"]
    ctx.put("pred_direction_mismatches", mismatched, text=", ".join(mismatched) or "none")
    if mismatched:
        ctx.note("prospective direction not reproduced by the model for: " + ", ".join(mismatched))
    save_table(
        ctx,
        "T3_predictions",
        rows,
        [
            "id",
            "statement",
            "assay",
            "compound",
            "conc_M",
            "readout",
            "predicted_direction",
            "observed_direction",
            "predicted_effect",
            "ci_low",
            "ci_high",
            "d_model_internal",
            "endpoint",
            "suggested_n_per_group",
            "live_protocol",
            "status",
            "live_result",
        ],
        "Prospective (not pre-registered) predictions H1-H7. `predicted_effect` and its CI are "
        "model-internal (teaching-EC50 jitter, drive jitter, RNG seed) and carry no "
        "biological variance; `suggested_n_per_group` therefore caps the standardised "
        "effect at d = 1.0 before the power calculation. `live_result` is null for "
        "every row and stays null until a real table is imported by hand.",
        md_fields=["id", "assay", "compound", "readout", "predicted_direction", "predicted_effect", "suggested_n_per_group", "status"],
    )


def step_mixtures(ctx: Ctx) -> None:
    """T5: binary mixtures against the published honeybee verdicts."""
    from flylab.pharm.mixtures import validate_against_published

    v = validate_against_published()
    s = v["summary"]
    ctx.put("mix_n_targets", int(s["n_targets"]))
    ctx.put("mix_n_matching", int(s["n_matching"]))
    ctx.put("mix_all_match", bool(s["all_match"]))
    ctx.put("mix_n_substituted", int(s["n_substituted"]))
    ctx.put("mix_conc_M", float(v["conc_M"]), text=_fmt_M(v["conc_M"]))
    ctx.put("mix_species", v["species"])
    rows = []
    for r in v["results"]:
        rows.append(
            {
                "partner": r["partner"],
                "components": "; ".join(f"{c['compound']}@{c['conc_M']:.0e}" for c in r["components"]),
                "published": r["published"],
                "published_as_flylab_verdict": r["published_as_flylab_verdict"],
                "model_verdict": r["model_verdict"],
                "null_model": r["null_model"],
                "matches": r["matches"],
                "shared_receptors": r.get("shared_receptors"),
                "observed_effect_fraction": (r.get("bliss") or {}).get("observed_effect_fraction"),
                "expected_effect_fraction": (r.get("bliss") or {}).get("expected_effect_fraction"),
                "combination_index": (r.get("loewe") or {}).get("combination_index"),
                "substituted": json.dumps(r.get("substituted"), default=str) if r.get("substituted") else "",
            }
        )
    save_table(
        ctx,
        "T5_mixtures",
        rows,
        [
            "partner",
            "components",
            "published",
            "published_as_flylab_verdict",
            "model_verdict",
            "null_model",
            "matches",
            "shared_receptors",
            "observed_effect_fraction",
            "expected_effect_fraction",
            "combination_index",
            "substituted",
        ],
        "Binary-mixture verdicts against Zhu et al. 2017 (*Apis mellifera*). The null "
        "model is chosen per pair: Loewe additivity where the components share a "
        "FlyLab receptor key, Bliss independence where they do not. Two partners are "
        "not in the library and are represented by a same-class substitute, which is "
        "recorded per row.",
        md_fields=["partner", "published", "model_verdict", "null_model", "matches", "shared_receptors"],
    )


def step_expression(ctx: Ctx) -> None:
    """Expression-weighted sensitivity layer: coverage and the motor-neuron gap."""
    from flylab.pharm.expression import coverage, expression_table, run_weighted_subgraph_assay

    cov = coverage()
    ctx.put("expr_coverage_overall", float(cov["fraction_known_overall"]))
    for receptor, d in cov["graphs"]["named"].items():
        if isinstance(d, dict) and "fraction_known" in d:
            ctx.put(f"expr_cov_named_{receptor}", float(d["fraction_known"]))
    tab = expression_table()
    ctx.put("expr_rows", len(tab["rows"]))
    gap = tab.get("motor_neuron_gap") or {}
    ctx.put("expr_motor_gap", bool(gap))
    ctx.put(
        "expr_motor_gap_note",
        (gap.get("statement") or gap.get("note") or "motor-neuron receptor expression is a confirmed gap"),
    )

    nb = run_weighted_subgraph_assay("imidacloprid", PAPER_CONC)
    r = nb["readouts"]
    ctx.put("expr_mn9_uniform_hz", float(r["uniform"]["mn9_hz"]), unit="Hz")
    ctx.put("expr_mn9_weighted_hz", float(r["mn9_hz"]), unit="Hz")
    ctx.put("expr_mn9_delta_hz", float(r["delta_vs_uniform"]["mn9_hz"]), unit="Hz")
    ctx.put("expr_mean_uniform_hz", float(r["uniform"]["mean_hz"]), unit="Hz")
    ctx.put("expr_mean_weighted_hz", float(r["mean_hz"]), unit="Hz")


def step_ablation(ctx: Ctx) -> None:
    """F13 + T13, T13b, T24, T25: the ablation ladder -- which layer carries it?

    Four models of the same prediction: receptor engagement alone, mechanism
    gains on the graph's transmitter composition alone, the real cut with one
    generic multiplier, and the full model. Each level's *ordering of the whole
    library* is correlated with the full model's, with **signed** rho.

    Two corrections from the second referee round change what this step can
    conclude. Level C's generic multiplier was depression-only, so it could not
    express disinhibition and no compound could come out positive; it is now
    direction-aware (generic magnitude, one bit of sign from the mechanism
    table), with the old rule kept and reported as a floor. And the
    composition-versus-full correlation is meaningless without a reference
    distribution, because both levels are functions of the same gain vector:
    T24 supplies one, and T25 asks whether the verdict survives the engine's
    row normalisation being changed.
    """
    from flylab.analysis.baselines import (
        C_FLOOR_LEVEL,
        DEFAULT_GENERIC_RULE,
        DIRECTION_AWARE_MULTIPLIER_RULE,
        GENERIC_MULTIPLIER_RULE,
        GENERIC_RULE_NOTES,
        LEVELS,
        REPRODUCES_RHO,
        ablation,
        ablation_table,
        composition_dominance_under_normalisations,
        composition_reference_distribution,
        glutamate_sign_reconciliation,
    )

    plt = _plt()
    concs = (1e-6,) if ctx.fast else (1e-8, 1e-7, 1e-6, 1e-5)
    tab = ablation_table(concs_M=concs)
    ctx.put("abl_compounds", len(tab["compounds"]))
    ctx.put("abl_concs", len(tab["concs_M"]))
    ctx.put("abl_levels", list(LEVELS), text=", ".join(LEVELS))
    ctx.put("abl_generic_rule", DIRECTION_AWARE_MULTIPLIER_RULE)
    ctx.put("abl_generic_rule_name", DEFAULT_GENERIC_RULE)
    ctx.put("abl_generic_floor_rule", GENERIC_MULTIPLIER_RULE)
    ctx.put("abl_generic_floor_note", GENERIC_RULE_NOTES["depressant_floor"])
    ctx.put("abl_reproduces_rho", float(REPRODUCES_RHO))

    per_level: dict[str, list[float]] = {lvl: [] for lvl in LEVELS}
    gain_rows: list[dict[str, Any]] = []
    for conc_key, gain in tab["information_gain"].items():
        cmp_block = gain.get("generic_rule_comparison") or {}
        for lvl in gain["levels"]:
            name = lvl["level"]
            rho = lvl.get("spearman_rho_vs_full")
            if rho is not None:
                per_level[name].append(float(rho))
            gain_rows.append(
                {
                    "conc_M": float(conc_key),
                    "level": name,
                    "unit": lvl.get("unit"),
                    "spearman_rho_vs_full": rho,
                    "pearson_r_vs_full": lvl.get("pearson_r_vs_full"),
                    "residual_rms_standardised": lvl.get("residual_rms_standardised"),
                    "residual_rms_hz": lvl.get("residual_rms_hz"),
                    "reproduces_full_ordering": lvl.get("reproduces_full_ordering"),
                    "information_added_vs_previous": lvl.get("information_added_vs_previous"),
                }
            )
        floor = cmp_block.get("floor_row") or {}
        if floor.get("spearman_rho_vs_full") is not None:
            gain_rows.append(
                {
                    "conc_M": float(conc_key),
                    "level": C_FLOOR_LEVEL,
                    "unit": floor.get("unit") or "Hz",
                    "spearman_rho_vs_full": floor.get("spearman_rho_vs_full"),
                    "pearson_r_vs_full": floor.get("pearson_r_vs_full"),
                    "residual_rms_standardised": floor.get("residual_rms_standardised"),
                    "residual_rms_hz": floor.get("residual_rms_hz"),
                    "reproduces_full_ordering": floor.get("reproduces_full_ordering"),
                    "information_added_vs_previous": None,
                }
            )
    # "A_receptor_only" -> "receptor": the key half of a results.json name, so
    # it must stay an identifier; the axis label spaces it out separately.
    short = {lvl: lvl.split("_", 1)[1].replace("_only", "") for lvl in LEVELS}
    for name, vals in per_level.items():
        key = short.get(name, name)
        if not vals:
            continue
        ctx.put(f"abl_rho_{key}_min", float(min(vals)))
        ctx.put(f"abl_rho_{key}_max", float(max(vals)))
    paper = tab["information_gain"][f"{1e-6:.3e}"]
    by_name = {lvl["level"]: lvl for lvl in paper["levels"]}
    ctx.put("abl_rho_topology_paper", _f(by_name["C_topology_only"].get("spearman_rho_vs_full")))
    ctx.put("abl_rms_topology_hz", _f(by_name["C_topology_only"].get("residual_rms_hz")), unit="Hz")
    ctx.put("abl_rho_composition_paper", _f(by_name["B_composition_only"].get("spearman_rho_vs_full")))
    ctx.put("abl_rho_receptor_paper", _f(by_name["A_receptor_only"].get("spearman_rho_vs_full")))
    paper_cmp = paper.get("generic_rule_comparison") or {}
    ctx.put(
        "abl_rho_topology_floor_paper",
        _f((paper_cmp.get("depressant_floor") or {}).get("spearman_rho_vs_full")),
    )
    ctx.put("abl_statements", tab["statements"], text=" ".join(tab["statements"]))

    # ---- the ladder's ranking is not constant across the dose range ----
    conc_dep = tab["concentration_dependence"]
    worst = {r["conc_M"]: r["worst_level"] for r in conc_dep["rows"]}
    worst_floor = {r["conc_M"]: r["worst_level_with_floor_rule"] for r in conc_dep["rows"]}
    ctx.put(
        "abl_worst_level_by_conc",
        {f"{k:.0e}": v for k, v in worst.items()},
        text="; ".join(
            f"{_fmt_M(r['conc_M'])}: {r['worst_level']} ({r['worst_rho']:.3f})"
            for r in conc_dep["rows"]
            if r["worst_rho"] is not None
        ),
    )
    ctx.put(
        "abl_worst_level_by_conc_floor_rule",
        {f"{k:.0e}": v for k, v in worst_floor.items()},
        text="; ".join(
            f"{_fmt_M(r['conc_M'])}: {r['worst_level_with_floor_rule']}"
            for r in conc_dep["rows"]
        ),
    )
    ctx.put("abl_worst_level_constant", len(set(worst.values())) == 1)
    ctx.put(
        "abl_worst_level_constant_floor_rule", len(set(worst_floor.values())) == 1
    )
    ctx.put("abl_concentration_statements", conc_dep["statements"], text=" ".join(conc_dep["statements"]))
    if len(set(worst_floor.values())) > 1:
        ctx.note(
            "the ablation ladder's ranking is concentration-dependent: with the "
            "historical depression-only level-C rule the weakest level is %s. Any "
            "sentence of the form 'the topology-only level is the worst of the four' "
            "is true of one column and false of the others."
            % ctx.text("abl_worst_level_by_conc_floor_rule")
        )

    # within one mechanism class, where the ordering problem is hardest
    one = ablation("imidacloprid", PAPER_CONC)
    nic = (one.get("information_gain_by_class") or {}).get("nicotinic agonist")
    if nic:
        nic_by = {lvl["level"]: lvl for lvl in nic["levels"]}
        ctx.put("abl_nicotinic_n", int(nic_by["B_composition_only"].get("n_compounds") or 0))
        ctx.put("abl_rho_composition_nicotinic", _f(nic_by["B_composition_only"].get("spearman_rho_vs_full")))
        ctx.put("abl_rho_receptor_nicotinic", _f(nic_by["A_receptor_only"].get("spearman_rho_vs_full")))

    if (ctx.get("abl_rho_topology_paper") or 0) >= (ctx.get("abl_rho_topology_floor_paper") or 0) + 0.3:
        ctx.note(
            "the level-C conclusion is withdrawn: with a direction-aware generic rule "
            "the topology-only baseline's rank correlation with the full model at 1 uM "
            "rises from %s (depression-only floor, which cannot express disinhibition "
            "and whose every entry is at most zero) to %s. 'The connectome without the "
            "pharmacology is not a cheap substitute' was a statement about a "
            "sign-broken baseline."
            % (
                ctx.text("abl_rho_topology_floor_paper"),
                ctx.text("abl_rho_topology_paper"),
            )
        )

    # ---- T24: what rho would have been unsurprising -------------------
    ref = composition_reference_distribution(
        conc_M=PAPER_CONC,
        n_draws=10 if ctx.fast else DEFAULT_EFFORT["abl_ref_n_draws"],
        seed=ctx.seed,
    )
    cond = ref["conditions"]
    sparse = cond["random_single_gain_vectors"]["spearman_rho"] or {}
    dense = cond["random_gain_vectors"]["spearman_rho"] or {}
    ctx.put("abl_ref_observed", _f(ref["observed"]["spearman_rho"]))
    ctx.put("abl_ref_n_draws", int(sparse.get("n") or 0))
    ctx.put("abl_ref_matched_median", _f(sparse.get("median")))
    ctx.put("abl_ref_matched_p05", _f(sparse.get("p05")))
    ctx.put("abl_ref_matched_p95", _f(sparse.get("p95")))
    ctx.put("abl_ref_dense_median", _f(dense.get("median")))
    ctx.put("abl_ref_observed_percentile", _f(ref["observed_percentile_of_matched_reference"]))
    ctx.put("abl_ref_shuffled", _f(cond["shuffled_compound_assignment"]["spearman_rho"]))
    ctx.put(
        "abl_ref_shuffled_identical",
        bool(cond["shuffled_compound_assignment"]["identical_to_observed"]),
    )
    ctx.put("abl_ref_no_floor", _f(cond["library_no_floor"]["spearman_rho"]))
    ctx.put("abl_ref_no_floor_n", int(cond["library_no_floor"]["n_compounds"]))
    ctx.put("abl_ref_statements", ref["statements"], text=" ".join(ref["statements"]))
    if cond["shuffled_compound_assignment"]["identical_to_observed"]:
        ctx.note(
            "the composition-versus-full correlation carries no compound-level "
            "information: shuffling the compound-to-gain assignment leaves it "
            "identically %s, and pharmacology-free pseudo-compounds whose gain vectors "
            "have the shape the mechanism rules produce already reach a median of %s. "
            "Both levels are functions of the same gain vector, so the number is a "
            "statement about the map from gains to readouts."
            % (ctx.text("abl_ref_shuffled"), ctx.text("abl_ref_matched_median"))
        )
    save_table(
        ctx,
        "T24_composition_reference",
        [
            {
                "condition": "observed (library gain vectors)",
                "spearman_rho": ref["observed"]["spearman_rho"],
                "n_compounds": ref["observed"]["n_compounds"],
                "draws": None,
                "note": "what the paper used to quote on its own",
            },
            {
                "condition": "library, floor-saturated compounds dropped",
                "spearman_rho": cond["library_no_floor"]["spearman_rho"],
                "n_compounds": cond["library_no_floor"]["n_compounds"],
                "draws": None,
                "note": cond["library_no_floor"]["note"],
            },
            {
                "condition": "matched reference: one gain moved per pseudo-compound",
                "spearman_rho": sparse.get("median"),
                "n_compounds": ref["n_compounds"],
                "draws": sparse.get("n"),
                "note": cond["random_single_gain_vectors"]["note"],
            },
            {
                "condition": "dense reference: every gain moved per pseudo-compound",
                "spearman_rho": dense.get("median"),
                "n_compounds": ref["n_compounds"],
                "draws": dense.get("n"),
                "note": cond["random_gain_vectors"]["note"],
            },
            {
                "condition": "compound labels shuffled",
                "spearman_rho": cond["shuffled_compound_assignment"]["spearman_rho"],
                "n_compounds": ref["n_compounds"],
                "draws": 1,
                "note": cond["shuffled_compound_assignment"]["note"],
            },
        ],
        ["condition", "spearman_rho", "n_compounds", "draws", "note"],
        "Reference distribution for the composition-versus-full rank correlation. The "
        "composition level and the full model are not independent models: both are "
        "functions of the same gain vector, so the correlation has a large structural "
        "floor and a bare value near 0.99 is not interpretable. The matched reference "
        "draws pseudo-compounds with the shape the shipped mechanism rules produce "
        "(one receptor, one transmitter, one gain moved) and no pharmacology at all; "
        "shuffling the compound labels leaves the observed value unchanged, which is a "
        "proof rather than a coincidence.",
        md_fields=["condition", "spearman_rho", "n_compounds", "draws"],
    )

    # ---- T25: the row normalisation, and the glutamate sign ------------
    norm = composition_dominance_under_normalisations(conc_M=PAPER_CONC)
    norm_rows: list[dict[str, Any]] = []
    for mode, block in norm["by_normalisation"].items():
        ctx.put(f"abl_rho_composition_{mode}", _f(block["spearman_rho_b_vs_d"]))
        ctx.put(f"abl_composition_survives_{mode}", bool(block["reproduces_full_ordering"]))
        norm_rows.append(
            {
                "conc_M": PAPER_CONC,
                "normalisation": mode,
                "spearman_rho_b_vs_d": block["spearman_rho_b_vs_d"],
                "pearson_r_b_vs_d": block["pearson_r_b_vs_d"],
                "reproduces_full_ordering": block["reproduces_full_ordering"],
            }
        )
    if True:
        low = composition_dominance_under_normalisations(conc_M=1e-8)
        for mode, block in low["by_normalisation"].items():
            ctx.put(f"abl_rho_composition_{mode}_1e8", _f(block["spearman_rho_b_vs_d"]))
            norm_rows.append(
                {
                    "conc_M": 1e-8,
                    "normalisation": mode,
                    "spearman_rho_b_vs_d": block["spearman_rho_b_vs_d"],
                    "pearson_r_b_vs_d": block["pearson_r_b_vs_d"],
                    "reproduces_full_ordering": block["reproduces_full_ordering"],
                }
            )
    survives = [m for m, b in norm["by_normalisation"].items() if b["reproduces_full_ordering"]]
    ctx.put(
        "abl_composition_survives_normalisations",
        survives,
        text=", ".join(survives) or "none",
    )
    ctx.put("abl_normalisation_statements", norm["statements"], text=" ".join(norm["statements"]))
    if "degree" not in survives:
        ctx.note(
            "the composition verdict is not robust to the engine's row normalisation: "
            "the B-versus-D rank correlation at 1 uM is %s under the shipped `row_abs` "
            "mode and %s under a degree-corrected one. The normalisation is an "
            "undocumented modelling choice that makes each cell's recurrent input a "
            "composition-weighted average of its presynaptic gains."
            % (
                ctx.text("abl_rho_composition_row_abs"),
                ctx.text("abl_rho_composition_degree"),
            )
        )

    glu = glutamate_sign_reconciliation(concs_M=concs)
    ctx.put("abl_glutamate_max_delta_rho", _f(glu["max_abs_delta_rho"]))
    ctx.put("abl_glutamate_statements", glu["statements"], text=" ".join(glu["statements"]))
    for r in glu["rows"]:
        norm_rows.append(
            {
                "conc_M": r["conc_M"],
                "normalisation": "glutamate sign: wholens (B excitatory)",
                "spearman_rho_b_vs_d": r["rho_wholens"],
                "pearson_r_b_vs_d": r["r_wholens"],
                "reproduces_full_ordering": None,
            }
        )
        norm_rows.append(
            {
                "conc_M": r["conc_M"],
                "normalisation": "glutamate sign: rate_engine (B inhibitory)",
                "spearman_rho_b_vs_d": r["rho_rate_engine"],
                "pearson_r_b_vs_d": r["r_rate_engine"],
                "reproduces_full_ordering": None,
            }
        )
    save_table(
        ctx,
        "T25_engine_normalisation",
        norm_rows,
        [
            "conc_M",
            "normalisation",
            "spearman_rho_b_vs_d",
            "pearson_r_b_vs_d",
            "reproduces_full_ordering",
        ],
        "The composition-versus-full rank correlation under each of the rate engine's "
        "row normalisations, and under each of the two glutamate sign conventions the "
        "two levels of the ablation ladder use. The shipped engine divides every row "
        "of the signed weight matrix by its own total absolute input, which makes each "
        "cell's recurrent input a composition-weighted average of its presynaptic "
        "gains; that is a modelling choice, it was documented nowhere, and it is what "
        "the composition-dominance verdict is measured on. Levels are not comparable "
        "across normalisations -- the unnormalised operator is supercritical and the "
        "r_max clip shapes its rates -- only the orderings are.",
        md_fields=["conc_M", "normalisation", "spearman_rho_b_vs_d", "reproduces_full_ordering"],
    )

    save_table(
        ctx,
        "T13_ablation",
        tab["rows"],
        ["compound", "class", "conc_M"] + list(LEVELS) + list(tab["extra_series"]),
        "Ablation ladder. Each level is a prediction for the same "
        "compound-concentration cell from a model that has been denied one layer of "
        "information: receptor engagement only (dimensionless), mechanism gains on the "
        "graph's transmitter composition only (excitation index), the real cut with one "
        "generic multiplier instead of mechanism-specific gains (Hz), and the full "
        "model (Hz). The generic multiplier is direction-aware: generic in magnitude, "
        "with one bit of sign taken from the mechanism table, which is the minimum a "
        "connectome-without-pharmacology model needs to order a library containing "
        "disinhibitors. `C_topology_only_floor` is the historical depression-only rule, "
        "kept as a floor: every one of its entries is at most zero, so it cannot "
        "express disinhibition and is not a competitive baseline. Levels have "
        "different units, so they are compared by ordering, never by value.",
        md_fields=["compound", "class", "conc_M"] + list(LEVELS),
    )
    save_table(
        ctx,
        "T13b_ablation_information",
        gain_rows,
        [
            "conc_M",
            "level",
            "unit",
            "spearman_rho_vs_full",
            "pearson_r_vs_full",
            "residual_rms_standardised",
            "residual_rms_hz",
            "reproduces_full_ordering",
            "information_added_vs_previous",
        ],
        "How much of the full model's ordering of the library each ablated level "
        "recovers, per concentration. `reproduces_full_ordering` is applied to the "
        "**signed** rank correlation: an ordering that is a perfect inversion of the "
        "full model's reproduces nothing, and the previous absolute-value test would "
        "have credited it.",
    )

    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    xs = list(range(len(LEVELS)))
    for i, conc_key in enumerate(sorted(tab["information_gain"])):
        gain = tab["information_gain"][conc_key]
        by = {lvl["level"]: lvl for lvl in gain["levels"]}
        ys = [by[lvl].get("spearman_rho_vs_full") for lvl in LEVELS]
        ax.plot(
            xs, [0.0 if y is None else y for y in ys], "o-",
            color=SEQ[min(i + 1, 4)], lw=1.8, ms=7, label=_fmt_M(float(conc_key)),
        )
        floor_rho = ((gain.get("generic_rule_comparison") or {}).get("depressant_floor") or {}).get(
            "spearman_rho_vs_full"
        )
        if floor_rho is not None:
            ax.plot(
                [xs[LEVELS.index("C_topology_only")]], [floor_rho], "x",
                color=SEQ[min(i + 1, 4)], ms=9, mew=2,
            )
    ax.axhline(0.0, color=TEXT2, lw=1)
    ax.axhline(REPRODUCES_RHO, color=GOOD, lw=1, ls="--")
    ax.text(
        0.02, REPRODUCES_RHO + 0.02, "reproduces the full ordering", fontsize=7,
        color=GOOD, transform=ax.get_yaxis_transform(),
    )
    ax.set_xticks(xs)
    ax.set_xticklabels([short.get(lvl, lvl).replace("_", " ") for lvl in LEVELS])
    ax.set_ylabel("signed Spearman rho of the library ordering vs the full model")
    ax.set_xlabel("information the level is allowed")
    ax.set_title(f"Ablation ladder over the {len(tab['compounds'])}-compound library")
    ax.legend(fontsize=7.5, title="concentration", title_fontsize=7.5)
    ax.grid(lw=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F13_ablation",
        f"Each level predicts the same {len(tab['compounds'])} compounds with one layer "
        "of information removed, and is scored by how well it reproduces the full "
        "model's *ordering* of the library (levels have different units, so values are "
        "not comparable; the correlation is signed, so a perfect inversion scores -1 "
        "rather than 1). Crosses mark the historical depression-only level-C rule, "
        "which applied a multiplier of at most 1.0 to every transmitter alike and so "
        "could not express disinhibition: its poor correlation was structural, not a "
        "finding about connectomes, and the conclusion drawn from it is withdrawn. "
        f"With a direction-aware rule level C rises to "
        f"{ctx.text('abl_rho_topology_paper')} at 1 uM from "
        f"{ctx.text('abl_rho_topology_floor_paper')}. The composition level's agreement "
        "with the full model must be read against the matched reference distribution "
        "of T24, not on its own.",
    )

def step_stability(ctx: Ctx) -> None:
    """F14 + T14, T15: specification robustness and threshold sensitivity (RQ3).

    Re-derives every principal conclusion of this paper under a prespecified
    family of admissible engagement-to-gain transformations, and recomputes the
    amplify/buffer split on the 3x3 threshold grid.
    """
    from flylab.analysis.robustness import (
        MECHANISM_FAMILY,
        default_spec,
        mechanism_spec,
        spec_by_name,
    )

    plt = _plt()
    payload = _background(ctx, "stability", _bg_stability, ctx.fast, ctx.seed)
    res, grid = payload["stability"], payload["threshold"]

    ctx.put("stab_family_size", int(res["family_size"]))
    ctx.put("stab_default_spec", str(res["default_spec"]))
    ctx.put("stab_n_conclusions", len(res["rows"]))
    ctx.put("stab_n_shuffles", int(res["n_shuffles"]))
    ctx.put("stab_shuffle_resolution", float(res["shuffle_resolution"]))
    ctx.put("stab_topology_alpha", float(res["topology_alpha"]))
    ctx.put("stab_n_structural_tests", int(res["n_structural_tests"]))
    ctx.put("stab_topology_engine", str(res["topology_engine"]))
    ctx.put("stab_fdr_can_reject", bool(res["topology_fdr"]["can_reject"]))
    ctx.put("stab_fdr_min_rejections", res["topology_fdr"].get("min_rejections"))
    ctx.put("stab_fdr_n_rejected", int(res["topology_fdr"]["n_rejected"]))
    fragile: list[str] = []
    for row in res["rows"]:
        key = row["conclusion"].split("_")[0]
        ctx.put(f"stab_{key}_claim", row["claim"])
        ctx.put(f"stab_{key}_readout", row.get("readout"))
        ctx.put(f"stab_{key}_kind", row.get("kind"))
        ctx.put(f"stab_{key}_retained", int(row["n_retained"]))
        ctx.put(f"stab_{key}_lost", int(row["n_lost"]))
        ctx.put(f"stab_{key}_undecidable", int(row["n_undecidable"]))
        ctx.put(f"stab_{key}_fraction", float(row["fraction_retained"]))
        ctx.put(
            f"stab_{key}_failing",
            row["failing_specs"],
            text=", ".join(row["failing_specs"]) or "none",
        )
        if row.get("n_retained_uncorrected") is not None:
            ctx.put(f"stab_{key}_retained_uncorrected", int(row["n_retained_uncorrected"]))
            ctx.put(f"stab_{key}_multiplicity", str(row["multiplicity"]))
            eq = row.get("equivalence_breakdown") or {}
            ctx.put(f"stab_{key}_equivalent", int(eq.get("equivalent_within_tolerance", 0)))
            ctx.put(f"stab_{key}_indeterminate", int(eq.get("indeterminate", 0)))
            ctx.put(f"stab_{key}_distinguishable", int(eq.get("distinguishable", 0)))
        if row["fragile"]:
            fragile.append(row["conclusion"])
    if (ctx.get("stab_C3_equivalent") or 0) or (ctx.get("stab_C3_indeterminate") or 0):
        ctx.note(
            "the C3 retention rate is no longer vacuous and is no longer uniform: at "
            "n_shuffles = %d (resolution %.4f, against the previous 6 at which p <= "
            "0.05 was arithmetically unattainable) the conclusion is retained by %d of "
            "%d specifications, but only %d of those reach equivalence within the "
            "prespecified margin; %d are merely indeterminate."
            % (
                int(res["n_shuffles"]),
                float(res["shuffle_resolution"]),
                ctx.get("stab_C3_retained") or 0,
                int(res["family_size"]),
                ctx.get("stab_C3_equivalent") or 0,
                ctx.get("stab_C3_indeterminate") or 0,
            )
        )
    ctx.put("stab_fragile", fragile, text=", ".join(fragile) or "none")
    ctx.put(
        "stab_specification_independent",
        [r["conclusion"] for r in res["rows"] if not r["fragile"]],
        text=", ".join(r["conclusion"] for r in res["rows"] if not r["fragile"]) or "none",
    )

    # the reversal, in its own numbers: the same drug at the same engagement
    c1 = next(r for r in res["rows"] if r["conclusion"].startswith("C1"))
    ctx.put("stab_C1_n_reversing", int(c1["n_lost"]))
    reversing = c1["failing_specs"]
    if reversing:
        # prefer a counter-specification at the *same* coefficient scale as the
        # shipped rule, so the contrast isolates the shape and not the size
        scale = str(res["default_spec"]).partition("@")[2]
        name = next((s for s in reversing if s.endswith(f"@{scale}")), reversing[0])
        ctx.put("stab_C1_example_spec", name)
        alt = spec_by_name(name, MECHANISM_FAMILY)
        from flylab.assays.subgraph import run_subgraph_assay

        for label, spec in (("default", default_spec()), ("monotone", alt)):
            with mechanism_spec(spec):
                nb = run_subgraph_assay("imidacloprid", PAPER_CONC)
            r = nb["readouts"]
            ctx.put(f"stab_{label}_treated_hz", _f(r.get("mean_hz")), unit="Hz")
            ctx.put(f"stab_{label}_vehicle_hz", _f((r.get("vehicle") or {}).get("mean_hz")), unit="Hz")
            ctx.put(f"stab_{label}_g_ach", _f((nb.get("gains") or {}).get("g_ach")))
        t, v = ctx.get("stab_monotone_treated_hz"), ctx.get("stab_monotone_vehicle_hz")
        if t and v:
            ctx.put("stab_monotone_pct_change", 100.0 * (t - v) / v, text=f"{100.0 * (t - v) / v:+.0f}%")
        ctx.note(
            "the suppression conclusion is specification-dependent: %d of %d "
            "prespecified specifications reverse it. Under %s the same compound at the "
            "same engagement EXCITES the network (%s Hz treated vs %s Hz vehicle) "
            "instead of suppressing it. Suppression is a property of the biphasic "
            "desensitisation term, not of the pharmacology-connectome integration."
            % (
                c1["n_lost"],
                c1["n_specs"],
                name,
                ctx.text("stab_monotone_treated_hz"),
                ctx.text("stab_monotone_vehicle_hz"),
            )
        )

    # per-specification diagnostics that the paper quotes
    diag = res["diagnostics"]
    gaps_default = diag.get(res["default_spec"], {})
    ctx.put("stab_mean_gap_nicotinic_default", _f(gaps_default.get("mean_gap_nicotinic")))
    # Under some specifications the circuit never crosses the threshold that
    # defines the index, so the gap is undefined there; average over the
    # reversing specifications that do define one, and say how many that was.
    rev_gaps = [
        float(diag[s0]["mean_gap_nicotinic"])
        for s0 in reversing
        if diag.get(s0, {}).get("mean_gap_nicotinic") is not None
    ]
    ctx.put(
        "stab_mean_gap_nicotinic_monotone",
        (sum(rev_gaps) / len(rev_gaps)) if rev_gaps else None,
    )
    ctx.put("stab_mean_gap_nicotinic_monotone_n", len(rev_gaps))
    # I11: the aggregate hides the spread, and "never less" is false. Report
    # the range over the reversing specifications, and whether any of them
    # buffers *less* than the shipped rule.
    ctx.put("stab_mean_gap_nicotinic_monotone_min", min(rev_gaps) if rev_gaps else None)
    ctx.put("stab_mean_gap_nicotinic_monotone_max", max(rev_gaps) if rev_gaps else None)
    default_gap = ctx.get("stab_mean_gap_nicotinic_default")
    weaker = (
        sorted(
            s0
            for s0 in reversing
            if diag.get(s0, {}).get("mean_gap_nicotinic") is not None
            and default_gap is not None
            and float(diag[s0]["mean_gap_nicotinic"]) > float(default_gap)
        )
        if rev_gaps
        else []
    )
    ctx.put(
        "stab_monotone_buffer_less",
        weaker,
        text=", ".join(weaker) or "none",
    )
    ctx.put("stab_monotone_buffer_less_n", len(weaker))
    if weaker:
        ctx.note(
            "the supplement's claim that the monotone rules 'buffer at least as much, "
            "never less' is false: %d of the %d reversing specifications with a defined "
            "index buffer LESS than the shipped rule (%s), and the aggregate of %s is "
            "carried by the extreme of a range running %s to %s."
            % (
                len(weaker),
                len(rev_gaps),
                ", ".join(weaker),
                ctx.text("stab_mean_gap_nicotinic_monotone"),
                ctx.text("stab_mean_gap_nicotinic_monotone_max"),
                ctx.text("stab_mean_gap_nicotinic_monotone_min"),
            )
        )
    save_table(
        ctx,
        "T14b_specification_diagnostics",
        [
            {
                "spec": name,
                "g_ach_imidacloprid": (d.get("gains_imidacloprid") or {}).get("g_ach"),
                "g_gaba_fipronil": (d.get("gains_fipronil") or {}).get("g_gaba"),
                "mean_si_gap_nicotinic": d.get("mean_gap_nicotinic"),
                "mean_si_gap_nav_ache": d.get("mean_gap_nav_ache"),
                "runtime_s": round(res["runtime_by_spec_s"].get(name, 0.0), 2),
            }
            for name, d in diag.items()
        ],
        ["spec", "g_ach_imidacloprid", "g_gaba_fipronil", "mean_si_gap_nicotinic", "mean_si_gap_nav_ache", "runtime_s"],
        "Per-specification diagnostics behind T14: the gains each specification "
        "produces for the two reference compounds at 1 uM, and the mean "
        "circuit-minus-receptor selectivity gap over the nicotinic and Nav/AChE "
        "sets. An empty gap means no compound of that set crossed the circuit "
        "threshold under that specification, so the index is undefined rather "
        "than zero.",
    )

    save_table(
        ctx,
        "T14_conclusion_stability",
        [
            {
                "conclusion": r["conclusion"],
                "claim": r["claim"],
                "readout": r["readout"],
                "n_specs": r["n_specs"],
                "n_retained": r["n_retained"],
                "n_lost": r["n_lost"],
                "n_undecidable": r["n_undecidable"],
                "fraction_retained": r["fraction_retained"],
                "fraction_retained_decidable": r["fraction_retained_decidable"],
                "n_retained_uncorrected": r.get("n_retained_uncorrected"),
                "multiplicity": r.get("multiplicity"),
                "n_equivalent_within_tolerance": (r.get("equivalence_breakdown") or {}).get(
                    "equivalent_within_tolerance"
                ),
                "n_indeterminate": (r.get("equivalence_breakdown") or {}).get("indeterminate"),
                "failing_specs": r["failing_specs"],
                "undecidable_specs": r["undecidable_specs"],
            }
            for r in res["rows"]
        ],
        [
            "conclusion",
            "claim",
            "readout",
            "n_specs",
            "n_retained",
            "n_lost",
            "n_undecidable",
            "fraction_retained",
            "fraction_retained_decidable",
            "n_retained_uncorrected",
            "multiplicity",
            "n_equivalent_within_tolerance",
            "n_indeterminate",
            "failing_specs",
            "undecidable_specs",
        ],
        f"Conclusion stability matrix: every principal conclusion re-derived under each "
        f"of {res['family_size']} prespecified engagement-to-gain specifications "
        "(biphasic, monotone-linear and saturating shapes at several coefficient "
        "scales). `n_undecidable` counts specifications under which the readout does "
        "not exist (for example a ratio whose denominator is silenced); the headline "
        "fraction counts those against the conclusion. The two topology rows are "
        f"decided on {res['n_shuffles']} permutations per specification (resolution "
        f"{res['shuffle_resolution']:.4f}) and on Benjamini-Hochberg adjusted "
        f"probabilities across the {res['n_structural_tests']} structural tests of the "
        "run, with the uncorrected count beside them. Until v0.6.1 the specification "
        "context rebound the gain function on the assay modules but not on the engine "
        "the permutation path resolves, so every specification fed the *default* gains "
        "to its nulls and the matrix's topology rows were empty; they are now computed "
        "per specification. For a conclusion shaped as a failure to reject, "
        "`n_equivalent_within_tolerance` is the part of the retention that is evidence "
        "of equivalence and `n_indeterminate` the part that is only a non-rejection.",
        md_fields=[
            "conclusion",
            "n_specs",
            "n_retained",
            "n_lost",
            "n_undecidable",
            "n_retained_uncorrected",
            "n_equivalent_within_tolerance",
            "n_indeterminate",
        ],
    )

    # ---- threshold grid ---------------------------------------------
    ctx.put("thr_stable", bool(grid["stable"]))
    ctx.put(
        "thr_compounds",
        list(grid["compounds"]),
        text=", ".join(grid["compounds"]),
    )
    ctx.put("thr_n_compounds", len(grid["compounds"]))
    # I12: T8 and T15 disagreed on how many compounds the circuit buffers at
    # the same default thresholds. They do not measure the same set: T8 scores
    # the whole library on both cuts, T15 only the compounds of the two
    # mechanism classes the conclusion is about. Name the difference rather
    # than leaving two tables to contradict each other.
    default_cell = [
        r
        for r in grid["per_compound"]
        if abs(r["circuit_frac"] - 0.50) < 1e-9 and abs(r["vert_limit"] - 0.20) < 1e-9
    ]
    t15_scored = {r["compound"] for r in default_cell if r["si_gap_circuit_minus_receptor"] is not None}
    t15_unscored = {r["compound"] for r in default_cell if r["si_gap_circuit_minus_receptor"] is None}
    t8_buffered = set(ctx.get("landscape_buffer_compounds") or [])
    t8_amplified = set(ctx.get("landscape_amplify_compounds") or [])
    t8_scored = t8_buffered | t8_amplified
    only_t8 = sorted(t8_scored - set(grid["compounds"]))
    both_disagree = sorted(t8_scored & t15_unscored)
    ctx.put("thr_t15_scored", sorted(t15_scored), text=", ".join(sorted(t15_scored)) or "none")
    ctx.put("thr_t15_unscored", sorted(t15_unscored), text=", ".join(sorted(t15_unscored)) or "none")
    ctx.put(
        "thr_outside_t15_scope",
        only_t8,
        text=", ".join(only_t8) or "none",
    )
    ctx.put("thr_n_outside_t15_scope", len(only_t8))
    ctx.put(
        "thr_scored_in_t8_unscored_in_t15",
        both_disagree,
        text=", ".join(both_disagree) or "none",
    )
    if only_t8 or both_disagree:
        ctx.note(
            "T8 and T15 have different scopes and must be read that way: T8 scores the "
            "whole library on both cuts, T15 only the %d compounds of the two mechanism "
            "classes the conclusion is about. %s %s a circuit index in T8 and %s "
            "outside T15's compound set entirely%s."
            % (
                len(grid["compounds"]),
                ", ".join(only_t8) or "no compound",
                "carries" if len(only_t8) == 1 else "carry",
                "is" if len(only_t8) == 1 else "are",
                (
                    "; " + ", ".join(both_disagree) + " are in T15's set but have no "
                    "defined gap at the default thresholds and are reported unscored "
                    "there"
                )
                if both_disagree
                else "",
            )
        )
    ctx.put("thr_circuit_fracs", grid["circuit_fracs"], text=", ".join(f"{f:.0%}" for f in grid["circuit_fracs"]))
    ctx.put("thr_vert_limits", grid["vert_limits"], text=", ".join(f"{v:.0%}" for v in grid["vert_limits"]))
    for check in grid["class_checks"]:
        k = check["class"]
        ctx.put(f"thr_{k}_matching", int(check["n_matching"]))
        ctx.put(f"thr_{k}_cells", int(check["n_cells"]))
        ctx.put(f"thr_{k}_stable", bool(check["stable"]))
        ctx.put(f"thr_{k}_failing", check["failing_cells"], text=", ".join(check["failing_cells"]) or "none")
    for k in ("nicotinic", "nav_ache"):
        vals = [c.get(f"mean_gap_{k}") for c in grid["grid"] if c.get(f"mean_gap_{k}") is not None]
        if vals:
            ctx.put(f"thr_{k}_gap_min", float(min(vals)))
            ctx.put(f"thr_{k}_gap_max", float(max(vals)))
    if not grid["stable"]:
        ctx.note(
            "the amplify/buffer split is NOT stable across the 3x3 threshold grid: "
            + "; ".join(
                f"{c['class']} fails at {', '.join(c['failing_cells'])}"
                for c in grid["class_checks"]
                if not c["stable"]
            )
            + ". The paper quotes the range over the grid, never the point estimate."
        )
    save_table(
        ctx,
        "T15_threshold_grid",
        grid["grid"],
        [
            "circuit_frac",
            "vert_limit",
            "n_amplify",
            "n_buffer",
            "n_unscored",
            "mean_gap_nicotinic",
            "verdict_nicotinic",
            "mean_gap_nav_ache",
            "verdict_nav_ache",
        ],
        "Threshold sensitivity of the amplify/buffer split. Both thresholds entering "
        "the circuit selectivity index are conventions: the relative circuit change "
        "that defines C_circuit and the vertebrate engagement that defines C_vert. The "
        "mechanism-level verdict is recomputed on every cell of the grid. The counts "
        "are over the compounds of the two mechanism classes the conclusion is about "
        "and not over the whole library, which is why they are smaller than T8's at "
        "the same thresholds; a compound with no circuit threshold inside the tested "
        "ladder is reported unscored here, never as buffered.",
    )

    # ---- F14: the stability matrix -----------------------------------
    import numpy as np

    specs = list(res["specs"])
    concl = [r["conclusion"] for r in res["rows"]]
    M = np.full((len(concl), len(specs)), np.nan)
    for i, name in enumerate(concl):
        for j, spec in enumerate(specs):
            v = res["matrix"][name].get(spec)
            M[i, j] = 1.0 if v is True else (0.0 if v is False else 0.5)
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    cmap = ListedColormap([CRIT, WARN, GOOD])
    fig, ax = plt.subplots(figsize=(max(7.0, 0.42 * len(specs)), 3.6))
    ax.imshow(M, cmap=cmap, norm=BoundaryNorm([-0.25, 0.25, 0.75, 1.25], cmap.N), aspect="auto")
    ax.set_xticks(range(len(specs)))
    ax.set_xticklabels(specs, rotation=90, fontsize=6.5)
    ax.set_yticks(range(len(concl)))
    ax.set_yticklabels(concl, fontsize=8)
    ax.set_title(
        f"Conclusion stability over {len(specs)} admissible gain specifications "
        f"(default: {res['default_spec']})",
        fontsize=10,
    )
    ax.legend(
        handles=[
            Patch(facecolor=GOOD, label="retained"),
            Patch(facecolor=WARN, label="undecidable"),
            Patch(facecolor=CRIT, label="reversed"),
        ],
        fontsize=7,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
    )
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F14_conclusion_stability",
        "Each row is a conclusion this paper could state; each column is one "
        "prespecified, admissible way of turning receptor engagement into synaptic "
        "gain. Green retains the conclusion, red reverses it, amber means the readout "
        "does not exist under that specification. The topology rows are decided on "
        f"{res['n_shuffles']} permutations per specification and on Benjamini-Hochberg "
        f"adjusted probabilities across the {res['n_structural_tests']} structural "
        "tests of the run; until v0.6.1 the specification never reached the engine "
        "those rows were computed on, so their uniformity was guaranteed by "
        "construction rather than measured. Retained by every specification: "
        + (", ".join(r["conclusion"] for r in res["rows"] if not r["fragile"]) or "none")
        + ". The nicotinic suppression conclusion is reversed by every monotone rule, "
        "which is how a result that depends on a modelling choice looks when it is "
        "tested rather than asserted; and for a conclusion shaped as a failure to "
        "reject, the fraction retained is only as strong as the equivalence behind it "
        "(T14's `n_equivalent_within_tolerance` against `n_indeterminate`).",
    )


def step_uncertainty(ctx: Ctx) -> None:
    """F15 + T16, T17: the global uncertainty budget and the value of information.

    A variance-based (Sobol') attribution over every assumption the model
    exposes, validated against the analytic Ishigami benchmark, and then
    converted into a ranking of the experiments that would remove the most
    model variance.
    """
    from flylab.analysis.uncertainty_global import (
        ISHIGAMI_REFERENCE,
        ishigami_unit,
        sobol_analysis,
    )
    from flylab.analysis.voi import voi

    plt = _plt()
    res = _background(ctx, "sobol", _bg_sobol, ctx.fast, ctx.seed)

    ctx.put("unc_n_base", int(res["n_base"]))
    ctx.put("unc_n_factors", int(res["n_factors"]))
    ctx.put("unc_evaluations", int(res["n_evaluations"]))
    ctx.put("unc_estimator", str(res["estimator"]))
    ctx.put("unc_output_mean", _f(res["output_mean"]), unit="Hz")
    ctx.put("unc_variance", _f(res["output_variance"]), unit="Hz^2")
    ctx.put("unc_interaction_share", _f(res["interaction_share"]))
    ctx.put("unc_interaction_share_clipped", _f(res["interaction_share_clipped"]))
    ctx.put("unc_sum_first_order", _f(res["sum_first_order"]))
    ctx.put("unc_sum_first_order_clipped", _f(res["sum_first_order_clipped"]))
    for r in res["rows"]:
        ctx.put(f"unc_S_{r['factor']}", _f(r["first_order"]))
        ctx.put(f"unc_ST_{r['factor']}", _f(r["total_order"]))
        ci = r.get("first_order_ci")
        if ci:
            ctx.put(
                f"unc_CI_{r['factor']}",
                [float(ci[0]), float(ci[1])],
                text=f"({ci[0]:+.3f}, {ci[1]:+.3f})",
            )
    ranked = sorted(res["rows"], key=lambda r: -r["first_order"])
    ctx.put("unc_top_factor", ranked[0]["factor"])
    ctx.put("unc_second_factor", ranked[1]["factor"])

    # I5: the noise floor is the largest |negative S1|, not the null factor's
    # own draw. A true first-order index cannot be negative, so that magnitude
    # is a measured lower bound on the estimator's error -- and on the shipped
    # run it is several times the null factor's.
    resolution = res["resolution"]
    ctx.put("unc_noise_floor", _f(resolution["noise_floor"]))
    ctx.put("unc_noise_floor_factor", str(resolution["noise_floor_factor"]))
    ctx.put("unc_noise_floor_definition", str(resolution["noise_floor_definition"]))
    noise = next((r for r in res["rows"] if r["factor"] == "lif_seed"), None)
    if noise:
        ctx.put("unc_null_factor_S", _f(noise["first_order"]))
    ctx.put(
        "unc_resolved",
        list(resolution["resolved"]),
        text=", ".join(resolution["resolved"]) or "none",
    )
    ctx.put("unc_n_resolved", len(resolution["resolved"]))
    ctx.put(
        "unc_at_noise_floor",
        list(resolution["unresolved"]),
        text=", ".join(resolution["unresolved"]) or "none",
    )
    ctx.put(
        "unc_null_control_factors",
        list(resolution["null_control"]),
        text=", ".join(resolution["null_control"]) or "none",
    )
    ctx.put("unc_resolution_statement", str(resolution["statement"]))
    if resolution["noise_floor_factor"] and resolution["noise_floor_factor"] != "lif_seed":
        ctx.note(
            "the estimator's noise floor is %s, set by the most negative first-order "
            "estimate (%s), not by the declared null factor lif_seed at %s. Only %s "
            "%s a first-order confidence interval that excludes zero; every other "
            "factor is unresolved at this sample size, and a negative point estimate "
            "is estimator error rather than a negative contribution."
            % (
                ctx.text("unc_noise_floor"),
                resolution["noise_floor_factor"],
                ctx.text("unc_null_factor_S"),
                ctx.text("unc_resolved"),
                "has" if len(resolution["resolved"]) == 1 else "have",
            )
        )
    conv = res["convergence"][-1]
    ctx.put("unc_max_delta_first", _f(conv.get("max_abs_delta_first")))

    # estimator validation against an analytic benchmark
    n_ish = 1024 if ctx.fast else DEFAULT_EFFORT["ish_n_base"]
    ish = sobol_analysis(
        model=ishigami_unit, factor_names=["x1", "x2", "x3"], n_base=n_ish, seed=ctx.seed, n_boot=0
    )
    got = [r["first_order"] for r in sorted(ish["rows"], key=lambda r: r["factor"])]
    want = ISHIGAMI_REFERENCE["first_order"]
    err = max(abs(g - w) for g, w in zip(got, want))
    ctx.put("ish_n_base", n_ish)
    ctx.put("ish_estimated", [round(g, 4) for g in got], text=", ".join(f"{g:.4f}" for g in got))
    ctx.put("ish_analytic", [round(w, 4) for w in want], text=", ".join(f"{w:.4f}" for w in want))
    ctx.put("ish_max_abs_error", float(err))

    save_table(
        ctx,
        "T16_uncertainty_budget",
        [
            {
                "source": b["source"],
                "kind": b.get("kind"),
                "share_of_variance_S1": b["share_of_variance"],
                "share_with_interactions_ST": b.get("share_with_interactions"),
                "variance_removed_hz2": b["variance_removed"],
                "ci_low": (b.get("ci") or [None, None])[0],
                "ci_high": (b.get("ci") or [None, None])[1],
            }
            for b in res["budget"]
        ],
        ["source", "kind", "share_of_variance_S1", "share_with_interactions_ST", "variance_removed_hz2", "ci_low", "ci_high"],
        f"Global uncertainty budget for the neighbourhood mean rate "
        f"({res['compound']} at {_fmt_M(res['conc_M'])}, {res['n_evaluations']} model "
        "evaluations, Jansen estimators on a Saltelli cross-sample). S1 is the variance "
        "share resolving that factor alone would remove; ST includes its interactions. "
        "Estimates are **unclipped**, so a factor whose true index is zero can come out "
        "negative; the largest such magnitude is the estimator's measured noise floor "
        "and is reported with the factor that set it. A factor counts as resolved only "
        "when its bootstrap interval excludes zero, which is a stronger test than "
        "exceeding that floor. The shares do not sum to one; the remainder is the "
        "interaction row, and it is reported both on the raw estimates and with the "
        "negative ones clipped, because leaving them in counts estimator noise as "
        "interaction. A null factor with no effect on the deterministic rate engine is "
        "included on purpose as a control.",
    )

    v = voi(result=res)
    ctx.put("voi_var_total", _f(v["output_variance"]), unit="Hz^2")
    for i, r in enumerate(v["rows"][:3], 1):
        ctx.put(f"voi_rank{i}_factor", r["factor"])
        ctx.put(f"voi_rank{i}_var", _f(r["voi_var"]), unit="Hz^2")
        ctx.put(f"voi_rank{i}_fraction", _f(r["voi_fraction"]))
        ctx.put(f"voi_rank{i}_experiment", r["experiment"])
        ctx.put(f"voi_rank{i}_cost", r["cost"])
        ctx.put(f"voi_rank{i}_state", str(r.get("state")))
    ctx.put(
        "voi_unresolved_factors",
        list(v["unresolved_factors"]),
        text=", ".join(v["unresolved_factors"]) or "none",
    )
    ctx.put(
        "voi_resolved_factors",
        list(v["resolved_factors"]),
        text=", ".join(v["resolved_factors"]) or "none",
    )
    ctx.put("voi_recommendation", v["recommendation"], text=str(v["recommendation"]))
    free = [r["factor"] for r in v["rows"] if r.get("cost") == "none (compute only)" and r["voi_fraction"] > 0]
    ctx.put("voi_no_experiment_needed", free, text=", ".join(free) or "none")
    save_table(
        ctx,
        "T17_value_of_information",
        [
            {
                "rank": r["rank"],
                "factor": r["factor"],
                "state": r.get("state"),
                "voi_fraction_of_var": r["voi_fraction"],
                "voi_fraction_raw": r.get("voi_fraction_raw"),
                "clipped_from_negative": r.get("clipped_from_negative"),
                "voi_var_hz2": r["voi_var"],
                "voi_upper_bound_var_hz2": r["voi_upper_bound_var"],
                "sd_reduction_hz": r["sd_reduction"],
                "experiment": r["experiment"],
                "resolves": r["resolves"],
                "cost": r["cost"],
                "blocking_gate": r["blocking_gate"],
            }
            for r in v["rows"]
        ],
        [
            "rank",
            "factor",
            "state",
            "voi_fraction_of_var",
            "voi_fraction_raw",
            "clipped_from_negative",
            "voi_var_hz2",
            "voi_upper_bound_var_hz2",
            "sd_reduction_hz",
            "experiment",
            "resolves",
            "cost",
            "blocking_gate",
        ],
        "Value of information. VOI_j = S_j x Var(Y) is the model variance that would "
        "disappear if assumption j were resolved exactly while everything else stayed "
        "as uncertain as it is; the upper bound uses the total-order index and is what "
        "resolving j *last* would buy. A value of information cannot be negative, so a "
        "negative first-order estimate is clipped to zero for the decision value and "
        "kept unclipped in `voi_fraction_raw`. `state` says whether the sample could "
        "resolve the factor at all: an unresolved factor carries no ranking claim, "
        "which is a statement about the sample size rather than about the factor. This "
        "is variance of a model output under assumed input ranges, not an expected gain "
        "in accuracy about a living fly.",
        md_fields=["rank", "factor", "state", "voi_fraction_of_var", "voi_var_hz2", "experiment", "cost"],
    )

    import numpy as np

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), gridspec_kw={"width_ratios": [1.25, 1.0]})
    ax = axes[0]
    order = sorted(res["rows"], key=lambda r: r["first_order"])
    y = np.arange(len(order))
    ax.barh(y + 0.19, [r["total_order"] for r in order], height=0.36, color=SEQ[1], label="total order (ST)")
    ax.barh(y - 0.19, [r["first_order"] for r in order], height=0.36, color=SEQ[3], label="first order (S1)")
    for r, yy in zip(order, y):
        ci = r.get("first_order_ci")
        if ci:
            ax.plot(
                [float(ci[0]), float(ci[1])], [yy - 0.19, yy - 0.19],
                color=TEXT1, lw=1.0, solid_capstyle="butt", zorder=3,
            )
    floor = float(resolution["noise_floor"])
    ax.axvline(floor, color=CRIT, lw=1, ls="--")
    ax.text(
        floor,
        len(order) - 0.4,
        f" estimator noise floor ({resolution['noise_floor_factor']})",
        fontsize=7,
        color=CRIT,
        va="top",
    )
    ax.axvline(0.0, color=TEXT2, lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [
            r["factor"] + ("" if resolution["by_factor"].get(r["factor"]) == "resolved" else " *")
            for r in order
        ],
        fontsize=8,
    )
    ax.set_xlabel("share of the variance of the neighbourhood mean rate")
    ax.set_title(f"Uncertainty budget ({res['n_evaluations']} evaluations)", fontsize=10)
    ax.legend(fontsize=7.5)
    ax.grid(axis="x", lw=0.5)
    ax.set_axisbelow(True)

    ax2 = axes[1]
    top = v["rows"][:6][::-1]
    y2 = np.arange(len(top))
    ax2.barh(y2, [r["voi_var"] for r in top], color=[GOOD if r.get("cost") == "none (compute only)" else SEQ[2] for r in top])
    ax2.set_yticks(y2)
    ax2.set_yticklabels([r["factor"] for r in top], fontsize=8)
    ax2.set_xlabel("model variance an exact answer would remove (Hz$^2$)")
    ax2.set_title(f"Value of information (Var(Y) = {res['output_variance']:.3g} Hz$^2$)", fontsize=10)
    ax2.grid(axis="x", lw=0.5)
    ax2.set_axisbelow(True)
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F15_uncertainty_voi",
        "Left: the share of the variance of the neighbourhood mean rate attributable to "
        "each assumption, first-order and total-order, from a Saltelli cross-sample "
        "with Jansen estimators, with bootstrap first-order confidence intervals. A "
        "factor is **resolved** only when its interval excludes zero; the starred rows "
        "are not, and a factor's index exceeding the noise floor does not resolve it. "
        "The dashed line is that floor, defined as the largest magnitude among the "
        "*negative* first-order estimates: a true first-order index cannot be negative, "
        "so that magnitude is a measured lower bound on the estimator's error, and on "
        "this run it is several times the declared null factor's own draw "
        f"({ctx.text('unc_noise_floor')} set by {ctx.get('unc_noise_floor_factor')}, "
        f"against lif_seed at {ctx.text('unc_null_factor_S')}). Right: the same shares "
        "scaled back into the readout's variance and mapped onto the experiment that "
        "would resolve each assumption; a negative estimate is clipped to zero for the "
        "decision value, because a value of information cannot be negative. Green marks "
        "a factor that needs no experiment at all, only a re-analysis of data already "
        "held.",
    )


def step_claims(ctx: Ctx) -> None:
    """T18: the provenance chain behind one result, generated rather than asserted.

    ``flylab.analysis.claims`` walks the dependency chain of a readout and
    labels every link by what it rests on, then splits the same result into
    facts, model inference and unknowns. It is the machine-readable form of the
    manuscript's separation between model-generated observations and encoded
    assumptions.
    """
    from flylab.analysis.claims import (
        CLAIMS_VERSION,
        claim_audit,
        fact_inference_unknown,
        label_counts,
    )

    audit = claim_audit(compound="imidacloprid", conc_M=PAPER_CONC, assay="subgraph")
    counts = label_counts(audit)
    split = fact_inference_unknown(audit=audit)

    ctx.put("claims_version", str(CLAIMS_VERSION))
    ctx.put("claims_chain_links", len(audit["chain"]))
    for label, n in counts.items():
        ctx.put(f"claims_{label.lower().replace('-', '_')}", int(n))
    ctx.put("claims_n_facts", len(split["facts"]))
    ctx.put("claims_n_inference", len(split["model_inference"]))
    ctx.put("claims_n_unknown", len(split["unknown"]))
    # the noun has to agree with the count, and the prose cannot know it
    ctx.put(
        "claims_facts_phrase",
        len(split["facts"]),
        text=_plural(len(split["facts"]), "fact"),
    )
    ctx.put(
        "claims_inference_phrase",
        len(split["model_inference"]),
        text=_plural(len(split["model_inference"]), "model inference"),
    )
    ctx.put(
        "claims_unknown_phrase",
        len(split["unknown"]),
        text=_plural(len(split["unknown"]), "unknown"),
    )
    ctx.put(
        "claims_observed_phrase",
        int(counts.get("OBSERVED", 0)),
        text=_plural(int(counts.get("OBSERVED", 0)), "link"),
    )
    observed = [link["step"] for link in audit["chain"] if link["label"] == "OBSERVED"]
    ctx.put("claims_observed_steps", observed, text=", ".join(observed) or "none")
    if len(observed) <= 1:
        ctx.note(
            "the provenance chain behind a circuit readout has %d links and exactly "
            "%d of them (%s) is a measurement of the system being simulated; the rest "
            "are literature reuse, asserted modelling choices and computation."
            % (len(audit["chain"]), len(observed), ", ".join(observed) or "none")
        )

    save_table(
        ctx,
        "T18_claim_provenance",
        [
            {
                "step": link["step"],
                "label": link["label"],
                "classification": link["classification"],
                "contributes": link.get("contributes") or link.get("statement"),
                "assumptions": "; ".join(link.get("assumptions") or []),
                "unknowns": "; ".join(link.get("unknowns") or []),
            }
            for link in audit["chain"]
        ],
        ["step", "label", "classification", "contributes", "assumptions", "unknowns"],
        "Provenance chain behind one circuit readout (imidacloprid at 1 uM on the "
        "`named` cut), generated by `flylab.analysis.claims`. `OBSERVED` is a "
        "measurement of the system being simulated; `LITERATURE-DERIVED` is a "
        "measurement of something else reused here; `MODEL-ASSUMPTION` is asserted by "
        "FlyLab; `COMPUTED` is arithmetic on the links above it. The same function "
        "splits the result into facts, model inference and unknowns.",
        md_fields=["step", "label", "classification", "contributes"],
    )



def step_scale(ctx: Ctx) -> None:
    """Where the dependence verdict might settle with the size of the cut.

    The verdict reversed between the two committed cuts, and nothing in the
    method says where -- or whether -- it stabilises.  A scaling study
    (:func:`flylab.analysis.scale.dependence_vs_scale`) answers that on a
    ladder of nested cuts of 1k to 50k cells built by one fixed recipe.  Only
    the 1k rung is small enough to commit; the rest are CI artifacts, so the
    study is run outside this pipeline and its result is **read** here from
    ``<outdir>/scale_study.json`` (or ``data/derived/scale_study.json``) when
    one exists.  Without it, this step still records what is on disk, what the
    study would cost at each rung (T26) and an explicit "not yet measured"
    statement, so the manuscript's open question is rendered from a key rather
    than typed by hand.
    """
    from flylab.analysis.dependence import cut_census
    from flylab.analysis.scale import (
        LADDER_SIZES,
        cut_path,
        feasibility_frontier,
        scale_report,
        verdict_stability,
    )

    ladder = [f"scale_{k // 1000}k" for k in LADDER_SIZES]
    ctx.put("scale_ladder_sizes", list(LADDER_SIZES), text=", ".join(f"{k // 1000}k" for k in LADDER_SIZES))
    # check existence by path only -- never JSON-load the multi-hundred-MB
    # ladder rungs (5k-50k) just to count them; the committed 1k rung is the
    # only one this step reads, for its census.
    on_disk = [name for name in ladder if cut_path(name) is not None]
    absent = [name for name in ladder if cut_path(name) is None]
    ctx.put("scale_rungs_on_disk", on_disk, text=", ".join(on_disk) or "none")
    ctx.put("scale_rungs_missing", absent, text=", ".join(absent) or "none")
    ctx.put("scale_n_rungs_on_disk", len(on_disk))

    # the committed rung's structure, beside the two cuts of T19
    small_path = cut_path("scale_1k")
    if small_path is not None:
        cen = cut_census(str(small_path))
        ctx.put("scale_1k_name", "scale_1k")
        ctx.put("scale_1k_nodes", int(cen.get("n_nodes") or 0))
        ctx.put("scale_1k_edges", int(cen.get("n_edges") or 0))
        ctx.put("scale_1k_mean_degree", _f(cen.get("mean_degree")), text=f"{float(cen.get('mean_degree') or 0):.1f}")
        share = cen.get("share_onto_seeds")
        ctx.put("scale_1k_share_onto_seeds", _f(share), text=(f"{100 * float(share):.0f}%" if share is not None else "n/a"))
        ctx.put("scale_1k_in_star", bool(cen.get("in_star")), text="yes" if cen.get("in_star") else "no")
    else:
        for key, txt in (
            ("scale_1k_name", "n/a"), ("scale_1k_nodes", "n/a"), ("scale_1k_edges", "n/a"),
            ("scale_1k_mean_degree", "n/a"), ("scale_1k_share_onto_seeds", "n/a"), ("scale_1k_in_star", "n/a"),
        ):
            ctx.put(key, None, text=txt)

    # what the study costs, rung by rung (measured constants, one core)
    fr = feasibility_frontier()
    rows = []
    for r in fr["rows"]:
        rows.append(
            {
                "scale": r["scale"],
                "n_nodes": r["n_nodes"],
                "n_edges": r["n_edges"],
                "single_run_s": round(float(r["single_run_s"]), 3),
                "dependence_profile_h": round(float(r["dependence_profile_h"]), 2),
                "dependence_profile_feasible": r["dependence_profile_feasible"],
                "landscape_days": round(float(r["landscape_days"]), 2),
                "landscape_feasible": r["landscape_feasible"],
                "lif_dense_gb": round(float(r["lif_dense_gb"]), 2),
                "lif_feasible": r["lif_feasible"],
            }
        )
    by_scale = {r["scale"]: r for r in fr["rows"]}
    top = by_scale.get(f"scale_{LADDER_SIZES[-1] // 1000}k")
    ctx.put("scale_frontier_n_paper", int(fr["n_paper"]))
    ctx.put("scale_largest_feasible_profile", fr["largest_cut_with_feasible_profile"], text=str(fr["largest_cut_with_feasible_profile"] or "none"))
    ctx.put("scale_largest_feasible_lif", fr["largest_cut_with_feasible_lif"], text=str(fr["largest_cut_with_feasible_lif"] or "none"))
    ctx.put("scale_top_rung", f"{LADDER_SIZES[-1] // 1000}k")
    ctx.put("scale_top_profile_hours", _f(top["dependence_profile_h"]) if top else None, text=(f"{top['dependence_profile_h']:.1f}" if top else "n/a"))
    ctx.put("scale_top_landscape_days", _f(top["landscape_days"]) if top else None, text=(f"{top['landscape_days']:.0f}" if top else "n/a"))
    whole = by_scale.get("whole_cns_w5")
    ctx.put("scale_whole_cns_profile_hours", _f(whole["dependence_profile_h"]) if whole else None, text=(f"{whole['dependence_profile_h']:.0f}" if whole else "n/a"))
    save_table(
        ctx,
        "T26_scale_frontier",
        rows,
        list(rows[0].keys()) if rows else ["scale"],
        "What a dependence analysis costs at each scale, from the measured per-shuffle "
        f"and per-edge constants of `flylab.analysis.scale` on one core: a profile at n = {fr['n_paper']} "
        f"permutations, an {fr['landscape_cells']}-cell landscape at n = {fr['n_landscape']}, and the dense "
        "matrix the LIF engine materialises. `feasible` means one result in under a day on one core, "
        "which is a generous bar for a single result and a hopeless one for a study repeated over "
        "the specification family. Node and edge counts are measured on MaleCNS v1.0.",
        md_fields=["scale", "n_nodes", "n_edges", "dependence_profile_h", "dependence_profile_feasible", "landscape_days", "lif_dense_gb", "lif_feasible"],
    )

    # the study itself, if it has been run
    study = None
    for cand in (ctx.outdir / "scale_study.json", REPO_ROOT / "data" / "derived" / "scale_study.json"):
        if cand.is_file():
            try:
                study = json.loads(cand.read_text())
                ctx.put("scale_study_source", str(cand.relative_to(REPO_ROOT)) if cand.is_relative_to(REPO_ROOT) else str(cand))
                break
            except (OSError, ValueError) as exc:
                ctx.note(f"scale study at {cand} could not be read: {exc}")
    if isinstance(study, dict) and study.get("rows"):
        rep = scale_report(study)
        stab = study.get("stability") or verdict_stability(study)
        clean = verdict_stability(study, recipe="scale_ladder")
        cuts_run = sorted({r["cut"] for r in study["rows"]}, key=lambda c: next((int(r["n_nodes"]) for r in study["rows"] if r["cut"] == c), 0))
        ctx.put("scale_study_status", "reported")
        ctx.put("scale_study_cuts", cuts_run, text=", ".join(cuts_run))
        ctx.put("scale_study_n_cuts", len(cuts_run))
        ctx.put("scale_study_compounds", list(study.get("compounds") or []), text=", ".join(study.get("compounds") or []))
        ctx.put("scale_study_settled", bool(stab.get("all_settled")), text="yes" if stab.get("all_settled") else "no")
        ctx.put("scale_study_settled_clean_ladder", bool(clean.get("all_settled")), text="yes" if clean.get("all_settled") else "no")
        ctx.put("scale_study_statement", rep["statements"], text=" ".join(rep["statements"]))
        seqs = []
        for comp, b in (clean.get("by_compound") or {}).items():
            seqs.append(f"{comp}: " + " -> ".join(f"{c} ({v})" for c, v in zip(b["cuts"], b["sequence"])))
        ctx.put("scale_study_sequences", seqs, text="; ".join(seqs) or "n/a")

        # --- the paper-quotable shape of the result -------------------
        # Rows, smallest cut first.  The key methodological point is that the
        # inversion does NOT track node count: `scale_1k` has fewer nodes than
        # `named` and a mean degree twenty times higher, and it already comes
        # out topology-dependent.  What separates them is recurrence.
        rows_by_size = sorted(study["rows"], key=lambda r: (int(r["n_nodes"]), int(r["n_edges"])))
        ns = [int(r["n"]) for r in rows_by_size]
        ctx.put("scale_study_n_range", [min(ns), max(ns)], text=f"{min(ns)}-{max(ns)}")
        res = [float(r["p_resolution"]) for r in rows_by_size if r.get("p_resolution")]
        ctx.put("scale_study_resolution_range", [min(res), max(res)], text=f"{min(res):.4f}-{max(res):.4f}")
        comp_cuts = sorted({r["cut"] for r in rows_by_size if r.get("class") == "composition-dominated"})
        ctx.put("scale_study_composition_dominated_cuts", comp_cuts, text=", ".join(comp_cuts) or "none")
        ctx.put("scale_study_n_composition_dominated", len(comp_cuts))
        every_null = sorted({r["cut"] for r in rows_by_size if r.get("necessary_information_level") == "real_connectome"})
        ctx.put("scale_study_every_null_cuts", every_null, text=", ".join(every_null) or "none")
        # Rungs whose permutation budget is so small that the smallest
        # attainable probability is within a factor of two of alpha: a
        # rejection there is the smallest the test can express, not a strong
        # one, and every such rejection sits exactly at the floor.
        alpha = float(study.get("alpha") or 0.05)
        coarse = [
            r for r in rows_by_size
            if r.get("p_resolution") and float(r["p_resolution"]) >= alpha / 2.0
        ]
        coarse_cuts = sorted({r["cut"] for r in coarse}, key=lambda c: next(int(x["n_nodes"]) for x in rows_by_size if x["cut"] == c))
        ctx.put("scale_study_coarse_cuts", coarse_cuts, text=", ".join(coarse_cuts) or "none")
        ctx.put(
            "scale_study_coarse_note",
            {"cuts": coarse_cuts, "alpha": alpha},
            text=(
                (
                    "at " + ", ".join(coarse_cuts) + " the permutation budget falls to n = "
                    + str(min(int(r["n"]) for r in coarse))
                    + ", where the smallest attainable probability is "
                    + f"{max(float(r['p_resolution']) for r in coarse):.3f}"
                    + f" against alpha = {alpha:g}: a rejection there is the smallest the test can express"
                )
                if coarse
                else "every rung was run at a resolution comfortably finer than alpha"
            ),
        )
        for comp in (study.get("compounds") or []):
            rs = [r for r in rows_by_size if r["compound"] == comp]
            if not rs:
                continue
            ctx.put(
                f"scale_study_{comp}_sequence",
                [{"cut": r["cut"], "n_nodes": r["n_nodes"], "class": r["class"]} for r in rs],
                text="; ".join(f"{r['cut']} ({r['n_nodes']} nodes, mean degree {float(r.get('mean_degree') or 0):.1f}): {r['class']}" for r in rs),
            )
            rel = [abs(float(r["relative_effect"])) for r in rs if r.get("relative_effect") is not None]
            ctx.put(
                f"scale_study_{comp}_relative_range",
                [min(rel), max(rel)] if rel else None,
                text=(f"{min(rel):.2f}-{max(rel):.2f}" if rel else "n/a"),
            )
        # the density-not-size sentence, generated rather than asserted
        instar = [r for r in rows_by_size if r.get("in_star")]
        dense_small = [
            r for r in rows_by_size
            if not r.get("in_star") and r.get("class") == "topology-dependent"
            and instar and int(r["n_nodes"]) <= min(int(x["n_nodes"]) for x in instar)
        ]
        if instar and dense_small:
            a, b = instar[0], dense_small[0]
            ctx.put(
                "scale_study_density_not_size",
                {"in_star": a["cut"], "dense": b["cut"]},
                text=(
                    f"the inversion does not track the number of cells: `{b['cut']}` has "
                    f"{int(b['n_nodes'])} nodes against `{a['cut']}`'s {int(a['n_nodes'])} and a mean "
                    f"degree of {float(b.get('mean_degree') or 0):.1f} against {float(a.get('mean_degree') or 0):.2f}, "
                    f"and it is already {b['class']}"
                ),
            )
        else:
            ctx.put("scale_study_density_not_size", None, text="n/a")
        save_table(
            ctx,
            "T27_scale_study",
            rep["rows"],
            list(rep["columns"]),
            "The same dependence profile across the cuts the scaling study could afford, "
            "smallest first: the cut's structure, the permutation count it was run at (which "
            "falls with the edge count, so the resolution coarsens up the ladder), the per-mode "
            "empirical probabilities, the necessary information level and the class. `named` "
            "and `taste_motor` are built by a different recipe from the `scale_*` rungs, so the "
            "clean scale axis is the `scale_*` rows alone. " + str(stab.get("caveat") or ""),
            md_fields=["cut", "n_nodes", "n_edges", "mean_degree", "compound", "n", "p_degree", "p_weight", "p_sign_wm", "necessary_information_level", "class"],
        )
        ctx.note(
            "scale study: %d cuts (%s); verdicts settled on the clean ladder: %s. %s"
            % (len(cuts_run), ", ".join(cuts_run), "yes" if clean.get("all_settled") else "no", " ".join(rep["statements"]))
        )
    else:
        ctx.put("scale_study_status", "not yet run")
        ctx.put("scale_study_source", None, text="none")
        for key in (
            "scale_study_cuts",
            "scale_study_compounds",
            "scale_study_sequences",
            "scale_study_n_range",
            "scale_study_resolution_range",
            "scale_study_composition_dominated_cuts",
            "scale_study_every_null_cuts",
            "scale_study_coarse_cuts",
            "scale_study_coarse_note",
            "scale_study_density_not_size",
        ):
            ctx.put(key, None, text="n/a")
        ctx.put("scale_study_n_composition_dominated", 0)
        for comp in ("imidacloprid", "fipronil"):
            ctx.put(f"scale_study_{comp}_sequence", None, text="n/a")
            ctx.put(f"scale_study_{comp}_relative_range", None, text="n/a")
        ctx.put("scale_study_n_cuts", 0)
        ctx.put("scale_study_settled", None, text="not measured")
        ctx.put("scale_study_settled_clean_ladder", None, text="not measured")
        ctx.put(
            "scale_study_statement",
            None,
            text=(
                "Where the verdict settles with the size of the cut has not been measured at "
                "this commit: the scale ladder's rungs on disk are "
                + (", ".join(on_disk) or "none")
                + " and the study over the full ladder had not been run."
            ),
        )
        ctx.note(
            "the scaling study has not been run at this commit (rungs on disk: %s; missing: %s); "
            "the manuscript reports the open question, not a verdict."
            % (", ".join(on_disk) or "none", ", ".join(absent) or "none")
        )


def step_architecture(ctx: Ctx) -> None:
    """F1: the architecture / dataflow schematic (PNG + SVG)."""
    plt = _plt()
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    fig, ax = plt.subplots(figsize=(9.4, 6.0))
    ax.set_xlim(-2, 102)
    ax.set_ylim(-11, 101)
    ax.axis("off")

    def box(x, y, w, h, title, body, fc, ec=None, fontsize=8):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.6,rounding_size=1.6",
                linewidth=1.1,
                facecolor=fc,
                edgecolor=ec or TEXT2,
            )
        )
        ax.text(x + w / 2, y + h - 3.2, title, ha="center", va="top", fontsize=fontsize + 0.8, weight="bold", color=TEXT1)
        ax.text(x + w / 2, y + h - 8.0, body, ha="center", va="top", fontsize=fontsize, color=TEXT2, linespacing=1.35)

    def arrow(x1, y1, x2, y2, label=None, colour=TEXT2):
        ax.add_patch(
            FancyArrowPatch(
                (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12, linewidth=1.2, color=colour,
                shrinkA=1, shrinkB=1,
            )
        )
        if label:
            ax.text((x1 + x2) / 2 + 1.2, (y1 + y2) / 2, label, fontsize=7.2, color=colour, ha="left", va="center")

    box(30, 88, 40, 11, "compound + concentration", "CLI / API / bench UI", "#ffffff")
    box(26, 68, 48, 15, "occupancy engine", "Hill occupancy, sourced EC50 + tier\ngenotype shift  |  mixtures  |  exposure C(t)", SEQ[0])
    box(2, 46, 42, 16, "insect scorecard", "nAChR, RDL, GluCl, AChE, Nav, OctR\n-> mechanism table -> gain patch", "#dbe9fb", INSECT)
    box(56, 46, 42, 16, "vertebrate scorecard", "a4b2, a7, GABA-A, GlyR, AChE, Nav1.x\noccupancy only: no vertebrate circuit", "#fbe2d6", VERTEBRATE)
    box(2, 23, 42, 19, "MaleCNS netlist (CC-BY)", "1-hop MN9/DNp01 cut  1126 / 1360\ntaste-motor cut  1841 / 19066\ncensus 165122 traced cells", "#e6f6ef", "#1baf7a")
    box(2, 4, 42, 15, "circuit runtime", "rate network (deterministic)\nLIF, Shiu-style + background Poisson", "#ffffff")
    box(56, 19, 42, 22, "analysis + validation", "null models (4 degradations)\nselectivity landscape\nrank validation  |  mixtures\nHill fit + bootstrap  |  predictions", "#ffffff")
    box(56, 1, 42, 15, "notebook JSON + provenance", "readouts, gains, warnings, live_lab\nversion / git sha / library sha256", SEQ[0])

    arrow(50, 88, 50, 83.5)
    arrow(38, 68, 26, 62.5, colour=INSECT)
    arrow(62, 68, 74, 62.5, colour=VERTEBRATE)
    arrow(23, 46, 23, 42.5, "gain patch", INSECT)
    arrow(23, 23, 23, 19.5, "netlist", "#1baf7a")
    arrow(44, 11.5, 56, 24, "readouts")
    arrow(77, 46, 77, 41.5, "same dose", VERTEBRATE)
    arrow(77, 19, 77, 16.5)
    ax.text(
        50,
        -7.5,
        "Netlist and patch stay in separate files: a reviewer can tell which number came from Janelia and which came from a paper EC50.",
        ha="center",
        fontsize=7.6,
        style="italic",
        color=TEXT2,
    )
    save_fig(
        ctx,
        fig,
        "F1_architecture",
        "FlyLab dataflow. A compound and a free concentration enter the occupancy "
        "engine, which reads the sourced library (optionally shifted by a genotype, "
        "combined with a mixture partner, or driven by an exposure profile). Occupancy "
        "splits into two scorecards. The insect side passes through the mechanism table "
        "to a gain patch on transmitter classes of a named MaleCNS cut; the vertebrate "
        "side is scored at the same dose and never touches a circuit. The circuit "
        "runtime (rate or LIF) produces named-cell readouts, which the analysis and "
        "validation layers turn into empirical permutation probabilities with their "
        "three-way verdicts, selectivity indices, rank correlations and prospective "
        "predictions. Everything lands in one notebook JSON with a provenance block.",
        svg=True,
    )


#: templates rendered by :func:`step_paper`, and the value key holding each
#: rendered document's body word count
PAPER_TEMPLATES: tuple[tuple[str, str], ...] = (
    ("IJRC_FlyLab_draft.md", "paper_words_body"),
    ("SUPPLEMENT.md", "supplement_words"),
)

#: body word count = everything before this heading, minus tables, code blocks,
#: figure/table captions and the front matter, so the number means what a venue
#: means by it
_BODY_END = "## References"


def _body_words(rendered: str) -> int:
    text = rendered.split(_BODY_END)[0]
    text = re.sub(r"```.*?```", " ", text, flags=re.S)  # code blocks
    text = re.sub(r"^\|.*$", " ", text, flags=re.M)  # tables
    text = re.sub(r"^>.*$", " ", text, flags=re.M)  # block quotes (standfirst)
    return len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'’\-]*", text))


def step_paper(ctx: Ctx) -> None:
    """Render the manuscript and its supplement from their templates.

    Every ``{{key}}`` is substituted from ``results.json``: the prose is written
    around the numbers, and a number that the pipeline did not produce renders
    as ``[[MISSING:key]]`` rather than as a plausible value.
    """
    word_re = re.compile(r"[A-Za-z0-9][A-Za-z0-9'’\-]*")
    missing: list[str] = []

    def sub(m: re.Match[str]) -> str:
        key = m.group(1).strip()
        row = ctx.values.get(key)
        if row is None:
            missing.append(key)
            return f"[[MISSING:{key}]]"
        return str(row["text"])

    for name, words_key in PAPER_TEMPLATES:
        template = ctx.outdir / f"{name}.in"
        if not template.exists():
            ctx.say(f"    skip: no template at {template}")
            continue
        # the template's leading HTML comment is authoring guidance, not manuscript
        text = re.sub(r"\A<!--.*?-->\s*", "", template.read_text(), flags=re.S)
        rendered = re.sub(r"\{\{([A-Za-z0-9_]+)\}\}", sub, text)
        out = ctx.outdir / name
        out.write_text(rendered)
        body = _body_words(rendered)
        ctx.put(words_key, body, unit="words")
        if name == "IJRC_FlyLab_draft.md":
            ctx.put("paper_words", len(word_re.findall(rendered)))
            ctx.put("paper_words_in_target", bool(5000 <= body <= 7000))
            if not 5000 <= body <= 7000:
                ctx.note(f"body word count {body} is outside the 5000-7000 target")
        ctx.say(f"    paper {out.name}: {body} body words")

    ctx.put("paper_missing_keys", sorted(set(missing)), text=", ".join(sorted(set(missing))) or "none")
    if missing:
        ctx.note(f"{len(set(missing))} template keys had no value: {', '.join(sorted(set(missing))[:8])}")


# --------------------------------------------------------------------------
# step registry
# --------------------------------------------------------------------------
Step = tuple[str, str, Callable[[Ctx], None]]

STEPS: list[Step] = [
    ("provenance", "versions, hashes, map ids, graph sizes", step_provenance),
    ("census", "MaleCNS census + labellar GRN type proof (T0)", step_census),
    ("scorecard", "insect vs vertebrate panel + library table (F2, T1)", step_scorecard),
    ("evidence", "typed evidence census, schema v3 (F11, T10)", step_evidence),
    ("curves", "engagement curves and selectivity windows (F3)", step_curves),
    ("mechanisms", "mechanism -> gain rules (T2)", step_mechanisms),
    ("graph", "taste-motor graph and GRN->MN9 paths (F4, T7)", step_graph),
    ("veto", "map-based bitter veto, rate vs LIF (F5)", step_veto),
    ("dependence", "connectome dependence + landscape (F6, F12, T6, T11, T12)", step_dependence),
    ("ablation", "ablation ladder over the library (F13, T13, T13b)", step_ablation),
    ("landscape", "receptor vs circuit selectivity (F7, T8)", step_landscape),
    ("dose", "model dose-response + sensitivity tornado (F8)", step_dose),
    ("genotype", "resistance-allele panels (F9, T9)", step_genotype),
    ("validation", "rank comparison vs published orders (F10, T4)", step_validation),
    ("predictions", "prospective predictions (T3)", step_predictions),
    ("mixtures", "binary mixtures vs published verdicts (T5)", step_mixtures),
    ("expression", "expression-weighted sensitivity coverage", step_expression),
    ("stability", "specification robustness + threshold grid (F14, T14, T14b, T15)", step_stability),
    ("uncertainty", "global uncertainty budget + VOI (F15, T16, T17)", step_uncertainty),
    ("claims", "claim provenance chain (T18)", step_claims),
    ("scale", "where the dependence verdict settles with the cut (T26, T27 when the study exists)", step_scale),
    ("architecture", "architecture schematic (F1)", step_architecture),
    ("paper", "render the draft and supplement from their templates", step_paper),
]

#: steps whose compute is started in a worker process when the run begins
BACKGROUND_STEPS: dict[str, tuple[str, Callable[..., dict[str, Any]]]] = {
    "uncertainty": ("sobol", _bg_sobol),
    "stability": ("stability", _bg_stability),
}
STEP_NAMES = [s[0] for s in STEPS]

#: Every key the manuscript or the supplement substitutes.  A test asserts
#: each one exists in the committed record after a full run, and another that
#: none has rotted into a wish list; regenerate this tuple from the two
#: templates whenever their prose changes.
PAPER_KEYS: tuple[str, ...] = (
    "abl_composition_survives_normalisations",
    "abl_compounds",
    "abl_concs",
    "abl_generic_floor_note",
    "abl_generic_floor_rule",
    "abl_generic_rule",
    "abl_glutamate_max_delta_rho",
    "abl_ref_matched_median",
    "abl_ref_matched_p05",
    "abl_ref_matched_p95",
    "abl_ref_n_draws",
    "abl_ref_no_floor",
    "abl_ref_shuffled",
    "abl_rho_composition_degree",
    "abl_rho_composition_degree_1e8",
    "abl_rho_composition_none",
    "abl_rho_composition_paper",
    "abl_rho_composition_row_abs",
    "abl_rho_receptor_max",
    "abl_rho_receptor_min",
    "abl_rho_topology_floor_paper",
    "abl_rho_topology_paper",
    "abl_worst_level_by_conc",
    "abl_worst_level_by_conc_floor_rule",
    "bal_matched_max_deviation",
    "bal_n",
    "bal_plain_ach_max",
    "bal_plain_ach_mean",
    "bal_plain_ach_min",
    "bal_plain_ach_percentile",
    "bal_plain_ach_sd",
    "bal_real_ach_share",
    "bal_tol",
    "claims_chain_links",
    "claims_computed",
    "claims_facts_phrase",
    "claims_inference_phrase",
    "claims_literature_derived",
    "claims_model_assumption",
    "claims_observed_phrase",
    "claims_unknown_phrase",
    "dep_land_below_relative_floor",
    "dep_land_below_relative_floor_cells",
    "dep_land_cells",
    "dep_land_composition_dominated",
    "dep_land_composition_dominated_absfloor",
    "dep_land_composition_dominated_raw",
    "dep_land_compounds",
    "dep_land_concs",
    "dep_land_fdr_can_reject",
    "dep_land_fipronil_class",
    "dep_land_fipronil_class_raw",
    "dep_land_fipronil_q",
    "dep_land_mixed",
    "dep_land_n_structural_tests",
    "dep_land_no_effect",
    "dep_land_no_effect_absfloor",
    "dep_land_no_effect_raw",
    "dep_land_non_monotone",
    "dep_land_non_monotone_cells",
    "dep_land_taste_cells",
    "dep_land_taste_composition_dominated",
    "dep_land_taste_imidacloprid_class",
    "dep_land_taste_imidacloprid_p_rewire_degree_preserving",
    "dep_land_taste_imidacloprid_q",
    "dep_land_taste_n_structural_tests",
    "dep_land_taste_no_effect",
    "dep_land_taste_non_monotone",
    "dep_land_taste_topology_dependent",
    "dep_land_topology_dependent",
    "dep_land_topology_dependent_absfloor",
    "dep_land_topology_dependent_raw",
    "dep_land_verdicts_distinguishable",
    "dep_land_verdicts_equivalent_within_tolerance",
    "dep_land_verdicts_indeterminate",
    "dep_n",
    "dep_n_landscape",
    "dep_n_landscape_taste",
    "dep_n_taste",
    "dep_named_fipronil_class",
    "dep_named_fipronil_effect",
    "dep_named_fipronil_level",
    "dep_named_fipronil_min_structural_p",
    "dep_named_fipronil_p_rewire_degree_preserving",
    "dep_named_fipronil_p_weight_permute",
    "dep_named_imidacloprid_class",
    "dep_named_imidacloprid_delta",
    "dep_named_imidacloprid_delta_frac",
    "dep_named_imidacloprid_effect",
    "dep_named_imidacloprid_equivalent_modes",
    "dep_named_imidacloprid_gap_sign_permute",
    "dep_named_imidacloprid_gap_sign_permute_weight_matched",
    "dep_named_imidacloprid_indeterminate_modes",
    "dep_named_imidacloprid_level",
    "dep_named_imidacloprid_level_verdict",
    "dep_named_imidacloprid_p_erdos_renyi",
    "dep_named_imidacloprid_p_rewire_degree_preserving",
    "dep_named_imidacloprid_p_sign_permute",
    "dep_named_imidacloprid_p_sign_permute_weight_matched",
    "dep_named_imidacloprid_p_weight_permute",
    "dep_named_imidacloprid_verdict_erdos_renyi",
    "dep_named_imidacloprid_verdict_rewire_degree_preserving",
    "dep_named_imidacloprid_verdict_sign_permute",
    "dep_named_imidacloprid_verdict_sign_permute_weight_matched",
    "dep_named_imidacloprid_verdict_weight_permute",
    "dep_p_floor",
    "dep_p_within_0p02_from_n",
    "dep_p_within_0p02_worst_n",
    "dep_taste_fipronil_equivalent_modes",
    "dep_taste_fipronil_indeterminate_modes",
    "dep_taste_fipronil_max_abs_z",
    "dep_verdict_stable_from_n",
    "dep_verdict_stable_to_n",
    "dia_insect_rdl_occ",
    "ev_binding_proxy_row_names",
    "ev_binding_row_names",
    "ev_compounds",
    "ev_dist_E0",
    "ev_dist_E1",
    "ev_dist_E2",
    "ev_dist_E3",
    "ev_dist_E4",
    "ev_distance_labels",
    "ev_library_problems",
    "ev_model_binding_engagement_proxy",
    "ev_model_binding_occupancy",
    "ev_model_functional_engagement",
    "ev_model_functional_engagement_proxy",
    "ev_model_not_modelled",
    "ev_param_EC50",
    "ev_param_IC50",
    "ev_param_Kd",
    "ev_param_unknown",
    "ev_receptor_keys",
    "ev_rows",
    "ev_rows_not_modelled",
    "expr_coverage_overall",
    "expr_mn9_delta_hz",
    "fip_ec50_ratio",
    "fip_insect_occ",
    "fip_vert_occ",
    "flylab_version",
    "geno_ddt_para_M918T_superkdr_fold",
    "geno_deltamethrin_para_M918T_superkdr_fold",
    "geno_fipronil_rdl_A301S_nlug_fold",
    "geno_fipronil_rdl_A302G_dsim_fold",
    "geno_fipronil_rdl_A302G_dsim_occ",
    "geno_gaba_rdl_A301S_nlug_fold",
    "git_sha",
    "ic50_M",
    "ic50_hi_M",
    "ic50_lo_M",
    "ic50_n_boot",
    "ic50_r2",
    "imi_ec50_ratio",
    "imi_insect_occ",
    "imi_vert_occ",
    "ish_analytic",
    "ish_estimated",
    "ish_max_abs_error",
    "ish_n_base",
    "landscape_amplify_compounds",
    "landscape_buffer_compounds",
    "landscape_unscored_compounds",
    "library_sha256",
    "map_id",
    "mechanism_rules",
    "mix_n_matching",
    "mix_n_substituted",
    "mix_n_targets",
    "mix_species",
    "named_cells_acetylcholine",
    "named_cells_gaba",
    "named_cells_glutamate",
    "named_cells_octopamine",
    "named_edges",
    "named_mean_degree",
    "named_n_seed_nodes",
    "named_nodes",
    "named_nodes_with_input",
    "named_share_onto_seeds",
    "named_share_recurrent_edges",
    "named_share_unclear_cells",
    "norm_note_degree",
    "norm_note_none",
    "norm_note_row_abs",
    "pred_H1_effect",
    "pred_H2_effect",
    "pred_H3_effect",
    "pred_H4_effect",
    "pred_H5_effect",
    "pred_H6_effect",
    "pred_H7_effect",
    "pred_direction_mismatches",
    "pred_max_n",
    "pred_n_hypotheses",
    "pred_n_rep",
    "rank_binomial_p",
    "rank_concordant",
    "rank_discordant",
    "rank_discrepancy_ids",
    "rank_entries_evaluable",
    "rank_evaluated",
    "rank_known_discrepancies",
    "rank_mean_rho_informative",
    "rank_n_degenerate",
    "rank_n_informative",
    "rank_shared_source",
    "rank_shared_source_be",
    "rank_shared_source_verb",
    "rank_skipped",
    "rank_source_disjoint",
    "rank_source_disjoint_be",
    "scale_1k_edges",
    "scale_1k_in_star",
    "scale_1k_mean_degree",
    "scale_1k_nodes",
    "scale_frontier_n_paper",
    "scale_ladder_sizes",
    "scale_study_settled_clean_ladder",
    "scale_study_statement",
    "scale_top_landscape_days",
    "scale_top_profile_hours",
    "scale_top_rung",
    "scale_whole_cns_profile_hours",
    "sign_asserted",
    "sign_asserted_transmitters",
    "stab_C1_claim",
    "stab_C1_failing",
    "stab_C1_lost",
    "stab_C1_retained",
    "stab_C1_undecidable",
    "stab_C2_claim",
    "stab_C2_lost",
    "stab_C2_retained",
    "stab_C2_undecidable",
    "stab_C3_claim",
    "stab_C3_equivalent",
    "stab_C3_indeterminate",
    "stab_C3_lost",
    "stab_C3_retained",
    "stab_C3_retained_uncorrected",
    "stab_C3_undecidable",
    "stab_C4_claim",
    "stab_C4_equivalent",
    "stab_C4_indeterminate",
    "stab_C4_lost",
    "stab_C4_retained",
    "stab_C4_retained_uncorrected",
    "stab_C4_undecidable",
    "stab_C5_claim",
    "stab_C5_lost",
    "stab_C5_retained",
    "stab_C5_undecidable",
    "stab_C6_claim",
    "stab_C6_lost",
    "stab_C6_retained",
    "stab_C6_undecidable",
    "stab_C7_claim",
    "stab_C7_lost",
    "stab_C7_retained",
    "stab_C7_undecidable",
    "stab_default_g_ach",
    "stab_default_spec",
    "stab_default_treated_hz",
    "stab_default_vehicle_hz",
    "stab_family_size",
    "stab_mean_gap_nicotinic_default",
    "stab_mean_gap_nicotinic_monotone",
    "stab_mean_gap_nicotinic_monotone_max",
    "stab_mean_gap_nicotinic_monotone_min",
    "stab_mean_gap_nicotinic_monotone_n",
    "stab_monotone_buffer_less",
    "stab_monotone_buffer_less_n",
    "stab_monotone_g_ach",
    "stab_monotone_pct_change",
    "stab_monotone_treated_hz",
    "stab_monotone_vehicle_hz",
    "stab_n_shuffles",
    "stab_n_structural_tests",
    "stab_shuffle_resolution",
    "stab_specification_independent",
    "stab_topology_alpha",
    "taste_edges",
    "taste_motor_cells_octopamine",
    "taste_motor_mean_degree",
    "taste_motor_share_onto_seeds",
    "taste_motor_share_recurrent_edges",
    "taste_nodes",
    "thr_circuit_fracs",
    "thr_n_compounds",
    "thr_nav_ache_cells",
    "thr_nav_ache_failing",
    "thr_nav_ache_gap_min",
    "thr_nav_ache_matching",
    "thr_nicotinic_cells",
    "thr_nicotinic_gap_max",
    "thr_nicotinic_gap_min",
    "thr_nicotinic_matching",
    "thr_outside_t15_scope",
    "thr_stable",
    "thr_vert_limits",
    "tornado_zero_span_params",
    "unc_CI_drive",
    "unc_ST_gain_transform",
    "unc_ST_weight_threshold",
    "unc_S_drive",
    "unc_S_gain_transform",
    "unc_S_weight_threshold",
    "unc_at_noise_floor",
    "unc_evaluations",
    "unc_interaction_share",
    "unc_interaction_share_clipped",
    "unc_n_base",
    "unc_n_factors",
    "unc_n_resolved",
    "unc_noise_floor",
    "unc_noise_floor_factor",
    "unc_null_factor_S",
    "unc_resolved",
    "unc_variance",
    "val_power_by_n",
    "val_power_by_strength",
    "val_power_control_hits",
    "val_power_control_runs",
    "val_power_fpr",
    "val_power_full_strengths",
    "val_power_ns",
    "val_power_replicates",
    "val_power_strengths",
    "val_power_weak_strengths",
    "val_recovery_detected",
    "val_recovery_edges",
    "val_recovery_n",
    "val_recovery_nodes",
    "val_recovery_positive",
    "val_recovery_statement",
    "veto_lif_fipronil_ratio",
    "veto_lif_vehicle_ratio",
    "veto_lif_vehicle_sugar",
    "veto_n_sweet_grn",
    "veto_rate_fipronil_ratio",
    "veto_rate_vehicle_ratio",
    "veto_rate_vehicle_sugar",
    "veto_sugar_hz",
    "voi_rank1_experiment",
    "voi_rank1_factor",
    "voi_rank1_var",
    "voi_rank2_factor",
    "voi_rank2_var",
    "voi_var_total",
    "window_log10_imidacloprid",
)


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
def _summary_table(ctx: Ctx, total_s: float) -> str:
    w = max(len(n) for n in ctx.steps) if ctx.steps else 10
    lines = [
        "",
        "=" * (w + 58),
        f"{'step'.ljust(w)}  {'status':8}  {'sec':>7}  {'figs':>4}  {'tabs':>4}  {'values':>6}",
        "-" * (w + 58),
    ]
    for name, s in ctx.steps.items():
        lines.append(
            f"{name.ljust(w)}  {s['status']:8}  {s['runtime_s']:7.1f}  "
            f"{s['figures']:4d}  {s['tables']:4d}  {s['values']:6d}"
        )
    lines.append("-" * (w + 58))
    lines.append(
        f"{'TOTAL'.ljust(w)}  {'':8}  {total_s:7.1f}  "
        f"{len(ctx.figures):4d}  {len(ctx.tables):4d}  {len(ctx.values):6d}"
    )
    lines.append("=" * (w + 58))
    if ctx.notes:
        lines.append("Findings recorded by this run (also in results.json.notes):")
        for n in ctx.notes:
            lines.append(f"  ! {n}")
    lines.append(f"mode: {'FAST (reduced replicates)' if ctx.fast else 'default'}   outdir: {ctx.outdir}")
    return "\n".join(lines)


def _write_results(ctx: Ctx, total_s: float) -> Path:
    prov = ctx.values.pop("_provenance", None)
    payload = {
        "schema": RESULTS_SCHEMA,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "fast": ctx.fast,
        "seed": ctx.seed,
        "jobs": ctx.jobs,
        "runtime_s": round(total_s, 2),
        "provenance": (prov or {}).get("value"),
        "steps": ctx.steps,
        "values": ctx.values,
        "figures": ctx.figures,
        "tables": ctx.tables,
        "notes": ctx.notes,
        "partial": ctx.partial,
        "honesty": [
            "Every number in this file is model_derived: it comes from the teaching "
            "EC50 library, the gain patch rules and one or two hops-limited MaleCNS "
            "cuts. Nothing here was measured in a living fly.",
            "live_lab is null in every notebook this pipeline produced.",
            "--fast reduces replicate and shuffle counts only; z-scores and CIs from a "
            "fast run are noisier than the paper's.",
        ],
    }
    if prov is not None:
        ctx.values["_provenance"] = prov
    out = ctx.outdir / "results.json"
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    return out


def _write_captions(ctx: Ctx) -> Path:
    lines = ["# Figure captions", "", "Generated by `scripts/reproduce_paper.py`. Do not edit by hand.", ""]
    seen: set[str] = set()
    for entry in sorted(ctx.figures, key=lambda e: str(e.get("figure") or e["name"])):
        name, caption = entry.get("figure"), entry.get("caption")
        if not name or not caption or name in seen:
            continue
        seen.add(name)
        lines.append(f"**{name}.** {caption}")
        lines.append("")
    p = ctx.figdir / "captions.md"
    ctx.ensure_dirs()
    p.write_text("\n".join(lines))
    return p


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="reproduce_paper",
        description="Regenerate every number, figure and table in the FlyLab IJRC paper.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="steps: " + ", ".join(STEP_NAMES),
    )
    ap.add_argument("--outdir", type=Path, default=REPO_ROOT / "papers", help="output directory (default: papers/)")
    ap.add_argument("--fast", action="store_true", help="fewer replicates/shuffles (see module docstring)")
    ap.add_argument(
        "--only",
        action="append",
        default=None,
        metavar="STEP",
        help="run only this step (repeatable); implies the steps it depends on are skipped",
    )
    ap.add_argument("--list", action="store_true", help="list the steps and exit")
    ap.add_argument("--seed", type=int, default=0, help="global RNG seed (default 0)")
    ap.add_argument(
        "--jobs",
        type=int,
        default=0,
        help=(
            "worker processes for the permutation analyses and for the two backgrounded "
            "steps (0 = choose from the core count). Results are identical at any value; "
            "only the wall clock changes."
        ),
    )
    ap.add_argument("--quiet", action="store_true", help="only print the summary table")
    return ap


GUIDANCE = """\
FlyLab paper reproduction is driven by the script, not by this shim, because it
takes options (and writes into papers/). Run one of:

  python scripts/reproduce_paper.py            full run, ~8-10 min, no downloads
  python scripts/reproduce_paper.py --fast     ~2 min, fewer permutations/replicates
  python scripts/reproduce_paper.py --only dependence --outdir /tmp/x
  python scripts/reproduce_paper.py --list     the named steps

Outputs: papers/figures/*.png(+svg), papers/tables/T*.csv|md, papers/results.json,
and papers/IJRC_FlyLab_draft.md + papers/SUPPLEMENT.md re-rendered from their
templates. --fast is a sanity check, not the paper: it changes statistical effort
only, and a test refuses a committed results.json that came from one.
"""


def main(argv: list[str] | None = None) -> int:
    """Run the pipeline.

    ``argv`` is the command line to parse. ``None`` means *no command line was
    supplied at all* — the ``flylab reproduce-paper`` shim calls this as a
    library function — and prints :data:`GUIDANCE` instead of running anything.
    Two reasons: inheriting the host process's ``sys.argv`` would make a test
    runner's own flags look like pipeline options, and a library call that
    silently spends minutes writing into ``papers/`` is a bad neighbour. Pass
    ``[]`` for an explicit default run; the module entry point passes
    ``sys.argv[1:]``.
    """
    if argv is None:
        print(GUIDANCE, end="")
        return 0
    raw = list(argv)
    # tolerate the subcommand name if a caller forwards it verbatim
    if raw and raw[0] == "reproduce-paper":
        raw = raw[1:]
    args = build_parser().parse_args(raw)

    if args.list:
        for name, desc, _ in STEPS:
            print(f"{name:14} {desc}")
        return 0

    only = set(args.only or [])
    unknown = only - set(STEP_NAMES)
    if unknown:
        print(f"unknown step(s): {', '.join(sorted(unknown))}\nknown: {', '.join(STEP_NAMES)}", file=sys.stderr)
        return 2

    cores = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    jobs = int(args.jobs) if int(args.jobs) > 0 else max(1, min(4, cores))
    ctx = Ctx(
        outdir=Path(args.outdir).resolve(),
        fast=bool(args.fast),
        seed=int(args.seed),
        quiet=bool(args.quiet),
        jobs=jobs,
    )
    ctx.ensure_dirs()

    # provenance is cheap and results.json is meaningless without it, so it
    # runs even when --only names other steps.
    if only and only != {"paper"}:
        only.add("provenance")

    # --only is incremental: keep what an earlier run in this outdir produced
    # for the steps we are NOT re-running, so a partial run cannot silently
    # truncate results.json (or captions.md) to the steps it happened to touch.
    merged = False
    prev_path = ctx.outdir / "results.json"
    if only and prev_path.exists():
        try:
            data = json.loads(prev_path.read_text())
        except (OSError, ValueError):
            data = None
        if isinstance(data, dict):
            merged = True
            ctx.values = {
                k: v
                for k, v in (data.get("values") or {}).items()
                if isinstance(v, dict) and v.get("step") not in only
            }
            ctx.figures = [f for f in (data.get("figures") or []) if f.get("step") not in only]
            ctx.tables = [t for t in (data.get("tables") or []) if t.get("step") not in only]
            ctx.steps = {k: v for k, v in (data.get("steps") or {}).items() if k not in only}
            ctx.notes = [
                n for n in (data.get("notes") or [])
                if not any(str(n).startswith(f"[{step}]") for step in only)
            ]
            if data.get("provenance") and "provenance" not in only:
                ctx.values["_provenance"] = {
                    "value": data["provenance"],
                    "text": "",
                    "unit": None,
                    "note": "carried over from the previous run in this outdir",
                    "step": "provenance",
                }

    ctx.say(
        f"FlyLab paper reproduction  |  mode {'fast' if ctx.fast else 'default'}  |  "
        f"seed {ctx.seed}  |  jobs {ctx.jobs}  |  outdir {ctx.outdir}"
        + ("  |  merging into the previous results.json" if merged else "")
    )
    t_all = time.perf_counter()

    # start the two long analyses now; they are joined at their own step
    pool = None
    wanted = [s for s in BACKGROUND_STEPS if (not only or s in only)]
    if wanted and ctx.jobs > 1:
        try:
            import multiprocessing as _mp
            from concurrent.futures import ProcessPoolExecutor

            pool = ProcessPoolExecutor(max_workers=len(wanted), mp_context=_mp.get_context("fork"))
            for step_name in wanted:
                key, fn = BACKGROUND_STEPS[step_name]
                ctx.background[key] = pool.submit(fn, ctx.fast, ctx.seed)
            ctx.say(f"  background: {', '.join(k for k, _ in (BACKGROUND_STEPS[s] for s in wanted))}")
        except Exception as exc:  # no fork, no /dev/shm, restricted sandbox...
            ctx.background.clear()
            pool = None
            ctx.say(f"  background pool unavailable ({type(exc).__name__}); running everything in process")

    failures = 0
    for name, desc, fn in STEPS:
        if only and name not in only:
            continue
        ctx._step = name
        n_fig, n_tab, n_val = len(ctx.figures), len(ctx.tables), len(ctx.values)
        ctx.say(f"  [{name}] {desc}")
        t0 = time.perf_counter()
        status = "ok"
        try:
            fn(ctx)
        except Exception as exc:  # a broken step must not lose the other steps
            failures += 1
            status = "FAILED"
            ctx.notes.append(f"[{name}] step failed: {type(exc).__name__}: {exc}")
            ctx.say(f"    ! FAILED: {type(exc).__name__}: {exc}")
        ctx.steps[name] = {
            "status": status,
            "runtime_s": round(time.perf_counter() - t0, 2),
            "description": desc,
            "figures": len(ctx.figures) - n_fig,
            "tables": len(ctx.tables) - n_tab,
            "values": len(ctx.values) - n_val,
        }

    if pool is not None:
        for fut in ctx.background.values():
            fut.cancel()
        pool.shutdown(wait=False, cancel_futures=True)

    if ctx.figures:
        _write_captions(ctx)
    # "partial" means the record is missing analysis; re-rendering the prose
    # from an unchanged record (--only paper) does not make it partial.
    ctx.partial = bool(only - {"paper", "provenance"})
    total = time.perf_counter() - t_all
    if only and not ctx.partial:
        # A pure re-render must not overwrite the analysis runtime with its own.
        # The per-step runtimes survive the merge, so their sum is the cost of
        # the analysis this record describes, whoever last re-rendered the prose.
        total = sum(float(st.get("runtime_s") or 0.0) for st in ctx.steps.values()) or total
    results = _write_results(ctx, total)
    print(_summary_table(ctx, total))
    print(f"results -> {results}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
