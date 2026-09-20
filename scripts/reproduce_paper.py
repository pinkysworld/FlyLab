#!/usr/bin/env python3
"""Regenerate every number, figure and table in the FlyLab IJRC paper.

One command, committed data only::

    python scripts/reproduce_paper.py            # full run  (~4-6 min)
    python scripts/reproduce_paper.py --fast     # reduced   (~1-2 min)
    python scripts/reproduce_paper.py --only nulls
    flylab reproduce-paper                       # same thing through the CLI

Nothing here downloads anything and nothing needs the 1.1 GB MaleCNS weight
matrix: every input is either ``data/derived/`` (two committed neighbourhood
cuts, the census fallback and the gustatory seed list), ``data/literature/``
(the six sourced YAML datasets) or ``flylab/pharm/library.yaml``.

Outputs (all under ``--outdir``, default ``papers/``):

===========================  ==================================================
``figures/F*.png`` / ``.svg``  300 dpi figures, colour-blind-safe, UI palette
``figures/captions.md``        one caption per figure
``tables/T*.csv`` / ``.md``    the six paper tables
``results.json``               every scalar quoted in the paper, keyed, with a
                               provenance block (flylab version, git sha,
                               library sha256, map ids, seeds)
``IJRC_FlyLab_draft.md``       rendered from ``IJRC_FlyLab_draft.md.in`` by
                               substituting ``{{key}}`` from ``results.json``
===========================  ==================================================

``--only`` is **incremental**: when ``<outdir>/results.json` already exists the
run keeps every value, figure, table and note produced by the steps it is not
re-running, so a partial run cannot silently truncate the record. The document
then carries ``"partial": true``. ``--only paper`` therefore re-renders the
draft from whatever the last full run stored.

``--fast`` changes only the *statistical* effort, never the model:

==================  ===============  ===============
knob                default          ``--fast``
==================  ===============  ===============
null shuffles       50 (subgraph)    10
                    10 (taste map)   3
selectivity ladder  13 half-decades  7 decades
IC50 bootstrap      200 resamples    40
IC50 replicates     4                2
prediction repl.    6                2
==================  ===============  ===============

A ``--fast`` z-score or CI is therefore noisier than the paper's; the point
estimates that do not depend on replicate count are identical.

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
        """
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
    ctx.put("dia_insect_rdl_occ", float(dia_rdl["occupancy"]))
    ctx.put("dia_vert_gabaa_occ", float(dia_gaba["occupancy"]))
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
                    "ec50_M": spec["ec50_M"],
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
            "ec50_M",
            "hill_n",
            "direction",
            "efficacy",
            "evidence_tier",
            "source",
        ],
        "The sourced compound library. Every row carries an evidence tier and a "
        "named source; `class_placeholder` rows carry no number (EC50 0.01 M, "
        "direction `none`) and are inert by construction.",
        md_fields=["compound", "class", "receptor", "ec50_M", "direction", "evidence_tier"],
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
                    "no sourced EC50",
                    xy=(i, 0.5),
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=TEXT2,
                    zorder=3,
                )
            ax.bar(i - 0.2, ins[i], 0.38, color=INSECT, hatch=hatch, edgecolor="white", linewidth=0.6, zorder=2)
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
    axes[0].set_ylabel("fractional receptor occupancy (0-1)")
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
    """T2: the mechanism -> gain patch rules, verbatim from the code."""
    from flylab.pharm.mechanisms import GAIN_KEYS, mechanism_table_rows

    rows = mechanism_table_rows()
    ctx.put("mechanism_rules", len(rows))
    ctx.put("gain_keys", list(GAIN_KEYS), text=", ".join(GAIN_KEYS))
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
    """F4: the MaleCNS taste-motor graph, GRN seeds, and GRN -> MN9 paths."""
    import numpy as np

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


def _f(x: Any) -> float | None:
    try:
        return None if x is None else float(x)
    except (TypeError, ValueError):
        return None


def step_nulls(ctx: Ctx) -> None:
    """F6 + T6: the drug effect against four connectome null models."""
    import numpy as np

    from flylab.analysis.nullmodels import (
        MODES,
        connectome_information_score,
        null_distribution,
    )

    plt = _plt()
    n = 10 if ctx.fast else 50
    n_taste = 3 if ctx.fast else 10
    ctx.put("null_n_shuffles", n)
    ctx.put("null_n_shuffles_taste", n_taste)

    # null_panel() drops the raw draws, and F6 needs them, so the arms are run
    # directly here; the row set is identical to null_panel()'s.
    arms = (("subgraph", "mean_hz", n), ("taste_map", "bitter_veto_ratio", n_taste))
    draws: dict[tuple[str, str], list[float]] = {}
    rows: list[dict[str, Any]] = []
    for compound in ("imidacloprid", "fipronil"):
        for assay, readout, n_i in arms:
            for mode in MODES:
                d = null_distribution(
                    assay, compound, PAPER_CONC, mode, n=n_i, seed=ctx.seed, readout=readout
                )
                rows.append(
                    {
                        "compound": compound,
                        "assay": assay,
                        "readout": readout,
                        "mode": mode,
                        "n": d["n"],
                        "n_ok": d["n_ok"],
                        "real_effect": d["real_effect"],
                        "null_mean": d["null_mean"],
                        "null_median": d["null_median"],
                        "null_sd": d["null_sd"],
                        "z": d["z"],
                        "p_two_sided": d["p_two_sided"],
                        "runtime_s": round(d["runtime_s"], 2),
                    }
                )
                if assay == "subgraph":
                    draws[(compound, mode)] = [v for v in d["null_effects"] if v is not None]
                    ctx.put(f"null_z_{compound}_{mode}", _f(d["z"]))
                    ctx.put(f"null_p_{compound}_{mode}", _f(d["p_two_sided"]))
                    ctx.put(f"null_real_effect_{compound}", _f(d["real_effect"]), unit="Hz")
                else:
                    ctx.put(f"nulltaste_z_{compound}_{mode}", _f(d["z"]))
        score = connectome_information_score(
            panel=[r for r in rows if r["compound"] == compound]
        )
        ctx.put(f"null_score_{compound}", _f(score["score"]))

    save_table(
        ctx,
        "T6_nullmodels",
        rows,
        [
            "compound",
            "assay",
            "readout",
            "mode",
            "n",
            "n_ok",
            "real_effect",
            "null_mean",
            "null_sd",
            "z",
            "p_two_sided",
            "runtime_s",
        ],
        "Connectome null models. `real_effect` is treated minus vehicle on the real "
        "MaleCNS cut; the null distribution is the same contrast on `n` degraded "
        "copies of that cut. `p` is the empirical two-sided permutation p, "
        "(k+1)/(n+1), so its floor is 1/(n+1).",
        md_fields=["compound", "assay", "readout", "mode", "n", "real_effect", "z", "p_two_sided"],
    )

    # honest read-out of what the panel says
    imi_top = [abs(ctx.get(f"null_z_imidacloprid_{m}") or 0.0) for m in MODES if m != "erdos_renyi"]
    if max(imi_top) < 2.0:
        ctx.note(
            "imidacloprid's neighbourhood mean-rate effect does NOT beat any "
            "structure-preserving null (max |z| = %.2f over sign_permute, "
            "weight_permute and rewire_degree_preserving); it beats only the "
            "Erdos-Renyi null. The mean-rate readout is an E/I-balance effect for "
            "nAChR agonists, not a wiring effect." % max(imi_top)
        )
    fip_top = [abs(ctx.get(f"null_z_fipronil_{m}") or 0.0) for m in ("weight_permute", "rewire_degree_preserving")]
    ctx.put("null_imi_max_structural_z", float(max(imi_top)))
    ctx.put("null_fip_max_structural_z", float(max(fip_top)))
    taste_z = [
        abs(ctx.get(f"nulltaste_z_fipronil_{m}") or 0.0) for m in MODES
    ]
    ctx.put("null_fip_taste_max_z", float(max(taste_z)))
    if max(taste_z) < 2.0:
        ctx.note(
            "the fipronil bitter-veto ratio on the taste-motor cut is NOT "
            "distinguishable from its shuffled controls either (max |z| = %.2f at "
            "n=%d shuffles); the veto ratio is reported as a model behaviour, not as "
            "evidence that the specific MaleCNS wiring carries it." % (max(taste_z), n_taste)
        )

    fig, axes = plt.subplots(2, len(MODES), figsize=(12.0, 5.4))
    for row, compound in enumerate(("imidacloprid", "fipronil")):
        sub = {
            r["mode"]: r
            for r in rows
            if r["assay"] == "subgraph" and r["compound"] == compound
        }
        for col, mode in enumerate(MODES):
            ax = axes[row, col]
            vals = np.asarray(draws.get((compound, mode), []), dtype=float)
            real = sub[mode]["real_effect"]
            if vals.size:
                ax.hist(vals, bins=min(18, max(5, vals.size // 2)), color=SEQ[1], edgecolor="white", linewidth=0.5)
            if real is not None:
                ax.axvline(float(real), color=CRIT, lw=2, label="real MaleCNS cut")
            z = sub[mode]["z"]
            ax.set_title(
                f"{mode}\nz = {'n/a' if z is None else f'{z:.2f}'}",
                fontsize=8.5,
                color=CRIT if (z is not None and abs(z) > 2) else TEXT2,
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
        f"Drug effect on the real MaleCNS neighbourhood vs four null models "
        f"(n = {n} shuffles, network mean rate, {_fmt_M(PAPER_CONC)})",
        y=1.02,
    )
    fig.tight_layout()
    save_fig(
        ctx,
        fig,
        "F6_null_models",
        "Histograms are the drug effect (treated minus vehicle on the network mean "
        f"rate) measured on {n} degraded copies of the 1-hop MN9/DNp01 cut; the red line "
        "is the same contrast on the real cut. Columns are the four degradations in "
        "increasing order of destruction: transmitter labels permuted (E/I histogram "
        "preserved), synapse weights permuted, degree-preserving double-edge rewiring, "
        "and Erdos-Renyi. Imidacloprid beats only the Erdos-Renyi null "
        f"(|z| = {abs(ctx.get('null_z_imidacloprid_erdos_renyi') or 0):.0f}) and not the "
        f"structure-preserving ones (max |z| = {ctx.get('null_imi_max_structural_z'):.2f}); "
        "fipronil beats the topology nulls "
        f"(|z| = {abs(ctx.get('null_z_fipronil_weight_permute') or 0):.2f} weight-permute, "
        f"{abs(ctx.get('null_z_fipronil_rewire_degree_preserving') or 0):.2f} "
        "degree-preserving rewiring). Read this as a negative control that worked: for "
        "an nAChR agonist the mean-rate readout is largely an E/I-balance effect.",
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
    n_boot = 40 if ctx.fast else 200
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
    ctx.put("rank_mean_rho_occupancy", float(v["assays"]["occupancy"]["mean_rho"]))
    ctx.put("rank_mean_rho_subgraph", float(v["assays"]["subgraph"]["mean_rho"]))
    ctx.put("rank_exact_occupancy", int(v["assays"]["occupancy"]["n_exact"]))
    ctx.put("rank_known_discrepancies", len(KNOWN_DISCREPANCIES))
    ctx.put(
        "rank_discrepancy_ids",
        [d["id"] for d in KNOWN_DISCREPANCIES],
        text=", ".join(d["id"] for d in KNOWN_DISCREPANCIES),
    )

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
            "n_compounds",
            "skipped",
            "spearman_rho",
            "kendall_tau",
            "exact_match",
            "receptor",
            "source",
        ],
        "Retrospective rank validation against `data/literature/published_rank_orders.yaml`. "
        "`skipped` entries are those the teaching library cannot cover (missing "
        "compounds or placeholder-only rows); they are reported, not dropped.",
        md_fields=["assay", "id", "species", "n_compounds", "skipped", "spearman_rho", "exact_match"],
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
        "teaching library can cover, with Spearman's rho per entry and the mean across "
        f"entries ({s['mean_rho']:.2f} over {s['n_evaluated']} evaluated entries; "
        f"{s['n_skipped']} entry-assay combinations are skipped because the library "
        "cannot cover them). Right: the two inversions the validation module records "
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
    """T3: the pre-registered predictions with effect sizes and suggested n."""
    from flylab.analysis.predictions import prediction_table

    n_rep = 2 if ctx.fast else 6
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
        ctx.note("pre-registered direction not reproduced by the model for: " + ", ".join(mismatched))
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
        "Pre-registered predictions H1-H7. `predicted_effect` and its CI are "
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
        "validation layers turn into null-model z-scores, selectivity indices, rank "
        "correlations and pre-registered predictions. Everything lands in one notebook "
        "JSON with a provenance block.",
        svg=True,
    )


def step_paper(ctx: Ctx) -> None:
    """Render papers/IJRC_FlyLab_draft.md from its template and results.json."""
    template = ctx.outdir / "IJRC_FlyLab_draft.md.in"
    if not template.exists():
        ctx.say(f"    skip: no template at {template}")
        return
    text = template.read_text()
    missing: list[str] = []

    def sub(m: re.Match[str]) -> str:
        key = m.group(1).strip()
        row = ctx.values.get(key)
        if row is None:
            missing.append(key)
            return f"[[MISSING:{key}]]"
        return str(row["text"])

    # the template's leading HTML comment is authoring guidance, not manuscript
    text = re.sub(r"\A<!--.*?-->\s*", "", text, flags=re.S)
    rendered = re.sub(r"\{\{([A-Za-z0-9_]+)\}\}", sub, text)
    out = ctx.outdir / "IJRC_FlyLab_draft.md"
    out.write_text(rendered)
    word_re = re.compile(r"[A-Za-z0-9][A-Za-z0-9'\-]*")
    words = len(word_re.findall(rendered))
    ctx.put("paper_words", words)
    ctx.put("paper_words_body", len(word_re.findall(rendered.split("## References")[0])))
    ctx.put("paper_missing_keys", sorted(set(missing)), text=", ".join(sorted(set(missing))) or "none")
    if missing:
        ctx.note(f"{len(set(missing))} template keys had no value: {', '.join(sorted(set(missing))[:8])}")
    ctx.say(f"    paper {out.name}: {words} words, {len(set(missing))} missing keys")


# --------------------------------------------------------------------------
# step registry
# --------------------------------------------------------------------------
Step = tuple[str, str, Callable[[Ctx], None]]

STEPS: list[Step] = [
    ("provenance", "versions, hashes, map ids, graph sizes", step_provenance),
    ("census", "MaleCNS census + labellar GRN type proof (T0)", step_census),
    ("scorecard", "dual scorecard + library table (F2, T1)", step_scorecard),
    ("curves", "occupancy curves and selectivity windows (F3)", step_curves),
    ("mechanisms", "mechanism -> gain rules (T2)", step_mechanisms),
    ("graph", "taste-motor graph and GRN->MN9 paths (F4, T7)", step_graph),
    ("veto", "map-based bitter veto, rate vs LIF (F5)", step_veto),
    ("nulls", "connectome null models (F6, T6)", step_nulls),
    ("landscape", "receptor vs circuit selectivity (F7, T8)", step_landscape),
    ("dose", "model IC50 + sensitivity tornado (F8)", step_dose),
    ("genotype", "resistance-allele panels (F9, T9)", step_genotype),
    ("validation", "rank validation vs published orders (F10, T4)", step_validation),
    ("predictions", "pre-registered predictions (T3)", step_predictions),
    ("mixtures", "binary mixtures vs published verdicts (T5)", step_mixtures),
    ("expression", "expression-weighted sensitivity coverage", step_expression),
    ("architecture", "architecture schematic (F1)", step_architecture),
    ("paper", "render the draft from its template", step_paper),
]
STEP_NAMES = [s[0] for s in STEPS]

#: keys the paper template is allowed to reference; tests assert they all exist
#: after a full run.  Keep in sync with papers/IJRC_FlyLab_draft.md.in.
PAPER_KEYS: tuple[str, ...] = (
    "flylab_version",
    "git_sha",
    "library_sha256",
    "map_id",
    "named_nodes",
    "named_edges",
    "taste_nodes",
    "taste_edges",
    "census_traced",
    "census_ach",
    "census_gaba",
    "census_glu",
    "census_gustatory",
    "lb_type_count",
    "sweet_grn_cells",
    "bitter_grn_cells",
    "sweet_grn_breakdown",
    "bitter_grn_breakdown",
    "imi_insect_occ",
    "imi_vert_occ",
    "imi_ec50_ratio",
    "imi_insect_ec50",
    "imi_vert_ec50",
    "fip_insect_occ",
    "fip_vert_occ",
    "fip_ec50_ratio",
    "dia_insect_rdl_occ",
    "dia_vert_gabaa_occ",
    "library_compounds",
    "library_rows",
    "library_rows_sourced",
    "library_rows_placeholder",
    "mechanism_rules",
    "paths_sweet",
    "paths_bitter",
    "veto_rate_vehicle_ratio",
    "veto_lif_vehicle_ratio",
    "veto_rate_fipronil_ratio",
    "veto_lif_fipronil_ratio",
    "veto_n_sweet_grn",
    "veto_n_bitter_grn",
    "null_n_shuffles",
    "null_z_imidacloprid_sign_permute",
    "null_z_imidacloprid_erdos_renyi",
    "null_z_fipronil_weight_permute",
    "null_z_fipronil_rewire_degree_preserving",
    "null_imi_max_structural_z",
    "null_fip_taste_max_z",
    "null_score_imidacloprid",
    "null_score_fipronil",
    "landscape_amplify_compounds",
    "landscape_buffer_compounds",
    "landscape_amplify_mean_gap",
    "landscape_buffer_mean_gap",
    "landscape_gap_split",
    "landscape_unscored_compounds",
    "ic50_M",
    "ic50_lo_M",
    "ic50_hi_M",
    "ic50_slope",
    "ic50_r2",
    "tornado_top_param",
    "tornado_zero_span_params",
    "rank_mean_rho",
    "rank_evaluated",
    "rank_skipped",
    "rank_known_discrepancies",
    "mix_n_targets",
    "mix_n_matching",
    "pred_n_hypotheses",
    "pred_max_n",
    "expr_coverage_overall",
    "expr_mn9_delta_hz",
    "n_genotypes",
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
    ap.add_argument("--quiet", action="store_true", help="only print the summary table")
    return ap


GUIDANCE = """\
FlyLab paper reproduction is driven by the script, not by this shim, because it
takes options (and writes into papers/). Run one of:

  python scripts/reproduce_paper.py            full run, ~3.5 min, no downloads
  python scripts/reproduce_paper.py --fast     ~1.5 min, fewer shuffles/replicates
  python scripts/reproduce_paper.py --only nulls --outdir /tmp/x
  python scripts/reproduce_paper.py --list     the named steps

Outputs: papers/figures/*.png(+svg), papers/tables/T*.csv|md, papers/results.json,
and papers/IJRC_FlyLab_draft.md re-rendered from its template.
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

    ctx = Ctx(outdir=Path(args.outdir).resolve(), fast=bool(args.fast), seed=int(args.seed), quiet=bool(args.quiet))
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
        f"seed {ctx.seed}  |  outdir {ctx.outdir}"
        + ("  |  merging into the previous results.json" if merged else "")
    )
    t_all = time.perf_counter()
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

    if ctx.figures:
        _write_captions(ctx)
    ctx.partial = bool(only)
    total = time.perf_counter() - t_all
    results = _write_results(ctx, total)
    print(_summary_table(ctx, total))
    print(f"results -> {results}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
