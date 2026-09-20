# FlyLab handoff (v0.5, 2026-09-19)

Read this first if you are Claude, GPT, or another agent continuing the repo.
Owner: https://github.com/pinkysworld/FlyLab
Package version: **0.5.0**.

## What the project is

A **virtual pharmacology bench** on the public **MaleCNS v1.0** adult male *Drosophila* CNS (brain + ventral nerve cord). One compound + one free concentration → insect occupancy panel + vertebrate occupancy panel at the same dose → gain patch on named MaleCNS cells → circuit readouts → one notebook JSON with provenance.

It is **not** a new connectome, not GLP tox, not a 166k-cell LIF product, and not live-animal results.

| Document | What it is for |
|---|---|
| `docs/RESEARCH_MAP.md` | the canonical plan, gates and what v0.6 is |
| `docs/NOVELTY.md` | what is actually new (read before accepting a PR) |
| `docs/ARCHITECTURE.md` | dataflow and layer ownership |
| `docs/STACK.md` | technology choices and the browser constraints they impose |
| `docs/PAGES.md` | the static (Pyodide) bench |
| `papers/IJRC_FlyLab_draft.md` | the manuscript (generated — edit the `.in` template) |
| `papers/STATUS.md`, `papers/SUBMISSION_CHECKLIST.md` | where the paper stands |

## Do this first

```bash
git pull
python -m pip install -e ".[dev,viz]"
python -m pytest -q                       # fast suite; `-m slow` adds the long checks
flylab occupancy imidacloprid --conc 1e-6
flylab assay-subgraph --compound imidacloprid --conc 1e-6
python scripts/reproduce_paper.py --fast  # ~1 min sanity check of the whole pipeline
flylab serve                              # http://127.0.0.1:8765
```

## What changed in v0.5

Everything below is new since the v0.4 handoff.

**Pharmacology**
- Library grew to **21 compounds / 99 receptor rows / 12 receptor classes**, schema v2, every row carrying `source` and `evidence_tier` (43 `literature_order`, 56 `class_placeholder`).
- `flylab/pharm/mechanisms.py` is now the **single source of truth** for receptor → gain. The four v0.4 rules are verbatim; four families were added (GluCl, AChE, Nav, OctR). No assay carries its own copy.
- New modules: `binding.py` (Gaddum/Schild/operational), `exposure.py` (one-compartment C(t)), `uncertainty.py` (log-EC50 Monte-Carlo), `genotype.py` (per-compound allele shifts), `mixtures.py` (Bliss/Loewe + published validation), `expression.py` (cross-atlas weighting, sensitivity only).
- **A library error was corrected**: fipronil's `vertebrate_GABA_A` was 1.0e-5 M while citing a paper reporting 1.1e-6 M. It is now the cited value, and fipronil's vertebrate occupancy at 1 µM is 0.48 — "vertebrate-safe" is no longer sayable.

**Circuit**
- `circuit/rate.py` (deterministic, ~50 ms per run) and `circuit/lif.py` (Shiu-style, dt 0.1 ms) extracted as shared runtimes.
- Second committed cut: `data/derived/malecns_taste_motor_neighborhood.json`, 1841 nodes / 19 066 edges, 1 hop + induced closure ≥ 10 synapses.
- `assays/taste_map.py` drives the **real** labellar GRN body IDs and reads MN9, on either engine.
- `assays/ensemble.py` (replicates, CIs, sensitivity tornado, model IC50) and `assays/experiment.py` (batch designs).

**Analysis and validation (new packages)**
- `analysis/nullmodels.py` — four connectome degradations + connectome-information score.
- `analysis/selectivity.py` — receptor SI vs circuit SI landscape.
- `analysis/predictions.py` — pre-registered H1–H7 with power.
- `analysis/{fit,impact,layout}.py` — Hill fits with bootstrap, per-edge/node/path impact, deterministic viewer layout.
- `validation/rank.py` — Spearman/Kendall against published orderings, with recorded discrepancies.
- `data/literature/` — six sourced YAML datasets, 55 citations with DOI/PMID, plus a list of ten places where the literature contradicts the library.

**Interfaces and reproduction**
- 34 API routes, eleven UI panels, an expanded CLI, and the same wheel running in the browser on Pyodide (`docs/PAGES.md`).
- **`scripts/reproduce_paper.py`** — one command regenerates 10 figures, 10 tables and `papers/results.json` from committed data in ~3.5 min (~1 min with `--fast`). The manuscript is rendered from a template by substituting those keys, so a number cannot drift between code and prose.

## Issues

- #1 Gate 0 occupancy CLI — **closed**.
- #2 Gate 1 map-extracted taste GRNs — **CLOSED in v0.5**. All 16 `LB*` types are present in the MaleCNS `type` column: 34 sweet (LB3b 11, LB3c 23) and 38 bitter (LB1a 11, LB1b 6, LB1c 16, LB1d 5) seed cells, driven by `flylab assay-taste-map`.
- #3 Gate 2 five-drug panel — **closed** (now 21 compounds).
- **Open (new): version strings disagree.** `flylab/__init__.py` reports `0.2.0`, `pyproject.toml` says `0.4.0`, and notebook provenance reports `0.5.0`. Pick one. Not fixed here because those files belong to other owners.
- **Open: expression coverage is 0.21** and adult motor-neuron receptor expression is a confirmed literature gap — MN9, the main readout, carries the `unknown` default weight.
- **Open: `insect_nAChR` is one key for several receptors** (α6 spinosad vs β1 neonicotinoids). A schema v3 change.

## Findings that constrain what you may claim

Read these before writing any new claim into the paper or the README.

1. **Imidacloprid's mean-rate effect is not a wiring result.** It beats an Erdős–Rényi null (|z| ≈ 46) and fails sign-permutation, weight-permutation and degree-preserving rewiring (max |z| ≈ 1.2). For nAChR agonists this readout tracks global E/I balance. Fipronil *does* beat the topology nulls (z ≈ 3). Never claim "the connectome matters" without naming the compound and the readout.
2. **The bitter veto ratio is not yet distinguishable from a shuffled graph** (max |z| ≈ 0.7 at 10 shuffles). H6 is a pharmacology prediction, not a wiring prediction.
3. **At 1 µM the model is exactly insensitive to the EC50 and the Hill coefficient** (zero sensitivity span). The gain-rule coefficient is the dominant parameter and is the correct target for the one allowed fit.
4. **Imidacloprid at 1 µM silences MN9 on the map path**, so its bitter-veto ratio is *undefined*, not 0.00. That is the agonist rule at its 0.05 floor.
5. **Eight compounds have no circuit selectivity index**: RDL/GluCl block cannot reach the 50 % threshold on these cuts.
6. **`data/literature/README.md` lists ten literature-vs-library contradictions.** Report them; do not silently edit the library to match.

## Hard constraints (do not break)

- Organism is **fly** (*Drosophila*); map is **MaleCNS v1.0** (MN9/DNp01 need the cord, so not FlyWire-only).
- **Patch rules may not change silently.** `insect nAChR agonist: g_ach = max(0.05, 1 + 0.4θ − 1.6θ²)`; `insect nAChR antagonist: g_ach = max(0.05, 1 − θ)`; `insect RDL antagonist: g_gaba = max(0.05, 1 − θ)`; `insect RDL agonist/PAM: g_gaba = max(0.05, 1 + 0.4θ)`. Changing one means editing `flylab/pharm/mechanisms.py`, its rationale docstring, and T2 — and saying so.
- **No invented EC50.** New rows need a `source` string and an `evidence_tier`. A missing number stays `class_placeholder`.
- **At most one fitted parameter later.** Never refit MaleCNS weights.
- **Notebooks keep their warnings** when the circuit is reduced or hops-limited.
- **`live_lab` stays `null`.** No code path may write a live-animal number. `protocols/` are unexecuted drafts.
- **Genotype shifts are per compound, never per receptor** (Rdl A301S moves GABA but not fipronil; para M918T moves deltamethrin but not DDT).
- **Never commit `*.feather`.** Do not add hops to a cut; raise the weight floor. Do not hand-edit `data/derived/*.json`.
- The science core must not import `pydantic`, `fastapi`, `typer`, `pandas` or `pyarrow` at module scope — the browser build has none of them.

## Big files (never commit)

From `https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome`:
`body-annotations-male-cns-v1.0-minconf-0.5.feather`, `body-neurotransmitters-male-cns-v1.0.feather`, `connectome-weights-male-cns-v1.0-minconf-0.5.feather` (~1.1 GB).
Local dir: `$FLYLAB_MALECNS` or `~/.flylab/malecns_v1`. `.gitignore` drops `*.feather`. The weight columns used are `body_pre`, `body_post`, `weight`; stream with `pyarrow.ipc.open_file` batches — a full pandas read OOMs. **None of this is needed to reproduce the paper.**

## Commands that matter

```bash
# pharmacology
flylab list-drugs
flylab occupancy imidacloprid --conc 1e-6
flylab compare imidacloprid fipronil --conc 1e-6

# circuits
flylab assay-subgraph  --compound imidacloprid --conc 1e-6     # rate, named cut
flylab assay-spiking   --compound fipronil     --conc 1e-6     # LIF
flylab assay-taste-map --compound fipronil     --conc 1e-6     # real GRNs -> MN9
flylab assay-cns       --compound imidacloprid --conc 1e-6     # census
flylab graph info --graph taste_motor

# analysis
flylab ic50 --compound imidacloprid
flylab sensitivity --compound imidacloprid --conc 1e-6
flylab null-panel --compound fipronil
flylab selectivity --conc 1e-6
flylab predictions

# paper
python scripts/reproduce_paper.py              # full, ~3.5 min
python scripts/reproduce_paper.py --fast       # ~1 min
python scripts/reproduce_paper.py --only nulls --outdir /tmp/x
python scripts/reproduce_paper.py --only paper # re-render prose from results.json
python scripts/reproduce_paper.py --list
flylab reproduce-paper                         # prints the recipe; does not run it

# data (optional, never committed)
flylab download-malecns           # ~55 MB atlas
flylab download-malecns --full    # + 1.1 GB weights
flylab extract-subgraph --types MN9,DNp01,LB1a,LB1b,LB1c,LB1d,LB3b,LB3c
```

## Files not to rewrite from scratch

`flylab/maps/extract.py` (streaming), `flylab/pharm/mechanisms.py` (the patch rules), `flylab/pharm/library.yaml` (cite sources if you change an EC50), `data/derived/*.json` (CI products), `data/literature/*.yaml` (sourced datasets), `papers/results.json` (generated).

## Next work, in order

1. **v0.6 gate 1 — expression coverage.** Either source adult motor-neuron receptor expression or declare the 0.21 coverage a permanent bound. A data problem, not a modelling one.
2. **v0.6 gate 2 — subunit-resolved nicotinic keys** (library schema v3). Until then every nicotinic selectivity number is capped.
3. **v0.6 gate 3 — one live assay.** Climbing first: a full published protocol and control statistics exist. Import by hand into `live_lab`; nothing else may write it.
4. **Editorial:** the manuscript body is ~8 000 words against a 5 000–7 000 target. Trim §4/§5 or move the mechanism rationale and null-model definitions to a supplement.
5. **Housekeeping:** reconcile the three version strings; mint the Zenodo DOI and replace the placeholder in the draft and `CITATION.cff`.

**Do not** start a 166k-cell LIF, add a second connectome, invent live PER/climbing numbers, or accept a PR that only adds a prettier viewer.
