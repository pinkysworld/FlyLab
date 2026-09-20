# Roadmap status

Full research plan: `docs/RESEARCH_MAP.md`. Numbers: `papers/results.json`.

| Gate | Usable when | Status |
|---|---|---|
| 0 Occupancy CLI | `flylab occupancy imidacloprid --conc 1e-6` | **done** |
| 1 Taste direction | bitter vetoes MN9 on the reduced circuit | **done** |
| 1b Named subgraph | MN9/DNp01 1-hop cut, 1126 / 1360 | **done** |
| 1c Gustatory types | LB1/LB3 seeds resolved on MaleCNS `type` | **done in v0.5** — all 16 `LB*` types proven present (was open issue #2) |
| 1d Map taste path | sugar/bitter drive on real GRN body IDs, MN9 read out | **done** (rate + LIF) |
| 2 Drug panel | insect vs vertebrate across the library | **done** — 21 compounds, 99 rows, evidence-tiered |
| 3 Local bench | assays + export in one UI | **done** — 34 routes, eleven panels |
| 4 Whole-CNS census | traced transmitter counts | **done** — committed fallback, no download |
| 4b 166k LIF | full matrix as a spike network | **out of scope** |
| 5a Mechanism table | one source of truth, rationale per rule | **done in v0.5** |
| 5b Second runtime | Shiu-style LIF beside the rate model | **done in v0.5** |
| 5c Uncertainty | MC on log-EC50, bootstrap Hill fit, tornado | **done in v0.5** |
| 5d Literature layer | genotype, mixtures, expression, rank validation | **done in v0.5** |
| 5e Null models | four connectome degradations + score | **done in v0.5** |
| 5f Selectivity landscape | receptor SI vs circuit SI | **done in v0.5** |
| 5g Pre-registration | H1–H7 with effect sizes and suggested *n* | **done in v0.5** |
| 5h Reproduction | `scripts/reproduce_paper.py` | **done in v0.5** |
| 5i Static bench | same wheel on Pyodide via GitHub Pages | **done in v0.5** (`docs/PAGES.md`) |
| 6 IJRC draft | methods / software article | **done in v0.5** — drafted, not submitted |
| 7 Live fly | one PER or climbing table in `live_lab` | **blocked — protocol only** |

## v0.6, in order

1. **Expression coverage** — currently 0.21 across five insect receptor keys, with adult motor neurons a confirmed literature gap. Either source it or declare the bound permanently.
2. **Subunit-resolved receptors** — split `insect_nAChR` into α6 / β1-dependent keys so spinosad and the neonicotinoids stop sharing a target. `library.yaml` schema v3.
3. **One live assay** — climbing first (a full published protocol and control statistics exist); PER's specific baseline could not be sourced and is recorded as `null`.

Not v0.6: a 166k LIF, a second connectome, a 3D viewer, or more than one fitted parameter.

## Routine maintenance

- Refresh a motor graph: Actions → `extract-malecns-subgraph` (never commit a feather).
- Rebuild the paper: `python scripts/reproduce_paper.py`, then commit `papers/`.
- Publish the static bench: push to `main` (Actions → `pages`).
