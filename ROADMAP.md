# Roadmap status

Full research plan: `docs/RESEARCH_MAP.md`. Numbers: `papers/results.json` (regenerate with `python scripts/reproduce_paper.py`).

| Gate | Usable when | Status |
|---|---|---|
| 0 Occupancy CLI | `flylab occupancy imidacloprid --conc 1e-6` | **done** |
| 1 Taste direction | bitter vetoes MN9 on the reduced circuit | **done** |
| 1b Named subgraph | MN9/DNp01 1-hop cut, 1126 / 1360 | **done** |
| 1c Gustatory types | LB1/LB3 seeds resolved on MaleCNS `type` | **done in v0.5** — all 16 `LB*` types proven present (was open issue #2) |
| 1d Map taste path | sugar/bitter drive on real GRN body IDs, MN9 read out | **done** (rate + LIF) |
| 2 Drug panel | insect vs vertebrate across the library | **done** — 21 compounds, 102 rows |
| 3 Local bench | assays + export in one UI | **done** |
| 4 Whole-CNS census | traced transmitter counts | **done** — committed fallback, no download |
| 4b 166k LIF | full matrix as a spike network | **out of scope** |
| 5a Mechanism table | one source of truth, rationale per rule | **done in v0.5** |
| 5b Second runtime | Shiu-style LIF beside the rate model | **done in v0.5** |
| 5c Uncertainty (local) | MC on log-potency, bootstrap Hill fit, tornado | **done in v0.5** |
| 5d Literature layer | genotype, mixtures, expression, rank comparison | **done in v0.5** |
| 5e Null models | four connectome degradations | **done in v0.5**, superseded by gate 6b |
| 5f Selectivity landscape | receptor SI vs circuit SI | **done in v0.5**, now reported over a threshold grid |
| 5g Prospective predictions | H1–H7 with effect sizes | **done in v0.5** — renamed from "pre-registered" in v0.6 |
| 5h Reproduction | `scripts/reproduce_paper.py` | **done in v0.5**, extended in v0.6 |
| 5i Static bench | same wheel on Pyodide via GitHub Pages | **done in v0.5** (`docs/PAGES.md`) |
| **6a Typed evidence** | schema v3: param type + source relation decide the transformation; unsupported rows return N/A | **done in v0.6** |
| **6b Connectome dependence** | permutation *p* first, dependence class, necessary information level, landscape over the library | **done in v0.6**, reworked in v0.6.1 |
| **6c Ablation ladder** | receptor-only / composition-only / topology-only / full | **done in v0.6** |
| **6d Specification robustness** | 25 prespecified gain specifications × 7 conclusions, plus the 3×3 threshold grid | **done in v0.6**; the topology rows were vacuous until v0.6.1 |
| **6e Global uncertainty + VOI** | Sobol' budget with a validated estimator, ranked experiments | **done in v0.6** |
| **6f Subunit-resolved nicotinic keys** | α6 (spinosad) and β1 (imidacloprid) split out; the rest stay on a documented aggregate | **done in v0.6** |
| 6g IJRC research article | Intro / Methods / Results / **Discussion** / Conclusion, inside the word target | **done in v0.6** — drafted, not submitted |
| **6h Second-round referee revision** | weight-matched transmitter null, three-way equivalence verdicts, FDR with confirmatory/exploratory labelling, a relative effect floor, a direction-aware ablation baseline with a matched reference distribution, the engine's row normalisation documented and swept, the transmitter sign table published, ground-truth validation of the dependence ladder, and the landscape repeated on the denser cut | **done in v0.6.1** — and it **reversed** the central RQ2 claim |
| 6i Scale ladder | the same profile across nested extracts of 1k-50k cells, with the permutation budget scaled to edge count | **saved result in v0.6.1** — T27; exact 5k-50k cuts are not present in this checkout, so those rungs cannot currently be rebuilt |
| 7 Expression evidence | map adult VNC motor-neuron labels to reprocessed counts and resolve MN9 separately | **open** — GSE141807 now has VNC-wide gene-detection fractions; cell classes and MN9 remain unresolved |
| 7b Scale power | rerun exact 5k-50k cuts at higher permutation resolution | **gated** — cuts missing here; profile/optimize exact edge-swap cost before raising n |
| 8 Live fly | one PER or climbing table in `live_lab` | **blocked — protocol only** |
| 9 Archival | Zenodo deposit; only then may "prospective" become "pre-registered" | **open** |

## What v0.6.1 changed about the plan

The dependence verdict turned out to be **substrate-dependent**. In the saved scale study the composition-dominated verdict occurs on exactly one of seven extracts — the in-star `named` cut — while `scale_1k`, with *fewer* nodes (1000 vs 1126) and twenty times the mean degree, is already topology-dependent. This is consistent with recurrence tracking the change, but the study does not isolate a causal feature. The 5k–50k source cuts are absent from this checkout, so T27 cannot be independently rebuilt here. The generalisation "most predictions do not need the connectome" is withdrawn. Consequences for the roadmap: (a) no dependence result may be quoted without the structural statistics of the extract it was measured on (`flylab.analysis.dependence.cut_census`, T19, T27), (b) the top rungs run at n = 20 (resolution 0.048) and confirm rather than establish a verdict, and (c) the paper is framed as a computing contribution with the fly pharmacology as its demonstration domain.

## v0.7, in order

0. **A better-powered top of the ladder.** T27 is an imported saved result; the exact 5k-50k cuts are absent from this checkout. Restore those exact cuts, profile the degree-preserving edge swaps, then increase n only after the optimized sampler preserves degree, simplicity, seeded reproducibility and the same permutation test. The existing 25k-50k results use n = 20 (resolution 0.048), where a rejection is the smallest the test can express.
1. **Expression evidence (gate 7).** Determine whether adult VNC motor-neuron labels can be joined to the reprocessed GEO GSE141807 barcodes and, if so, report class-specific detection with pinned annotations, sample/sex composition, detection rule and dropout treatment. The existing counts are VNC-wide and do not identify MN9: it is a proboscis motor neuron in the brain subesophageal zone and needs a separate class-matched atlas or targeted assay. Keep uniform gains until evidence supports a change.
2. **The one allowed fit.** The uncertainty budget now names it without ambiguity: the engagement→gain transformation (S1 0.410, VOI 1.29 of 3.149 Hz²), not a potency value, which the model cannot see at a saturating dose.
3. **The free win.** `weight_threshold` carries S1 0.294 and needs **no experiment**: re-run the dependence and selectivity analyses against synapse-confidence strata of data already held.
4. **One live assay (gate 8).** Climbing first — a full published protocol and control statistics exist; the PER baseline could not be sourced and is recorded as `null`.
5. **Archive a tagged release (gate 9)**, then and only then call the predictions pre-registered.

Not v0.7: a 166k LIF, a second connectome, a 3D viewer, or more than one fitted parameter.

## Routine maintenance

- Refresh a motor graph: Actions → `extract-malecns-subgraph` (never commit a feather).
- Rebuild the paper: `python scripts/reproduce_paper.py` (~35–45 min; `--fast` ~4 min is a sanity check, not the paper), then commit `papers/`. A test refuses a committed record that came from `--fast` **or** that was produced at any statistical effort other than the shipped defaults.
- Publish the static bench: push to `main` (Actions → `pages`).
