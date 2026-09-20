# FlyLab research map

The research plan, not a results paper. Software that exists today is marked **done**. Live-animal work is marked **blocked**. Numbers quoted here come from `papers/results.json` (regenerate with `python scripts/reproduce_paper.py`).

## Claim the project is allowed to make

Same compound, same concentration, two panels, **typed**:

1. Insect: engagement of insect nAChR (α6 / β1 / aggregate) / RDL / GluCl / AChE / Nav / OctR → gain patch on *named* MaleCNS objects.
2. Vertebrate: engagement of α4β2 / α7 / GABA-A / GlyR / AChE / Nav1.x — scored, never simulated.

Say **engagement**, not occupancy, unless the row carries a measured `Kd`/`Ki` (exactly one does: imidacloprid at `insect_nAChR_beta1`). A row with no sourced value returns **N/A**, enforced by the type system, and is excluded from every numeric result.

One notebook JSON. Map version and library hash recorded. No invented potency value. No fake live *n*.

And, since v0.6, no prediction is stated without the answer to *what does it depend on*: every headline claim carries its connectome-dependence class and the weakest graph model that reproduces it.

## Map objects (MaleCNS v1.0)

| Object | What it is | Status |
|---|---|---|
| Atlas census | 165 122 traced bodies, consensus transmitters | **done** (committed fallback, no download) |
| Named motor seeds | type `MN9`, `DNp01` | **done** (real body IDs) |
| 1-hop neighbourhood (`named`) | 1126 nodes / 1360 edges, weight ≥ 5 | **done** (Actions + `data/derived/`) |
| Labellar GRN types | LB1a–d bitter, LB3b/c sweet | **done — PROVEN present in MaleCNS v1.0** (see Gate 1c) |
| Taste-motor cut (`taste_motor`) | 1841 nodes / 19 066 edges, 1 hop + induced closure ≥ 10 | **done** |
| Full 25M-edge LIF | every traced cell | **out of scope** |

Public feathers (`body-annotations…`, `body-neurotransmitters…`, `connectome-weights…` ≈ 1.1 GB) are never committed. Cut on GitHub Actions; keep each derived JSON ≤ ~1 MB by raising the weight floor, never by adding hops.

## Gates

| Gate | Usable when | Status |
|---|---|---|
| 0 Occupancy CLI | `flylab occupancy imidacloprid --conc 1e-6` | **done** |
| 1 Taste direction | bitter vetoes MN9 on the reduced circuit | **done** |
| 1b Named subgraph | MN9/DNp01 1-hop cut from the real weight matrix | **done** |
| 1c Gustatory types | LB1/LB3 seeds resolved from MaleCNS `type` | **done — this was open issue #2 and it is now closed** |
| 1d Map taste path | sugar/bitter drive on extracted GRNs, read MN9 | **done** (both engines) |
| 2 Drug panel | ≥ 5 compounds, insect vs vertebrate | **done** (21 compounds, 101 typed rows, 14 receptor keys) |
| 3 Local bench | assays + export in one UI | **done** |
| 4 Whole-CNS census | traced transmitter counts | **done** |
| 4b 166k LIF | full matrix as a spike network | **out of scope** |
| 5 Live fly | PER / climbing *n* | **blocked — protocol only** |
| 6 IJRC | research article with a Discussion, or a live table | **v0.6 draft complete, not submitted** |
| 6a Typed evidence | parameter type + source relation decide the transformation | **done in v0.6** |
| 6b Connectome dependence | permutation *p*, dependence class, necessary information level | **done in v0.6** |
| 6c Ablation ladder | which layer carries the ordering information | **done in v0.6** |
| 6d Specification robustness | 25 gain specifications x 7 conclusions + threshold grid | **done in v0.6** |
| 6e Global uncertainty + VOI | Sobol' budget, validated estimator, ranked experiments | **done in v0.6** |
| 7 Expression coverage | adult motor-neuron receptor expression, or a declared bound | **open — 0.205, MN9 unresolved** |
| 8 Archival | Zenodo deposit; only then may "prospective" become "pre-registered" | **open** |

### Gate 1c is closed (former issue #2)

All 16 `LB*` types named by the Cell 2026 gustatory nomenclature appear in the MaleCNS v1.0 `type` column. The six FlyLab drives are **LB3b (11) + LB3c (23) = 34 sweet** and **LB1a (11) + LB1b (6) + LB1c (16) + LB1d (5) = 38 bitter** cells, out of 1428 traced gustatory neurons and 163 labellar-bristle neurons. The seed list in `data/derived/malecns_gustatory_seeds.json` is real body IDs, not a wish list, and `flylab assay-taste-map` drives them.

Two caveats travel with it, permanently: the sweet/bitter assignment of the types is a **working hypothesis** from the published nomenclature and not a MaleCNS annotation; and some of these cells (all of LB1b) carry `unclear` as their predicted transmitter and contribute no sign to the model.

## v0.5 gates (all done)

| Gate | What it added |
|---|---|
| 5a Mechanism table | one source of truth for receptor → gain, with a rationale per rule (T2) |
| 5b Second runtime | Shiu-style LIF beside the rate model, with a stated calibration and a background Poisson drive |
| 5c Uncertainty | log-EC50 Monte-Carlo, replicate ensembles, bootstrap Hill fits, sensitivity tornado |
| 5d Literature layer | six sourced YAML datasets; genotype shifts, mixtures, expression weighting, rank validation |
| 5e Null models | four connectome degradations (superseded in v0.6 by connectome-dependence analysis) |
| 5f Selectivity | receptor SI vs circuit SI landscape over the library and both cuts |
| 5g Prospective predictions | H1–H7 with model-internal effect sizes and an illustrative planning *n* |
| 5h Reproduction | `scripts/reproduce_paper.py`: every figure, table and number from committed data |
| 5i Static bench | the same wheel on Pyodide, published from GitHub Pages (`docs/PAGES.md`) |

## Pharmacological hypotheses (prospective, not yet tested in vivo)

`flylab predictions` prints these with model-internal CIs and an **illustrative planning minimum** *n* per group (capped at d = 1.0; not a power calculation). They are *prospective*, not pre-registered: a mutable repository is not a registry. Full table: `papers/tables/T3_predictions.csv`.

H1 dual panel gap at 1 µM · H2 imidacloprid lowers the neighbourhood mean rate · H3 fipronil raises it · H4 diazepam moves nothing · H5 bitter vetoes sugar-driven MN9 on the map path · H6 fipronil relieves that veto · H7 picrotoxin perturbs the circuit with ~zero receptor selectivity.

**Software check: done for all seven. Live check: blocked for all seven.** `live_result` is `null` in every row.

## What v0.6 learned that changes the plan

Every number here comes from `papers/results.json`.

1. **Most predictions do not need the connectome; a chemically coherent minority do.** Over 21 compounds × 4 concentrations: 20 cells topology-dependent, 51 composition-dominated, 13 no-effect, 0 mixed, and the topology-dependent set is exactly the chloride-channel blockers. Imidacloprid at 1 µM is composition-dominated (p = 0.275 / 0.586 / 0.472 against the three structure-preserving nulls at n = 1000; necessary level `degree_sequence`); fipronil is topology-dependent (p = 0.0060 weight, 0.0070 degree; necessary level `wiring_without_transmitter_identity`). Never claim "the connectome matters" without naming the compound, the concentration and the readout.
2. **The bitter-veto arm is a properly powered negative.** At n = 300 on the taste-motor cut, fipronil's veto ratio gives p = 0.42 / 0.88 / 0.20 / 0.24 with all |z| < 0.25, replacing the ten-shuffle hand-wave of v0.5. H6 is a pharmacology prediction, not a wiring prediction, and that is now measured rather than assumed. Imidacloprid's veto ratio there is **undefined** (MN9 silenced), not zero.
3. **Nearly all the ordering information is in the mechanism rules.** Composition-only reproduces the full model's ordering (ρ 0.906–0.989; ρ 0.953 within the nicotinic agonists); receptor-only orders it backwards (ρ −0.05 to −0.51); topology-only with a generic multiplier is the worst level (ρ 0.24, 11.7 Hz rms). The connectome without the pharmacology is not a cheap substitute for the connectome with it.
4. **The suppression result is specification-dependent.** "A nicotinic agonist suppresses circuit activity" is retained by 15 of 25 prespecified gain specifications and reversed by all ten monotone ones, where the same drug at the same engagement *excites* the network. Nicotinic buffering (18/25) is reversed by none; Nav/AChE amplification is 19/25; RDL disinhibition, both topology conclusions and the map bitter-veto direction are 25/25. Quote the specification-independent ones freely; caveat the rest.
5. **The uncertainty is in two places only.** Sobol' at n_base 1024 (11 264 evaluations, Var(Y) = 3.149 Hz²): `gain_transform` S1 0.410 (ST 0.626), `weight_threshold` S1 0.294 (ST 0.542), drive 0.019, everything else at the ±0.007 noise floor measured by a deliberately null factor. Potency and Hill *n* are invisible at a saturating dose. Interactions carry 0.344.
6. **The next experiment is named by the model, not by taste.** VOI: `gain_transform` 1.29 Hz² (a synaptic-gain calibration at a known synapse) and `weight_threshold` 0.925 Hz² — which needs **no experiment at all**, only a re-run against synapse-confidence strata.
7. **The amplify/buffer split is threshold-sensitive.** Nicotinic buffering holds 9/9 on the 3×3 grid (gap −0.362 to −1.224); Nav/AChE amplification holds 8/9 and flips at the strictest corner (gap −0.075). Quote ranges, never the point estimate.
8. **Typing the evidence removed numbers that never existed.** 101 rows: 34 EC50, 10 IC50, 1 Kd, 56 not modelled. Diazepam at insect RDL is N/A, not 1e-4.

## v0.6 gates (all done)

| Gate | What it added |
|---|---|
| 6a Typed evidence | schema v3: parameter type + source relation decide the transformation; unsupported rows return N/A (T10, F11) |
| 6b Connectome dependence | permutation *p* first, dependence class, necessary information level, a landscape over the library, a permutation-count sweep (T6, T11, T12, F6, F12) |
| 6c Ablation ladder | receptor-only / composition-only / topology-only / full (T13, F13) |
| 6d Specification robustness | 25 gain specifications × 7 conclusions, plus the 3×3 threshold grid (T14, T15, F14) |
| 6e Global uncertainty + VOI | Sobol' budget with an Ishigami-validated estimator, ranked experiments (T16, T17, F15) |
| 6f Subunit-resolved nicotinic keys | α6 (spinosad), β1 (imidacloprid); the rest stay on a documented aggregate because their sources used hybrid constructs |
| 6g IJRC research article | Introduction / Methods / Results / Discussion / Conclusion plus a supplement, inside the word target |

## v0.7 (the next things, in order)

1. **Expression coverage.** 0.205 across the insect receptor keys, and motor neurons — MN9's class — are a confirmed gap. Either source it or declare the bound permanent. A data problem, not a modelling one.
2. **The one allowed fit, now named by the model.** The engagement→gain transformation. Not a potency: at a saturating dose the model cannot see one.
3. **The free win.** `weight_threshold` needs no bench time: re-run against synapse-confidence strata.
4. **One live assay.** Climbing first; the PER baseline could not be sourced and is `null`.
5. **Archive a tagged release** before the word "pre-registered" may be used anywhere in this repository.

Not v0.7: a 166k-cell LIF, a second connectome, a 3D viewer, or fitting more than one parameter.

## Parameter policy

Allowed to fit later, at most one coefficient per paper: the engagement→gain transformation or its coefficient — the variance budget says so — and never a potency value.
Not allowed: refitting MaleCNS weights, inventing transmitters, hiding a reduced-circuit warning, or editing a library EC50 to improve a rank correlation.

## Assays

| Assay | Command | Reads | Writes |
|---|---|---|---|
| Occupancy | `flylab occupancy` | `library.yaml` | receptor table + selectivity |
| Taste (reduced control) | `flylab assay` | 5 rate units | MN9 Hz, gains |
| Taste (map) | `flylab assay-taste-map` | `taste_motor` cut | MN9 Hz, veto ratio |
| Neighbourhood (rate) | `flylab assay-subgraph` | `named` cut | named-cell Hz, mean Hz, gains |
| Neighbourhood (LIF) | `flylab assay-spiking` | `named` cut | spikes, rates |
| Census | `flylab assay-cns` | census fallback or atlas | transmitter counts, excitation index |
| Live PER / climbing | `protocols/*.md` | flies | `live_lab` only, by hand |

## Paper path (IJRC)

`papers/IJRC_FlyLab_draft.md` is a research article (Introduction / Methods / Results / Discussion / Conclusion) with `papers/SUPPLEMENT.md`, rendered from `papers/IJRC_FlyLab_draft.md.in`. Submit when either the journal accepts a software article with no live table, or one PER/climbing table exists in a notebook. Do not submit H2–H7 as empirical biology without that table. Checklist: `papers/SUBMISSION_CHECKLIST.md`.

## What this map is not

Not a new connectome. Not GLP toxicology. Not a 166k-cell product. Not a claim that the teaching EC50s are the fly's EC50s.
