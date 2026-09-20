# FlyLab research map

The research plan, not a results paper. Software that exists today is marked **done**. Live-animal work is marked **blocked**. Numbers quoted here come from `papers/results.json` (regenerate with `python scripts/reproduce_paper.py`).

## Claim the project is allowed to make

Same compound, same concentration, two scorecards:

1. Insect: occupancy on insect nAChR / RDL / GluCl / AChE / Nav / OctR → gain patch on *named* MaleCNS objects.
2. Vertebrate: occupancy on α4β2 / α7 / GABA-A / GlyR / AChE / Nav1.x — scored, never simulated.

One notebook JSON. Map version and library hash recorded. No invented EC50. No fake live *n*.

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
| 2 Drug panel | ≥ 5 compounds, insect vs vertebrate | **done** (21 compounds, teaching tiers) |
| 3 Local bench | assays + export in one UI | **done** |
| 4 Whole-CNS census | traced transmitter counts | **done** |
| 4b 166k LIF | full matrix as a spike network | **out of scope** |
| 5 Live fly | PER / climbing *n* | **blocked — protocol only** |
| 6 IJRC | software note or live table | **draft complete, not submitted** |

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
| 5e Null models | four connectome degradations and a connectome-information score |
| 5f Selectivity | receptor SI vs circuit SI landscape over the library and both cuts |
| 5g Pre-registration | H1–H7 with model-internal effect sizes and a suggested *n* |
| 5h Reproduction | `scripts/reproduce_paper.py`: every figure, table and number from committed data |
| 5i Static bench | the same wheel on Pyodide, published from GitHub Pages (`docs/PAGES.md`) |

## Pharmacological hypotheses (pre-registered, not yet tested in vivo)

`flylab predictions` prints these with model-internal CIs and a suggested *n* per group. Full table: `papers/tables/T3_predictions.csv`.

H1 dual scorecard gap at 1 µM · H2 imidacloprid lowers the neighbourhood mean rate · H3 fipronil raises it · H4 diazepam moves nothing · H5 bitter vetoes sugar-driven MN9 on the map path · H6 fipronil relieves that veto · H7 picrotoxin perturbs the circuit with ~zero receptor selectivity.

**Software check: done for all seven. Live check: blocked for all seven.** `live_result` is `null` in every row.

## What v0.5 learned that changes the plan

1. **The mean-rate readout is not a wiring result for nicotinic agonists.** Imidacloprid beats an Erdős–Rényi null decisively and fails every structure-preserving null (max |z| ≈ 1.2). Fipronil beats the topology nulls (z ≈ 3). Any future claim that "the connectome matters" must name the compound and the readout.
2. **The bitter veto ratio is not yet distinguishable from a shuffled graph** (max |z| ≈ 0.7 at 10 shuffles). More shuffles and a readout less dominated by the drive are needed before H6 can be called a wiring prediction.
3. **At a saturating dose the model cannot see its own EC50** (zero sensitivity span). Dose ladders, not single doses, constrain the receptor layer; the gain-rule coefficient is the parameter that matters.
4. **RDL/GluCl blockers have no circuit selectivity index** on these cuts: disinhibition is bounded by the GABA share and the 0.05 gain floor.
5. **Correcting one library row removed a claim.** Fipronil's vertebrate GABA-A EC50 contradicted its own citation by ~9×; after correction its vertebrate occupancy at 1 µM is 0.48, and "vertebrate-safe" is no longer sayable.

## v0.6 (the next three things, in order)

1. **Expression coverage.** Coverage is 0.21 across the five insect receptor keys, and motor neurons — the class MN9 belongs to — are a confirmed gap. Either find per-cell-type receptor expression for adult leg/labellar motor neurons, or state the gap as a permanent bound on the weighted layer. Not a new model; a data problem.
2. **Subunit-resolved receptors.** `insect_nAChR` is one key where α1/α5/α6/α7 and β1/β2 are different, partly overlapping cell populations. Spinosad (α6) and the neonicotinoids (β1) must not share a key. This is a `library.yaml` schema change (v3) plus per-subunit expression weights, and it caps how much the selectivity landscape can mean until it is done.
3. **One live assay.** Climbing is the better first target: a full published protocol and control statistics exist (Martelli 2020; RING), whereas the specific PER baseline FlyLab needs was not found in the literature and is recorded as `null`. One table in `live_lab` converts this project from a methods article to a results article.

Not v0.6: a 166k-cell LIF, a second connectome, a 3D viewer, or fitting more than one parameter.

## Parameter policy

Allowed to fit later, at most one coefficient per paper: a gain-rule coefficient (now known to be the sensitive one) or one EC50.
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

`papers/IJRC_FlyLab_draft.md` is a methods / software article, rendered from `papers/IJRC_FlyLab_draft.md.in`. Submit when either the journal accepts a software article with no live table, or one PER/climbing table exists in a notebook. Do not submit H2–H7 as empirical biology without that table. Checklist: `papers/SUBMISSION_CHECKLIST.md`.

## What this map is not

Not a new connectome. Not GLP toxicology. Not a 166k-cell product. Not a claim that the teaching EC50s are the fly's EC50s.
