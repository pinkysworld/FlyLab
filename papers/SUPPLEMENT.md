# Supplementary material

**Supplement to:** *Auditing Input Dependence in Literature-Parameterised Network Simulations: A Connectome Case Study*.

FlyLab 0.6.0, library SHA-256 `359096467b33`, map `male-cns:v1.0`. An internal code revision is recorded in the provenance manifest; a public revision identifier is omitted from blind-review text. Reference numbers `[n]` use the 76-entry source-paper bibliography. The standard script regenerates core model values from the current parameter library; the larger scale-study values are imported from a saved record whose 5k–50k cuts are absent here. No value was measured in a living animal. The library source-value audit remains incomplete.

## S1. Mechanism rules: the gain patch and its rationale

`flylab/pharm/mechanisms.py` is the single source of truth mapping engagement rows onto circuit gains; no assay carries its own copy, which is what makes Table T2 a truthful description of the software. Gains are dimensionless multipliers on the synaptic weight of a transmitter class (or, for `g_nav`, on node excitability); 1.0 is vehicle, and every gain is floored at 0.05 because a network with a gain of exactly zero is not informative. The 15 rules are listed in full in T2. The four core rules are:

* **Insect nAChR, activating:** `g_ach = max(0.05, 1 + 0.4θ − 1.6θ²)`. Rises at low engagement and falls at high engagement, encoding excitation followed by desensitisation and depolarisation block. At θ = 0.99 it is at the floor, which is why a saturating dose of a nicotinic agonist silences cholinergic transmission in the model rather than merely reducing it — and, as §3.4 of the manuscript shows, why the suppression conclusion is a property of *this* rule.
* **Insect nAChR, blocking:** `g_ach = max(0.05, 1 − θ)`.
* **Insect RDL, blocking:** `g_gaba = max(0.05, 1 − θ)`. Blocking the GABA-gated chloride channel removes inhibition.
* **Insect RDL, activating:** `g_gaba = max(0.05, 1 + 0.4θ)`.

Four additional families use stated modelling rationales: GluCl activation uses a larger coefficient (`1 + 0.8θ`) than RDL because avermectin opening is effectively irreversible [37]; AChE inhibition raises an ambient cholinergic tone (`ach_tone = 1 + 2.0θ`) multiplying `g_ach` together with the agonist curve, reproducing the organophosphate excitation-then-block picture; Nav modulation scales node excitability globally (`g_nav = 1 + 1.5θ`), since pyrethroids and DDT act on every cell regardless of transmitter [22], [66]; octopamine agonism uses the smallest coefficient (`1 + 0.5θ`) because it is modulatory.

**The level-C generic rule.** The ablation ladder gives its topology-only level one signed multiplier per compound rather than a mechanism-specific vector. The direction-aware rule is `g_all = max(0.05, 1 - theta_max) if the compound's dominant mechanism lowers net synaptic drive, else 1 + theta_max`; its magnitude matches the depression-only comparison `g_all = max(0.05, 1 - theta_max)`. The historical rule: max(0.05, 1 - theta_max) on every transmitter. Structurally incapable of a positive effect, so it cannot order a library containing RDL/GluCl blockers. Reported as a FLOOR (what a depression-only connectome model achieves), never as a competitor. Both are reported for 21 compounds at 4 concentrations in T13 and T13b. Their difference assesses sensitivity to the generic rule and is not a biological contrast.

**These coefficients and shapes are modelling choices, not measurements.** The manuscript evaluates 25 gain-rule specifications (§2.4), but the five asserted transmitter sign magnitudes (glutamate -0.4, histamine -0.5, dopamine +0.2, serotonin +0.2, octopamine +0.2; T21) and transmitter-defined gain assignment remain fixed. The label-permutation null changes which nodes receive labels; it does not vary the sign magnitudes. Therefore the family does not bound uncertainty from these assumptions.

### S1.1 The specification family

Five shapes — monotone-linear, saturating, weakly biphasic, the shipped biphasic rule and a strongly biphasic rule — at five coefficient scales (0.50×, 0.75×, 1.00×, 1.25×, 1.50×), giving 25 specifications, of which `flylab_biphasic@1.00x` is the shipped one. Installation is process-global, so specifications are evaluated serially. Full matrix: T14; per-specification diagnostics: T14b.

| Conclusion | Claim | Retained | Reversed | Undecidable |
|---|---|---|---|---|
| C1 | a nicotinic agonist (imidacloprid, 1 uM) suppresses network activity | 15 | 10 | 0 |
| C2 | an RDL antagonist (fipronil, 1 uM) disinhibits the network | 25 | 0 | 0 |
| C3 | imidacloprid's mean-rate effect is not distinguishable from any structure-preserving degradation of the cut (not topology-dependent) | 25 | 0 | 0 |
| C4 | fipronil's mean-rate effect is topology-dependent where imidacloprid's is not, on the same cut, shuffles and criterion | 25 | 0 | 0 |
| C5 | the neighborhood buffers nicotinic receptor selectivity | 18 | 0 | 7 |
| C6 | the neighborhood amplifies Nav / AChE receptor selectivity | 18 | 7 | 0 |
| C7 | on the map path, bitter drive lowers MN9 under fipronil | 25 | 0 | 0 |

**Specification-test scope.** The topology predicates use 100 shuffles per specification, allowing empirical probabilities to reach α = 0.050. The specification is passed explicitly into the null engine. The C1–C7 rows of T14 are predicates about model outputs under this fixed implementation and selected specification family; they are not biological validation.

The topology predicates are retained in 25/25 and 25/25 cases after Benjamini-Hochberg adjustment of the 100 structural component tests (uncorrected counts 25 and 25). This adjustment does not give a cell-level or predicate-level FDR guarantee. Among C3 cases, 17 pass the descriptive point-gap criterion and 8 remain indeterminate; the C4 split is 0 / 0. These counts are conditional on the fixed sign table, assignment rule and common shuffle stream, not independent replications.

Specifications reversing C1: linear@0.50x, linear@0.75x, linear@1.00x, linear@1.25x, linear@1.50x, saturating@0.50x, saturating@0.75x, saturating@1.00x, saturating@1.25x, saturating@1.50x. The mean selectivity-index gap in the nicotinic set under the shipped gain rule is -0.884 log₁₀ units; over 6 reversing specifications with a defined index it is -0.893. The per-specification range is -0.199 to -1.734; 3 buffer less than the shipped rule (linear@1.50x, saturating@1.25x, saturating@1.50x). An aggregate gap should not be interpreted as a uniform response across specifications.

## S2. Null-model definitions, and the conventions around them

`flylab/analysis/nullmodels.py` and `flylab/analysis/dependence.py` degrade the same committed cut, in increasing order of destruction:

* **`sign_permute_weight_matched`** — a constrained shuffle of consensus transmitter labels: random pairs of differently-labelled nodes are swapped, and a swap is kept only if every tracked transmitter's share of *total outgoing synaptic weight* stays within 0.002 of the real graph's. The node count of each label is preserved exactly, the weighted share to within the tolerance. This is the ladder's rank-3 rung.
* **`sign_permute`** — the plain permutation of the node label array. It is a **joint target-set-and-sign null**, not a rung: see below.
* **`weight_permute`** — synapse counts are permuted across edges. Topology and transmitters are kept; the pairing between wiring and strength is destroyed.
* **`rewire_degree_preserving`** — double-edge swaps, roughly 10 × |E| accepted swaps. Every node keeps its in-degree, out-degree and transmitter; the wiring pattern is destroyed.
* **`erdos_renyi`** — the same node and edge counts, endpoints drawn uniformly, weights drawn from the empirical weight list. Only the size of the graph and its census survive.

**Transmitter-label null.** Plain label permutation preserves class counts but changes their weighted outgoing shares because node out-strength varies. In the `named` cut, the cholinergic outgoing-weight share is 0.619 on the real graph and 0.544 ± 0.027 across 1000 permutations (range 0.461–0.616), with the observed share at the 100th percentile of that null. It changes both the compound target set and the signed matrix. The matched variant constrains each tracked weighted share to within 0.002. The imidacloprid median point gap changes from 2.652 to 0.236 rate-model units against margin 0.332. These nulls still hold the transmitter sign magnitudes fixed.

Every mode preserves the node set, the node order and the seed block, and never creates self-loops or duplicate edges. The effect is always `treated − vehicle` with **both arms on the same graph**. The reported probability is the empirical two-sided permutation *p*, `(k+1)/(n+1)`, with resolution `1/(n+1)`; on a landscape the same `n` shuffled graphs are reused by every cell (a paired design), and shuffle *i* depends only on `(seed, i)`, so each cell's numbers are identical to running its profile alone.

**Why `z` is secondary.** A null distribution built on degraded graphs is not normal, and for the Erdős–Rényi mode its variance can collapse, producing numerically enormous *z* values that mean only "decisively different".

### S2.1 The three-way verdict and the equivalence margin

Each graph mode receives `distinguishable`, `equivalent_within_tolerance` or `indeterminate` under a margin δ set at 5% of the real graph’s vehicle readout. `equivalent_within_tolerance` requires non-rejection plus an observed point gap below δ; it is a descriptive rule, not a formal equivalence test based on a confidence interval. A zero or undefined baseline yields no usable margin and non-rejected modes are indeterminate.

**The confirmatory profile in full.** For imidacloprid at 1 µM on `named`, at 1000 permutations against a margin of 0.332 rate-model units, the per-mode probabilities and verdicts are: weight permutation *p* = 0.586, equivalent_within_tolerance; degree-preserving rewiring *p* = 0.472, equivalent_within_tolerance; the weight-matched transmitter null *p* = 0.389, equivalent_within_tolerance; the joint label-permutation null *p* = 0.275, indeterminate; and the Erdős–Rényi control *p* = 9.99e-04, distinguishable. The resulting necessary level is `degree_sequence` with verdict equivalent_within_tolerance.

The **necessary information level** is the weakest rung the test could not distinguish from the real cut, reported with its verdict. Ranks 2 and 3 are not strictly nested — `weight_permute` keeps transmitter identity and loses the weight-topology pairing, the transmitter nulls the reverse — so a richer model can be distinguishable where a weaker one is not. Such cells are **non-monotone** and are flagged: 0 of 84 on the `named` cut (none) and 1 of 84 on `taste_motor`. Their level must be read as "the cheapest graph model this test cannot tell apart from the real one", never as "everything above it is indistinguishable too".

### S2.2 Effect floor, multiplicity, and confirmatory versus exploratory

**Effect floor.** The absolute floor of 1e-6 rate-model units is below relevant baseline magnitudes, so the effective floor is `max(absolute, 1 % × |vehicle readout|)`. Landscapes report both. On `named`, the relative floor moves 15 cells to `no-effect` (acetylcholine 10 nM (0.000761 rate-model units), acetylcholine 100 nM (0.0118 rate-model units), ddt 10 nM (0.00309 rate-model units), ddt 100 nM (0.0309 rate-model units), gaba 100 nM (-2.74e-05 rate-model units), gaba 1 µM (-0.000864 rate-model units), gaba 10 µM (-0.025 rate-model units), ivermectin 10 nM (-0.0125 rate-model units), ivermectin 100 nM (-0.0312 rate-model units), ivermectin 1 µM (-0.0367 rate-model units), ivermectin 10 µM (-0.0373 rate-model units), picrotoxin 10 nM (0.0076 rate-model units), thiamethoxam 10 nM (0.000303 rate-model units), thiamethoxam 100 nM (0.00302 rate-model units), thiamethoxam 1 µM (0.0289 rate-model units)), changing counts from 0 / 59 / 25 to 0 / 44 / 40.

**Multiplicity.** The landscape is exploratory. Benjamini-Hochberg adjusts its declared structural component-test probabilities (168 on `named`, 168 on `taste_motor`), and classes are recomputed from those adjusted values. This does not imply an FDR guarantee for the derived cell labels. The report also states the empirical permutation resolution and whether any adjusted component test can cross α at that resolution (`fdr_can_reject` = yes).

**Planned versus exploratory comparisons.** The imidacloprid and fipronil profiles at *n* ≥ 1000 are planned comparisons reported without multiplicity adjustment. Landscape screening is exploratory. For fipronil at 1 µM, the planned component-test *p* is 0.006; in the landscape the adjusted component-test value is 0.062 and the derived class is composition-dominated versus topology-dependent before adjustment.

### S2.3 Does the instrument work? Ground-truth recovery and power

The planted-cycle synthetic test examines whether the procedure recovers one known wiring-dependent contrast. A synthetic cut of 250 nodes and about 1500 edges is built with log-normal weights and a transmitter composition close to the real cut’s. A directed cholinergic cycle through the driven seeds is planted at known strength; strength 0 is the negative control. This test does not estimate power for arbitrary circuit effects or validate the pharmacological parameter library.

the ladder recovered 2 of 3 planted topology-dependent effects at n = 200 and returned no false positive on the unplanted control. In a finite grid of 6 independently generated graphs per cell, planted strengths 0, 0.25, 0.5, 1 and permutation counts 50, 200, 1000, detection by strength is strength 0.25: 0.00-0.17; strength 0.5: 0.67-0.83; strength 1: 1.00-1.00 and by count n=50: 0.61, n=200: 0.61, n=1000: 0.67. The unplanted control has 0 detections in 18 runs (0.000). These frequencies are specific to this planted mechanism and grid.

**Permutation counts in this document.** 1000 per mode for the `named` profiles (resolution 0.0010), 300 for the `taste_motor` path arm, 1000 for each of the 84 `named` landscape cells and 300 for each of the 84 `taste_motor` cells. The sweep (T11) shows the verdict at α = 0.05 settling from *n* ≈ 50–400, and a mid-range *p* reaching ±0.02 of its final value between *n* = 25 and 400.

### S2.4 Where does the verdict settle with the size of the cut?

The descriptive dependence labels differ between the two committed cuts. A saved scale-study record contains a nested-cut ladder of 1k, 5k, 10k, 25k, 50k cells grown by a fixed ranked-BFS recipe. The 5k–50k input cuts are absent from this package, so their recorded results cannot be independently rebuilt here. The two committed cuts use another construction rule, and contrasts between them confound size, density, composition and extraction.

Cost is the constraint, not correctness. From the per-shuffle and per-edge constants measured on this machine (T26), a profile at *n* = 1000 on the 50k rung is 21.1 hours on one core and the landscape over the specification family is 7 days per specification; a whole-CNS profile at the weight-5 floor is 59 hours and the dense LIF matrix stops fitting in memory well before it. The permutation count must therefore fall with the cut, and `permutation_budget` refuses to go below the resolution floor `1/α − 1`, so every rung records the *n* it was run at.

The saved ladder's smallest rung, `scale_1k`, provides a descriptive comparison with the committed cuts of T19: 1000 nodes, 22 857 edges, mean degree 22.9, in-star no — fewer nodes than `named`, but a denser cut under a different construction rule.

**Imported result** (T27, from saved `scale_study.json`). Across 7 recorded cuts (scale_1k, named, taste_motor, scale_5k, scale_10k, scale_25k, scale_50k), the composition-dominated label occurs on 1 (named). Per-rung imidacloprid labels are scale_1k (1000 nodes, mean degree 22.9): topology-dependent; named (1126 nodes, mean degree 1.2): composition-dominated; taste_motor (1841 nodes, mean degree 10.4): topology-dependent; scale_5k (5000 nodes, mean degree 56.0): topology-dependent; scale_10k (10000 nodes, mean degree 62.2): topology-dependent; scale_25k (25000 nodes, mean degree 54.6): topology-dependent; scale_50k (50000 nodes, mean degree 44.3): topology-dependent; the record labels fipronil topology-dependent on each rung. The record summaries are imidacloprid: the verdict is 'topology-dependent' on the largest 5 of 7 cuts tested, first reached at 'taste_motor'. fipronil: the verdict is 'topology-dependent' on the largest 7 of 7 cuts tested, first reached at 'scale_1k'. and the inversion does not track the number of cells: `scale_1k` has 1000 nodes against `named`'s 1126 and a mean degree of 22.9 against 1.21, and it is already topology-dependent. It reports scale_5k, scale_10k cuts at which every tested null is distinguishable; the imidacloprid relative effect lies within 0.52-0.93 of the vehicle readout, and the clean-ladder summary is yes. These are conditional observations from the saved run; the absent larger input cuts prevent independent regeneration of all rungs.

**Scale-study limits.** Permutation counts decrease as edge count grows (20-1000 draws; resolutions 0.0010-0.0476), and at scale_10k, scale_25k, scale_50k the permutation budget falls to n = 20, where the smallest attainable probability is 0.048 against alpha = 0.05: a rejection there is the smallest the test can express. The `named` and `taste_motor` cuts are hops-limited whereas the `scale_*` ladder is ranked-BFS. No comparison here identifies recurrence, density or scale as an isolated cause of the class change.

## S3. Runtimes and their calibration

**Rate runtime.** `flylab/circuit/rate.py` constructs `W[j,i] = sign(nt_i) × synapses(i→j)` and row-normalises it as `W_signed = W / max(Σ_i |W[j,i]|, 1)` before applying presynaptic gains: `W_eff = W_signed · diag(g)`. The clipped iteration is `r ← clip((1 − α)·r + α·g_nav·(d + W_eff·r), 0, r_max)` with α = 0.3 and 80 steps. Its numerical readout is in **rate-model units**, not calibrated firing hertz; drive and clipping parameters are set on the model scale. This deterministic runtime takes about 50 ms for one 1126-node evaluation with single-threaded BLAS.

The normalisation is an **asserted modelling choice**, not a measurement, and it is consequential for RQ2: a readout in which each cell averages its presynaptic gains builds part of "composition-dominated" into itself. Three modes are implemented and swept in §3.3 (T25):

| mode | denominator | note |
|---|---|---|
| `row_abs` (shipped) | `max(Σ\|W\|, 1)` | denom[j] = max(sum_i |W[j, i]|, 1): each cell's recurrent input is a composition-weighted average of its presynaptic gains. Shipped default; every published FlyLab number uses it. |
| `none` | 1 | denom[j] = 1: raw signed synapse counts. The recurrent operator is supercritical (spectral radius ~150 on the named cut) and the r_max clip, not the wiring, bounds the rates. |
| `degree` | `max(in-degree, 1)` | denom[j] = max(in-degree[j], 1): degree-corrected. Relative input strengths survive (unlike row_abs) but hubs are discounted by fan-in. |

The composition-versus-full rank correlation at 1 µM is 0.921 under `row_abs`, 0.935 under `none` and 0.551 under `degree` (-0.583 at 10 nM): the composition verdict survives removing the normalisation but not replacing it with a degree-corrected one. Levels are not comparable across modes — the unnormalised operator is supercritical on this cut, so the `r_max` clip rather than the wiring bounds the rates — only the orderings are.

**LIF runtime.** `flylab/circuit/lif.py` is a current-based leaky integrate-and-fire network with the parameters of Shiu *et al.* [7]: τ_m = 20 ms, V_rest = −52 mV, V_th = −45 mV, V_reset = −52 mV, refractory 2.2 ms, exponential synaptic current with τ_s = 5 ms, forward Euler at dt = 0.1 ms. A presynaptic spike on edge *i*→*j* injects `w_scale × synapses × sign(nt_i) × gain(nt_i)` mV into *j*.

Two LIF choices are modelling decisions, not measured values. The synaptic scale `w_scale` = 0.05 mV per synapse and drive weight 0.5 mV per event were selected by a one-time scan, without fitting to recordings. Background Poisson drive at 0.65 × seed rate is supplied to non-seed cells to stand in for inputs omitted by the 1-hop cut. It sets the LIF operating point. The rate-model output is 0.293 rate-model units; the LIF output is 105 Hz under the standard 150 Hz seed drive. The values are not directly calibrated to one another.

## S4. Supporting capabilities

These are implemented, tested and used, but are supporting capabilities rather than contributions of this paper.

**Exposure.** `flylab/pharm/exposure.py` maps a dose and route onto *C(t)* and hence engagement(t) in one compartment. Its constants are class placeholders chosen so that the three routes differ in the direction the literature describes, and are fitted to no *Drosophila* ADME dataset because that dataset does not exist: of the seven quantities such a model needs, `data/literature/fly_pharmacokinetics.yaml` records only metabolite identity as well established [51], [53], [54], [57]. Adult haemolymph volume, adult body mass, cuticular penetration rate and any *Drosophila* elimination half-life are `null`; what is recorded is 25 nL recoverable per adult [56], 174 ± 81 nL in third-instar larvae [55], a water-partitioning measurement [52], an ethanol entry that could not supply a rate [64], and an imidacloprid half-life of 4.5–5 h in *honeybee* [57] — the only verified insect half-life for a library compound. Two structural warnings travel with every profile: haemolymph is the lowest-radioactivity compartment for imidacloprid, so a one-compartment haemolymph model under-estimates CNS exposure [58]; and a substantial metabolic route runs through the gut microbiome, so a single-enzyme clearance term is wrong for this compound [53].

**Genotype.** `flylab/pharm/genotype.py` applies published target-site resistance alleles as multiplicative potency shifts from `data/literature/resistance_alleles.yaml` [12]–[24], enforcing one rule: **a shift applies per compound, never per receptor.** Two published findings make the naive implementation wrong — *Rdl* A301S shifts GABA potency 33.9-fold but fipronil only 1.1-fold on the expressed planthopper receptor [14], and *para* M918T shifts deltamethrin 100-fold while leaving DDT at 1-fold [21]. An allele/compound pair with no sourced number is left unchanged and listed under `unshifted`, never extrapolated. One row is deliberately ugly: *Rdl* A301G in *D. simulans* carries a 20000-fold shift derived from a whole-animal resistance ratio, which drives fipronil engagement to 4.63e-04 and is almost certainly an overstatement of the target-site change [18]. The notebook says so. Panel: T9, F9.

**Mixtures.** `flylab/pharm/mixtures.py` implements Bliss independence [44] and Loewe additivity [47], [48] as explicit null models, following the response-surface framing of Greco *et al.* [45], and selects between them per pair: Loewe where the components share a receptor key, Bliss where their targets are disjoint. The three-case comparison against the honeybee dataset [50] is a qualitative concordance check and not validation. Table T5.

**Expression weighting.** `flylab/pharm/expression.py` maps class-level atlas statements [38]–[43] to MaleCNS superclasses by hand and interpolates gains. This is a sensitivity layer: the atlas labels do not match MaleCNS superclasses, every mapping is tagged `cross_atlas_inference`, coverage is 0.205 of cells, and identified MN9 motor-neuron expression remains unknown. The rate-model MN9 output shifts by -2.622 rate-model units at 1 µM imidacloprid. This is a model sensitivity, not an experimental correction.

**Model dose-response.** `flylab/assays/ensemble.py` builds the model's own concentration-response curve, fits a four-parameter Hill function by Nelder–Mead in log-concentration space and bootstraps it. For imidacloprid on the `named` cut the circuit midpoint is 34.7 nM with a bootstrap 95 % CI of 35.2 nM to 71.5 nM (r² = 0.9999). The interval does not bracket the point estimate, and that is reported rather than hidden: the bootstrap resamples replicates within each concentration, the lower plateau of the ladder is short, and the midpoint is therefore poorly identified. This is in any case the **model's** circuit midpoint, not a receptor or an animal IC50. The one-at-a-time tornado at 1 µM finds exactly zero span for ec50, hill_n — the engagement is saturated — which the global analysis (§3.5) confirms and quantifies.

## S5. Prospective predictions

The following 7 hypotheses are written before any experiment and labelled `software_prediction`; `live_result` is `null` in every row and stays null until a real table is imported by hand. They are **prospective, not pre-registered**: a mutable repository is not a registry, and the wording will change only if a tagged release is externally archived with a timestamp. Full table: T3.

* **H1** (receptor engagement). At 1 µM, imidacloprid insect nicotinic engagement ≥ 0.9 while vertebrate α4β2 stays ≤ 0.2. Predicted gap 0.854. Checkable against published expressed-receptor work rather than against flies.
* **H2** (PER, `protocols/per.md`). At the stated dose, the model predicts a change of -6.145 rate-model units in the MN9/DNp01 neighbourhood. The sign is specification-dependent: the shipped rule gives 6.631 to 0.458 rate-model units, while a monotone member gives 6.631 to 10.1; 15 of 25 specifications retain suppression. This is a model-output range for prospective testing, not a validated fly effect. The nearest sourced behavioural control concerns food-deprivation time [60], which differs from the model readout.
* **H3** (climbing, `protocols/climbing.md`). The fipronil arm predicts 0.896 rate-model units in the neighbourhood. Its comparison with H2 inherits H2’s specification dependence; 25 of 25 gain-rule specifications retain the fipronil predicate. Negative geotaxis provides a prospective assay [59], [63], [65], not current validation.
* **H4** (negative control). Diazepam engages vertebrate GABA-A and not insect RDL in this library, so it moves no fly readout: predicted effect 0.000. A nonzero biological effect would challenge the combined model assumptions; identifying its cause would require further experiments.
* **H5** (within-animal control). The extracted sensory path predicts a change of -0.296 rate-model units when bitter input accompanies sugar drive. Its direction is qualitatively concordant with the published sensory-veto result [7], [62]; the present cut and parameter choices are not an independent validation of that experiment.
* **H6**. Fipronil relieves the bitter veto on the map path: predicted rise in the veto ratio 0.413. The dependence analysis shows this readout is not distinguishable from its shuffled controls. A matching biological result would provide concordance with the predicted direction, but these null tests cannot determine whether the wiring contributes.
* **H7** (selectivity negative control). The model predicts a fly-neighbourhood change of 0.400 rate-model units for picrotoxin while its selectivity index is near 0. This checks the software’s selectivity interpretation; receptor-level cross-species data [30] do not validate the circuit output.

Rows where the model's own output disagrees with the prospective direction: none.

**On sample sizes.** Model-internal intervals (6 replicates of library jitter, drive jitter and the RNG seed) contain no biological variance, so effect sizes derived from them cannot support a biological power calculation. T3 therefore caps *d* at 1.0 before computing *n*, which yields 23 per group for every row. **This is an illustrative calculation, not a study power analysis**: it assumes a large standardised effect and independent flies, whereas vial-based protocols may violate independence; it must not be used to size a study.

## S6. Evidence typing: the full rules

The permitted transformation is a function of two facts — the **parameter type** the source measured and the **evidence distance** between that source and this compound at this receptor in this species. The distance grades are E0 (on-target); E1 (cross-species); E2 (related receptor); E3 (class extrapolation); E4 (unsupported). The full (parameter type × distance) → model table is generated as T20; its shape is:

| Parameter type | E0 | E1–E2 | E3 | E4 |
|---|---|---|---|---|
| `Kd`, `Ki` | `binding_occupancy` | `binding_engagement_proxy` | `functional_engagement_proxy` | `not_modelled` |
| `EC50`, `IC50`, `Kb` | `functional_engagement` | `functional_engagement_proxy` | `functional_engagement_proxy` | `not_modelled` |
| `relative_potency`, `class_order`, `unknown` | `not_modelled` | `not_modelled` | `not_modelled` | `not_modelled` |

A proxy computes the same Hill expression as its non-proxy sibling but reports itself as a transfer and carries a provenance warning naming what was transferred — the species, the preparation or the chemical class. Asking for a binding occupancy from an EC50, or for any numeric engagement from a placeholder row, raises `EvidenceTypeError`; the library validator rejects a placeholder that carries a value and a non-placeholder that does not. The library currently reports none validation problems.

**The admission invariant.** *No numeric engagement reaching an enumerated circuit-readout path may originate from a row whose (parameter type, evidence distance) pair does not admit one.* `tests/test_pharm_evidence.py` checks every shipped row and tested entry point, verifies that dropping inadmissible rows leaves the tested gain vectors unchanged, and detects a rule-table mutation. This tests software admission logic; it does not verify that every admitted number matches its cited assay.

**Census.** 102 rows over 21 compounds and 15 receptor keys (Table T10, Figure F11). Parameter types: 29 EC50, 8 IC50, 2 Kd, 63 unknown. Evidence distances: 7 E0, 15 E1, 8 E2, 9 E3, 63 E4. Engagement models: 1 `binding_occupancy`, 1 `binding_engagement_proxy`, 6 `functional_engagement`, 31 `functional_engagement_proxy`, 63 `not_modelled`. Rows entitled to a binding-occupancy model: imidacloprid at insect_nAChR_native_dmel (Kd 2e-09 M) [76]. Rows demoted to a binding proxy by distance: imidacloprid at insect_nAChR_beta1 (Kd 8.3e-11 M).

## S7. Reproduction details

| Knob | value in the run that produced this document | `--fast` |
|---|---|---|
| dependence permutations (`named` profiles) | 1000 | 20 |
| dependence permutations (`taste_motor` path) | 300 | 20 |
| landscape permutations per cell (`named`) | 1000 | 20 |
| landscape permutations per cell (`taste_motor`) | 300 | 20 |
| ladder recovery permutations | 200 | 20 |
| ladder power replicates per cell | 6 | 2 |
| transmitter-balance draws | 1000 | 20 |
| specification family | 25 | 9 (subsample) |
| shuffles per specification | 100 | 25 |
| composition reference draws | 50 | 10 |
| Sobol' base sample | 1024 (11 264 evaluations) | 64 |
| Ishigami validation base sample | 16 384 | 1024 |
| dose-response bootstrap | 200 resamples | 40 |
| prediction replicates | 6 | 2 |

The middle column is read from `results.json`, so it describes the run that rendered this document rather than an aspiration. `--fast` changes statistical effort only, never the model: point estimates that do not depend on replicate or permutation count are identical, and a `--fast` *p* is coarser. Two tests guard the committed record: it must not come from a `--fast` run, and every recorded statistical parameter must match the shipped default.

Two of the analyses are started in worker processes when the run begins and joined at their own step, because they dominate the wall clock and depend on nothing the other steps produce. `--jobs` controls that and the permutation workers; results are identical at any value.

## S8. Model-versus-implementation warnings carried in every notebook

Every assay returns a notebook (schema 0.3) whose `warnings` list ships with the artefact rather than in a console message: hops-limited neighbourhood, predicted (not measured) transmitter signs, teaching-tier potency values, reduced-circuit substitution where it applies, engagement-is-not-occupancy, a provenance warning naming the transfer for every proxy row, and the not-modelled reason for every N/A row. The `provenance` block records the FlyLab version, git commit, library SHA-256, map name and citation, RNG seed and platform. The library hash is what makes the honesty rules auditable: if a number changes, the hash changes, and a notebook produced before the change cannot be confused with one produced after.

## S9. Model-specific caveats in full

These are summarised in §4 and given here at length.

**(i) The transmitter sign table is mostly asserted.** Acetylcholine at +1.0 is the normalisation every other entry is scaled against, and GABA at −1.0 is the convention that inhibition is as strong per synapse as excitation. The remaining five — glutamate -0.4, histamine -0.5, dopamine +0.2, serotonin +0.2, octopamine +0.2 — are asserted magnitudes with no measurement behind them (T21). Glutamate in particular is signed inhibitory at 40 % of the stated strength, which is a hedge against the cells whose glutamatergic output is excitatory rather than a measurement of anything; any such cell is mis-signed in this model. None of these five is a member of the specification family, and the transmitter nulls do not bound them, because they permute *assignments* and leave the table untouched. Separately, the ablation ladder's level B inherits a whole-CNS excitation index that scores glutamate as partly excitatory while level D signs it inhibitory; the two levels of one ladder therefore disagree about a transmitter sign, and reconciling them moves the composition-versus-full rank correlation by at most 0.086.

**(ii) The row normalisation is load-bearing** (§S3, §3.3).

**(iii) Receptor effects are applied across transmitter-defined synapses.** A gain on `g_ach` scales every cholinergic synapse, whether or not that synapse expresses the receptor the compound binds. This is the assumption the expression layer relaxes, and it is reported as a bound rather than an improvement: coverage is 0.205 of cells, every mapping row is a cross-atlas inference, and receptor expression in identified adult motor neurons — MN9's class — is a confirmed gap in the atlas literature, so the paper's principal readout carries the default weight. Weighting shifts the rate-model MN9 readout by -2.622 rate-model units at 1 µM imidacloprid.

**(iv) Hops-limited cuts.** A 1-hop neighbourhood with a weight floor is a fraction of MN9's real input; missing inhibition can flip the sign of any readout, and the disinhibition readout is bounded by the GABA share that survived the cut. Growing the cut is deliberately constrained: past about 1 MB of JSON the weight floor is raised rather than hops added, because a 2-hop cut of a motor neuron in a 25-million-edge CNS is most of the animal. §3.2 is the quantitative demonstration of how much the choice of cut matters — the dependence verdict itself reverses between the two committed cuts.

**(v) The two cuts differ in what they can express.** The `named` cut carries 612 cholinergic, 295 GABAergic and 56 glutamatergic cells, and contains 0 octopaminergic cells and 0 on `taste_motor`, so a compound whose only insect mechanism acts through `g_oct` cannot move either readout for substrate reasons rather than pharmacological ones; 14% of `named` cells carry no consensus transmitter prediction and are presynaptically inert. Full structural census: T19.

**(vi) The rate-versus-LIF discrepancy** (§3.6) is agreement of direction only, and follows from the LIF background drive of §S3.

## S10. Domain readouts, selectivity and literature concordance

These are the domain-side results the method is demonstrated on. They are supporting material: none of them is a contribution of this paper, and §3 evaluates the instruments rather than these numbers.

**The evidence census in full.** Of 102 rows over 21 compounds and 15 receptor keys, 29 carry an EC50, 8 an IC50 and 2 a binding constant. By evidence distance the split is 7 / 15 / 8 / 9 / 63 rows at E0–E4. At the 1 µM working concentration the two headline compounds engage their insect targets at 0.991 and 0.985 against 0.091 and 0.476 at the vertebrate counterparts (500-fold and 37-fold; 2.70 log₁₀ units of window for the first). Placeholder rows contribute nothing: the negative-control row reports "not modelled (N/A)".

**Sensory-veto readout.** Driving 34 sweet GRN IDs at 150 Hz produces the published *direction* of bitter veto [7], [62] in both runtimes (vehicle ratios 0.000 rate model and 0.000 LIF; fipronil ratios 0.413 and 0.159). This is qualitative concordance only. Sugar-driven MN9 is 0.293 rate-model units versus 105 Hz in the LIF network; these numbers have different calibration status and are not cross-model validation.

**Path-null result.** At 300 permutations, the extracted-path veto ratio is not distinguishable from the shuffled controls (maximum |z| 0.220); none modes meet the descriptive point-gap criterion and erdos_renyi, rewire_degree_preserving, sign_permute, sign_permute_weight_matched, weight_permute are indeterminate. The planted-cycle study of §S2.3 does not establish power for this endpoint, so the result cannot rule out a wiring effect. The comparator ratio is undefined at this dose, not zero.

**Selectivity landscape.** The split is by mechanism rather than potency: ddt, deltamethrin are amplified by the circuit and acetamiprid, acetylcholine, clothianidin, imidacloprid, nicotine, nitenpyram buffered, subject to the threshold caveat of §3.4. Compounds never reaching the circuit threshold on either cut (chlordimeform, dieldrin, fipronil, gaba, ivermectin, permethrin, picrotoxin, thiamethoxam) are reported **unscored** rather than buffered, and one of them is unscored for a substrate reason rather than a pharmacological one: chlordimeform's only insect mechanism acts through octopamine, and neither cut contains an octopaminergic cell.

**Literature concordance** is reported per entry with its *n* (T4) and never as a bare mean: of 8 evaluable entry–assay combinations, 4 order exactly two compounds, where Spearman's ρ can only be ±1 whatever the model does, and over the 4 entries with more than two the mean ρ is 0.756. As a concordance count the result is 8 concordant against 0 discordant orderings (two-sided binomial *p* = 0.008); 16 combinations are skipped for want of coverage and 2 inversions are recorded rather than repaired (nitenpyram_rank_inverted, clothianidin_vs_imidacloprid_inverted). Of the 6 distinct published orderings the library can cover, 1 comes from a publication the library already cites for a compound in that ordering and is therefore *not* out-of-sample, while 5 show no detected overlap among recorded DOI/PMID identifiers; missing identifiers leave source independence unresolved. None is an independent test of the *circuit* model. The mixture arm is a **limited qualitative concordance check** (3 of 3 verdicts reproduced, reference organism Apis mellifera, 2 partners same-class substitutes, §S4).
