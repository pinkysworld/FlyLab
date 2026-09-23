# Current manuscript: Round 1 independent statistical review

**Snapshot reviewed:** repository `b82b8ed`; manuscript/results provenance `f95fd88`. This is a fresh review of the current source and artifacts. Earlier findings that the revision has fixed are not repeated.

**Recommendation: major revision.** The planted-effect calibration, relative effect floor, FDR layer, weight-matched transmitter null, and explicit scale study are substantive improvements. The current headline still overreaches in three places: the code's within-margin rule is not a valid equivalence test; the scale study does not identify recurrence as the cause of the verdict; and the synthetic power surface is too small to support the general design rules claimed for it.

## Blocking

### B1. `equivalent_within_tolerance` is a point-estimate heuristic, not statistical evidence of equivalence

**Repository evidence.** `flylab/analysis/dependence.py:697-759` assigns `equivalent_within_tolerance` when two conditions hold: an ordinary difference-style permutation test does not reject, and the point gap `|real_effect - median(null)|` is below δ. No uncertainty interval is formed for the null median or for the gap. With this rule, even a very small shuffle sample can be called equivalent whenever its sample median happens to lie near the real effect. In landscapes, `_decision_p` uses the BH-adjusted probability (`dependence.py:790-804`), so multiplicity adjustment can change a raw rejection into non-rejection and thereby help trigger an “equivalent” label. FDR correction of a difference test cannot supply evidence for equivalence.

The rule also makes difference and equivalence mutually exclusive, although they answer different questions and can both hold. T6 supplies a live example: fipronil versus degree-preserving rewiring has `p = 0.00699`, gap `0.243 Hz`, and δ `0.332 Hz`; it is within the stated practical margin but the code labels it `distinguishable`. The same occurs for weight permutation (`p = 0.00599`, gap `0.310 < 0.332`). If δ defines practical sameness, a small difference-test p should not automatically overrule it.

**Implication.** The abstract's “prespecified equivalence margin,” §2.3's “same effect to within the margin,” 281 landscape equivalence verdicts, and the 17/25 C3 equivalence count state more than the procedure establishes. The 5% margin is also a modelling convention justified by another model threshold, not an empirically grounded smallest important difference; calling it “prespecified” should mean fixed before examining these results, not merely encoded before the latest run.

**Required fix.** Choose the estimand first: equivalence of the real effect to the null ensemble's median, or a distributional statement about null realizations. For median equivalence, require a confidence interval for the null median (bootstrap or order-statistic interval) to lie wholly inside `[real_effect - δ, real_effect + δ]`; an equivalent two-one-sided construction is also acceptable. For a distributional claim, prespecify the required mass inside the margin and estimate it with uncertainty. Report difference evidence and equivalence evidence on separate axes, and correct the equivalence family directly if many equivalence claims are made. Until then rename the current output `within_margin_point_estimate` and do not call it equivalence.

### B2. The data rule out a simple node-count explanation; they do not establish “recurrence, not size”

**Repository evidence.** The clean `scale_*` ladder changes node set, edge count, mean degree, transmitter composition, seed share, vehicle rate, and recurrence together. Every clean rung is topology-dependent, so there is no within-recipe verdict variation from which to identify a structural predictor. The only composition-dominated cut is `named`, built by a different extraction recipe. `flylab/analysis/scale.py:602-695` correctly describes the structure comparison as observational and says that a minority class of one cannot distinguish candidate statistics; lines 615-619 explicitly require a controlled cut with the same node set and altered edge structure. `papers/scale_study.json` contains no such intervention.

The manuscript's key size counterexample is especially non-diagnostic. At `scale_1k`, imidacloprid is called topology-dependent because **weight permutation** rejects (`p = 0.00670`), while degree-preserving rewiring, the mode that destroys who-connects-to-whom while retaining degrees, does not (`p = 0.493`) and is labelled within tolerance. Its necessary level is paradoxically `size_and_composition` because the ladder is non-monotone. Thus `scale_1k` does not show that destroying recurrence changes the answer; it shows sensitivity to the weight-topology pairing.

**Implication.** The repeated headline “recurrence, not size, predicts which way” in the abstract, §3.2, Discussion, and Conclusion is causal/selective language unsupported by the study. The evidence shows that node count alone is insufficient and that the verdict differs across structurally and procedurally different extracts. It does not select recurrence over density, composition, degree distribution, weight assignment, seed coverage, or extraction recipe.

**Required fix.** Replace the headline with: “the verdict changes across extracts and does not follow node count monotonically.” Describe recurrence as a candidate associated feature. To retain the stronger claim, run the controlled experiment already specified by `scale.py`: same node set, transmitter labels, degree/weight distributions, readout, and permutation budget, with recurrence selectively broken or varied. A factorial synthetic study varying size and recurrence independently would also support the methodological claim.

### B3. The top-of-ladder verdict has not statistically “settled”

**Repository evidence.** `scale_25k` and `scale_50k` each use only 20 permutations. For imidacloprid, both structural probabilities equal the smallest possible value, `1/21 = 0.047619`, and every mode has `stabilised = false` in `papers/scale_study.json`. Nevertheless `verdict_stability` calls the sequence settled whenever the last two categorical labels agree (`scale.py:475-580`). Its `under_powered_cuts` check only detects a resolution *coarser* than alpha; a floor barely below 0.05 passes. Agreement between two thresholded results at the only rejectable value is not evidence that a trend has stabilized.

**Required fix.** Call the upper rungs exploratory and resolution-limited, and replace “settled” with “the observed class was unchanged across the tested clean rungs.” Do not use the two `n = 20` rungs to establish stability. A settled claim needs a prespecified precision or convergence criterion met on at least two upper rungs, ideally with a common adequate permutation budget or sequential Monte Carlo bounds around the p-values.

## Major

### M1. The synthetic experiment is a useful positive control, but the claimed power and false-positive performance are too precise

**Repository evidence.** T22 uses six independently generated graphs per strength-by-*n* cell. Consequently every detection estimate moves in steps of 0.167. Illustrative 95% Wilson intervals are broad: `4/6 = 0.667` gives `[0.300, 0.903]`, `5/6 = 0.833` gives `[0.436, 0.970]`, and `6/6` gives `[0.610, 1.000]`. The control has 0 detections in 18 runs, whose 95% Wilson upper bound is 0.176; it does not estimate the false-positive rate as zero with useful precision. Reusing the same six graph seeds across permutation counts is a sensible paired design, but those cells are not independent replications.

The control is also not identical apart from loop strength: `synthetic_cut` adds ten cycle edges and relabels the ring nodes cholinergic only when strength is positive (`dependence.py:1577-1597`). That is a valid planted positive control, but it varies degree, composition, weight allocation, and recurrence together. It tests one motif and one gain mechanism, not general instrument validity.

**Required fix.** Report binomial counts and confidence intervals in T22/F16. Replace “validated,” “full power,” “effect size dominates,” and “permutations past a few hundred are wasted” with “positive-control calibration” and the narrower observations from this grid. Increase graph replicates substantially before claiming a power surface or false-positive rate. Use a matched control containing the same nodes, labels, degrees, and cycle-edge weights but with the loop selectively broken; then vary graph size and loop strength independently.

### M2. BH validity under the observed dependence has not been demonstrated

**Repository evidence.** The supplement asserts that sharing one shuffle stream makes the structural tests positively dependent and therefore licenses Benjamini-Hochberg. A common random stream does not by itself establish positive regression dependence: compound effects have different signs and nonlinear gain rules, so test statistics can be positively or negatively associated. No dependence diagnostic or proof is supplied. The same concern applies to the 100 structural tests in the specification run.

**Required fix.** Either establish the dependence condition required by BH for these statistics, use Benjamini-Yekutieli for arbitrary dependence, or calibrate family error/FDR empirically with joint null simulations. State the exact family once: landscapes currently adjust two structural modes per cell, while other reported ladder modes remain descriptive.

### M3. The scale result is re-renderable but not independently recomputable from the repository

**Repository evidence.** Only `malecns_scale_1k.json` is committed. The 5k-50k paths recorded in `papers/scale_study.json` point to `/root/.flylab/...` and are absent in this checkout. The study JSON has no code SHA, cut hashes, library hash, generation timestamp, or per-cut seed/provenance block. `scripts/reproduce_paper.py:4622-4635, 4717-4728` explicitly reads this cached result rather than rerunning the study. Current manuscript/results provenance is `f95fd88`, while the reviewed source is `b82b8ed`.

**Implication.** One command can re-render T27 from a committed cache; it cannot regenerate the scale analysis from its source graphs. “Every figure, table and number” and RQ4 need to distinguish artifact replay from analysis reproduction.

**Required fix.** Add a schema/provenance block to `scale_study.json` with code, library, map, cut hashes, analysis parameters, and seeds; publish the large cuts as durable artifacts or provide a verified workflow that recreates them. Have the paper pipeline validate the cache against those hashes and the checked-out release. Rephrase the current claim as “all ordinary analyses rerun; the scale table is deterministically re-rendered from a committed, provenance-checked study record” until the external cuts are available.

## Statistical language that can remain

The following current claims are supported if kept scoped: the verdict differs between `named` and `taste_motor`; node count alone does not order the observed verdicts; the scale ladder stays topology-dependent over the tested rungs; the instrument detects the strong planted loop in this synthetic design and is insensitive to its weak version; and a negative dependence result must be reported with the extract and runtime on which it was measured.
