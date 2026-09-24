# Research positioning

This document is the claim guardrail for FlyLab.

Its purpose is not to market the project. It exists to keep the paper, README, website and future PRs aligned with what the software actually contributes and with adjacent published work.

## The narrow contribution

**Frame the work as a computing contribution.** The problem FlyLab addresses is general: a simulation built on a measured network and parameterised from published literature cannot say which of its inputs its predictions depend on. The four instruments — typed evidence propagation, input-dependence testing, the ablation ladder with its matched reference, and specification-family robustness with variance attribution — are the contribution, and none is specific to flies. The fly pharmacology is the domain the method is demonstrated on, and it supplies the hard case rather than the claim.

FlyLab should not be presented as "connectomics plus pharmacology is new."

Executable Drosophila circuits already exist. Connectome-based LIF models already exist. Receptor-informed pharmacological perturbation of whole-brain network models already exists.

FlyLab's contribution is narrower:

1. **Typed pharmacological evidence**  
   The software records what kind of parameter a source reported, where the source sits relative to the modelled receptor and species, and whether a numeric transformation is supported at all.

2. **Compound- and concentration-specific perturbation of named, synapse-resolution Drosophila circuits**  
   The insect circuit and vertebrate receptor scorecard are evaluated at the same supplied free concentration, while only the insect side can alter a circuit.

3. **Connectome-dependence analysis**  
   A real-graph drug contrast is compared with distributions produced by degraded graph models that retain different levels of network information.

4. **Model-auditing analyses**  
   Ablation, specification robustness, global uncertainty and claim provenance make it possible to determine whether a conclusion comes from measured structure, sourced pharmacology, an asserted model choice, or their interaction.

5. **Reproducible computational artifacts**  
   The scientific core, notebooks and paper pipeline retain version, library, map and random-seed provenance.

The transferable idea is not simply "simulate a drug on a connectome." It is:

> **Make every connectome-based pharmacology prediction explain what evidence and network information it actually depends on.**

## Adjacent work that must be acknowledged

### Executable fly circuits

FlyBrainLab already provides an interactive environment for constructing and comparing executable Drosophila circuits from connectome data.

Do not claim FlyLab is the first software to simulate named fly circuits.

### Connectome LIF models

Published whole-brain Drosophila LIF work demonstrates connectome-driven sensorimotor simulation and provides an important directional control for FlyLab.

Reproducing an existing sensorimotor direction is not a FlyLab novelty claim.

### Receptor-informed whole-brain pharmacology

Human whole-brain models have already used receptor maps to simulate pharmacological interventions.

Do not write that FlyLab is the first project to combine network structure and receptor pharmacology.

The relevant distinction is scale, evidence typing, named synapse-resolution circuits, cross-species scorecards, and explicit dependence testing.

## How to describe connectome dependence

The dependence framework uses degraded graph ensembles and empirical permutation probabilities.

flylab/analysis/dependence.py reports one of three verdicts per null mode, and the safe interpretation follows them:

- a small empirical p means the real-graph effect is **distinguishable** from that null ensemble at the chosen permutation effort,
- a non-rejection whose gap from the null ensemble's median is below the prespecified margin is reported as **equivalent within tolerance**; the margin is declared in advance, not fitted, and defaults to a small fraction of the vehicle readout,
- any other non-rejection is **indeterminate**: consistent with a difference this study cannot resolve, and not evidence of equality or biological equivalence.

Accordingly, prefer:

> "The real-graph effect was not distinguishable from the degree-preserving null ensemble."

over:

> "The degree-preserving graph reproduces the effect."

Only the equivalent-within-tolerance verdict licenses the word "reproduces", and even then it means "within the prespecified margin", never "the shuffle gave the same effect". Check the verdict against a run rather than assuming: for imidacloprid at the headline dose and n = 1000, the weight, degree-preserving and weight-matched transmitter modes fall inside the margin and report **equivalent_within_tolerance**, while the plain transmitter permutation reports **indeterminate** (gap 2.65 Hz against a margin of 0.33 Hz) and the Erdos-Renyi control is **distinguishable**.

### The verdict depends on the extract; the saved scale ladder points to recurrence

The single most important thing to know before writing any dependence sentence: **the same analysis, with the same instrument and no parameter changed, gives opposite answers on different extracts of one connectome.** On the sparse 1-hop `named` cut — an in-star, most of whose edges terminate on four seed cells and most of whose nodes receive no input — most cells come back composition-dominated. On the denser `taste_motor` cut, none does.

The saved scaling study (`flylab.analysis.scale.dependence_vs_scale`, T27)
contains seven extracts. The composition-dominated verdict occurs only on the
sparse `named` in-star; `scale_1k` has fewer nodes than `named` (1000 against
1126), a mean degree of 22.9 instead of 1.21, and is already topology-dependent.
This is consistent with recurrence tracking the difference, but it does not
isolate a causal feature: the extracts differ in more than node count. The
larger `scale_5k`–`scale_50k` inputs are absent from the current checkout, so
T27 cannot be independently rebuilt here.

At `scale_5k` and `scale_10k`, the imidacloprid effect is distinguishable from
all three ranked edge/wiring nulls (`weight_permute`,
`rewire_degree_preserving`, `erdos_renyi`). The plain transmitter-label null is
indeterminate at `scale_10k`; it is an orthogonal check outside the ranked
necessary-information ladder. At `scale_25k` and `scale_50k`, the three ranked
edge/wiring nulls remain distinguishable, while the two transmitter-label
checks are not both distinguishable. The 25k and 50k runs use n = 20, where the
minimum p is 0.048 against α = 0.05, so those rejections are the smallest the
test can express and confirm rather than establish the pattern.

1. **Never state a dependence class without the extract it was measured on**, and quote `flylab.analysis.dependence.cut_census` (mean degree, share of edges onto seeds, share of nodes with any input, recurrence budget) beside it.
2. **"Most predictions do not need the connectome" is withdrawn**, along with "the topology-dependent set is exactly the chloride-channel blockers". What replaces them is a methodological claim: the dependence verdict is a joint property of the prediction and the extract; the saved extracts associate the inversion more with recurrence than node count, but do not show that recurrence alone causes it.
3. **Do not quote the top of the ladder as strong evidence.** The 25k and 50k rungs run at n = 20. Raising n requires restoring the exact cuts and validating any faster degree-preserving sampler against the same graph invariants and permutation procedure.

The planted-cycle checks provide limited evidence that the ladder recovers that one synthetic mechanism; they do not establish general power for biological endpoints. The `named` versus `taste_motor` landscape uses the same n = 1000 effort, so its contrast is conditional on those two extracts and settings. The T27 scale rungs use different, often smaller permutation counts; that trend is descriptive and cannot rule out power differences.

## Current headline comparisons

### Imidacloprid

At the current headline concentration and readout, the imidacloprid real-graph effect is not distinguishable from several structure-preserving null ensembles at n = 1000, while it differs strongly from the Erdős-Rényi control. Those non-rejections carry either **equivalent_within_tolerance** or **indeterminate** depending on how far the effect sits from the null median, so quote the verdict the run returned rather than assuming the weaker one. Both are weaker than equality, and neither is evidence about the network the extract was cut from: this verdict is measured on a 1126-cell in-star and reverses on extracts with recurrence.

The current interpretation is that detailed wiring is not supported as necessary for this particular model output, which is weaker than saying it is unnecessary.

### Fipronil

At the same headline concentration and readout, the fipronil effect differs from the weight-permuted and degree-preserving rewired ensembles at n = 1000.

That is evidence that this output is more sensitive to wiring structure than the imidacloprid comparison.

### Landscape

The 21-compound by 4-concentration landscape is useful as an exploratory screen, but it contains many comparisons.

Do **not** summarize it as:

> "Exactly the chloride-channel blockers are topology-dependent."

That statement is biologically inaccurate and too categorical for the current exploratory classification. The topology-sensitive set is chlorpyrifos_oxon, dieldrin, fipronil, gaba, ivermectin and picrotoxin, and three of those are not chloride-channel blockers: gaba is an RDL agonist, ivermectin a GluCl opener, and chlorpyrifos_oxon an AChE inhibitor that appears at one concentration only.

Multiplicity is handled rather than pending. The landscape applies Benjamini-Hochberg across its structural tests, records p_adjusted, q_value and fdr_alpha per cell, classifies on the adjusted values while keeping the raw ones beside them, and labels the run **confirmatory** or **exploratory** according to whether it was prespecified and whether its permutation resolution can support a rejection at the adjusted threshold at all. Any count taken from the landscape must be quoted with that label.

Describe the landscape as:

> "Topology-sensitive classifications are concentrated among perturbations affecting inhibitory signalling, with mechanism and concentration-specific exceptions."

The exact counts and classifications should be taken from the generated `papers/results.json` and identified with the shipped defaults. T27 is an imported scale-study record; its 5k–50k input cuts are not in this checkout.

## Specification robustness

The gain-rule family is one of FlyLab's most important self-audits.

The current analysis demonstrates that some qualitative conclusions are sensitive to the engagement-to-gain functional form. In particular, a nicotinic suppression conclusion reverses under monotone alternatives.

That is a useful result because it identifies a model choice that should not be mistaken for a connectome-derived biological conclusion.

However, distinguish two things:

1. **Qualitative circuit conclusions re-derived under alternative gain rules**
2. **Statistical dependence conclusions derived from permutation nulls**

If a robustness predicate uses a small null sample or z-score shortcut, do not present it as equivalent to the full n = 1000 permutation analysis. The manuscript and README should reserve strong statistical language for analyses that actually use the corresponding permutation evidence.

## Withdrawn ablation claims

Two ablation readings were withdrawn after review. Both still appear elsewhere in the repository and in earlier manuscript text, so neither may be restated without the correction below.

### "The connectome without the pharmacology is not a cheap substitute for the connectome with it"

Withdrawn. The claim rested on the topology-only level of the ablation ladder, which applied a depression-only generic multiplier to every transmitter alike. A baseline that can only reduce activity cannot order a library containing compounds whose full-model effect is positive, so its poor score was a property of the rule, not of the connectome. Under the direction-aware rule now shipped as the default (GENERIC_RULES in flylab/analysis/baselines.py), which keeps a generic magnitude and takes one sign bit from the mechanism table, that level's rank correlation with the full model rises from about 0.24 to about 0.96. The old rule survives only as a floor, reported beside the fair rule, and is never a competitor.

### A bare rank correlation of about 0.99 between composition-only and the full model

Withdrawn as stated. The two levels are functions of the same gain vector: composition-only applies it to a fixed census, the full model applies it per cell on a row-normalised matrix in which each cell's recurrent input is already a composition-weighted average of its presynaptic gains. A high correlation is therefore expected for algebraic reasons. composition_reference_distribution shows that pharmacology-free pseudo-compounds, drawn as arbitrary gain vectors, already reach a median rank correlation around 0.8, so the observed value is the structural floor plus a small excess. composition_dominance_under_normalisations shows that the agreement does not survive a degree-corrected engine: it falls to about 0.67 at 1e-6 M and inverts, to about -0.37, at 1e-8 M.

Never quote the correlation without its matched reference distribution and the engine normalisation it was computed under.

## Typed evidence

FlyLab's evidence model is a methodological contribution because unsupported rows remain unsupported.

Important distinctions:

- EC50 and IC50 are functional potency parameters.
- Kd and Ki are binding parameters.
- source relation remains relevant even when the parameter type is strong.
- a binding parameter measured in another species is not automatically a direct measurement of the modelled Drosophila receptor state.
- missing evidence must not become a small numerical effect.

Public descriptions should usually use the umbrella term **engagement** and allow detailed outputs to expose the more specific evidence type.

## Uncertainty and value of information

Sobol and VOI outputs describe variance of a **model output under declared input ranges**.

They are not estimates of biological variability.

Finite-sample Sobol first-order estimates can be slightly negative around zero. Negative raw estimates may be retained for statistical transparency, but they must not be described as negative physical "variance removed" or negative experimental value.

The useful qualitative question is:

> Which uncertain assumption is currently responsible for the largest resolvable share of model variance?

## Literature concordance

Use **literature concordance**, not **validation**, when the comparison:

- reuses a source represented in the pharmacology library,
- evaluates receptor-level ordering rather than circuit output,
- uses another species,
- or uses a coarse qualitative endpoint.

Reserve "independent out-of-sample validation" for data that did not enter the library, model construction or calibration.

The current circuit model has no such independent biological validation.

## Prospective predictions

The software predictions are **prospective**, not pre-registered, while they exist only in a mutable repository.

An archived and timestamped release can establish a fixed prediction record for future experiments, but repository history alone should not be described as a formal registry.

## Supporting capabilities, not primary novelty

The following are useful capabilities but should not carry equal novelty weight:

- exposure C(t),
- mixtures,
- genotype shifts,
- expression weighting,
- the graph viewer,
- a Hill fit to the model's own dose-response,
- reproduction of a previously published bitter-veto direction.

## Claims FlyLab should not make

Do not describe FlyLab as:

- a digital twin of a fly,
- a validated predictor of in vivo drug effects,
- a vertebrate safety model,
- a replacement for regulatory toxicology,
- proof that a given compound "needs the connectome",
- proof that non-significant null comparisons are equivalent,
- the first connectome pharmacology simulator,
- the first executable Drosophila circuit environment,
- a demonstration that the connectome without the pharmacology is not a cheap substitute for the connectome with it,
- a composition-versus-topology verdict resting on a bare rank correlation, without its reference distribution and the engine normalisation used.

## Claims FlyLab can make

FlyLab can accurately be described as:

- a provenance-first computational pharmacology workbench,
- a framework for applying typed receptor evidence to named MaleCNS-derived circuits,
- a system for testing sensitivity of model predictions to degraded connectome information,
- a reproducible model-auditing framework,
- a tool for computational hypothesis triage,
- a research object that records when its own conclusions are assumption-dependent.

## One preferred novelty sentence

For the paper or a project summary:

> FlyLab combines typed pharmacological evidence with compound-specific perturbation of named, synapse-resolution Drosophila circuits and adds a validated model-auditing layer that tests which network information and modelling assumptions each computed conclusion actually depends on — and which shows that the answer depends on the cut of the connectome it is asked about.

That sentence is deliberately narrower than the software's full feature list.
