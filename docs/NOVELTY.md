# Research positioning

This document is the claim guardrail for FlyLab.

Its purpose is not to market the project. It exists to keep the paper, README, website and future PRs aligned with what the software actually contributes and with adjacent published work.

## The narrow contribution

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

The safe interpretation is:

- a small empirical p means the real-graph effect is distinguishable from that null ensemble at the chosen permutation effort,
- a large p means the analysis did not distinguish the real effect from that null ensemble,
- a large p is not proof of equality or biological equivalence.

Accordingly, prefer:

> "The real-graph effect was not distinguishable from the degree-preserving null ensemble."

over:

> "The degree-preserving graph reproduces the effect."

The latter is stronger than a non-significant permutation comparison establishes unless an explicit equivalence criterion is added.

## Current headline comparisons

### Imidacloprid

At the current headline concentration and readout, the imidacloprid real-graph effect is not distinguishable from several structure-preserving null ensembles at n = 1000, while it differs strongly from the Erdős-Rényi control.

The current interpretation is that detailed wiring is not supported as necessary for this particular model output.

### Fipronil

At the same headline concentration and readout, the fipronil effect differs from the weight-permuted and degree-preserving rewired ensembles at n = 1000.

That is evidence that this output is more sensitive to wiring structure than the imidacloprid comparison.

### Landscape

The 21-compound by 4-concentration landscape is useful as an exploratory screen, but it contains many comparisons.

Do **not** summarize it as:

> "Exactly the chloride-channel blockers are topology-dependent."

That statement is biologically inaccurate and too categorical for the current exploratory classification. The current set includes perturbations with different pharmacological actions and at least one low-concentration AChE case.

Until multiplicity handling or a clearly declared exploratory framing is finalized, describe the landscape as:

> "Topology-sensitive classifications are concentrated among perturbations affecting inhibitory signalling, with mechanism and concentration-specific exceptions."

The exact counts and classifications should be taken from papers/results.json, not copied by hand into durable project descriptions.

## Specification robustness

The gain-rule family is one of FlyLab's most important self-audits.

The current analysis demonstrates that some qualitative conclusions are sensitive to the engagement-to-gain functional form. In particular, a nicotinic suppression conclusion reverses under monotone alternatives.

That is a useful result because it identifies a model choice that should not be mistaken for a connectome-derived biological conclusion.

However, distinguish two things:

1. **Qualitative circuit conclusions re-derived under alternative gain rules**
2. **Statistical dependence conclusions derived from permutation nulls**

If a robustness predicate uses a small null sample or z-score shortcut, do not present it as equivalent to the full n = 1000 permutation analysis. The manuscript and README should reserve strong statistical language for analyses that actually use the corresponding permutation evidence.

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
- the first executable Drosophila circuit environment.

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

> FlyLab combines typed pharmacological evidence with compound-specific perturbation of named, synapse-resolution Drosophila circuits and adds a model-auditing layer that tests which network information and modelling assumptions each computed conclusion actually depends on.

That sentence is deliberately narrower than the software's full feature list.
