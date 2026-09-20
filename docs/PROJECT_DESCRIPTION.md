# Project description and repository metadata

This file keeps the short public descriptions of FlyLab in one place so the repository, paper, website and package metadata can stay aligned.

## Recommended GitHub About description

**Provenance-first computational pharmacology on the Drosophila MaleCNS connectome, with typed evidence, circuit simulation and connectome-dependence testing.**

Shorter alternative:

**Typed pharmacology + Drosophila connectomics + model-dependence testing.**

## Recommended homepage

https://pinkysworld.github.io/FlyLab/

## Recommended repository topics

- drosophila
- connectomics
- computational-neuroscience
- computational-pharmacology
- scientific-software
- reproducible-research
- systems-biology
- pyodide
- fastapi
- uncertainty-quantification
- sensitivity-analysis
- null-models

## One-sentence description

FlyLab is a provenance-first computational pharmacology workbench that maps typed receptor evidence onto named circuits in the adult Drosophila MaleCNS connectome and tests how strongly each model prediction depends on wiring, assumptions and evidence quality.

## Short project summary

FlyLab connects pharmacology and connectomics without hiding the assumptions between them. A compound and free concentration are translated through a typed, sourced receptor library into mechanism-level perturbations of named MaleCNS circuits. The resulting model outputs can be inspected with deterministic and spiking runtimes, then challenged with degraded-graph null models, ablation, specification robustness and global uncertainty analysis.

The project is designed to make negative conclusions visible. A prediction may turn out not to need detailed wiring, a qualitative effect may reverse under another admissible gain rule, or an uncertainty analysis may show that a modelling assumption matters more than a potency value. Those outcomes are treated as useful results rather than failures of the software.

FlyLab is a computational methods project. It does not generate live-animal data, does not model vertebrate circuits, and does not claim biological efficacy or safety from model output alone.

## Academic description

FlyLab is an open-source framework for provenance-preserving pharmacological perturbation of synapse-resolution Drosophila circuits. It combines a typed receptor-evidence layer, compound- and concentration-specific mechanism transforms, MaleCNS-derived circuit models, and a suite of model-auditing analyses. These analyses include graph-degradation null models for connectome dependence, layer ablation, specification robustness, variance-based global sensitivity analysis and model-oriented value of information.

The framework separates literature-derived inputs, measured connectome structure, modelling assumptions and computed outputs. This separation is carried into exported notebook provenance and the paper reproduction pipeline. The current repository evaluates the computational method and model behaviour; circuit-level biological validation remains future work.

## What FlyLab should not be described as

Avoid descriptions such as:

- "the first virtual Drosophila pharmacology lab"
- "the first connectome pharmacology simulator"
- "a digital twin of Drosophila"
- "a predictor of insecticide safety"
- "a validated virtual fly"
- "a replacement for animal experiments"

Those descriptions overstate either novelty, biological realism or validation.

## Preferred public language

Good:

> FlyLab tests pharmacological perturbations on named Drosophila circuits and audits whether the resulting model predictions actually depend on detailed wiring.

Good:

> FlyLab combines typed receptor evidence with MaleCNS-derived circuits and makes model assumptions inspectable.

Good:

> FlyLab is a computational hypothesis-triage and model-auditing tool.

Too strong:

> FlyLab predicts what a drug will do to a fly.

Too strong:

> FlyLab proves which drugs require the connectome.

## Versioning note

The public project description should describe the stable architecture, not a single paper result. Numerical findings belong in papers/results.json, the manuscript, and research documentation rather than the GitHub About field.
