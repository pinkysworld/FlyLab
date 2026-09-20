# FlyLab documentation

This directory contains the technical and research documentation for FlyLab v0.6.

If you are new to the project, start with the main [README](../README.md). This page explains where the deeper documentation lives and which files are authoritative.

## Start here

| Document | Use it for |
|---|---|
| [README](../README.md) | public project overview, installation and first commands |
| [Architecture](ARCHITECTURE.md) | scientific dataflow, package boundaries and provenance model |
| [Research positioning](NOVELTY.md) | what FlyLab claims as new, what adjacent work already exists, and claim-language guardrails |
| [Research map](RESEARCH_MAP.md) | current analyses, generated results and open research gaps |
| [Browser build](PAGES.md) | Pyodide bridge, static site and GitHub Pages details |
| [Software stack](STACK.md) | dependencies, engineering choices and platform constraints |
| [Paper status](../papers/STATUS.md) | current manuscript state and what remains before submission |
| [Submission checklist](../papers/SUBMISSION_CHECKLIST.md) | paper-specific mechanical and scientific checks |
| [Contributor handoff](../HANDOFF.md) | internal constraints for contributors and coding agents |

## What is authoritative

FlyLab has evolved quickly, so not every design document has equal status.

### Current

These describe the current v0.6 project:

1. [README](../README.md)
2. [ARCHITECTURE.md](ARCHITECTURE.md)
3. [RESEARCH_MAP.md](RESEARCH_MAP.md)
4. [NOVELTY.md](NOVELTY.md)
5. [papers/STATUS.md](../papers/STATUS.md)
6. the source code and generated papers/results.json record

### Historical

[DESIGN_v0.5.md](DESIGN_v0.5.md) is intentionally retained as the v0.5 interface contract. It contains useful historical context and endpoint names, but it is not authoritative for the current typed-evidence model or the current analysis layer.

When a historical document disagrees with the current code or architecture document, the current code and v0.6 documentation win.

## Scientific model in one page

FlyLab separates four kinds of information that are easy to blur together in simulation work:

1. **Published evidence**  
   Receptor-level parameters and literature-derived supporting datasets.

2. **Model assumptions**  
   Engagement-to-gain transformations, circuit cuts, transmitter signs, external drive and other modelling choices.

3. **Measured connectome structure**  
   MaleCNS-derived edges and annotations used by the committed research cuts.

4. **Computed results**  
   Circuit readouts, dependence profiles, ablation results, specification robustness, uncertainty budgets and notebooks.

That separation is the core of the project. A result should be traceable from a reported number back through every layer that produced it.

## Main research instruments

### Typed evidence

The library records the type and provenance of receptor evidence rather than treating every number as an EC50.

The current schema distinguishes binding constants, functional potency parameters and unsupported rows. Missing evidence stays missing.

Relevant files:

- flylab/pharm/library.yaml
- flylab/pharm/evidence.py
- flylab/pharm/occupancy.py

### Connectome-dependence analysis

The dependence analysis compares a real-graph effect with distributions produced by degraded graph models that retain different amounts of network information.

Its safest interpretation is comparative:

- a small empirical permutation p means the real-graph effect differs from that null ensemble at the chosen permutation effort,
- a large p means the analysis did not distinguish the real effect from that null ensemble,
- a large p is not, by itself, proof of mathematical or biological equivalence.

Relevant files:

- flylab/analysis/dependence.py
- flylab/analysis/nullmodels.py

### Ablation

The ablation ladder asks which layer carries the ordering information across compounds.

Relevant file:

- flylab/analysis/baselines.py

### Specification robustness

The robustness layer re-runs qualitative conclusions across a prespecified family of engagement-to-gain transformations.

This is a model-robustness analysis, not a posterior probability over model specifications.

Relevant file:

- flylab/analysis/robustness.py

### Global uncertainty and value of information

The uncertainty layer attributes variance in a selected model output to declared uncertain inputs and modelling choices.

The VOI layer converts first-order model-variance contributions into a prioritization aid. These are statements about **model variance under assumed input ranges**, not biological variability or expected clinical/toxicological value.

Relevant files:

- flylab/analysis/uncertainty_global.py
- flylab/analysis/voi.py

## Reproducibility

The native publication pipeline is:

~~~bash
python scripts/reproduce_paper.py
~~~

It regenerates figures, tables, papers/results.json, and rendered manuscript text from committed artifacts.

The browser build runs the same Python science core through Pyodide. Browser/server parity tests cover representative route families. That is different from claiming the browser regenerates the complete publication pipeline.

See [PAGES.md](PAGES.md) for details.

## Data documentation

The main data locations are:

| Path | Meaning |
|---|---|
| flylab/pharm/library.yaml | typed pharmacology library |
| data/literature/ | supporting sourced datasets |
| data/derived/ | committed MaleCNS-derived research cuts |
| papers/results.json | machine-readable publication result record |
| papers/tables/ | generated paper tables |
| papers/figures/ | generated paper figures |

Large upstream connectome files are not committed to the repository.

## Claim-language guardrails

Public documentation and manuscripts should distinguish the following carefully.

Prefer:

- "not distinguishable from the null ensemble at n permutations"
- "model-derived"
- "functional engagement" for EC50/IC50-derived Hill responses
- "binding-derived occupancy" only when the underlying evidence and preparation support that wording
- "prospective software prediction"
- "literature concordance"

Avoid:

- "the null reproduces the effect" when only non-significance was shown
- "validated" for source-overlapping or model-internal comparisons
- "pre-registered" for predictions that exist only in a mutable repository
- "safe" based on the vertebrate receptor scorecard
- "the connectome matters" without naming the compound, concentration, readout and null comparison
- biological conclusions from circuit outputs that have no independent biological validation

See [NOVELTY.md](NOVELTY.md) for the research-positioning version of these rules.

## Development checks

Fast contract suite:

~~~bash
python -m pytest -q
~~~

Longer end-to-end checks:

~~~bash
pytest -m slow
~~~

Paper regeneration:

~~~bash
python scripts/reproduce_paper.py
~~~

Before changing the scientific core, read [HANDOFF.md](../HANDOFF.md).

## Project description

A set of suggested GitHub About text, longer project summaries and repository topics is maintained in [PROJECT_DESCRIPTION.md](PROJECT_DESCRIPTION.md).
