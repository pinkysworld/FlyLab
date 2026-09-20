# Architecture

FlyLab is built around a strict separation:

- **pharmacological evidence** is a sourced input,
- **connectome structure** is measured upstream data,
- **gain rules and runtime settings** are modelling assumptions,
- **readouts and analyses** are computed outputs.

That separation is what makes provenance and model auditing possible.

## Scientific dataflow

~~~text
compound + free concentration
        |
        v
typed receptor evidence
        |
        |  parameter type
        |  source relation
        |  species
        |  evidence tier
        v
receptor engagement representation
        |
        +-------------------------> vertebrate receptor scorecard
        |                           scored only
        |                           never patches a circuit
        v
mechanism -> gain transformation
        |
        v
MaleCNS-derived circuit
        |
        +--> deterministic rate runtime
        +--> LIF runtime
        |
        v
named-cell and network readouts
        |
        +--> dependence analysis
        +--> ablation
        +--> specification robustness
        +--> global uncertainty
        +--> value of information
        +--> provenance audit
        |
        v
notebook + generated research artifacts
~~~

## Evidence layer

The pharmacology library is schema v3.

Each receptor row records:

- parameter type,
- parameter value when supported,
- Hill coefficient,
- direction,
- species or preparation,
- source relation,
- source text,
- evidence tier.

The central rule is simple:

> Missing evidence remains missing.

Unsupported rows return no numeric engagement and are excluded from numeric aggregation paths that require a sourced value.

The evidence layer also distinguishes a binding parameter such as Kd/Ki from functional potency such as EC50/IC50. Source relation remains visible because a strong parameter measured in another species or preparation is still not a direct measurement of the modelled Drosophila system.

Parameter type is only half of the rule. Each source relation is also graded as an *evidence distance*, E0 to E4, where E0 is this compound at this receptor in this species and E4 is a row nothing supports. The pair (parameter type, evidence distance) decides which engagement representation may be computed, and the representation is carried in the result rather than inferred by the reader:

- binding occupancy, for a binding parameter measured on target,
- binding engagement proxy, for a binding parameter transferred across species or onto a related preparation,
- functional engagement, for a functional potency measured on target,
- functional engagement proxy, for a functional potency transferred, or for a binding parameter carried as far as a class statement,
- not modelled, for anything at E4 and for row types that carry no potency at all.

A proxy is a labelled extrapolation, not a measurement, and the label travels with the number into the notebook.

Primary files:

- flylab/pharm/library.yaml
- flylab/pharm/evidence.py
- flylab/pharm/occupancy.py

## Mechanism layer

The mechanism layer transforms receptor engagement into a gain patch.

The circuit runtime does not compute pharmacology itself.

Primary source of mechanism rules:

- flylab/pharm/mechanisms.py

Typical outputs include:

- g_ach
- g_gaba
- g_glu
- g_oct
- g_nav
- ach_tone

These transforms are model assumptions. Some are intentionally stress-tested by the robustness layer because their functional form is not measured directly.

## Connectome layer

FlyLab uses small committed circuit cuts derived from MaleCNS v1.0 rather than requiring the full source matrix for routine reproduction.

Current research artifacts include:

- named neighbourhood cut,
- taste-motor cut,
- compact cell census fallback.

These derived files are research inputs and must not be hand-edited.

Primary locations:

- data/derived/
- flylab/maps/

## Circuit runtimes

### Rate model

The deterministic rate runtime is the workhorse for large permutation and uncertainty analyses.

It is designed to make repeated graph perturbation computationally practical.

### LIF model

The LIF runtime provides a spiking implementation with explicit stochastic drive.

Absolute rate agreement between the rate and LIF engines should not be assumed. Where the engines agree only on direction, documentation must say so.

Primary files:

- flylab/circuit/rate.py
- flylab/circuit/lif.py

## Assays

Assays combine:

1. typed pharmacology,
2. a circuit or reduced model,
3. a runtime,
4. readouts,
5. warnings,
6. notebook provenance.

Primary package:

- flylab/assays/

The assay layer should not duplicate mechanism rules.

## Analysis layer

The analysis package is deliberately allowed to produce results that weaken the model's own claims.

### Connectome dependence

flylab/analysis/dependence.py compares a real-graph drug contrast with degraded-graph null ensembles.

The key statistical distinction is:

- **distinguishable from a null ensemble**
- **not distinguishable from a null ensemble**

A non-significant comparison does not, by itself, establish equivalence. Durable documentation should therefore avoid treating "not distinguishable" as identical to "reproduced" unless an explicit tolerance or equivalence rule is added.

### Null models

flylab/analysis/nullmodels.py generates graph degradations that preserve different information.

Current families include models that alter:

- edge placement,
- degree-preserving wiring,
- weight placement,
- transmitter identity.

### Ablation

flylab/analysis/baselines.py removes layers of information to determine which layer carries ordering information across compounds.

### Specification robustness

flylab/analysis/robustness.py re-runs qualitative claims across alternative engagement-to-gain specifications.

This is a robustness analysis over a declared model family, not a probability distribution over models.

### Global uncertainty

flylab/analysis/uncertainty_global.py attributes model-output variance to declared uncertain factors.

Interpretation is conditional on the chosen ranges.

### Value of information

flylab/analysis/voi.py translates first-order model-variance contributions into a prioritization aid for measurements or re-analysis.

This is value with respect to model uncertainty, not biological or clinical value.

### Claim provenance

The claims layer classifies dependency links behind a result so a reader can distinguish:

- observed structure,
- literature-derived input,
- model assumption,
- computed transformation,
- remaining unknowns.

## Validation layer

The validation package should use the word **validation** narrowly.

Published compound rankings or qualitative mixture outcomes are useful concordance checks, but they are not independent circuit-level validation when they share sources with the library or evaluate a different organism or endpoint.

Primary package:

- flylab/validation/

## Interfaces

FlyLab exposes the science core through three front ends:

| Interface | Implementation |
|---|---|
| CLI | flylab/cli.py |
| local web API | flylab/server.py |
| browser bridge | flylab/browser/bridge.py |

The interfaces should delegate to the same scientific implementation rather than carrying duplicate model logic.

See [PAGES.md](PAGES.md) for browser-specific constraints.

## Reproduction pipeline

scripts/reproduce_paper.py is the native publication pipeline.

It generates:

- paper figures,
- paper tables,
- papers/results.json,
- rendered manuscript text,
- rendered supplement text.

Generated numerical paper content should not be hand-edited.

The manuscript templates are the editable sources for prose containing generated values.

## Notebook provenance

A FlyLab notebook records enough context to identify the computational environment that produced it.

Typical provenance includes:

- FlyLab version,
- source revision where available,
- pharmacology library hash,
- map identity,
- RNG seed,
- platform,
- warnings.

The browser may not have a git SHA available inside Python. Static build metadata can record the build revision separately.

## Layer ownership

| Layer | Package or path | Responsibility |
|---|---|---|
| evidence | flylab/pharm/ | sourced receptor evidence and transformations |
| maps | flylab/maps/, data/derived/ | MaleCNS-derived inputs |
| runtime | flylab/circuit/ | numerical circuit execution |
| assays | flylab/assays/ | compose evidence, circuit and runtime |
| analysis | flylab/analysis/ | model auditing and inference about model dependence |
| validation | flylab/validation/ | literature comparisons and independence checks |
| interfaces | CLI, server, browser, static | transport and presentation |
| reproduction | scripts/reproduce_paper.py | publication artifact generation |

## Non-goals

The architecture intentionally does not provide:

- a vertebrate neural circuit,
- fly pharmacokinetics that predict receptor concentration from a dose,
- automatically generated live-animal observations,
- regulatory safety conclusions,
- a second JavaScript implementation of the scientific model,
- silent numeric defaults for unsupported receptor evidence.

## Design principle

The most important architectural principle is:

> A result should be inspectable as a chain from source evidence and measured structure, through explicit modelling assumptions, to computed output.

A feature that makes the UI richer but makes that chain harder to audit is usually a regression for FlyLab.
