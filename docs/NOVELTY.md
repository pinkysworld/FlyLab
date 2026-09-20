# Novelty (keep this short in the paper)

Read this before accepting a PR. If a PR only adds a prettier fly viewer, reject it.

## What this is *not* novel against

Say these out loud before claiming anything:

- **Executable fly circuits from connectome data exist.** FlyBrainLab (Lazar *et al.*, *eLife* 10:e62362, 2021, DOI 10.7554/eLife.62362) builds and compares executable circuits from fly brain data interactively. FlyLab is not the first thing to simulate named *Drosophila* circuits.
- **Connectome LIF runtimes exist.** Shiu *et al.* (*Nature* 634, 2024) reproduce sugar-driven motor output and the bitter veto on the whole adult brain without fitting. FlyLab uses that result as a directional control, not as a finding.
- **Receptor-informed pharmacological perturbation of a connectome exists.** On human structural connectomes, whole-brain models are perturbed by scaling regional synaptic responses with PET receptor-density maps (Mindlin *et al.*, *Commun. Biol.* 7:1176, 2024, DOI 10.1038/s42003-024-06852-9; Deco *et al.*, *Curr. Biol.* 28:3065, 2018). **That work also found the effect dominated by overall receptor presence rather than receptor placement** — the improvement correlated with the mean density of activated receptors across the brain. Our composition-dominated result is the same phenomenon at synapse resolution in an insect. Cite it: the convergence makes FlyLab more credible, not less novel.

Never write "no tool has combined connectomes and pharmacology". It is false and a reviewer will catch it.

## What *is* new, narrowly

1. **Typed pharmacological evidence.** Each row records the parameter its source measured (`Kd`/`Ki` vs `EC50`/`IC50`/`Kb` vs nothing) and the source's distance from this compound/receptor/species. Those two facts decide which transformation is permitted: binding occupancy, functional engagement, or **not modelled**. Unsupported rows return N/A and are excluded from numeric results — enforced by the type system (`EvidenceTypeError`), not by convention. In the shipped library: 101 rows, 34 EC50, 10 IC50, **exactly 1 Kd** (imidacloprid at `insect_nAChR_beta1`, Bass 2011 saturation binding), 56 not modelled.
2. **Compound- and concentration-specific perturbation of named, synapse-resolution circuits**, with the insect and vertebrate panels scored at the same free concentration and only the insect side allowed to patch a circuit.
3. **Connectome-dependence analysis.** For each `compound × concentration × readout`, the empirical permutation *p* against four degraded graph models, the dependence class, and the **weakest graph model that already reproduces the effect**. This is the transferable method.
4. **Specification robustness and an uncertainty budget** that tell you which of your own conclusions are properties of a modelling choice.
5. **Provenance that survives**: library SHA-256 in every notebook, per-compound (never per-receptor) genotype shifts, ten recorded literature-vs-library contradictions, and one caught library error (fipronil's vertebrate GABA-A EC50 contradicted its own citation ~9-fold; after correction "vertebrate-safe" is no longer sayable).

## The findings that justify 3 and 4

**Dependence (n = 1000 permutations, `named` cut, mean rate).**

- Imidacloprid, −6.17 Hz: **composition-dominated**. p = 0.275 (sign), 0.586 (weight), 0.472 (degree), 0.0010 (ER, at the resolution floor). Necessary level: `degree_sequence`. A graph that knows only each node's degree and transmitter reproduces the drug effect.
- Fipronil, +0.90 Hz: **topology-dependent**. p = 0.0060 (weight), 0.0070 (degree), 0.113 (sign). Necessary level: `wiring_without_transmitter_identity`.
- Taste-motor arm, fipronil veto ratio at n = 300: p = 0.42 / 0.88 / 0.20 / 0.24, all |z| < 0.25 — **a properly powered negative**, not the ten-shuffle hand-wave of v0.5. Imidacloprid's veto ratio is **undefined** there (MN9 silenced), not zero.
- Landscape, 21 compounds × 4 concentrations: **20 topology-dependent, 51 composition-dominated, 13 no-effect, 0 mixed**. The topology-dependent set is exactly the chloride-channel blockers (dieldrin, fipronil, ivermectin, picrotoxin, gaba ≥ 1e-7, chlorpyrifos-oxon at 1e-8).
- Permutation-count sweep: verdicts settle from n ≈ 25–50; quoting a mid-range *p* to ±0.02 needs 400–1000.

**Ablation ladder (21 compounds).** Receptor engagement alone orders the library **backwards** (ρ −0.05 to −0.51). Composition-only (mechanism gains on the cut's transmitter proportions, no edges) **reproduces the full model's ordering** (ρ 0.906–0.989; ρ 0.953 within the nine nicotinic agonists). Topology-only (real cut, generic multiplier) is the **worst** level (ρ 0.24, 11.7 Hz rms). Nearly all the information is in the mechanism rules; the connectome without the pharmacology is not a cheap substitute for the connectome with it.

**Specification robustness (25 prespecified gain specifications).** "Nicotinic agonist suppresses circuit activity" retained 15/25 and **reversed by all ten monotone specifications** — under a monotone rule the same drug at the same engagement *excites* the network. "Nicotinic buffering" 18/25, reversed by none (7 undecidable). "Nav/AChE amplification" 19/25. RDL disinhibition, both topology conclusions and the map bitter-veto direction: **25/25**. So the topology results are specification-independent and the suppression result is not.

**Uncertainty budget (Sobol', n_base 1024, 11264 evaluations, Var(Y) = 3.149 Hz²).** gain_transform S1 0.410 (ST 0.626), weight_threshold S1 0.294 (ST 0.542), drive 0.019, everything else at the ±0.007 noise floor measured by a null factor. VOI: gain_transform 1.29 Hz² (needs a synaptic-gain calibration), weight_threshold 0.925 Hz² (**needs no experiment**, only synapse-confidence strata).

## Not novelty

- Reproducing the bitter veto. Shiu *et al.* did that; FlyLab uses it as a directional control.
- A Hill curve. Textbook. (And on an EC50 it is *engagement*, not occupancy.)
- Exposure C(t), mixtures, genotype shifts, expression weighting — supporting capabilities, documented in the supplement, not equal-weight novelty claims.
- A bigger graph, a 3D viewer, more compounds, or a spiking model of the whole CNS.
- Any claim about a living fly. There are none in this repository.
