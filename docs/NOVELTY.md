# Novelty (keep this short in the paper)

Existing virtual rats (RatCVS, Virtual Rat Web, PharmVR) model organs without neurons.

Existing fly connectome runtimes (FlyWire/MaleCNS LIF ports, flypoke, meme sims) model neurons without ligand–receptor occupancy.

FlyLab is the missing object: **the same compound, two scorecards, one notebook** — insect circuit plus vertebrate receptor panel, at the same dose, with the provenance of every number attached and a live-assay slot left empty.

If a PR only adds a prettier fly viewer, reject it.

## What is actually new, sharpened by v0.5

Three of these did not exist in v0.4 and are what the paper argues for.

### 1. The dual scorecard as one artefact

Not "we also computed a vertebrate number". The vertebrate panel is evaluated at the *same free concentration* as the insect panel, in the same call, and is structurally forbidden from patching a circuit. Every row on both sides carries an evidence tier, so a reader can see which half of the claim rests on a measured constant and which on a placeholder. Nothing else writes that object.

### 2. Null models for a connectome drug effect

This is the contribution most likely to outlive the rest of the tool, and it is a **method for disbelieving your own simulation**. Four degradations of the same cut — transmitter-label permutation (E/I histogram preserved), weight permutation, degree-preserving rewiring, Erdős–Rényi — with the drug effect measured on each, seeds held fixed so the same named cells are driven and read.

The finding is the point:

- **Imidacloprid's effect on the neighbourhood mean rate beats only the Erdős–Rényi null** (|z| ≈ 46) and fails all three structure-preserving nulls (max |z| ≈ 1.2). Because `sign_permute` preserves the transmitter histogram exactly, that says the effect follows the cut's global excitation/inhibition balance, not the identity of the cholinergic cells. An nAChR agonist that floors `g_ach` removes a fixed fraction of excitatory weight wherever it sits.
- **Fipronil does beat the topology nulls** (z ≈ 3.4 weight permutation, z ≈ 3.0 degree-preserving rewiring). Disinhibition depends on *where* the inhibition is, and that is a wiring fact.
- **The bitter veto ratio is not yet distinguishable from its shuffles** (max |z| ≈ 0.7 at 10 shuffles) and is reported as model behaviour, not as evidence about wiring.

Any connectome-simulation paper that reports only an Erdős–Rényi comparison is reporting that its graph has the right size. FlyLab ships the harder controls and publishes the ones it fails.

### 3. The selectivity landscape: does the circuit amplify or buffer the receptor margin?

Receptor selectivity index = log10(EC50_vert / EC50_insect), a property of two numbers in a YAML file. Circuit selectivity index = log10 of the window between the concentration that moves the simulated network by 50 % and the concentration at which the most potent vertebrate target reaches 20 % occupancy. Their difference is a quantity no receptor table and no connectome can produce alone.

The split is **by mechanism, not by potency**:

- **amplify** (+0.29 log10 mean): DDT, deltamethrin, permethrin, chlorpyrifos-oxon — Nav modulators and the AChE inhibitor, whose gains act globally.
- **buffer** (−0.88 log10 mean): imidacloprid, clothianidin, acetamiprid, nitenpyram, spinosad, nicotine, acetylcholine — the nicotinic set, because the agonist curve rises before it falls, so low occupancy barely moves the network.

A separation of about 1.2 log units, reproducible on both committed cuts, generated from sourced EC50s and a public netlist. Also reported: eight compounds get **no** circuit index at all, because RDL/GluCl block cannot reach the 50 % threshold on these cuts — a limit of the readout, stated rather than engineered around.

### 4. Provenance discipline that pays

Evidence tiers, a library SHA-256 in every notebook, per-compound (never per-receptor) genotype shifts, and a literature dataset that lists ten places where the sources contradict the library. It caught a real error: fipronil's vertebrate GABA-A EC50 was 1.0 × 10⁻⁵ M while citing a paper that reports 1.1 × 10⁻⁶ M. After correction the vertebrate occupancy at 1 µM is 0.48 and the "vertebrate-safe" framing is gone. A source string that can be checked eventually is.

## Not novelty

- Reproducing the bitter veto. Shiu *et al.* did that; FlyLab uses it as a directional control.
- Hill occupancy. Textbook.
- A bigger graph, a 3D viewer, more compounds, or a spiking model of the whole CNS.
- Any claim about a living fly. There are none in this repository.
