# FlyLab: A Dual-Scorecard Virtual Pharmacology Bench Grounded in the Drosophila Connectome

**Working draft for the International Journal of Research in Computing (IJRC)**  
Software: https://github.com/pinkysworld/FlyLab  
Status: architecture and prototype. No live-animal dataset is claimed.

## Abstract

Virtual animals used in pharmacology teaching (RatCVS and its web ports) implement dose-response at the organ level. They contain no neurons. Public adult Drosophila connectomes (FlyWire; MaleCNS v1.0) and spike-network ports implement neurons and predicted transmitter sign, but not ligand-receptor occupancy. FlyLab joins those objects: a compound and concentration map to occupancy on insect and vertebrate receptor classes; occupancy patches gain on a named fly circuit, on a 1-hop MaleCNS neighborhood of MN9 and DNp01, and on a whole-CNS cell census; the same dose is scored on a vertebrate receptor panel; the run is exported as notebook JSON. We do not claim a new EM map, a fitted in vivo EC50, or replacement of GLP toxicology.

**Index terms:** connectomics, in-silico pharmacology, Drosophila, receptor occupancy, scientific software, 3Rs

## 1. Introduction

Two software traditions sit next to each other and do not touch. RatCVS models a pithed-rat cardiovascular preparation [1], [2]. Neurons are absent by design. FlyWire and MaleCNS v1.0 are public synaptic wiring diagrams [5], [6]. Point-neuron simulations recover qualitative controls such as sugar-driven MN9 and bitter veto [7]. Those ports are netlists plus LIF. They do not contain binding curves or a vertebrate comparison at the same dose.

Pharmacology is a statement about receptors. Imidacloprid is far more potent at insect nAChRs than at many vertebrate subtypes [8], [9]. FlyLab evaluates the same compound on (i) an insect circuit / MaleCNS neighborhood / census and (ii) a reduced vertebrate receptor panel.

Map choice: **MaleCNS v1.0** (brain + ventral nerve cord) so motor and descending neurons (MN9, DNp01) exist in the same atlas.

## 2. Related work

Virtual Rat Web [2] and PharmVR [3] have no named CNS cells. Fly connectome sims [5]-[7], [11] have no occupancy. Insect nAChR hybrid-receptor EC50s [8], [9] supply teaching priors, not synapse-resolved occupancy. Gap: no public tool that writes one notebook with insect circuit plus vertebrate panel.

## 3. Design goals

Reproducible JSON notebook. Two scorecards. Named cells. Fail loud if the circuit is reduced. Cross-platform `flylab serve`. No invented kinetics. Pairable with a fly room later. Large connectome files stay off git; GitHub Actions cut a named neighborhood.

## 4. Methods

Occupancy is the Hill function. Library YAML carries EC50, n, direction, source. Agonist nAChR gain: g_ACh = max(0.05, 1 + 0.4 theta - 1.6 theta^2). Antagonist RDL gain: g_GABA = max(0.05, 1 - theta). Taste assay is reduced_taste_v0 (five rate units) constrained by the Shiu sugar/bitter direction [7]. Neighborhood assay loads the committed 1-hop graph around MN9 and DNp01 extracted from the public MaleCNS weight matrix (1126 nodes, 1360 edges at synapse weight >= 5). Whole-CNS assay loads MaleCNS traced neurons (~165k) and consensus transmitters, then applies occupancy to an excitation index weighted by real ACh/GABA/Glu cell counts. Vertebrate panel is occupancy only, not a mouse connectome. Software: Python 3.11, FastAPI, pytest, GitHub Actions. Maps are CC-BY and downloaded, not committed.

## 5. Prototype behaviour (software, not animals)

Imidacloprid at 1e-6 M: insect nAChR occupancy ~0.99 vs vertebrate a4b2 ~0.09 (teaching library). Bitter vetoes MN9 on reduced_taste_v0 by >75%. High-dose imidacloprid lowers sugar-driven MN9 and the whole-CNS excitation index versus vehicle. Fipronil occupancy is high at insect RDL and lower at vertebrate GABA-A in the same library. Census and neighborhood report MN9 and DNp01 body IDs from the real map. Notebooks warn when the circuit is reduced or hops-limited.

## 6. Discussion

A 166k-cell LIF on the 1.1 GB weight matrix is the next engineering gate (`flylab download-malecns --full`), not a change of thesis. A scientific result later requires one live Drosophila assay and at most one fitted parameter. Limits: predicted transmitters, Hill not full kinetics, teaching EC50s, no live data in this draft. IJRC fit is computing infrastructure.

## 7. Conclusion

FlyLab is a dual-scorecard bench on the Drosophila connectome programme. The shipped prototype computes occupancy, contrasts insect and vertebrate receptors, runs a reduced taste circuit, a MaleCNS named-cell neighborhood, and a whole-CNS census, and exports a notebook. It does not integrate the full Janelia weight matrix as a spike network and does not report live flies.

## Statements

AI used for drafting code and text. No generated live-animal experiment is included. Code: https://github.com/pinkysworld/FlyLab. Connectomes: official Janelia releases, CC-BY. Competing interests: none. Ethics: no new animal experiments in this draft.

## References

[1] J. Dempster, RatCVS, Univ. of Strathclyde.
[2] S. Okabe et al., Eur. J. Pharmacol., 997, 177618, 2025.
[3] PharmVR preprint, 2025.
[4] J. Educ. Eval. Health Prof., 2020.
[5] S. Dorkenwald et al., Nature, 634, 2024.
[6] HHMI Janelia and Google Research, MaleCNS v1.0, 2026.
[7] P. Shiu et al., Nature, 2024.
[8] C. Bass et al., Insect Biochem. Mol. Biol., 36, 86-96, 2006.
[9] H. Dederer et al., Insect Biochem. Mol. Biol., 2011.
[10] FlyConnectome neurotransmitter ground-truth repository.
[11] P. Mineault, Deconstructing viral fly sims, 2026.
[12] Chemoconnectomics of Drosophila, Neuron, 2019.
