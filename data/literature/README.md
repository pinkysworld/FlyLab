# `data/literature/` — literature-grounded datasets for FlyLab

Six YAML datasets compiled from the primary literature for the FlyLab virtual
pharmacology bench. They are **inputs and validation targets for the model**,
not model output, and not live-animal data.

Three rules govern every file:

1. **Nothing is invented.** Every number is traceable to a named publication
   (author, year, journal, and DOI or PMID).
2. **A missing number stays missing.** Where a value could not be found the
   field is `null` and the reason is stated in `notes`. Each file ends with a
   `not_sourced:` list. Absence of a value is itself a finding and is reported
   as one.
3. **Species is always named.** *Drosophila melanogaster* data are preferred.
   Where the only quantitative work is in another arthropod (aphid, planthopper,
   moth, mite, honeybee, house fly) the `species` field says so, and the entry
   must be flagged as a cross-species transfer in any notebook that uses it.

Nothing in this directory may be written into a notebook's `live_lab` field.
No Python file was touched to produce these datasets.

## Evidence tiers

The same three tiers used by `flylab/pharm/library.yaml`, applied here per
numeric field rather than per receptor row:

| tier | meaning |
|---|---|
| `literature_value` | The number is stated in the cited publication for this allele/compound/class/assay. A value derived by simple arithmetic from two published numbers is also tagged `literature_value`, and the derivation is written out in `notes` (e.g. a synergism ratio computed from two published resistance ratios). |
| `literature_order` | Only an order of magnitude, a direction, or an ordinal statement ("considerably more cells", "approximately 10-fold") is supported by the cited publication. Not a measured constant. |
| `class_placeholder` | No usable number exists. The field is `null` and `notes` says why. Never a filler value. |

A fourth tier, `measured_fit`, is reserved by the library for values fitted to
user-imported data. Nothing here uses it.

## The files

**`resistance_alleles.yaml` — 9 alleles.** Genotype toggles expressed as
EC50/IC50 shifts on the five FlyLab insect receptor keys. Covers *Rdl* A302S
(*D. melanogaster*), A301G+T350M (*D. simulans*) and A301S (*N. lugens*);
nAChR β1 R81T (*M. persicae*); *para* kdr L1014F and super-kdr M918T; the four
*Drosophila* *Ace* point mutations; and the GluCl abamectin-resistance
substitutions. Every shift carries an `assay_type` and a `use_as_ec50_shift`
flag, because FlyLab applies an allele as a multiplicative EC50 shift and that
is only defensible for expressed-receptor electrophysiology or binding.
Whole-animal resistance ratios (which contain penetration and metabolism as
well as target-site insensitivity) are marked `use_as_ec50_shift: false`. Two
entries are there specifically because they break the obvious assumption:
*Rdl* A301S does **not** shift fipronil potency on the expressed planthopper
receptor, and *para* M918T leaves DDT almost untouched while profoundly
reducing pyrethroid sensitivity — so a toggle that shifts all ligands at a
receptor together would be wrong in both cases. No GluCl "G36A" allele exists
in the literature; A309V, G315E, G323D and G326E are the published ones.

**`published_rank_orders.yaml` — 12 entries.** Retrospective validation
targets: ordered compound lists from named published assays, each with a
`use_in_flylab` field of `rank-order-only`, `fold-shift` or `curve`. Includes
neonicotinoid potency in *Drosophila* larval and adult bioassays, chronic
lifespan ordering, resistance ratios in *M. persicae*, the rat α7 partial
agonist EC50 set, the fipronil insect-versus-vertebrate GABA-A contrast, and
the three requested selectivity contrasts: benzodiazepines inactive at insect
RDL (flunitrazepam, 0.1–100 µM, no effect), fipronil far weaker at native
vertebrate GABA-A, and ivermectin's ~1000-fold GluCl selectivity over
mammalian GlyR and GABA-A. Two entries are deliberate counterweights: fipronil
is nanomolar-potent at the **human β3 homopentamer**, so "vertebrate-safe" is
an over-claim that depends on subunit composition; and no single *Drosophila*
RDL preparation was found that ranks fipronil, picrotoxin and dieldrin
together, so no such ordering is asserted.

**`receptor_expression_by_class.yaml` — 6 published atlas values, 10 MaleCNS
superclass rows, 4 transmitter-class rows.** Coarse, class-level expression of
nAChR subunits, Rdl, GluClα, octopamine receptors and Ace. **This is an
explicitly interpretive cross-atlas mapping**: none of the cited single-cell
studies is annotated with MaleCNS `superclass` labels, so every superclass row
is tagged `mapping: cross_atlas_inference` and its level is `unknown` unless a
cited paper says something directly about that class. Most rows are `unknown`.
The file never records `absent` from a single-nucleus dropout. Three facts are
genuinely quantitative and directly useful: the VNC is 40% cholinergic / 38%
GABAergic / 18% glutamatergic (Allen 2020); all seven nAChR α subunits are
expressed in the midbrain but α1/α5/α6/α7 in far more cells than α2/α3/α4, with
β1 in more than twice as many cells as β2 (Croset 2018) — meaning
`insect_nAChR` is not one receptor and subunit-selective compounds reach only
partly overlapping cell sets; and GluClα is broadly expressed in the VNC.
**Motor neurons are a confirmed gap**: no cited source reports Rdl or nAChR
expression for identified adult leg or labellar motor neurons, which is
recorded in a `motor_neuron_gap` block with a suggested next step.

**`mixtures.yaml` — 4 reference models, 7 empirical entries.** Bliss
independence and Loewe additivity are written out as equations with variable
definitions, applicability notes and their primary citations, plus the Greco
1995 response-surface review, and a recommendation that FlyLab use Loewe when
two compounds share a receptor key and Bliss when they do not. The empirical
entries are neonicotinoid + PBO in *M. persicae* (synergism ratios 7.2 and 8.7,
derived from published resistance ratios), PBO analogues (up to 290-fold —
explicitly *not* PBO itself), fipronil + PBO in *C. suppressalis*, and
imidacloprid's binary mixtures in honeybee. The honeybee set is the most useful
validation target in the file: neonicotinoid + pyrethroid and
neonicotinoid + organophosphate were only **additive**, and
imidacloprid + clothianidin — two agonists at the same target — showed **no
interaction**, which is what Loewe additivity predicts and is a direct test for
`flylab/pharm/binding.py`. Ivermectin + PBO and a numeric neonicotinoid +
pyrethroid synergism ratio could not be sourced and are `null`. Note the
sign reversal in the fipronil + PBO data: PBO synergises in the resistant
strain and **antagonises** in the susceptible one, so any FlyLab metabolism
module must be able to represent bioactivation, not only clearance.

**`fly_pharmacokinetics.yaml` — 15 entries, most of them null.** This is the
thinnest file and that is the finding. Of the seven quantities an exposure
model needs for the adult fly — haemolymph volume, body mass, oral
bioavailability, cuticular penetration rate, distribution volume, elimination
half-life, metabolite fractions — **only metabolite identity is well
established**. Adult total haemolymph volume, adult body mass, cuticular
penetration rate, any *Drosophila* elimination half-life, nicotine kinetics, and
ethanol clearance are all `null`. The widely repeated "~80 nL haemolymph"
figure could not be traced to a primary measurement and is deliberately **not**
recorded. What is recorded: 25 nL recoverable per adult by nanolitre sampling,
174 ± 81 nL in third-instar larvae at 1.39 ± 0.15 mg, an imidacloprid
half-life of 4.5–5 h in *honeybee* (the only verified insect half-life for a
library compound), the 4 ppm chronic *Drosophila* dose giving 50% mortality in
20–25 days, and the CYP6G1 / gut-microbe metabolism results. Two structural
warnings sit in this file: haemolymph is the **lowest**-radioactivity
compartment for imidacloprid, so a one-compartment haemolymph model will
under-estimate CNS exposure; and a substantial metabolic route runs through the
microbiome, so a single-enzyme clearance term is wrong for this compound.
`flylab/pharm/exposure.py` should be labelled `class_placeholder` in notebooks
until a real fly PK study is imported.

**`behavioral_assays.yaml` — 5 PER entries, 4 climbing entries.** Published
control statistics for power analysis, recorded so that a power calculation can
be done before any fly is touched. Climbing is well served: the Martelli 2020
imidacloprid study gives a full protocol (5 female flies per vial, 7 cm,
30 s cutoff, n = 100 per group per timepoint) and a large effect (6% → 44% →
68% failing at days 1, 10, 20 on 4 ppm imidacloprid), and the automated RING
protocol fixes n and the scoring window. **PER is not well served**: the
specific quantity requested — a control PER rate to 100 mM sucrose with its
variance — was **not found** and is `null`, largely because 100 mM sucrose is
commonly used as a pre-screen in which only 100%-responding flies are kept, so
no baseline rate is reported. The best available substitute is Çevik & Erden
2012 (600 mM sucrose, 40% responders at 1 h food deprivation rising to 68% at
2 h, n = 463 wild-type), which also shows that food-deprivation time is a
larger nuisance variable than most drug effects FlyLab would try to detect.
Quantitative bitter-suppression effect sizes are `null`, although French et al.
2015 supplies the qualitative result that `flylab/circuit/reduced_taste.py`
encodes. Binomial power figures in the `notes` are arithmetic derived here, are
labelled as such, and assume independent flies — which vial-based protocols
violate, so every n is a lower bound.

## Where the literature contradicts the current teaching EC50s

Reported, not edited. `flylab/pharm/library.yaml` was not modified.

1. **Fipronil at `vertebrate_GABA_A` is ~9× too weak.** The library uses
   `1.0e-5 M` and cites Ratra & Casida 2001, but that paper's vertebrate
   average IC50 for fipronil is **1103 nM ≈ 1.1e-6 M**. The library value is
   not the cited value.
2. **Fipronil at `vertebrate_GABA_A` is ~4000× too weak against the β3
   homopentamer.** Ratra, Kamita & Casida 2001 report 0.5–2.4 nM for fipronil,
   lindane and α-endosulfan on the human β3 homo-oligomer. The library's
   "vertebrate-safe" framing holds only for native heteromeric receptors and
   should say so.
3. **Fipronil at `insect_RDL` sits in the middle of a ~1500-fold literature
   spread.** The library uses `3.0e-8 M`. House-fly membrane binding gives
   3–12 nM (Ratra & Casida 2001) and ~2.4–6.3 nM (Cole 1993); *Drosophila* RDL
   in oocytes gives 0.3–3 µM (Lees 2014, already noted in the library's own
   source string). The chosen value is defensible as a midpoint but should not
   be quoted as a measured constant.
4. **The fipronil selectivity *ratio* is fine even though both absolute values
   are off.** Library ratio 1e-5 / 3e-8 = 333×; Ratra & Casida give
   1103 / (3–12) = 92–368×. The scorecard's selectivity column survives; the
   absolute occupancy curves do not.
5. **Nitenpyram's rank is inverted relative to the *Drosophila* bioassay.** The
   library puts nitenpyram at `5.0e-8 M` on `insect_nAChR`, among the most
   potent neonicotinoids. In *Drosophila* bioassays (Perry et al. 2012)
   nitenpyram is the **least** potent of the set tested — 8× weaker than
   clothianidin in larvae and 25× weaker in adults. Receptor potency and
   whole-animal potency legitimately differ (nitenpyram has poor persistence),
   but any FlyLab claim to reproduce a published *ordering* will fail on this
   compound unless the difference is stated.
6. **Clothianidin versus imidacloprid ordering is inverted.** The library makes
   imidacloprid (`2.0e-8`) slightly more potent than clothianidin (`3.0e-8`);
   both the acute *Drosophila* bioassays and the chronic-lifespan ordering put
   clothianidin first. The gap is within the uncertainty of either number, but
   the sign is opposite.
7. **`gaba` at `insect_RDL` is species-dependent by ~13×.** The library uses
   `5.0e-5 M` from *Drosophila* RDL (Lees 2014, pEC50 4.2–4.3). *N. lugens* RDL
   wild type gives pEC50 5.43 ≈ 3.7 µM (Garrood 2017). Not an error, but the
   RDL agonist EC50 is not a cross-species constant.
8. **A genotype toggle must not shift a whole receptor uniformly.** Two
   findings in `resistance_alleles.yaml` contradict the natural implementation:
   *Rdl* A301S shifts ethiprole and GABA but **not** fipronil on the expressed
   planthopper receptor; *para* M918T shifts deltamethrin and permethrin 100×
   but leaves DDT essentially unchanged. Shifts must be per-compound.
9. **`insect_nAChR` is not one receptor.** Croset et al. 2018 show α1/α5/α6/α7
   and β1 dominate the expressing population while α2/α3/α4 are sparse. A
   single `insect_nAChR` key cannot represent the fact that spinosad (α6) and
   the neonicotinoids/sulfoximines (β1-dependent) act on different, partly
   overlapping cell sets. This is a structural limitation of the library
   schema, not a wrong number.
10. **Super-kdr reduces efficacy, not only affinity.** Vais et al. 2000 report
    that the mutation cuts the number of deltamethrin binding sites per channel
    from two to one. A pure EC50 shift under-represents the allele; the maximal
    Nav effect should also be capped if efficacy is modelled.

Confirmed rather than contradicted: the library's `deltamethrin` /
`permethrin` / `ddt` insect Nav values (4.3e-8, 4.0e-7, 6.5e-5 M) match Burton
et al. 2011 exactly; picrotoxin's lack of insect-over-vertebrate selectivity
matches Ratra & Casida 2001; and the `class_placeholder` for diazepam at
`insect_RDL` is directly supported by Hosie & Sattelle 1996.

## Citations

Every source used across the six files, once each.

### Resistance alleles and target-site pharmacology

- Bass C, Puinean AM, Andrews M, Cutler P, Daniels M, Elias J, Paul VL, Crossthwaite AJ, Denholm I, Field LM, Foster SP, Lind R, Williamson MS, Slater R (2011) Mutation of a nicotinic acetylcholine receptor β subunit is associated with resistance to neonicotinoid insecticides in the aphid *Myzus persicae*. *BMC Neurosci* 12:51. DOI 10.1186/1471-2202-12-51. PMID 21627790.
- ffrench-Constant RH, Rocheleau TA, Steichen JC, Chalmers AE (1993) A point mutation in a *Drosophila* GABA receptor confers insecticide resistance. *Nature* 363:449–451. DOI 10.1038/363449a0.
- Garrood WT, Zimmer CT, Gutbrod O, Lüke B, Williamson MS, Bass C, Nauen R, Davies TGE (2017) Influence of the RDL A301S mutation in the brown planthopper *Nilaparvata lugens* on the activity of phenylpyrazole insecticides. *Pestic Biochem Physiol* 142:1–8. DOI 10.1016/j.pestbp.2017.01.007. PMID 29107231.
- Hosie AM, Baylis HA, Buckingham SD, Sattelle DB (1995) Actions of the insecticide fipronil, on dieldrin-sensitive and -resistant GABA receptors of *Drosophila melanogaster*. *Br J Pharmacol* 115(6):909–912. DOI 10.1111/j.1476-5381.1995.tb15896.x. PMID 7582519.
- Hosie AM, Sattelle DB (1996) Allosteric modulation of an expressed homo-oligomeric GABA-gated chloride channel of *Drosophila melanogaster*. *Br J Pharmacol* 117(6):1229–1237. DOI 10.1111/j.1476-5381.1996.tb16720.x. PMID 8882620.
- Kwon DH, Yoon KS, Clark JM, Lee SH (2010) A point mutation in a glutamate-gated chloride channel confers abamectin resistance in the two-spotted spider mite, *Tetranychus urticae* Koch. *Insect Mol Biol* 19:583–591. DOI 10.1111/j.1365-2583.2010.01017.x. PMID 20522121.
- Le Goff G, Hamon A, Bergé JB, Amichot M (2005) Resistance to fipronil in *Drosophila simulans*: influence of two point mutations in the RDL GABA receptor subunit. *J Neurochem* 92:1295–1305. DOI 10.1111/j.1471-4159.2004.02922.x. PMID 15748149.
- Menozzi P, Shi MA, Lougarre A, Tang ZH, Fournier D (2004) Mutations of acetylcholinesterase which confer insecticide resistance in *Drosophila melanogaster* populations. *BMC Evol Biol* 4:4. DOI 10.1186/1471-2148-4-4. PMID 15018651. PMCID PMC362867.
- Mutero A, Pralavorio M, Bride JM, Fournier D (1994) Resistance-associated point mutations in insecticide-insensitive acetylcholinesterase. *Proc Natl Acad Sci USA* 91(13):5922–5926. DOI 10.1073/pnas.91.13.5922. PMID 8016090.
- Usherwood PNR, Vais H, Khambay BPS, Davies TGE, Williamson MS (2005) Sensitivity of the *Drosophila para* sodium channel to DDT is not lowered by the super-kdr mutation M918T on the IIS4–S5 linker that profoundly reduces sensitivity to permethrin and deltamethrin. *FEBS Lett* 579(28):6317–6325. DOI 10.1016/j.febslet.2005.09.096. PMID 16263118.
- Vais H, Williamson MS, Goodson SJ, Devonshire AL, Warmke JW, Usherwood PNR, Cohen CJ (2000) Activation of *Drosophila* sodium channels promotes modification by deltamethrin: reductions in affinity caused by knock-down resistance mutations. *J Gen Physiol* 115(3):305–318. DOI 10.1085/jgp.115.3.305. PMID 10694259.
- Wang X, Puinean AM, O'Reilly AO, Williamson MS, Smelt CLC, Millar NS, Wu Y (2017) Mutations on M3 helix of *Plutella xylostella* glutamate-gated chloride channel confer unequal resistance to abamectin by two different mechanisms. *Insect Biochem Mol Biol* 86:50–57. DOI 10.1016/j.ibmb.2017.05.006. PMID 28576654.
- Zhang HG, ffrench-Constant RH, Jackson MB (1994) A unique amino acid of the *Drosophila* GABA receptor with influence on drug sensitivity by two mechanisms. *J Physiol* 479(Pt 1):65–75. DOI 10.1113/jphysiol.1994.sp020278. PMID 7527461.

### Potency rank orders and selectivity

- Cartereau A, Martin C, Thany SH (2018) Neonicotinoid insecticides differently modulate acetycholine-induced currents on mammalian α7 nicotinic acetylcholine receptors. *Br J Pharmacol* 175(11):1987–1998. DOI 10.1111/bph.14018.
- Cole LM, Nicholson RA, Casida JE (1993) Action of phenylpyrazole insecticides at the GABA-gated chloride channel. *Pestic Biochem Physiol* 46(1):47–54. DOI 10.1006/pest.1993.1035.
- Dawson GR, Wafford KA, Smith A, Marshall GR, Bayley PJ, Schaeffer JM, Meinke PT, McKernan RM (2000) Anticonvulsant and adverse effects of avermectin analogs in mice are mediated through the γ-aminobutyric acid A receptor. *J Pharmacol Exp Ther* 295(3):1051–1060. PMID 11082440.
- Lees K, Musgaard M, Suwanmanee S, Buckingham SD, Biggin P, Sattelle D (2014) Actions of agonists, fipronil and ivermectin on the predominant in vivo splice and edit variant (RDLbd, I/V) of the *Drosophila* GABA receptor expressed in *Xenopus laevis* oocytes. *PLoS ONE* 9(5):e97468. DOI 10.1371/journal.pone.0097468.
- Perry T, Heckel DG, McKenzie JA, Batterham P (2012) Effects of mutations in *Drosophila* nicotinic acetylcholine receptor subunits on sensitivity to insecticides targeting nicotinic acetylcholine receptors. *Pestic Biochem Physiol* 102(1):56–60. DOI 10.1016/j.pestbp.2011.10.010.
- Ratra GS, Casida JE (2001) GABA receptor subunit composition relative to insecticide potency and selectivity. *Toxicol Lett* 122(3):215–222. DOI 10.1016/s0378-4274(01)00366-6. PMID 11489356.
- Ratra GS, Kamita SG, Casida JE (2001) Role of human GABA_A receptor β3 subunit in insecticide toxicity. *Toxicol Appl Pharmacol* 172(3):233–240. DOI 10.1006/taap.2001.9154. PMID 11312652.
- Shan Q, Haddrill JL, Lynch JW (2001) Ivermectin, an unconventional agonist of the glycine receptor chloride channel. *J Biol Chem* 276(16):12556–12564. DOI 10.1074/jbc.M011264200.
- Tasman K, Rands SA, Hodge JJL (2021) The power of *Drosophila melanogaster* for modeling neonicotinoid effects on pollinators and identifying novel mechanisms. *Front Physiol* 12:659440. DOI 10.3389/fphys.2021.659440.
- Tomizawa M, Casida JE (2003) Selective toxicity of neonicotinoids attributable to specificity of insect and mammalian nicotinic receptors. *Annu Rev Entomol* 48:339–364. DOI 10.1146/annurev.ento.48.091801.112731. PMID 12208819.
- Tomizawa M, Casida JE (2005) Neonicotinoid insecticide toxicology: mechanisms of selective action. *Annu Rev Pharmacol Toxicol* 45:247–268. DOI 10.1146/annurev.pharmtox.45.120403.095930.
- Tomizawa M, Lee DL, Casida JE (2000) Neonicotinoid insecticides: molecular features conferring selectivity for insect versus mammalian nicotinic receptors. *J Agric Food Chem* 48(12):6016–6024. DOI 10.1021/jf000873c. PMID 11312774.
- Wolstenholme AJ (2012) Glutamate-gated chloride channels. *J Biol Chem* 287(48):40232–40238. DOI 10.1074/jbc.R112.406280.

### Single-cell atlases and receptor distribution

- Allen AM, Neville MC, Birtles S, Croset V, Treiber CD, Waddell S, Goodwin SF (2020) A single-cell transcriptomic atlas of the adult *Drosophila* ventral nerve cord. *eLife* 9:e54074. DOI 10.7554/eLife.54074. PMID 32314735. PMCID PMC7173974.
- Croset V, Treiber CD, Waddell S (2018) Cellular diversity in the *Drosophila* midbrain revealed by single-cell transcriptomics. *eLife* 7:e34550. DOI 10.7554/eLife.34550. PMID 29671739. PMCID PMC5927767.
- Davie K, Janssens J, Koldere D, De Waegeneer M, Pech U, Kreft Ł, Aibar S, Makhzami S, Christiaens V, Bravo González-Blas C, Poovathingal S, Hulselmans G, Spanier KI, Moerman T, Vanspauwen B, Geurs S, Voet T, Lammertyn J, Thienpont B, Liu S, Konstantinides N, Fiers M, Verstreken P, Aerts S (2018) A single-cell transcriptome atlas of the aging *Drosophila* brain. *Cell* 174(4):982–998.e20. DOI 10.1016/j.cell.2018.05.057. PMID 29909982. PMCID PMC6086935.
- Enell L, Hamasaka Y, Kolodziejczyk A, Nässel DR (2007) γ-Aminobutyric acid (GABA) signaling components in *Drosophila*: immunocytochemical localization of GABA_B receptors in relation to the GABA_A receptor subunit RDL and a vesicular GABA transporter. *J Comp Neurol* 505(1):18–31. DOI 10.1002/cne.21472. PMID 17729251.
- Konstantinides N, Kapuralin K, Fadil C, Barboza L, Satija R, Desplan C (2018) Phenotypic convergence: distinct transcription factors regulate common terminal features. *Cell* 174(3):622–635.e13. DOI 10.1016/j.cell.2018.05.021. PMID 29909983.
- Li H, Janssens J, De Waegeneer M, et al. (2022) Fly Cell Atlas: a single-nucleus transcriptomic atlas of the adult fruit fly. *Science* 375(6584):eabk2432. DOI 10.1126/science.abk2432. PMID 35239393. PMCID PMC8944923.

### Mixtures and synergy

- Bliss CI (1939) The toxicity of poisons applied jointly. *Ann Appl Biol* 26(3):585–615. DOI 10.1111/j.1744-7348.1939.tb06990.x.
- Greco WR, Bravo G, Parsons JC (1995) The search for synergy: a critical review from a response surface perspective. *Pharmacol Rev* 47(2):331–385. PMID 7568331. *(no DOI assigned)*
- Huang Q, Deng Y, Zhan T, He Y (2010) Synergistic and antagonistic effects of piperonyl butoxide in fipronil-susceptible and resistant rice stem borers, *Chilo suppressalis*. *J Insect Sci* 10:182. DOI 10.1673/031.010.14142. PMID 21062143.
- Loewe S, Muischnek H (1926) Über Kombinationswirkungen: 1. Mitteilung: Hilfsmittel der Fragestellung. *Naunyn-Schmiedebergs Arch Exp Pathol Pharmakol* 114:313–326. DOI 10.1007/BF01952257.
- Loewe S (1953) The problem of synergism and antagonism of combined drugs. *Arzneimittelforschung* 3(6):285–290. PMID 13081480.
- Philippou D, Borzatta V, Capparella E, Moroni L, Field L, Moores G (2016) The use of substituted alkynyl phenoxy derivatives of piperonyl butoxide to control insecticide-resistant pests. *Pest Manag Sci* 72(10):1946–1950. DOI 10.1002/ps.4234. PMID 26800141.
- Zhu YC, Yao J, Adamczyk J, Luttrell R (2017) Synergistic toxicity and physiological impact of imidacloprid alone and binary mixtures with seven representative pesticides on honey bee (*Apis mellifera*). *PLoS ONE* 12(5):e0176837. DOI 10.1371/journal.pone.0176837. PMID 28467462.

### Pharmacokinetics and metabolism

- Daborn PJ, Yen JL, Bogwitz MR, Le Goff G, Feil E, Jeffers S, Tijet N, Perry T, Heckel D, Batterham P, Feyereisen R, Wilson TG, ffrench-Constant RH (2002) A single P450 allele associated with insecticide resistance in *Drosophila*. *Science* 297(5590):2253–2256. DOI 10.1126/science.1074170. PMID 12351787.
- Folk DG, Han C, Bradley TJ (2001) Water acquisition and partitioning in *Drosophila melanogaster*: effects of selection for desiccation-resistance. *J Exp Biol* 204(19):3323–3331. DOI 10.1242/jeb.204.19.3323. PMID 11606606.
- Fusetto R, Denecke S, Perry T, O'Hair RAJ, Batterham P (2017) Partitioning the roles of CYP6G1 and gut microbes in the metabolism of the insecticide imidacloprid in *Drosophila melanogaster*. *Sci Rep* 7:11339. DOI 10.1038/s41598-017-09800-2. PMID 28900131. PMCID PMC5595926.
- Joußen N, Heckel DG, Haas M, Schuphan I, Schmidt B (2008) Metabolism of imidacloprid and DDT by P450 CYP6G1 expressed in cell cultures of *Nicotiana tabacum* suggests detoxification of these insecticides in Cyp6g1-overexpressing strains of *Drosophila melanogaster*, leading to resistance. *Pest Manag Sci* 64(1):65–73. DOI 10.1002/ps.1472. PMID 17912692.
- Piyankarage SC, Augustin H, Grosjean Y, Featherstone DE, Shippy SA (2008) Hemolymph amino acid analysis of individual *Drosophila* larvae. *Anal Chem* 80(4):1201–1207. DOI 10.1021/ac701785z.
- Piyankarage SC, Featherstone DE, Shippy SA (2012) Nanoliter hemolymph sampling and analysis of individual adult *Drosophila melanogaster*. *Anal Chem* 84(10):4460–4466. DOI 10.1021/ac3002319.
- Suchail S, Debrauwer L, Belzunces LP (2004) Metabolism of imidacloprid in *Apis mellifera*. *Pest Manag Sci* 60(3):291–296. DOI 10.1002/ps.772. PMID 15025241.
- Suchail S, De Sousa G, Rahmani R, Belzunces LP (2004) In vivo distribution and metabolisation of ^14C-imidacloprid in different compartments of *Apis mellifera* L. *Pest Manag Sci* 60(11):1056–1062. DOI 10.1002/ps.895.

### Behavioural assays

- Cao W, Song L, Cheng J, Yi N, Cai L, Huang FD, Ho MS (2017) An automated rapid iterative negative geotaxis assay for analyzing adult climbing behavior in a *Drosophila* model of neurodegeneration. *J Vis Exp* (127):56507. DOI 10.3791/56507. PMID 28931001.
- Çevik MÖ, Erden A (2012) The course of habituation of the proboscis extension reflex can be predicted by sucrose responsiveness in *Drosophila*. *PLoS ONE* 7(6):e39863. DOI 10.1371/journal.pone.0039863. PMID 22761915.
- Chen W, Gu X, Yang YT, Batterham P, Perry T (2022) Dual nicotinic acetylcholine receptor subunit gene knockouts reveal limits to functional redundancy. *Pestic Biochem Physiol* 184:105118. DOI 10.1016/j.pestbp.2022.105118. PMID 35715057.
- French AS, Sellier MJ, Ali Agha M, Guigue A, Chabaud MA, Reeb PD, Mitra A, Grau Y, Soustelle L, Marion-Poll F (2015) Dual mechanism for bitter avoidance in *Drosophila*. *J Neurosci* 35(9):3990–4004. DOI 10.1523/JNEUROSCI.1312-14.2015. PMID 25740527.
- Gargano JW, Martin I, Bhandari P, Grotewiel MS (2005) Rapid iterative negative geotaxis (RING): a new method for assessing age-related locomotor decline in *Drosophila*. *Exp Gerontol* 40(5):386–395. DOI 10.1016/j.exger.2005.02.005. PMID 15919590.
- Lee HY, Zhou L, Ghanta S, Asplund M, Zhu C, et al. (2014) Mechanisms of naturally evolved ethanol resistance in *Drosophila melanogaster*. *J Exp Biol*. DOI 10.1242/jeb.110510. *(cited only as the source of a null in `fly_pharmacokinetics.yaml`)*
- Martelli F, Zhongyuan Z, Wang J, Wong CO, Karagas NE, Roessner U, Rupasinghe T, Venkatachalam K, Perry T, Bellen HJ, Batterham P (2020) Low doses of the neonicotinoid insecticide imidacloprid induce ROS triggering neurological and metabolic impairments in *Drosophila*. *Proc Natl Acad Sci USA* 117(41):25840–25850. DOI 10.1073/pnas.2011828117. PMID 32989137.

### Cited in `flylab/pharm/library.yaml` and re-checked here

- Burton MJ, Mellor IR, Duce IR, Davies TGE, Field LM, Williamson MS (2011) Differential resistance of insect sodium channels with kdr mutations to deltamethrin, permethrin and DDT. *Insect Biochem Mol Biol* 41(9):723–732. DOI 10.1016/j.ibmb.2011.05.004. *(confirms the library's deltamethrin 4.3e-8, permethrin 4.0e-7 and DDT 6.5e-5 M values at insect Nav)*
