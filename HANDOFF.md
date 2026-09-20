# FlyLab handoff (v0.6.1, 2026-09-20)

Read this first if you are Claude, GPT, or another agent continuing the repo.
Owner: https://github.com/pinkysworld/FlyLab
Package version: **0.6.0** in `CITATION.cff` and in the paper. **`flylab/__init__.py` and `pyproject.toml` both still say `0.5.0`** — see Issues; those files are not mine to edit.

## What the project is

A **virtual pharmacology bench** on the public **MaleCNS v1.0** adult male *Drosophila* CNS (brain + ventral nerve cord). One compound + one free concentration → a **typed** insect engagement panel + a vertebrate panel at the same dose → gain patch on named MaleCNS cells → circuit readouts → one notebook JSON with provenance.

Since v0.6 it is also, and mainly, **a method for disbelieving its own predictions**: connectome-dependence analysis, an ablation ladder, a conclusion-stability matrix over alternative gain specifications, and a variance-based uncertainty budget that names the next experiment.

**The manuscript is framed for a computing venue.** IJRC is a computing research journal, so the paper leads with the general problem (a network simulation cannot say which of its inputs a prediction depends on), presents four transferable instruments, evaluates *those instruments*, and treats the fly pharmacology as the demonstration domain with the domain detail in the supplement. Do not re-centre the paper on the biology.

**v0.6.1 pointed those instruments at themselves and the central empirical claim reversed.** Two of them were broken in ways that flattered the conclusions (a transmitter null that moved the weighted excitation/inhibition balance; a specification family whose alternative gain rules never reached the engine its topology predicates ran on), the dependence ladder is now validated against planted ground truth, and repeating the landscape on the denser `taste_motor` cut gives the opposite verdict to the sparse `named` cut. Read "Findings that constrain what you may claim" below before writing any new claim.

It is **not** a new connectome, not GLP tox, not a 166k-cell LIF product, and not live-animal results.

| Document | What it is for |
|---|---|
| `docs/RESEARCH_MAP.md` | the canonical plan, gates and what v0.7 is |
| `docs/NOVELTY.md` | what is actually new — and what is **not** (read before accepting a PR) |
| `docs/ARCHITECTURE.md` | dataflow and layer ownership |
| `docs/STACK.md` | technology choices and the browser constraints they impose |
| `docs/PAGES.md` | the static (Pyodide) bench |
| `papers/IJRC_FlyLab_draft.md` | the manuscript (generated — edit the `.in` template) |
| `papers/SUPPLEMENT.md` | the supplement (also generated from a `.in` template) |
| `papers/STATUS.md`, `papers/SUBMISSION_CHECKLIST.md` | where the paper stands, and the point-by-point reviewer response |

## Do this first

```bash
git pull
python -m pip install -e ".[dev,viz]"
python -m pytest -q                        # fast suite; `-m slow` adds the long checks
flylab occupancy imidacloprid --conc 1e-6
python scripts/reproduce_paper.py --fast   # ~4 min sanity check of the whole pipeline (a default run is ~35-45 min)
flylab serve                               # http://127.0.0.1:8765
```

## What changed in v0.6 and v0.6.1

**Pharmacology — the evidence is typed (schema v3)**
- Every row records `param_type` (what the source measured: `Kd`/`Ki`/`EC50`/`IC50`/`Kb`/`relative_potency`/`class_order`/`unknown`) and `relation` (how far the source is from this compound/receptor/species). The relation is graded as an ordered **evidence distance** (`E0` exact compound+receptor+species, `E1` other species, `E2` related receptor, `E3` class extrapolation, `E4` unsupported), and the two together decide the **engagement model**: `binding_occupancy` (Kd/Ki at E0 only), `binding_engagement_proxy` (Kd/Ki at E1/E2), `functional_engagement` (a potency at E0), `functional_engagement_proxy` (any transferred number), or `not_modelled`. Proxy rows carry a `provenance_warning` naming what was transferred; `flylab.pharm.occupancy.evidence_distance_table()` is the census by distance.
- **Unsupported rows return N/A**, never a small number, and `flylab.pharm.evidence.check_transformation` *raises* if asked anyway. Diazepam at insect RDL is N/A, not 1e-4.
- Say **engagement**, not occupancy, unless the row carries a Kd measured on *this* compound, receptor **and** species. Exactly one does, and it is not the aphid row: Tomizawa, Latli and Casida 1996 on *Drosophila* head membranes, on the preparation-resolved key `insect_nAChR_native_dmel`. The aphid `insect_nAChR_beta1` Kd is a **cross-species proxy** — same number, different name, provenance warning attached.
- Library: 102 rows / 21 compounds / 15 receptor keys. 34 EC50, 10 IC50, 2 Kd, 56 not modelled. Evidence distance E0–E4: 10 / 15 / 10 / 11 / 56. Exactly one row earns `binding_occupancy` (imidacloprid at `insect_nAChR_native_dmel`, Kd 2e-09 M, a *Drosophila* head-membrane measurement); the aphid `insect_nAChR_beta1` Kd is demoted to a `binding_engagement_proxy`.
- Nicotinic key split by subunit: spinosad → α6, imidacloprid → β1; the rest stay on a documented aggregate because their sources used hybrid Dα/chicken-β2 constructs. Source correction: **Dederer 2011 is a cat-flea α1/α2 + chicken β2 hybrid, not a *Drosophila* receptor**, and the row now says so.

**Analysis — five new modules**
- `analysis/dependence.py` — connectome-dependence analysis: permutation *p* first, the profile, the dependence class, and the **necessary information level** (the weakest degraded graph the test could **not** distinguish from the real cut — never "reproduces"). v0.6.1 adds a three-way verdict against a prespecified equivalence margin, a weight-matched transmitter null on the ladder's rank-3 rung, a relative effect floor, Benjamini-Hochberg across a landscape's structural tests with confirmatory/exploratory labelling, `cut_census` for the substrate's own statistics, and `synthetic_cut` / `ladder_recovery` / `ladder_power` for ground-truth validation of the instrument.
- `analysis/baselines.py` — the ablation ladder: receptor-only / composition-only / topology-only / full, compared by ordering.
- `analysis/robustness.py` — 25 prespecified gain specifications × 7 conclusions, and the 3×3 amplify/buffer threshold grid.
- `analysis/uncertainty_global.py` — Sobol' budget with Jansen estimators, a null factor that measures the noise floor, and Ishigami validation.
- `analysis/voi.py` — `VOI_j = S_j · Var(Y)` mapped onto concrete experiments.
- `analysis/claims.py` (added by the dashboard agent, used by the pipeline) — `claim_audit()` walks the nine-link dependency chain behind a readout and labels each link OBSERVED / LITERATURE-DERIVED / MODEL-ASSUMPTION / COMPUTED; exactly one link (the MaleCNS edges) is a measurement of the system being simulated. `fact_inference_unknown()` gives the facts / model-inference / unknown split. Table T18.

**Paper and pipeline**
- The manuscript is now a **research article** (Introduction / Methods / Results / **Discussion** / Conclusion) organised around four research questions, with `papers/SUPPLEMENT.md` holding the mechanism rationale, the null-model definitions, the LIF calibration, the supporting capabilities and the prospective predictions. Body is inside the 5000–7000 target and a test enforces it.
- `scripts/reproduce_paper.py` now produces 16 figures and T0–T25, including T12b (the landscape on `taste_motor`), T19 (both cuts' structure), T20 (the transformation rule), T21 (the transmitter signs), T22 (instrument validation), T23 (transmitter-null balance), T24 (the composition reference distribution) and T25 (the normalisation sweep). Default ~35–45 min, `--fast` ~4 min. A test asserts the committed record was produced at the shipped statistical effort, not merely that it is self-consistent.
- "Pre-registered" is gone everywhere: the predictions are **prospective** until a tagged release is externally archived.

## Findings that constrain what you may claim

Read these before writing any new claim into the paper, the README or a PR description.

1. **The dependence verdict is a joint property of the prediction and the extract, and recurrence predicts it. That is the headline.** The same landscape, same instrument, no parameter changed: on the 1-hop `named` cut (mean degree 1.208, 85% of edges onto four seed cells, 184 of 1126 nodes with any input) 11 of 84 cells are topology-dependent, 42 composition-dominated and 31 no-effect at n = 1000 with Benjamini-Hochberg and a 1 % relative effect floor. On `taste_motor` (mean degree 10.4) the same cells give 54 / 0 / 30. The **scaling study** (`papers/scale_study.json`, T27, `flylab.analysis.scale`) then settles what the inversion tracks: over seven extracts the composition-dominated verdict survives on **exactly one**, the in-star, and `scale_1k` — *fewer* nodes than `named` (1000 vs 1126) but mean degree 22.9 — is already topology-dependent. **Recurrence, not size.** At 5k and 10k imidacloprid is distinguishable from *every* null; the top rungs run at n = 20 (resolution 0.048) and so confirm rather than establish. "Most predictions do not need the connectome" and "the topology-dependent set is exactly the chloride-channel blockers" are **withdrawn**. Never state a dependence result without the extract, its structure, the compound, the concentration, the readout, the *p* and the correction.
   1a. **Non-rejection is not equivalence.** Each mode carries a three-way verdict against a prespecified margin (0.332 Hz, 5% of the vehicle readout on `named`). For imidacloprid at n = 1000: equivalent within tolerance for rewire_degree_preserving, sign_permute_weight_matched, weight_permute, indeterminate for sign_permute. Never write "reproduces".
   1b. **The transmitter null was biased.** Plain label permutation moves the cholinergic share of out-weight from 0.619 to 0.544 ± 0.027, putting the real graph at the 100th percentile of its own null. It is now reported as a joint target-set-and-sign null; the ladder's rank-3 rung is the weight-matched shuffle.
   1c. **The instrument is validated.** Planted ground truth on a synthetic cut: the ladder recovered 2 of 3 planted topology-dependent effects at n = 200 and returned no false positive on the unplanted control. False-positive rate 0 of 18; detection by permutation count n=50: 0.61, n=200: 0.61, n=1000: 0.67, so effect size dominates.
2. **The bitter-veto arm is a properly powered negative** (n = 300: p = 0.42 / 0.88 / 0.20 / 0.24, all |z| < 0.25). H6 is a pharmacology prediction, not a wiring prediction. Imidacloprid's veto ratio there is **undefined** (MN9 silenced), not zero.
3. **The suppression result is a property of the gain rule — and the old stability matrix was vacuous.** "A nicotinic agonist suppresses circuit activity" survives 15/25 specifications and is reversed by all 10 monotone ones, where the same drug at the same engagement *excites* the network. The v0.6 topology rows were computed at 6 shuffles (p ≤ 0.05 unattainable) **and** under a context that never rebound the gain function on the engine the permutation path resolves, so every specification fed default gains to its nulls. Recomputed at 100 shuffles with the specification pushed into the null engine and FDR across 100 tests, C3 holds 25/25 and C4 25/25 — but only 17 of C3's retentions reach equivalence and 8 are indeterminate. Retained by every specification: C2_rdl_disinhibition, C3_imidacloprid_topology_not_distinguishable, C4_fipronil_topology_exceeds, C7_map_bitter_veto.
4. **The ablation conclusion is withdrawn, and the composition correlation needs its reference.** Level C's old generic multiplier could only depress, so it could not express disinhibition; with a direction-aware rule its rank correlation at 1 µM rises from 0.242 to 0.961, and "the connectome without the pharmacology is not a cheap substitute" is retracted. The composition-versus-full correlation is 0.988, but pharmacology-free pseudo-compounds already reach a median of 0.847 and shuffling compound labels leaves it identically 0.988 — it carries no compound-level information. It is also not robust to the engine: 0.666 at 1 µM and -0.372 at 10 nM under a degree-corrected denominator. Never quote it without its reference distribution and the normalisation it was computed under.
5. **The uncertainty lives in two assumptions**: gain_transform, weight_threshold are the only factors whose bootstrap first-order interval excludes zero (S1 0.410 and 0.294). The estimator's noise floor is 0.046, set by gain_coef and **not** by the declared null factor `lif_seed` (-0.007); `drive` sits inside it and is no longer a ranked finding. Potency and Hill *n* are invisible at a saturating dose. Interactions carry 0.344 raw and 0.272 with negative estimates clipped.
6. **The amplify/buffer split is threshold-sensitive** (`stable = False`): quote the range over the 3×3 grid, never the point estimate.
7. **The rate and LIF engines agree on direction only** (sugar-driven MN9 ≈ 0.29 Hz vs ≈ 105 Hz). That is not cross-model validation and must not be called one.
8. **There is no independent out-of-sample validation of the circuit model.** Rank comparisons are literature *concordance*, and the pipeline flags which of them share a source with the library.
9. **`data/literature/README.md` lists ten literature-vs-library contradictions.** Report them; do not silently edit the library to match.
10. **Five transmitter sign magnitudes are asserted** (glutamate -0.4, histamine -0.5, dopamine +0.2, serotonin +0.2, octopamine +0.2), vary in no specification family, and are not bounded by the label-permuting null. They are published as T21; quote them as assertions.
11. **The rate engine row-normalises its weight matrix** (`row_abs`), so each cell's recurrent input is a composition-weighted average of its presynaptic gains. Any composition-versus-topology statement must name the normalisation it was computed under.
12. **The literature-concordance mean is not quotable bare**: 4 of 8 comparisons order two compounds. Use the mean over the 4 informative ones (0.756) or the concordance count (8 vs 0, binomial p = 0.008).

## Issues

- #1 Gate 0 occupancy CLI — **closed**. #2 Gate 1 map-extracted taste GRNs — **closed in v0.5**. #3 Gate 2 five-drug panel — **closed**.
- **Open: version strings disagree.** `flylab/__init__.py` and `pyproject.toml` both report `0.5.0`; `CITATION.cff` and the manuscript say `0.6.0`. Pick one and bump the two package files. Not fixed here because those files belong to other owners.
- **Open: expression coverage is 0.205** and adult motor-neuron receptor expression is a confirmed literature gap — MN9, the main readout, carries the `unknown` default weight.
- **Open (new): `flylab/analysis/nullmodels.py` still exposes `connectome_information_score`**, a v0.5 summary superseded by the dependence class and the necessary information level. The paper no longer quotes it. Either delete it or document it as deprecated.
- **Open (new): docstrings in `flylab/` still say "occupancy" in places where the value is a functional engagement** (e.g. `flylab/pharm/occupancy.py` module name and several assay docstrings, `hill_occupancy` kept as a deprecated alias). The behaviour is correct and typed; the naming lags the paper.
- **Open (new): `flylab/analysis/predictions.py` still labels H1–H7 "pre-registered"** in its docstring and output fields. The paper now says *prospective*. Rename when convenient.

## Hard constraints (do not break)

- Organism is **fly** (*Drosophila*); map is **MaleCNS v1.0** (MN9/DNp01 need the cord, so not FlyWire-only).
- **Patch rules may not change silently.** `insect nAChR agonist: g_ach = max(0.05, 1 + 0.4θ − 1.6θ²)`; `insect nAChR antagonist: g_ach = max(0.05, 1 − θ)`; `insect RDL antagonist: g_gaba = max(0.05, 1 − θ)`; `insect RDL agonist/PAM: g_gaba = max(0.05, 1 + 0.4θ)`. Changing one means editing `flylab/pharm/mechanisms.py`, its rationale docstring, T2 **and the conclusion-stability family**, and saying so.
- **No invented potency values.** New rows need a `source`, an `evidence_tier`, a `param_type` and a `relation`. A missing number stays `unknown`/`unsupported` and returns N/A.
- **At most one fitted parameter later**, and the variance budget says which. Never refit MaleCNS weights.
- **Notebooks keep their warnings** when the circuit is reduced or hops-limited.
- **`live_lab` stays `null`.** No code path may write a live-animal number. `protocols/` are unexecuted drafts.
- **Genotype shifts are per compound, never per receptor.**
- **Never commit `*.feather`.** Do not add hops to a cut; raise the weight floor. Do not hand-edit `data/derived/*.json`.
- The science core must not import `pydantic`, `fastapi`, `typer`, `pandas` or `pyarrow` at module scope — the browser build has none of them.
- **Never write "pre-registered"** unless a tagged release is externally archived. A test enforces this on the manuscript template.

## Big files (never commit)

From `https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome`:
`body-annotations-male-cns-v1.0-minconf-0.5.feather`, `body-neurotransmitters-male-cns-v1.0.feather`, `connectome-weights-male-cns-v1.0-minconf-0.5.feather` (~1.1 GB).
Local dir: `$FLYLAB_MALECNS` or `~/.flylab/malecns_v1`. `.gitignore` drops `*.feather`. Stream with `pyarrow.ipc.open_file` batches — a full pandas read OOMs. **None of this is needed to reproduce the paper.**

## Commands that matter

```bash
# pharmacology
flylab list-drugs
flylab occupancy imidacloprid --conc 1e-6
flylab compare imidacloprid fipronil --conc 1e-6

# circuits
flylab assay-subgraph  --compound imidacloprid --conc 1e-6     # rate, named cut
flylab assay-spiking   --compound fipronil     --conc 1e-6     # LIF
flylab assay-taste-map --compound fipronil     --conc 1e-6     # real GRNs -> MN9
flylab assay-cns       --compound imidacloprid --conc 1e-6     # census

# analysis (v0.6 layers are in flylab.analysis; see each module's docstring)
python -c "from flylab.analysis.dependence import dependence_profile as d; print(d('fipronil',1e-6,n=200,n_jobs=4)['class'])"
python -c "from flylab.analysis.baselines import ablation_table as a; print(a()['statements'])"
python -c "from flylab.analysis.robustness import threshold_sensitivity as t; print(t()['stable'])"
flylab selectivity --conc 1e-6
flylab predictions

# paper
python scripts/reproduce_paper.py              # full, ~35-45 min, 16 figures / T0-T25
python scripts/reproduce_paper.py --fast       # ~2 min sanity check, NOT the paper
python scripts/reproduce_paper.py --only dependence --outdir /tmp/x
python scripts/reproduce_paper.py --only paper # re-render prose from results.json
python scripts/reproduce_paper.py --list
```

## Files not to rewrite from scratch

`flylab/maps/extract.py` (streaming), `flylab/pharm/mechanisms.py` (the patch rules), `flylab/pharm/evidence.py` (the type system the paper's RQ1 rests on), `flylab/pharm/library.yaml` (cite sources if you change a value), `data/derived/*.json` (CI products), `data/literature/*.yaml` (sourced datasets), `papers/results.json` (generated), `papers/*.md` (generated — edit the `.in` templates).

## Next work, in order

1. **v0.7 gate 7 — expression coverage.** Either source adult motor-neuron receptor expression or declare the 0.205 coverage a permanent bound. A data problem, not a modelling one.
2. **The one allowed fit.** The uncertainty budget names it: the engagement→gain transformation. Calibrate it against evoked postsynaptic responses at a known synapse; do not fit a potency value.
3. **The free win.** Re-run the dependence and selectivity analyses against synapse-confidence strata: `weight_threshold` carries S1 0.294 and costs no bench time.
4. **One live assay.** Climbing first: a full published protocol and control statistics exist. Import by hand into `live_lab`; nothing else may write it.
5. **Housekeeping:** reconcile the version strings, mint the Zenodo DOI, and only then consider the predictions pre-registered.

**Do not** start a 166k-cell LIF, add a second connectome, invent live PER/climbing numbers, claim that no tool has combined connectomes and pharmacology (FlyBrainLab and the receptor-map whole-brain models exist and are cited), or accept a PR that only adds a prettier viewer.
