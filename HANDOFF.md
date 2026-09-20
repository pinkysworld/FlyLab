# FlyLab handoff (v0.6, 2026-09-19)

Read this first if you are Claude, GPT, or another agent continuing the repo.
Owner: https://github.com/pinkysworld/FlyLab
Package version: **0.6.0** in `CITATION.cff` and in the paper. **`flylab/__init__.py` and `pyproject.toml` both still say `0.5.0`** — see Issues; those files are not mine to edit.

## What the project is

A **virtual pharmacology bench** on the public **MaleCNS v1.0** adult male *Drosophila* CNS (brain + ventral nerve cord). One compound + one free concentration → a **typed** insect engagement panel + a vertebrate panel at the same dose → gain patch on named MaleCNS cells → circuit readouts → one notebook JSON with provenance.

Since v0.6 it is also, and mainly, **a method for disbelieving its own predictions**: connectome-dependence analysis, an ablation ladder, a conclusion-stability matrix over alternative gain specifications, and a variance-based uncertainty budget that names the next experiment.

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
python scripts/reproduce_paper.py --fast   # ~2 min sanity check of the whole pipeline
flylab serve                               # http://127.0.0.1:8765
```

## What changed in v0.6

**Pharmacology — the evidence is typed (schema v3)**
- Every row records `param_type` (what the source measured: `Kd`/`Ki`/`EC50`/`IC50`/`Kb`/`relative_potency`/`class_order`/`unknown`) and `relation` (how far the source is from this compound/receptor/species). The relation is graded as an ordered **evidence distance** (`E0` exact compound+receptor+species, `E1` other species, `E2` related receptor, `E3` class extrapolation, `E4` unsupported), and the two together decide the **engagement model**: `binding_occupancy` (Kd/Ki at E0 only), `binding_engagement_proxy` (Kd/Ki at E1/E2), `functional_engagement` (a potency at E0), `functional_engagement_proxy` (any transferred number), or `not_modelled`. Proxy rows carry a `provenance_warning` naming what was transferred; `flylab.pharm.occupancy.evidence_distance_table()` is the census by distance.
- **Unsupported rows return N/A**, never a small number, and `flylab.pharm.evidence.check_transformation` *raises* if asked anyway. Diazepam at insect RDL is N/A, not 1e-4.
- Say **engagement**, not occupancy, unless the row carries a measured Kd. Exactly one does: imidacloprid at `insect_nAChR_beta1` (Kd 8.3e-11 M, Bass 2011 saturation binding).
- Library: 101 rows / 21 compounds / 14 receptor keys. 34 EC50, 10 IC50, 1 Kd, 56 not modelled. Relations: 9 exact/exact/exact, 15 other-species, 10 related-receptor, 11 class-extrapolation, 56 unsupported.
- Nicotinic key split by subunit: spinosad → α6, imidacloprid → β1; the rest stay on a documented aggregate because their sources used hybrid Dα/chicken-β2 constructs. Source correction: **Dederer 2011 is a cat-flea α1/α2 + chicken β2 hybrid, not a *Drosophila* receptor**, and the row now says so.

**Analysis — five new modules**
- `analysis/dependence.py` — connectome-dependence analysis: permutation *p* first, the four-mode profile, the dependence class, and the **necessary information level** (the weakest degraded graph that still reproduces the effect), plus a landscape over the whole library.
- `analysis/baselines.py` — the ablation ladder: receptor-only / composition-only / topology-only / full, compared by ordering.
- `analysis/robustness.py` — 25 prespecified gain specifications × 7 conclusions, and the 3×3 amplify/buffer threshold grid.
- `analysis/uncertainty_global.py` — Sobol' budget with Jansen estimators, a null factor that measures the noise floor, and Ishigami validation.
- `analysis/voi.py` — `VOI_j = S_j · Var(Y)` mapped onto concrete experiments.
- `analysis/claims.py` (added by the dashboard agent, used by the pipeline) — `claim_audit()` walks the nine-link dependency chain behind a readout and labels each link OBSERVED / LITERATURE-DERIVED / MODEL-ASSUMPTION / COMPUTED; exactly one link (the MaleCNS edges) is a measurement of the system being simulated. `fact_inference_unknown()` gives the facts / model-inference / unknown split. Table T18.

**Paper and pipeline**
- The manuscript is now a **research article** (Introduction / Methods / Results / **Discussion** / Conclusion) organised around four research questions, with `papers/SUPPLEMENT.md` holding the mechanism rationale, the null-model definitions, the LIF calibration, the supporting capabilities and the prospective predictions. Body is inside the 5000–7000 target and a test enforces it.
- `scripts/reproduce_paper.py` gained `evidence`, `dependence` (replacing `nulls`), `ablation`, `stability`, `uncertainty` and `claims` steps: 15 figures, 21 tables. Default ~8–10 min, `--fast` ~2 min. `--jobs` backgrounds the two dominant analyses; results are identical at any value.
- "Pre-registered" is gone everywhere: the predictions are **prospective** until a tagged release is externally archived.

## Findings that constrain what you may claim

Read these before writing any new claim into the paper, the README or a PR description.

1. **Most predictions do not need the connectome.** Landscape over 21 compounds × 4 concentrations: 20 topology-dependent, 51 composition-dominated, 13 no-effect, 0 mixed. The topology-dependent set is exactly the chloride-channel blockers. Imidacloprid at 1 µM is composition-dominated (p = 0.275 / 0.586 / 0.472 at n = 1000; necessary level `degree_sequence`); fipronil is topology-dependent (p = 0.0060 / 0.0070; necessary level `wiring_without_transmitter_identity`). Never say "the connectome matters" without the compound, the concentration, the readout and the *p*.
2. **The bitter-veto arm is a properly powered negative** (n = 300: p = 0.42 / 0.88 / 0.20 / 0.24, all |z| < 0.25). H6 is a pharmacology prediction, not a wiring prediction. Imidacloprid's veto ratio there is **undefined** (MN9 silenced), not zero.
3. **The suppression result is a property of the gain rule.** "A nicotinic agonist suppresses circuit activity" survives 15/25 specifications and is reversed by all ten monotone ones, where the same drug at the same engagement *excites* the network. Do not lead with it. RDL disinhibition, both topology conclusions and the map bitter-veto direction are 25/25 and may be quoted freely.
4. **Nearly all ordering information is in the mechanism rules**, not the wiring: composition-only reproduces the full model's ordering (ρ 0.906–0.989), topology-only with a generic multiplier is the worst level (ρ 0.24, 11.7 Hz rms).
5. **The uncertainty lives in two assumptions**: `gain_transform` (S1 0.410) and `weight_threshold` (S1 0.294). Potency and Hill *n* are invisible at a saturating dose. The one allowed fit is therefore the gain transformation, and `weight_threshold` can be reduced with **no experiment at all**.
6. **The amplify/buffer split is threshold-sensitive** (`stable = False`): quote the range over the 3×3 grid, never the point estimate.
7. **The rate and LIF engines agree on direction only** (sugar-driven MN9 ≈ 0.29 Hz vs ≈ 105 Hz). That is not cross-model validation and must not be called one.
8. **There is no independent out-of-sample validation of the circuit model.** Rank comparisons are literature *concordance*, and the pipeline flags which of them share a source with the library.
9. **`data/literature/README.md` lists ten literature-vs-library contradictions.** Report them; do not silently edit the library to match.

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
python scripts/reproduce_paper.py              # full, ~8-10 min, 15 figures / 21 tables
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
