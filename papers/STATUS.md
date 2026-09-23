# Manuscript status

**Files:** `papers/IJRC_FlyLab_draft.md` and `papers/SUPPLEMENT.md` — **both generated**. Edit `papers/IJRC_FlyLab_draft.md.in` / `papers/SUPPLEMENT.md.in` and re-render; never hand-edit a number into a rendered document.

**Venue:** IJRC (ijrcom.org), as a **Research Article** (Introduction, Methods, Results, Discussion, Conclusion, References).
**State:** additional internal-review and readability revision of main `38cad67`. Not submitted. Review decisions: `papers/reviews/current_revision_decisions.md` and `papers/reviews/round4_editorial_decisions.md`. Author details are supplied; an editable Word preparation copy is in `papers/submission/`.

**Current internal-review update:** the denser-cut study, synthetic recovery cases, corrected evidence propagation and specification fixes are retained. The current revision narrows the remaining claims: extract sensitivity is observed, recurrence is not isolated causally, the legacy equivalence label is only a median-gap diagnostic, and the small synthetic detection grid does not establish endpoint-specific power. Review reports and decisions are in `papers/reviews/`.

**This review's verification:** the existing default-effort numerical record is retained. The attempted full rerun was interrupted during dependence analysis. `scripts/paper_review_metadata.py` refreshes interpretation and captions while asserting that numerical value payloads and table CSV bytes are unchanged; `results.json.review_metadata` distinguishes this replay from a completed rerun. See `papers/reviews/current_verification.md` for completed checks.

**Reproduction boundary:** the standard pipeline recomputes core committed-cut analyses but imports the previously computed `papers/scale_study.json`. Only `scale_1k` is committed; large-rung results require upstream cuts to recompute and the saved study lacks source/input hashes. Re-rendering its table is not an independent rerun.

```bash
python scripts/reproduce_paper.py                # recompute core analyses and replay the cached scale table
python scripts/reproduce_paper.py --only paper   # re-render both documents from results.json
python -m pytest -q                              # must be green
```

## What the draft contains

| Part | State |
|---|---|
| Abstract (architecture + one topology result + one reproducibility result + the central limitation; individual numbers stripped) | done |
| §1 Introduction: three research questions, related work, and the narrow novelty claim with a comparison table | done; reproduction is supporting evidence rather than RQ4 |
| §2 Methods — typed evidence, the circuit substrate, connectome-dependence analysis, the ablation ladder and specification family, the evaluation taxonomy, global uncertainty and VOI | done (6 subsections) |
| §3 Results — the instruments evaluated first (ground truth and power, then the inversion and the scale ladder), with the domain readouts as supporting material; every number from `results.json` | done (7 subsections) |
| §4 **Discussion** — what was learned computationally (including the substrate dependence of the method's own verdict, and two instruments that were broken in ways that flattered the conclusions); model-generated vs encoded; adjacent work; **what the software cannot support**; five model-specific caveats, at length in §S9; prospective-not-pre-registered | done |
| §5 Conclusion | done |
| Statements — precise AI declaration, animal research, data/code availability, competing interests, funding | done |
| References — 76 entries with DOIs/PMIDs, including FlyBrainLab [67], the receptor-map whole-brain models [68], [69] and the *Drosophila* binding constant [76] | done |
| Supplement — mechanism rationale and the level-C rule (S1), the specification family and why the previous matrix was vacuous (S1.1), null-model definitions with the weight-matched transmitter null (S2), the three-way verdict and the equivalence margin (S2.1), effect floor and multiplicity (S2.2), ground-truth recovery and power (S2.3), the scale ladder and where the verdict settles (S2.4), runtimes and the row normalisation (S3), supporting capabilities (S4), prospective predictions with the H2/H3 caveat (S5), evidence-typing rules and the soundness invariant (S6), reproduction knobs (S7), notebook warnings (S8), model-specific caveats in full (S9), and the domain readouts, selectivity and literature concordance moved out of the body (S10) | done |
| 16 figures (F16 is the instrument validation) + T0–T25, regenerated at 300 dpi | done |

**Body word count is generated, not estimated**: `papers/results.json → values.paper_words_body`. The editorial test guard is 5000–8000. The concise revision is below the journal's typical Research Articles 6000–8000 range but within its general 5000–10000 guidance. Confirm the length rather than pad the text. Final format and blind-review packaging remain author checks.

## Implemented capabilities and committed artifacts

- Typed evidence: 102 rows / 21 compounds / 15 receptor keys; 34 EC50, 10 IC50, **2 Kd**, 56 not modelled and returning N/A, enforced by a tagged-record dispatch discipline and regression checks across all shipped rows and enumerated pipeline paths. Evidence distance E0–E4: 10 / 15 / 10 / 11 / 56.
- Two committed MaleCNS v1.0 cuts (1126/1360 and 1841/19 066) plus a 165 122-cell census fallback; all 16 labellar `LB*` types proven present.
- Connectome-dependence analysis at n = 1000 with three-way descriptive tolerance diagnostics, a weight-matched transmitter null, a permutation-count sweep, synthetic positive-control checks of the ladder, and a landscape classified using BH-adjusted component tests over 21 compounds × 4 concentrations **on both committed cuts**.
- Ablation ladder, conclusion-stability matrix over 25 gain specifications, a 3×3 threshold grid, a Sobol' budget (11 264 evaluations) with Ishigami validation, and a VOI ranking.
- Map-extracted GRN → MN9 path with a complete bitter veto in vehicle and partial relief under fipronil, on both engines.
- Core committed-cut analyses regenerated by the standard command; scale profiles are imported from the saved study, and a test asserts the committed record was produced at the shipped statistical effort.

## Interpretation limits

- Dependence labels change between extracts. Current evidence does not select recurrence over density, extraction recipe or correlated structural properties.
- `equivalent_within_tolerance` is a legacy name for non-rejection plus a small point gap from the null median, not a confidence-bounded equivalence test.
- BH-adjusted values are reported; shared shuffles alone do not prove the required dependence assumptions.
- The taste-motor veto-ratio non-rejection has no endpoint-specific power analysis.
- The six-network-per-setting synthetic grid is an initial recovery check with imprecise detection rates.
- Top scale rungs use 20 permutations and report probabilities at the resolution floor. Matching labels are not proof of convergence.
- High ablation correlation depends on shared assumptions; its observed value also exceeds the reported reference 95th percentile, so the reference does not explain it away entirely.
- Nicotinic suppression reverses under monotone gain rules.
- Rate-model units and LIF spikes/s are not interchangeable measured firing rates.
- No independent circuit-level or live-animal validation exists.

## Still missing for a **results** paper

1. One live PER or climbing table in `live_lab` (climbing first — the control statistics exist; the PER baseline does not).
2. The one allowed fit — and the uncertainty budget now names it without argument: the engagement→gain transformation, not a potency value.
3. Expression coverage above 0.205, with adult motor neurons resolved.
4. A source-pinned archival release and, for preregistration, a frozen experimental protocol dated before data collection.

## Submit when

Either the journal accepts a research article whose evaluation is of the method rather than of a biological result, **or** one PER/climbing table exists in a notebook. Do not submit H1–H7 as empirical biology without that table.
