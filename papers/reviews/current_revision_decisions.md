# Review decisions for the current manuscript

Base: main `38cad67`. This revision first compared the newer Fable/Opus changes with the earlier `round2_referee_report.md` and the Grok review supplied by the author. Earlier counts and defects were not assumed to describe the current manuscript. Review is internal and AI-assisted, not journal peer review.

## Review rounds

1. Independent statistical review of the updated manuscript, analysis code and saved results, plus a separate reference audit. Reports: `current_round1_statistical.md` and `current_reference_audit.md`.
2. Independent rereview of the revised templates against implementation details. Report: `current_round2_statistical.md`. This found residual component-test versus cell-level multiplicity claims, asymmetric interpretation of non-rejection, stale confirmatory labels, a contradictory C4 warning, and overstated invariant-test coverage.
3. Integrating review of manuscript, supplement, generator, project claim guidance and submission checklist. Accepted changes were checked against their implementation; remaining limitations are recorded below. This round also repaired prose and checked that the earlier reviews were applied to the current evidence rather than their superseded snapshots.

The statistical and implementation reviews used GPT-5.6 Sol with high reasoning; the focused reference audit used GPT-5.6 Luna with high reasoning. The integrating reviewer retained responsibility for accepting changes. Some earlier agent turns stopped at a usage limit and were resumed; an interrupted review was not counted as completed evidence.

## Earlier referee report

| Issue | Current disposition |
|---|---|
| Row normalisation and sparse-cut generalisation | Newer upstream work already documents normalisation, compares the denser taste-motor cut, and supplies scale/normalisation analyses. Retained. Further narrowed the claim: recurrence is not isolated from density, composition and extraction recipe. |
| Biased transmitter-label null | Upstream added a transmitter-outweight-matched null and documented imbalance of plain label permutation. Retained the distinction; no claim that plain label permutation isolates sign alone. |
| Vacuous specification analysis | Upstream fixed gain propagation and increased permutation effort. Retained. Further clarified that adjusted non-rejection does not strengthen C3, whereas C4 is a thresholded rejection predicate. |
| Non-rejection described as reproduction | Removed remaining affirmative equivalence/reproduction interpretations. The legacy within-tolerance identifier is a descriptive point-gap check, not a confidence-bounded equivalence test. |
| Direction-blind ablation and absolute rank correlation | Upstream corrected both. Retained the direction-aware comparator and signed rank correlation. |
| Composition/full correlation and reference | Retained the matched reference. The observed correlation exceeds its reference 95th percentile; joint permutation of both vectors is a mathematical invariance, not proof of no compound information. |
| Uncertainty, interaction residual and drive ranking | Retained upstream corrections. The largest negative estimate is a diagnostic, not a calibrated error bound or general noise floor. |
| Effect floor, nonmonotone ladder and library/substrate composition | Retained explicit reporting conventions, nonmonotonicity caveat, transmitter census and restricted interpretation. No claim that a null effect establishes biological inactivity. |
| Evidence provenance and pharmacological interpretation | Retained the newer evidence-distance rules and native Drosophila binding row. Narrowed invariant claims to tested rows and enumerated paths. |

## Grok review

| Comment | Decision |
|---|---|
| M1: false chloride-blocker class | Already corrected upstream. Keep mechanisms distinct; no chemical-class discovery claim. |
| M2: submitted header and missing metadata | Working draft, not submitted. Author identities, affiliations, ORCIDs and corresponding contact require human confirmation. No identity or DOI fabricated. |
| M3: exploratory landscape census | Earlier 20/84 count is obsolete. Retain current counts from generated results, described as exploratory classes based on BH-adjusted component tests. This is neither demonstrated cell-level FDR control nor validated arbitrary-dependence control. |
| M4: type-system inflation and explicit rules | Describe validated tagged records, dispatch/refusal rules and finite regression coverage. Add admissible/rejected examples in main Methods. This is evidence gating and labelling, not a new binding model. |
| M5: vertebrate comparison | Table I identifies a vertebrate receptor panel; no second nervous-system simulation or organism safety claim. |
| M6: incomparable Hz | Distinguish rate-model units, LIF spikes/s and input spike rates in prose, captions and plotted labels. Machine-readable legacy field names need not imply physical calibration. |
| Minor: LB3/LB1 modality | Explicit literature-based working assignment, not MaleCNS ground truth. |
| Minor: expression coverage | Abstract includes coverage and unresolved MN9. Coverage is of cell/receptor mapping, not the fraction of experimentally measured cells. |
| Minor: VOI | Caption and text identify variance under assumed model input ranges, not value of an experiment on a fly. |
| Minor: 3Rs | Removed unsupported index/keyword claim. No demonstrated replacement of animal experiments. |
| Minor: stale commit header | Removed prose commit pin; actual numerical provenance remains in generated records. An archival tagged release remains pending. |
| Minor: references | Corrected Virtual Rat citation, explicitly identified PharmVR as a preprint concept, replaced an unverifiable technical note, and labelled newsletter commentary as informal. See source audit. |
| Minor: abstract ablation result | Included the high composition-only/full-model agreement under shared assumptions. |

## New findings accepted

- Extract sensitivity is observed; recurrence is a candidate associated feature, not an isolated cause or validated predictor.
- Matching labels on the two largest scale rungs at 20 permutations do not establish convergence. Their probability estimates sit at the resolution floor.
- Six synthetic networks per setting give imprecise detection rates. The zero-strength control also changes edges and transmitter labels. This grid does not establish general calibration or power for the taste-motor ratio endpoint.
- BH is applied to structural component tests. A cell is classified by an OR over those tests, so component-test adjustment does not automatically control cell-level FDR. Shared random draws do not establish the usual BH dependence condition.
- High permutation effort and compound identity do not make a retrospectively selected profile confirmatory.
- Missing permutation probability cannot yield a within-tolerance verdict. The helper now returns indeterminate rather than accepting the point gap or failing while formatting a missing value.
- The recovery network has 250 nodes and is not size- or density-matched to the named cut; the main text and captions now say so.
- The standard pipeline imports the scale-study cache. Large source cuts and cache source/input hashes are missing, so replay is distinguished from independent recomputation.

## Remaining limits and author actions

The manuscript remains a computational-methods working draft. No animal outcome, clinical safety inference, independent circuit validation or formal statistical equivalence has been established. Formal equivalence, controlled recurrence interventions, stronger null calibration and a fully source-pinned scale archive would require additional research, not prose edits. The existing scale-study JSON remains a historical computed record; its legacy labels are explained rather than silently rewriting numerical evidence.

Before submission, confirm author metadata and AI-use disclosure, create a genuine archival deposit, confirm venue formatting and prepare submission materials. A DOI alone would not constitute preregistration; that requires a frozen protocol before data collection. These are not fabricated to make the checklist appear complete.

Validation details are recorded in `current_verification.md`. The full numerical rerun was interrupted; the delivered update preserves the committed default-effort estimates and explicitly records editorial replay. It does not claim a completed independent numerical reproduction.
