# Current manuscript: round-2 statistical review

Scope: current `papers/IJRC_FlyLab_draft.md.in` and `papers/SUPPLEMENT.md.in`, checked against the current analysis code and `papers/results.json`. The revised templates now appropriately qualify the median-gap diagnostic, the small synthetic grid, recurrence, and cached scale study. The remaining actionable issues are below.

## 1. Major — the multiplicity target does not match the reported cell-level claim

**Evidence.** `_apply_fdr` runs Benjamini–Hochberg over the two structural **mode tests** per cell (`flylab/analysis/dependence.py:2300-2329`), stores the minimum adjusted mode probability as the cell's `q_value` (`:2337-2347`), and then reclassifies the cell. A topology-dependent cell is an OR decision: at least one structural mode rejects. The specification analysis does the same explicitly (`flylab/analysis/robustness.py:1166-1174`). Nevertheless, the generated statements call **11/84** and **54/84 cells** “FDR-controlled” (`papers/results.json:3169-3174`, `:3609-3614`; originating at `flylab/analysis/dependence.py:2550-2554`).

**Why it matters.** BH control for 168 component hypotheses does not automatically become FDR control for 84 compound-level OR decisions. For example, a true-positive cell may contribute two true rejected modes while a null cell contributes one false mode; the test-level false-discovery proportion can then be smaller than the cell-level proportion. In addition, reused shuffle draws establish dependence, not the positive-dependence condition needed for ordinary BH; the supplement acknowledges this only for the specification run (`papers/SUPPLEMENT.md.in:48`), not for the landscapes.

**Fix.** Either (a) describe these as “classes based on BH-adjusted component tests” and remove “FDR-controlled” at the cell/count level, explicitly treating BH as a sensitivity analysis under unverified dependence, or (b) form one valid cell-level global-null probability (for example, a within-cell multiplicity-adjusted minimum or joint permutation maximum statistic), then adjust those 84 cell probabilities. Do not present the current count as cell-level FDR-controlled without such a change.

## 2. Major — “C3 survives correction” reverses the meaning of multiplicity adjustment

**Evidence.** The manuscript says both topology conclusions “still hold” after BH (`papers/IJRC_FlyLab_draft.md.in:167`), and the supplement says both “survive the correction” (`papers/SUPPLEMENT.md.in:48`). C3 is retained precisely when neither structural test rejects (`flylab/analysis/robustness.py:1027-1045`, `:1330-1331`). Its record is 25/25 both adjusted and unadjusted, with only 17 median-gap passes and 8 indeterminate cases (`papers/results.json`, keys `stab_C3_retained`, `stab_C3_retained_uncorrected`, `stab_C3_equivalent`, `stab_C3_indeterminate`).

**Why it matters.** Increasing component probabilities through multiplicity adjustment can only make a failure-to-reject predicate easier to retain. Thus “survives correction” is meaningful for positive C4 rejections, but it cannot strengthen C3 or serve as robustness evidence for absence of topology dependence. The following caveat about non-rejection is correct but does not repair that asymmetric logic.

**Fix.** Report the two directions separately: “C4 remains rejected in 25/25 specifications after adjustment; C3 has no structural rejection in 25/25, unchanged from the nominal analysis, with 17 descriptive median-gap passes and 8 indeterminate.” Avoid “both survive correction.”

## 3. Major — C4's audit warning describes a statistic the implementation no longer uses

**Evidence.** The run warning says “C4 is an effect-size contrast with no threshold” (`flylab/analysis/robustness.py:1257-1262`). In fact, C4 is defined as `class == topology-dependent`, requiring at least one structural permutation test to reject (`:1047-1060`), and its final stability value is overwritten from BH-adjusted rejection decisions (`:1311-1331`). The result record correctly calls it a “permutation classification (rejection), FDR-corrected,” so the warning contradicts both code and reported evidence.

**Why it matters.** This is exactly the sort of self-description mismatch the paper presents its auditability machinery as preventing. It can also mislead readers into treating 25/25 as a threshold-free effect-size comparison.

**Fix.** Replace the warning and stale C4 doc/name language with the actual predicate: a thresholded OR over the two structural permutation tests after the stated adjustment. Add a test that the emitted warning/readout agrees with `CONCLUSIONS[C4]` and the post-BH decision path. This wording fix need not change the numeric result.

## 4. Major — retrospective reference profiles remain labelled “confirmatory (prespecified)” in artifacts and API output

**Evidence.** The supplement correctly says the imidacloprid and fipronil profiles “were refined during iterative development and are not independent preregistered confirmations” (`papers/SUPPLEMENT.md.in:82`). The implementation nevertheless declares those two compounds confirmatory (`flylab/analysis/dependence.py:79-80`, `:365-369`) and automatically sets the flag solely from compound identity and `n >= PAPER_N` (`:1280-1283`), producing `design = "confirmatory (prespecified)"` (`:1465-1467`). Those contradictory labels are committed in `papers/results.json:1660-1672` and `:2065-2077`.

**Why it matters.** Computational effort and nominal single-profile testing do not make a retrospectively selected profile confirmatory. The metadata can be consumed independently of the narrowed prose and overstate the inferential status.

**Fix.** Rename these outputs to “high-effort reference profile” (or equivalent), make prospective/confirmatory status an explicit caller-supplied field backed by a dated protocol rather than inferred from compound and sample size, and regenerate the affected record. Nominal probabilities can still be reported as such.

## 5. Moderate — the soundness-invariant novelty claim overstates the mutation and exhaustiveness evidence

**Evidence.** The manuscript claims an invariant checked “over every row and every entry point” plus “a mutation check” (`papers/IJRC_FlyLab_draft.md.in:77`); the supplement says mutating the rule table makes the test fail with a row-specific message (`papers/SUPPLEMENT.md.in:164`). The main test does cover all shipped rows through a hand-maintained set of five route groups (`tests/test_pharm_evidence.py:428-499`). However, `test_the_invariant_would_catch_a_violation` neither mutates the rule table nor reruns the pipeline invariant: it makes two direct invalid calls to `check_transformation` and asserts they raise (`:502-516`). A maintained enumeration is also not proof that all present and future entry points are covered.

**Why it matters.** The invariant is presented as a primary computing contribution. The tests support a finite regression audit of the shipped library and enumerated paths, but not the stronger mutation-testing and whole-pipeline exhaustiveness claims.

**Fix.** Narrow the prose to “all shipped rows across the enumerated comparison, curve, uncertainty, assay, and dashboard paths, plus deliberate invalid-transformation tests.” Alternatively, add a real mutation test that alters the transformation rule/row, runs the end-to-end audit, and demonstrates the claimed diagnostic; add an architectural enforcement point if exhaustive path coverage is intended.
