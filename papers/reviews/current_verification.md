# Verification of the internal-review revision

Date: 2026-09-21. Base manuscript: main `38cad67`. Reviewed source: `3bd92cc2bb759774884887c7c0fd6741fe260ef9`. This note describes completed checks, not journal peer review or submission readiness.

## Completed

- Full default test suite: **1221 passed, 2 skipped, 12 deselected**, with eight dependency/fork deprecation warnings, in 364.84 seconds. The repository excludes slow tests by default. Command: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -m pytest -q`.
- Focused paper/replay checks: **21 passed, 1 deselected**. These cover rendering contracts, both cached-scale rendering branches, numerical/provenance preservation and repeatable editorial replay.
- Commit-attribution policy passed in an isolated checkout of this branch. The working repository retains pre-rewrite local recovery refs, so its all-refs scan includes superseded history; the proposed branch does not.
- Manuscript and supplement rendered from their templates. Main body: **6872 words**. All quoted keys resolve. Bibliography: 76 consecutive entries, with no missing citation numbers. No em dashes in the rendered documents.
- `git diff --check` passed. Updated receptor-panel and rate/LIF plots were visually inspected for labels and clipping.
- The architecture, receptor scorecard, curves and genotype plots were regenerated at ordinary effort. Their value payloads exactly match the saved numerical record. The rate/LIF plot was also regenerated before the full attempt was interrupted.
- The stored record is not `--fast`. Its declared effort remains 1000 permutations for the named profiles and named landscape, 300 for the taste-motor landscape, 100 per stability specification, and Sobol base size 1024. These are inherited recorded settings, not a claim that this review independently reran them all.

## Numerical preservation and reproduction limit

The full default-effort rerun was interrupted during the dependence step. It reached the reference-profile, transmitter-balance and synthetic-recovery outputs, but did not complete the complete landscape and paper pipeline. No completed full numerical reproduction is claimed.

The delivered `results.json` preserves the numerical values from the previously committed default-effort record and its original recorded source provenance (`f95fd88a2165e66d5203db6c6fa8c6a020c5a8e0`). That historical source identifier is retained as recorded, not replaced with a new hash to imply a rerun. Comparing all value payloads with main shows changes only to design-status/interpretation fields and document word counts. All table CSV files are byte-identical to main.

`scripts/paper_review_metadata.py` performs an explicit editorial replay. It corrects captions, notes and retrospective/confirmatory labels while checking that other value payloads and table CSV bytes are unchanged. The same caption/note corrections are applied by future ordinary pipeline runs. `results.json.review_metadata` records the distinction and the reviewed source SHA. A paper-only render does not itself establish analysis reproducibility.

The scale table remains a replay of `scale_study.json`. The uncommitted large source cuts and missing cache input/source hashes remain limitations. Formal equivalence, cell-level FDR control, endpoint-specific power and an isolated causal recurrence effect have not been established by these edits.

The submission checklist still requires verified author metadata, author review of sources and AI disclosure, a genuine archival deposit, venue formatting and submission materials. No author identity, DOI, animal result or completed rerun was fabricated.
