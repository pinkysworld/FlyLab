# IJRC submission checklist

Target: **International Journal of Research in Computing** (ijrcom.org), as a **Research Article** (Introduction, Methods, Results, Discussion, Conclusion, References).
Manuscript: `papers/IJRC_FlyLab_draft.md` (rendered from `IJRC_FlyLab_draft.md.in`), with `papers/SUPPLEMENT.md`.

Regenerate everything before checking anything:

```bash
python scripts/reproduce_paper.py      # ~8-10 min, no downloads
python -m pytest -q                    # must be green
```

## R. Response to the Major Revision review

Ten comments, what changed, and where. Every number cited below is generated into `papers/results.json` by the step named in its row.

| # | Reviewer comment | What changed | Where |
|---|---|---|---|
| 1 | A Hill curve on an EC50 is not "receptor occupancy". | The evidence is **typed**. `param_type` (what the source measured) and `relation` (its distance from this compound/receptor/species) jointly decide the transformation: `binding_occupancy` for a measured Kd/Ki, `functional_engagement` for EC50/IC50/Kb, `not_modelled` otherwise. The manuscript says **engagement** throughout and reserves *occupancy* for the single row that earns it (imidacloprid at `insect_nAChR_beta1`, Kd 8.3e-11 M). | §2.1, §3.1; F11, T10; supplement S6; step `evidence` |
| 2 | Missing evidence must not become a small quantitative response. | Unsupported rows return **N/A** and are excluded from every numeric result and from the selectivity block. This is enforced by the type system — `EvidenceTypeError` — not by convention, and the library validator rejects a placeholder that carries a value. Diazepam at insect RDL is now N/A, not 1e-4. 56 of 101 rows are not modelled. | §2.1, §3.1; supplement S6 |
| 3 | The gain rule may be doing the work attributed to the integration. | A **conclusion-stability matrix** over 25 prespecified specifications is now a first-class result, and the paper states plainly that the suppression conclusion is specification-dependent: retained 15/25 and reversed by all ten monotone rules, where the same compound at the same engagement *excites* the network. The topology conclusions, RDL disinhibition and the map bitter-veto direction are 25/25. | §2.4, §3.4, §4; F14, T14; supplement S1.1; step `stability` |
| 4 | "Validation" is used for things that are not validation. | Evaluation is split into **implementation consistency**, **literature concordance** and **independent out-of-sample validation**, and the paper states that it has none of the third at circuit level. The pipeline computes, from DOIs/PMIDs, which rank comparisons share a source with the library and therefore cannot be out-of-sample. The honeybee mixture test is reframed as a limited qualitative concordance check (non-*Drosophila* reference organism, two substituted compounds, three-valued endpoint). | §2.5, §3.6; T4; step `validation` |
| 5 | Null-model strength: lead with the empirical p and report a count sweep. | The dependence analysis reports the empirical two-sided permutation *p* first with its resolution `1/(n+1)`; *z* is explicitly secondary. n = 1000 on the `named` cut, 300 on `taste_motor`, 100 per landscape cell. The permutation-count sweep is computed inside the same draws and reported: verdicts settle early, a quotable mid-range *p* does not. | §2.3, §3.2; F6, T6, T11; supplement S2; step `dependence` |
| 6 | The selectivity thresholds are conventions. | The amplify/buffer split is recomputed on the full 3×3 grid (25/40/50 % circuit change × 10/20/30 % vertebrate engagement). Overall `stable = False`; nicotinic buffering holds 9/9 (gap −0.362 to −1.224), Nav/AChE amplification 8/9 and flips at the strictest corner (gap −0.075). The paper quotes ranges, never the point estimate. | §2.4, §3.4; T15; step `stability` |
| 7 | Novelty positioning overstates the gap. | "Two bodies of software sit next to each other and do not touch" and every "no tool has" claim are **deleted** (a test fails if they return). FlyBrainLab [68] and receptor-informed whole-brain models [69], [70] are cited; the paper discusses that Mindlin *et al.* [69] likewise found the effect dominated by overall receptor presence rather than placement, and treats that convergence as support. Novelty is stated narrowly with a comparison table; exposure, mixtures and genotypes are demoted to supporting capabilities in the supplement. | §1.1, Table I, §4; supplement S4 |
| 8 | A mutable GitHub repository is not a registry. | "Pre-registered" is replaced by **prospective** everywhere, with an explicit statement of what would change that (an externally archived, timestamped release). A test fails the manuscript if the word returns. | §4, supplement S5; `tests/test_reproduce.py` |
| 9 | The capped-d sample sizes are not a power analysis. | Removed from the main text and demoted to the supplement as an **illustrative planning minimum**, with the reasons it must not size a study (no biological variance; vial-based protocols violate independence). | supplement S5 |
| 10 | A Research Article needs a Discussion. | §4 Discussion added: what was learned computationally; model-generated observations separated from encoded assumptions; comparison to adjacent work; intended use; **what the software cannot support**; and the four model-specific caveats (glutamate signed inhibitory, receptor effects applied across transmitter-defined synapses, 20.5 % expression coverage with MN9 unresolved, the rate-vs-LIF discrepancy as direction-only agreement) discussed rather than listed. | §4 |
| + | Word count (8350 body words against a 5000–7000 target). | Mechanism rationale, null-model definitions, runtime calibration, supporting capabilities and the prospective predictions moved to `papers/SUPPLEMENT.md`. The body count is **generated** (`values.paper_words_body`) and a test fails outside the target. | §all; supplement |
| + | Organise around research questions. | The paper is organised around RQ1 (evidence without losing provenance or semantics), RQ2 (which effects need MaleCNS topology), RQ3 (which conclusions survive, which assumptions dominate), RQ4 (identical reproduction natively and in-browser). | §1, §3 |
| + | Separate model-generated observations from encoded assumptions. | Done in prose in §4 **and generated**: `flylab.analysis.claims` audits the dependency chain behind a readout (Table T18) and labels each link OBSERVED / LITERATURE-DERIVED / MODEL-ASSUMPTION / COMPUTED, with a companion fact / model-inference / unknown split. Exactly one link of nine is a measurement of the system being simulated. | §4; T18; step `claims` |
| + | AI declaration. | Replaced with a precise statement separating AI-assisted implementation and drafting from the author's scientific responsibility: the author defined the questions and methods, verified the source literature, executed and reviewed the analyses, and accepts responsibility for the scientific content. | Statements |

## A. What the journal needs

| Item | State | Where |
|---|---|---|
| Title, author block, affiliations | **missing** — authors are "FlyLab contributors" | needs real names, ORCIDs, affiliations, a corresponding author and an email |
| Abstract | done — numbers stripped, architecture + one topology result + one reproducibility result + the central limitation | draft §Abstract |
| Index terms / keywords | done (10) | draft §Index terms |
| Numbered IEEE-style references with DOIs | done — 76 entries | draft §References |
| In-text citation of every reference | done | check with the grep in §D |
| Figures, 300 dpi, captions | done — 15 figures, PNG (+SVG for F1) | `papers/figures/`, `figures/captions.md` |
| Tables | done — 21 tables as CSV and Markdown | `papers/tables/` |
| Word count | **done** — inside the 5000–7000 target, generated and test-enforced | `values.paper_words_body` |
| Required section structure (Intro/Methods/Results/Discussion/Conclusion/References) | done | draft §1–§5 |
| Supplementary material | done | `papers/SUPPLEMENT.md` |
| Data availability statement | done | draft §Statements |
| Code availability statement | done | draft §Statements |
| Ethics / animal statement | done — no animal experiments, `live_lab` null everywhere | draft §Statements |
| AI-assistance disclosure | done — precise, responsibility-assigning | draft §Statements |
| Competing interests, funding | done (none / none) | draft §Statements |
| Archival DOI | **placeholder** — Zenodo deposit not created | draft §3.7, `CITATION.cff` |
| Submission format (LaTeX/Word template) | **not done** — the draft is Markdown | convert once the venue's template is confirmed |
| Cover letter | **not written** | — |
| Suggested reviewers | **not chosen** | — |

## B. What is still missing for a **results** paper

1. **One live *Drosophila* assay.** Climbing first: a complete published protocol and control statistics exist. Until a real table is imported by hand into `live_lab`, H1–H7 stay prospective software predictions.
2. **The one fitted parameter.** The variance budget names it: the engagement→gain transformation (S1 0.410, VOI 1.29 of 3.149 Hz²), not a potency value, which the model cannot see at a saturating dose.
3. **Expression coverage above 0.205**, with adult motor-neuron receptor expression sourced or the bound declared permanent.
4. **An externally archived release**, without which "prospective" may not become "pre-registered".

## C. Known weaknesses a reviewer will find first

Have an answer ready for each. All are already stated in the manuscript.

1. **Most predictions do not need the connectome.** Answer: that is the finding, it is measured per prediction rather than argued, and the topology-dependent minority is chemically coherent (the chloride-channel blockers). An independent parallel exists in human whole-brain receptor-map models [69].
2. **The suppression result depends on the gain rule.** Answer: stated as a first-class result (§3.4), quantified over 25 specifications, and nominated by the variance budget as the experiment to do next. The conclusions we do quote are 25/25 specification-independent.
3. **The gain rules are asserted coefficients.** Answer: yes — and the paper says which of its own claims that costs it, rather than defending them.
4. **Every effect size is model-internal.** Answer: §2.5 and §4; the planning minima are demoted to the supplement and labelled illustrative.
5. **No independent validation.** Answer: stated explicitly in §2.5 and §3.6, with the shared-source rank comparisons flagged by the pipeline.
6. **The LIF background Poisson drive is a modelling choice.** Answer: supplement S3, stated plainly; the engines agree on direction only and the paper says so.
7. **Two rank inversions.** Answer: recorded, explained as receptor-parameter versus whole-animal potency, deliberately not fixed.
8. **`--fast` is not the paper.** Answer: a test refuses a committed `results.json` that came from a fast run.

## D. Pre-submission mechanical checks

```bash
# every figure and table referenced by the draft exists
grep -o 'F[0-9]\+\b' papers/IJRC_FlyLab_draft.md | sort -u
ls papers/figures/*.png
grep -o '\bT[0-9]\+\b' papers/IJRC_FlyLab_draft.md | sort -u
ls papers/tables/*.csv

# no unresolved template placeholders, no missing keys, in either document
grep -c '{{\|\[\[MISSING' papers/IJRC_FlyLab_draft.md papers/SUPPLEMENT.md   # must be 0
python -c "import json;print(json.load(open('papers/results.json'))['values']['paper_missing_keys']['value'])"

# the body word count the venue will check
python -c "import json;print(json.load(open('papers/results.json'))['values']['paper_words_body'])"

# the retracted framings must not have come back
grep -in 'pre-registered\|preregistered\|do not touch\|no public tool' papers/IJRC_FlyLab_draft.md

# every reference is cited at least once
python - <<'EOF'
import re
t = open('papers/IJRC_FlyLab_draft.md').read()
body, refs = t.split('## References')
declared = {int(n) for n in re.findall(r'^\[(\d+)\]', refs, re.M)}
cited = set()
for m in re.findall(r'\[(\d+)\](?:\s*[,–-]\s*\[(\d+)\])?', body):
    a = int(m[0]); cited.add(a)
    if m[1]:
        cited |= set(range(a, int(m[1]) + 1))
print('declared but never cited:', sorted(declared - cited))
print('cited but not declared :', sorted(cited - declared))
EOF

# the record came from a full (not --fast, not partial) run
python -c "import json;d=json.load(open('papers/results.json'));print('fast',d['fast'],'partial',d.get('partial'))"
```

## E. Final pass before sending

- [ ] Real authors, affiliations, ORCIDs, corresponding author.
- [ ] Zenodo deposit created; DOI substituted for the placeholder in the draft and in `CITATION.cff`; only then may the predictions be called pre-registered.
- [ ] Tag the release; confirm `git_sha` in `papers/results.json` matches the tag.
- [ ] Re-run `scripts/reproduce_paper.py` on a clean checkout and confirm the values are unchanged.
- [ ] `python -m pytest -q` green; `pytest -m slow` green.
- [ ] GitHub Pages bench live, and the link works from a private window.
- [ ] Body word count inside the venue's limit (generated value, not an estimate).
- [ ] Confirm that no sentence in the manuscript asserts a live-animal result, or that a prediction needs the connectome without naming its permutation *p*.
