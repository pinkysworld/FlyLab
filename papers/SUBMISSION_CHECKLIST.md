# IJRC submission checklist

Target: **International Journal of Research in Computing** (ijrcom.org), as a **Research Article** (Introduction, Methods, Results, Discussion, Conclusion, References).
Manuscript: `papers/IJRC_FlyLab_draft.md` (rendered from `IJRC_FlyLab_draft.md.in`), with `papers/SUPPLEMENT.md`.

Regenerate everything before checking anything:

```bash
python scripts/reproduce_paper.py      # ~35-45 min, no downloads
python -m pytest -q                    # must be green
```

## Framing

IJRC is a computing research journal and the manuscript is written as one. The lead is the general problem — *a simulation built on a measured network and literature parameters cannot say which of its inputs its predictions depend on* — and the contribution is four transferable instruments: typed evidence propagation under a checked soundness invariant; input-dependence testing against an ordered ladder with an equivalence margin and multiplicity control; an ablation ladder scored against a matched reference; and specification-family robustness with variance attribution. §3 evaluates **the instruments**. The fly pharmacology is the demonstration domain, with its census, readouts and literature concordance as supporting material in §S6 and §S10.

The headline is methodological: **a dependence verdict obtained on one extract of a network can invert on another extract of the same network, and the extract's recurrence — not its size — predicts which way** (§3.2, T27).

Internal review history lives in `papers/reviews/`. It is archived critique, not journal correspondence: the corrections it forced are in the manuscript and the code, and a formal point-by-point response belongs to an actual resubmission against real referees' numbered comments.

## A. What the journal needs

| Item | State | Where |
|---|---|---|
| Title, author block, affiliations | **missing** — authors are "FlyLab contributors" | needs real names, ORCIDs, affiliations, a corresponding author and an email |
| Abstract | done — computing problem first, the methodological headline, the central limitation | draft §Abstract |
| Index terms / keywords | done | draft §Index terms |
| Numbered IEEE-style references with DOIs | done — 77 entries | draft §References |
| In-text citation of every reference | check with the grep in §D | — |
| Figures, 300 dpi, captions | done — 16 figures (F16 is the instrument validation) | `papers/figures/`, `figures/captions.md` |
| Tables | done — T0–T27 as CSV and Markdown | `papers/tables/` |
| Word count | done — inside the 5000–7000 target, generated and test-enforced | `values.paper_words_body` |
| Required section structure | done | draft §1–§5 |
| Supplementary material | done — S1–S10 | `papers/SUPPLEMENT.md` |
| Data / code / ethics / AI statements | done | draft §Statements |
| Archival DOI | **placeholder** — Zenodo deposit not created | draft §3.7, `CITATION.cff` |
| Submission format (LaTeX/Word template) | **not done** — the draft is Markdown | convert once the venue's template is confirmed |
| Cover letter | **not written** | — |
| Suggested reviewers | **not chosen** | — |

## B. What is still missing for a **results** paper

1. **One live *Drosophila* assay.** Climbing first: a complete published protocol and control statistics exist. Until a real table is imported by hand into `live_lab`, H1–H7 stay prospective software predictions.
2. **The one fitted parameter.** The variance budget names it: the engagement→gain transformation, not a potency value, which the model cannot see at a saturating dose.
3. **Better power at the top of the scale ladder.** The ladder settles the verdict, but its 25k/50k rungs run at n = 20 (resolution 0.048), where a rejection is the smallest the test can express. More permutations there, or a cheaper null, would turn confirmation into evidence.
4. **Expression coverage above its current value**, with adult motor-neuron receptor expression sourced or the bound declared permanent.
5. **An externally archived release**, without which "prospective" may not become "pre-registered".

## C. Known weaknesses a reviewer will find first

1. **The central RQ2 result reverses between the two cuts.** Answer: that *is* the result, it is measured rather than argued, the instrument is validated against planted ground truth (F16, T22) so the reversal is not an artefact of power, the structural difference is reported in advance (T19), and the scaling study over seven extracts (T27) identifies recurrence rather than size as what predicts it.
2. **The top rungs of the scale ladder run at n = 20.** Answer: stated in §3.2 and §S2.4 with the resolution (0.048 against α = 0.05) — those rungs confirm a verdict already settled at better-powered rungs, and the paper says they do not establish it.
3. **The composition-versus-full correlation is nearly an algebraic identity.** Answer: stated as such, with a matched reference distribution (T24) and a normalisation sweep (T25); the generalisation drawn from it is withdrawn.
4. **The headline circuit result's direction depends on the gain rule.** Answer: a first-class result (§3.4), quantified over the prespecified specification family, and nominated by the variance budget as the measurement to make next.
5. **The gain rules and five transmitter signs are asserted coefficients.** Answer: published (T2, T21) and named as such; the paper says which of its own claims that costs it.
6. **No independent out-of-sample validation.** Answer: stated in §2.5 and §S10, with shared-source rank comparisons flagged by the pipeline. The ground-truth experiment validates the instrument, not the domain model, and the paper says so.
7. **Two of the instruments were broken.** Answer: found by us, published in §3.4, §S1.1 and §S2, and both fixed with the corrected numbers in the text. Neither failure was detectable by a determinism check, which is itself reported as a finding.
8. **`--fast` is not the paper.** Answer: one test refuses a committed record from a fast run, another asserts the record's statistical knobs equal the shipped defaults.

## D. Pre-submission mechanical checks

```bash
# every figure and table referenced by the draft exists
grep -o 'F[0-9]\+\b' papers/IJRC_FlyLab_draft.md | sort -u
ls papers/figures/*.png
grep -o '\bT[0-9]\+[a-z]\?\b' papers/IJRC_FlyLab_draft.md | sort -u
ls papers/tables/*.csv

# no unresolved template placeholders, no missing keys, in either document
grep -c '{{\|\[\[MISSING' papers/IJRC_FlyLab_draft.md papers/SUPPLEMENT.md   # must be 0
python -c "import json;print(json.load(open('papers/results.json'))['values']['paper_missing_keys']['value'])"

# the body word count the venue will check
python -c "import json;print(json.load(open('papers/results.json'))['values']['paper_words_body'])"

# the retracted framings must not have come back
grep -in 'pre-registered\|preregistered\|do not touch\|no public tool' papers/IJRC_FlyLab_draft.md
grep -in 'exactly the chloride-channel\|reproduces the effect\|preserved exactly' papers/IJRC_FlyLab_draft.md papers/SUPPLEMENT.md

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

# the record came from a full (not --fast, not partial) run at the shipped effort
python -c "import json;d=json.load(open('papers/results.json'));print('fast',d['fast'],'partial',d.get('partial'))"
python -m pytest -q tests/test_reproduce.py
```

## E. Final pass before sending

- [ ] Real authors, affiliations, ORCIDs, corresponding author.
- [ ] Zenodo deposit created; DOI substituted for the placeholder in the draft and in `CITATION.cff`; only then may the predictions be called pre-registered.
- [ ] Tag the release; confirm `git_sha` in `papers/results.json` matches the tag.
- [ ] Re-run `scripts/reproduce_paper.py` on a clean checkout and confirm the values are unchanged.
- [ ] `python -m pytest -q` green; `pytest -m slow` green.
- [ ] GitHub Pages bench live, and the link works from a private window.
- [ ] Body word count inside the venue's limit (generated value, not an estimate).
- [ ] Confirm that no sentence asserts a live-animal result, that no dependence claim is stated without the cut it was measured on, and that no non-rejection is described as a reproduction.
