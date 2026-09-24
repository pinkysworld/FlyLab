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

The result is methodological: two structurally different extracts yield different descriptive dependence-label counts. Density, composition, extraction rule and recurrent paths vary together; the comparison does not identify recurrence or network size as the cause. The larger-scale study is an imported record whose 5k–50k input cuts are absent from this review package.

Internal review history lives in `papers/reviews/`. It is archived critique, not journal correspondence: the corrections it forced are in the manuscript and the code, and a formal point-by-point response belongs to an actual resubmission against real referees' numbered comments.

## A. What the journal needs

| Item | State | Where |
|---|---|---|
| Title, author block, affiliations | prepared on a separate editorial title page; withheld from the blind manuscript | Michél Nguyen, University of the People, ORCID and correspondence details |
| Abstract | done — computing problem first, the methodological headline, the central limitation | draft §Abstract |
| Index terms / keywords | done | draft §Index terms |
| Numbered references | 76 source-paper entries; DOI and PMID identities checked, assay-level support remains incomplete | draft §References and review audit |
| In-text citation of every reference | check with the grep in §D | — |
| Figures, 300 dpi, captions | done — 16 figures (F16 is the instrument validation) | `papers/figures/`, `figures/captions.md` |
| Tables | done — T0–T27 as CSV and Markdown | `papers/tables/` |
| Word count | done — inside the 5000–7000 target, generated and test-enforced | `values.paper_words_body` |
| Required section structure | done | draft §1–§5 |
| Supplementary material | S1–S11 prepared, with cited figures and tables in a companion archive | `papers/SUPPLEMENT.md` and separate review files |
| Funding / competing interests | author declares no funding and no competing interests | draft §Statements and editorial title page |
| Author contributions / AI / ethics / SDG mapping | supplied; SDG 9.5 is a relevance mapping, not an impact result | draft §Statements and editorial title page |
| Data and code | numerical code at a commit-specific GitHub URL; masked review archive supplied | draft §Statements and review package |
| DOI | IJRC assigns the article DOI after publication; a separate code archive DOI is optional and currently absent | draft §Statements |
| Submission format | blinded Word manuscript, separate title page and formatted supplement prepared | local review package |
| Cover letter | **not written** | — |
| Suggested reviewers | **not chosen** | — |

## B. Work not claimed as completed

1. **Biological validation.** No live *Drosophila* assay was performed. H1–H7 remain prospective software predictions, not empirical effects.
2. **Engagement-to-gain calibration.** The transformation and five transmitter-sign magnitudes remain assumed model inputs; the reported specification family does not vary all of them.
3. **Independent large-scale reproduction.** The 5k–50k cuts and their inputs are absent, so ST27 is an imported exploratory record. Its low permutation counts also limit resolution.
4. **Assay-level library provenance.** Bibliographic identity does not establish the potency value used in every numeric receptor row.
5. **Preregistration.** Neither a public repository nor a later code DOI can retrospectively preregister H1–H7; a dated, frozen experimental protocol would be needed before a future experiment.

## C. Known weaknesses a reviewer will find first

1. **The two cuts yield different labels.** The manuscript reports descriptive differences and a limited planted-cycle positive control; it does not attribute the contrast causally to recurrence or claim general power for the paper endpoint.
2. **The top rungs of the saved scale study run at n = 20.** Their minimum attainable probability is 0.048, and the input cuts are unavailable here. The study is exploratory, not independent confirmation.
3. **The composition-versus-full correlation is nearly an algebraic identity.** Answer: stated as such, with a matched reference distribution (T24) and a normalisation sweep (T25); the generalisation drawn from it is withdrawn.
4. **One model-output direction depends on the gain rule.** The tested family is reported in §3.4 and §S1.1; fixed sign magnitudes and transmitter assignment remain outside it.
5. **The gain rules and five transmitter signs are asserted coefficients.** Answer: published (T2, T21) and named as such; the paper says which of its own claims that costs it.
6. **No independent out-of-sample validation.** This is stated in §2.5 and §S10, with detectable shared-source rank comparisons flagged by the pipeline. The planted-cycle test is a limited synthetic positive control; it does not validate the domain model or estimate general detection power.
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

- [x] Author, affiliation, ORCID and correspondence details are on the separate editorial title page.
- [x] Funding, competing interests, AI use, contributions, ethics, data/code and SDG statements are drafted; check editorial metadata entry at submission.
- [x] The article DOI is left to IJRC; the code uses a commit-specific URL and hash. No code DOI or preregistration is claimed.
- [ ] Confirm that the blinded manuscript, formatted supplement and masked code archive are accessible to reviewers without exposing the title page.
- [ ] Re-run `scripts/reproduce_paper.py` on a clean checkout and confirm the values are unchanged.
- [ ] `python -m pytest -q` green; `pytest -m slow` green.
- [ ] GitHub Pages bench live, and the link works from a private window.
- [ ] Body word count inside the venue's limit (generated value, not an estimate).
- [ ] Confirm that no sentence asserts a live-animal result, that no dependence claim is stated without the cut it was measured on, and that no non-rejection is described as a reproduction.
