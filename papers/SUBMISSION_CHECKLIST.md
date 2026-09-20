# IJRC submission checklist

Target: **International Journal of Research in Computing** (ijrcom.org), as a **methods / software article**.
Manuscript: `papers/IJRC_FlyLab_draft.md` (rendered from `IJRC_FlyLab_draft.md.in`).

Regenerate everything before checking anything:

```bash
python scripts/reproduce_paper.py      # ~3.5 min, no downloads
python -m pytest -q                    # must be green
```

## A. What the journal needs

| Item | State | Where |
|---|---|---|
| Title, author block, affiliations | **missing** — authors are "FlyLab contributors" | needs real names, ORCIDs, affiliations, a corresponding author and an email |
| Abstract | done | draft §Abstract |
| Index terms / keywords | done (8) | draft §Index terms |
| Numbered IEEE-style references with DOIs | done — 67 entries, [13]–[67] carry DOI or PMID | draft §References |
| In-text citation of every reference | done | check with the grep in §D below |
| Figures, 300 dpi, captions | done — 10 figures, PNG (+SVG for F1) | `papers/figures/`, `figures/captions.md` |
| Tables | done — 10 tables as CSV and Markdown | `papers/tables/` |
| Word count | **over** — body ≈ 8 000 words against a 5 000–7 000 target | see §C |
| Data availability statement | done | draft §Statements |
| Code availability statement | done | draft §Statements |
| Ethics / animal statement | done — no animal experiments, `live_lab` null everywhere | draft §Statements |
| AI-assistance disclosure | done | draft §Statements |
| Competing interests, funding | done (none / none) | draft §Statements |
| Archival DOI | **placeholder** — Zenodo deposit not created | draft §8, `CITATION.cff` |
| Licence for reuse of figures | MIT code, CC-BY 4.0 connectome data credited | draft §8 |
| Submission format (LaTeX/Word template) | **not done** — the draft is Markdown | convert once the venue's template is confirmed |
| Cover letter | **not written** | — |
| Suggested reviewers | **not chosen** | — |

## B. What is still missing for a **results** paper

This manuscript is deliberately a methods article. To publish the same work as a results paper, three things are needed, in order:

1. **One live *Drosophila* assay.** Climbing is the better first target: a complete published protocol and control statistics exist (RING; Martelli 2020, imidacloprid). PER is worse served — the specific control rate FlyLab needs (response to 100 mM sucrose with its variance) could not be sourced and is recorded as `null` in `data/literature/behavioral_assays.yaml`. Until a real table is imported by hand into `live_lab`, H2–H7 are predictions and must be labelled as such.
2. **At most one fitted parameter.** The sensitivity analysis says which one: a gain-rule coefficient, not an EC50 (at a saturating dose the model is *exactly* insensitive to the EC50). Fitting more than one would break the project's stated parameter policy.
3. **Subunit-resolved nicotinic receptors.** `insect_nAChR` is one key where α6 (spinosad) and β1-dependent (neonicotinoid, sulfoximine) receptors occupy different, partly overlapping cell sets. Any quantitative selectivity claim about the nicotinic set is capped by this until the library schema splits them.

Nice to have, not blocking: expression coverage above 0.21; more null-model shuffles so the bitter-veto arm can be tested properly; a second connectome (FlyWire) as an independent netlist.

## C. Known weaknesses a reviewer will find first

Have an answer ready for each. All are already stated in the manuscript.

1. **The central circuit result does not beat a degree-preserving null for imidacloprid** (max |z| ≈ 1.2). Answer: that is reported as the paper's main negative result, with the mechanistic reason (E/I-balance effect), and fipronil does beat those nulls.
2. **The bitter veto ratio is not distinguishable from its shuffles** (max |z| ≈ 0.7 at 10 shuffles). Answer: stated in §5.5 and in the limitations; H6 is offered as a pharmacology prediction, not a wiring prediction.
3. **Every effect size is model-internal.** Answer: §6 says so explicitly and caps d at 1.0 before any power calculation.
4. **The gain rules are invented coefficients.** Answer: yes, and §5.7 shows they are the dominant parameter; they are in one file with a rationale per rule, and they are the intended target of the one allowed fit.
5. **The LIF background Poisson drive is a modelling choice** that sets the operating point. Answer: §4.5, stated plainly; no result's direction depends on it.
6. **The word "veto" is borrowed from a published result.** Answer: `reduced_taste_v0` is retained explicitly as the directional control and the published work is cited as the source of the direction.
7. **Body word count exceeds the target.** Answer: trim §4 and §5 further, or move the mechanism-rule rationale and the null-model definitions into a supplementary section, keeping T2 and F6 in the main text.
8. **Two rank inversions.** Answer: recorded, explained as receptor-EC50 versus whole-animal potency, and deliberately not fixed.

## D. Pre-submission mechanical checks

```bash
# every figure and table referenced by the draft exists
grep -o 'F[0-9]\+\b' papers/IJRC_FlyLab_draft.md | sort -u
ls papers/figures/*.png
grep -o '\bT[0-9]\+\b' papers/IJRC_FlyLab_draft.md | sort -u
ls papers/tables/*.csv

# no unresolved template placeholders, no missing keys
grep -c '{{\|\[\[MISSING' papers/IJRC_FlyLab_draft.md          # must be 0
python -c "import json;print(json.load(open('papers/results.json'))['values']['paper_missing_keys']['value'])"

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
- [ ] Zenodo deposit created; DOI substituted for the placeholder in the draft and in `CITATION.cff`.
- [ ] Tag the release; confirm `git_sha` in `papers/results.json` matches the tag.
- [ ] Re-run `scripts/reproduce_paper.py` on a clean checkout and confirm the values are unchanged.
- [ ] `python -m pytest -q` green; `pytest -m slow` green.
- [ ] GitHub Pages bench live, and the link works from a private window.
- [ ] Word count inside the venue's limit, or a stated reason why not.
- [ ] Confirm that no sentence in the manuscript asserts a live-animal result.
