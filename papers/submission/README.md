# Word manuscript preparation

`FlyLab_manuscript.docx` is the editable, author-identified preparation copy of the reviewed manuscript. It is not a submitted paper and no submission has been made.

Author details supplied by the author on 22 September 2026:

- Michél Nguyen
- University of the People
- Corresponding email: michel_ng@icloud.com
- ORCID: https://orcid.org/0000-0001-6834-4422

## Rebuild

Run `python scripts/export_paper_word.py` with Python packages `python-docx` and the `pandoc` executable available. It reads the rendered Markdown, existing figure images and captions. It does not run a simulation or change numerical results. Render and visually inspect the exported document before delivery.

The Word manuscript embeds four main result figures, numbered 1 to 4 in reading order. They correspond to archive identifiers F6, F13, F14 and F11. Other figure and table references use supplementary SF and ST identifiers and point to the existing `papers/figures/`, `papers/tables/` and `papers/SUPPLEMENT.md` artifacts. Those supporting files must accompany a later submission; the Word manuscript alone is not the full evidence package.

## Venue checks before submission

Checked against https://ijrcom.org/index.php/ijrc/about/submissions on 22 September 2026. A4, 20 mm margins, single spacing and no page numbering are used. The Research Articles section suggests 6000 to 8000 words excluding references and figures; a separate general section gives 5000 to 10000. The concise revision falls below the article-specific typical range and within the general range. Confirm its length with the editor; text was not added merely to reach 6000 words. The automated 5000 to 8000 guard is an editorial check, not certification of venue compliance.

The journal states double-blind review. Prepare an anonymised manuscript and separate author information as required by the submission portal. The current file intentionally retains the author information requested by the author.

The reference instructions mix IEEE style with Word endnote and superscript instructions. This preparation copy retains the manuscript's numbered IEEE-style citations and reference list. Confirm the editor's current template and citation requirements before converting references or declaring final format compliance. Postal address and telephone information, where required by the portal, have not been invented.

An archival release and DOI, complete numerical rerun, final author verification of references and declarations, cover letter, and final supplementary-file packaging remain outstanding. The inherited numerical record and the interrupted-rerun disclosure are preserved; Word export is not numerical reproduction.
