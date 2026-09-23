"""Refresh review wording on saved paper artifacts without rerunning analyses.

The numerical-source provenance is retained. This is explicitly a metadata replay,
not numerical reproduction. All non-editorial value payloads are checked unchanged.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def review_caption(name: str, caption: str) -> str:
    """Correct interpretation of legacy captions without replacing data values."""
    text = caption.replace('—', ';')
    text = text.replace('prespecified equivalence margin', 'configured median-gap tolerance')
    text = text.replace('prespecified, admissible', 'fixed, admissible')
    text = text.replace('prespecified engagement-to-gain', 'fixed engagement-to-gain')
    text = text.replace('the part of the retention that is evidence of equivalence',
                        'the part passing a descriptive median-gap check, not statistical equivalence')
    text = text.replace('only as strong as the equivalence behind it',
                        'a non-rejection, not evidence of equivalence')
    text = text.replace('checked exhaustively over every row and every entry point',
                        'audited over shipped rows and enumerated pipeline paths')
    text = text.replace('The type system, not a convention, enforces this:',
                        'Runtime validation enforces this:')
    text = text.replace('the transformation the type system then permits',
                        'the transformation the evidence dispatch rule permits')
    if name.startswith(('F1_', 'F2_', 'F3_', 'F7_', 'F8_', 'F9_')):
        text = text.replace('occupancy', 'engagement').replace('Occupancy', 'Engagement')
    if name.startswith('F2_'):
        text = text.replace('Hatched bars rest on a class placeholder: the compound has no sourced EC50 at that pair and the row is inert.',
                            'Hatched entries indicate unsupported parameters and are interpreted as N/A.')
    if name.startswith('F12_'):
        text = text.replace('classified by the weakest graph model the permutation test could **not** distinguish from the real cut -- never by a model that reproduces the effect.',
                            'classified from structural component-test decisions and the configured effect floor.')
        text = re.sub(r'(\d+ of \d+) cells need the wiring pattern', r'\1 cells are classified topology-dependent', text)
        text = text.replace('do not move the readout at all', 'fall below the configured effect floor')
    if name.startswith(('F16_', 'T22_')):
        text = ('Synthetic recovery checks with an imposed cholinergic cycle and a gain patch. '
                'The synthetic network is not size- or density-matched to the named connectome cut. '
                'The zero-strength case omits cycle edges and associated transmitter relabelling. '
                'Detection frequencies in this small grid have broad uncertainty; networks are reused '
                'across permutation counts. These checks do not establish general calibration or '
                'power for the real taste-motor ratio endpoint. Counts and settings are shown in the artifact.')
    if name.startswith(('F6_', 'T6_', 'F14_', 'T14_')):
        note = 'The legacy within-tolerance label is a point-gap diagnostic, not a formal equivalence test.'
        if note not in text:
            text += ' ' + note
    if name.startswith(('F12_', 'T12', 'F14_', 'T14_')):
        note = ('BH adjusts structural component tests; cell/conclusion-level FDR control is not established, '
                'and the dependence assumptions remain unverified. Adjustment cannot strengthen non-rejection.')
        if note not in text:
            text += ' ' + note
    if name.startswith(('F15_', 'T16_')):
        text = text.replace('noise floor', 'negative-estimate diagnostic')
        note = ('The largest negative estimate is an observed error diagnostic, not a calibrated uncertainty bound. '
                'VOI concerns output variance under assumed model input ranges, not the benefit of an experiment on a fly.')
        if note not in text:
            text += ' ' + note
    if name.startswith('T27_'):
        note = ('This table replays a saved study; large-cut profiles are not recomputed. Matching labels at '
                'the largest rungs do not establish convergence or isolate recurrence causally.')
        if note not in text:
            text += ' ' + note
    if name.startswith('F5_'):
        text = text.replace('at 150 Hz', 'at a nominal drive of 150 rate-model input units or LIF input spikes/s')
        note = 'Left-axis values are rate-model units; right-axis values are simulated LIF spikes/s, not measured fly firing rates.'
    elif name.startswith(('F6_', 'F8_', 'F13_', 'F15_', 'F16_', 'T6_', 'T13', 'T16_', 'T17_', 'T22_')):
        note = 'Rate-engine quantities labelled Hz in legacy fields or axes are rate-model units, not calibrated physiological firing rates.'
    else:
        note = ''
    if note and note not in text:
        text += ' ' + note
    return text


def review_note(note: str) -> str:
    text = note.replace(' Hz', ' rate-model units')
    text = text.replace('is the substrate to generalise from', 'is a denser comparison, not evidence of whole-network generality')
    text = text.replace('matters to the equivalence claim and not to the probability', 'changes both the descriptive median-gap check and the ensemble probabilities')
    text = text.replace('changes the descriptive point-gap diagnostic, not the probability', 'changes both the descriptive point-gap diagnostic and the ensemble probabilities')
    text = text.replace('Effect size dominates the permutation count: full power at planted strength', 'All sampled networks were detected at planted strength')
    text = text.replace('carries no compound-level information', 'must be interpreted against shared assumptions and its reference distribution')
    text = text.replace('reach equivalence within the prespecified margin', 'pass the descriptive median-gap tolerance')
    text = text.replace('is no longer vacuous and is no longer uniform', 'uses an attainable rejection threshold but is a non-rejection criterion')
    text = text.replace("the estimator's noise floor", 'the largest observed negative-estimate magnitude')
    text = text.replace('verdicts settled on the clean ladder', 'labels match on the largest tested clean-ladder rungs')
    qualifiers = {
        '[dependence]': 'These diagnostics do not establish statistical equivalence, cell-level FDR control or endpoint-specific power.',
        '[ablation]': 'Jointly permuting paired scores preserves correlation by construction; the observed correlation exceeds the matched reference 95th percentile.',
        '[stability]': 'Adjusted non-rejection is not strengthened evidence for absence of topology dependence.',
        '[scale]': 'This is a cached descriptive comparison, not a new scale-study rerun, convergence proof or isolated recurrence effect.',
    }
    for prefix, qualification in qualifiers.items():
        if text.startswith(prefix) and qualification not in text:
            text += ' ' + qualification
    return text


def refresh(outdir: Path = ROOT / 'papers') -> dict:
    path = outdir / 'results.json'
    data = json.loads(path.read_text())
    before = copy.deepcopy(data['values'])
    allowed = set()
    for key, row in data['values'].items():
        if key.startswith('dep_') and key.endswith('_confirmatory'):
            row['value'] = False
            row['text'] = 'False'
            allowed.add(key)
        elif key.startswith('dep_') and key.endswith('_design'):
            row['value'] = row['text'] = ('exploratory high-effort reference profile'
                if '_named_' in key else 'exploratory')
            allowed.add(key)
        elif key in ('dep_land_statement', 'dep_land_taste_statement'):
            row['value'] = row['text'] = str(row['value']).replace(
                'FDR-controlled', 'classified from BH-adjusted component tests; no cell-level FDR guarantee')
            allowed.add(key)
    for row in data['figures']:
        if row.get('caption'):
            row['caption'] = review_caption(row.get('figure', row['name']), row['caption'])
        file = outdir / row['path']
        if file.exists():
            row['bytes'] = file.stat().st_size
    data['notes'] = [review_note(note) for note in data.get('notes', [])]
    csv_before = {p: p.read_bytes() for p in (outdir / 'tables').glob('*.csv')}
    for row in data['tables']:
        file = outdir / row['path']
        if file.suffix == '.md' and file.exists():
            first, rest = file.read_text().split('\n', 1)
            prefix, caption = first.split('** ', 1)
            file.write_text(prefix + '** ' + review_caption(file.stem, caption) + '\n' + rest)
        if file.exists():
            row['bytes'] = file.stat().st_size
    assert all(p.read_bytes() == content for p, content in csv_before.items())
    for key in before:
        if key not in allowed:
            assert data['values'][key] == before[key], key
    data['review_metadata'] = {
        'mode': 'editorial replay; not a completed numerical rerun',
        'numerical_source_git_sha': data.get('provenance', {}).get('git_sha'),
        'review_source_git_sha': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'changed_design_fields': sorted(allowed),
        'numerical_value_payloads_preserved': True,
        'table_csv_bytes_preserved': True,
        'full_rerun_status': 'interrupted during dependence analysis; no completed rerun claimed',
    }
    path.write_text(json.dumps(data, indent=2) + '\n')
    return data['review_metadata']


if __name__ == '__main__':
    print(json.dumps(refresh(), indent=2))
