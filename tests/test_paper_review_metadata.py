"""Editorial replay must not modify scientific estimates or their source."""
import copy
import json
from pathlib import Path
from scripts.paper_review_metadata import refresh, review_caption


def test_review_replay_preserves_estimates_provenance_and_table_data(tmp_path):
    source = Path(__file__).resolve().parents[1] / 'papers' / 'results.json'
    data = json.loads(source.read_text())
    before = copy.deepcopy(data)
    (tmp_path / 'results.json').write_text(json.dumps(data))
    tables = tmp_path / 'tables'
    tables.mkdir()
    csv = tables / 'control.csv'
    csv.write_bytes(b'compound,effect,p\nfipronil,0.9,0.001\n')
    original_csv = csv.read_bytes()
    metadata = refresh(tmp_path)
    after = json.loads((tmp_path / 'results.json').read_text())
    assert after['provenance'] == before['provenance']
    allowed = set(metadata['changed_design_fields'])
    for key, row in before['values'].items():
        if key not in allowed:
            assert after['values'][key] == row
    assert csv.read_bytes() == original_csv
    once = copy.deepcopy(after)
    refresh(tmp_path)
    assert json.loads((tmp_path / 'results.json').read_text()) == once


def test_legacy_plot_labels_are_qualified_without_changing_estimates():
    caption = 'Effect = -6.173 Hz; configured median-gap tolerance = 0.332 Hz.'
    revised = review_caption('F6_dependence_profile', caption)
    assert caption in revised
    assert 'not a formal equivalence test' in revised
    assert 'not calibrated physiological firing rates' in revised
    assert review_caption('F6_dependence_profile', revised) == revised
