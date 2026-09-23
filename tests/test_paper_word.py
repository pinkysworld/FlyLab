"""Check the committed Word copy without requiring the authoring dependencies."""
import hashlib
import json
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
WORD = ROOT / "papers/submission/FlyLab_manuscript.docx"
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def test_word_export_matches_inputs_and_payload():
    record = json.loads(WORD.with_suffix(".json").read_text())
    for path, expected in record["inputs_sha256"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
    assert hashlib.sha256(WORD.read_bytes()).hexdigest() == record["output_sha256"]
    assert record["main_figure_archive_ids"] == ["F6", "F13", "F14", "F11"]


def test_word_author_figures_and_references():
    with ZipFile(WORD) as archive:
        document = ET.fromstring(archive.read("word/document.xml"))
    text = "\n".join("".join(p.itertext()) for p in document.iter(W + "t"))
    for item in ("Michél Nguyen", "University of the People", "michel_ng@icloud.com", "0000-0001-6834-4422"):
        assert item in text
    assert len(list(document.iter(W + "drawing"))) == 4
    for number in range(1, 77):
        assert f"[{number}]" in text
    assert "[[MISSING:" not in text
    assert "Author name, affiliation" not in text
    assert "—" not in text
