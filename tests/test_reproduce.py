"""Cheap guards on the paper reproduction pipeline.

The expensive part of ``scripts/reproduce_paper.py`` is exercised once, in the
``slow`` test at the bottom. Everything else here runs two fast steps into a
temporary directory and then checks the *committed* ``papers/results.json``
against the paper template, which is what actually keeps the manuscript honest:
a number can only appear in the draft if the pipeline produced it.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.reproduce_paper import (  # noqa: E402
    PAPER_KEYS,
    RESULTS_SCHEMA,
    STEP_NAMES,
    main,
)

PAPERS = REPO_ROOT / "papers"
TEMPLATE = PAPERS / "IJRC_FlyLab_draft.md.in"
COMMITTED_RESULTS = PAPERS / "results.json"

#: two steps that touch no circuit runtime, so the fast run stays well under a
#: second: one writes a table, the other writes a figure.
QUICK_STEPS = ("census", "architecture")

_VALUE_KEY = re.compile(r"^[A-Za-z0-9_]+$")


# --------------------------------------------------------------------------
# the schema
# --------------------------------------------------------------------------
def validate_results(doc: object) -> list[str]:
    """Validate a ``results.json`` document. Returns a list of problems.

    Deliberately hand-rolled and small: the point is that the contract is
    readable here, not that it is expressed in a schema language FlyLab would
    then have to depend on.
    """
    problems: list[str] = []

    def need(cond: bool, msg: str) -> None:
        if not cond:
            problems.append(msg)

    need(isinstance(doc, dict), "top level is not an object")
    if not isinstance(doc, dict):
        return problems

    need(doc.get("schema") == RESULTS_SCHEMA, f"schema != {RESULTS_SCHEMA!r}")
    need(isinstance(doc.get("generated_utc"), str), "generated_utc missing or not a string")
    need(isinstance(doc.get("fast"), bool), "fast missing or not a bool")
    need(isinstance(doc.get("seed"), int), "seed missing or not an int")
    need(isinstance(doc.get("runtime_s"), (int, float)), "runtime_s missing or not a number")
    for key in ("steps", "values"):
        need(isinstance(doc.get(key), dict), f"{key} missing or not an object")
    for key in ("figures", "tables", "notes", "honesty"):
        need(isinstance(doc.get(key), list), f"{key} missing or not a list")

    prov = doc.get("provenance")
    need(isinstance(prov, dict), "provenance missing or not an object")
    if isinstance(prov, dict):
        for key in (
            "flylab_version",
            "library_sha256",
            "map_id",
            "graphs",
            "seeds",
            "platform",
        ):
            need(key in prov, f"provenance.{key} missing")
        need(isinstance(prov.get("graphs"), dict), "provenance.graphs is not an object")
        sha = prov.get("library_sha256")
        need(isinstance(sha, str) and len(sha) == 64, "provenance.library_sha256 is not a sha256")

    for name, step in (doc.get("steps") or {}).items():
        need(name in STEP_NAMES, f"steps.{name} is not a known step")
        need(isinstance(step, dict), f"steps.{name} is not an object")
        if isinstance(step, dict):
            need(step.get("status") in ("ok", "FAILED"), f"steps.{name}.status is not ok/FAILED")
            need(
                isinstance(step.get("runtime_s"), (int, float)),
                f"steps.{name}.runtime_s is not a number",
            )

    for key, row in (doc.get("values") or {}).items():
        need(bool(_VALUE_KEY.match(key)), f"values key {key!r} is not an identifier")
        need(isinstance(row, dict), f"values.{key} is not an object")
        if isinstance(row, dict):
            for field in ("value", "text", "unit", "note", "step"):
                need(field in row, f"values.{key}.{field} missing")
            need(isinstance(row.get("text"), str), f"values.{key}.text is not a string")

    for entry in (doc.get("figures") or []) + (doc.get("tables") or []):
        need(isinstance(entry, dict), "figure/table entry is not an object")
        if isinstance(entry, dict):
            for field in ("name", "path", "bytes"):
                need(field in entry, f"figure/table entry missing {field}")
    return problems


def template_keys() -> set[str]:
    """Every ``{{key}}`` the paper template substitutes."""
    if not TEMPLATE.exists():
        return set()
    text = re.sub(r"\A<!--.*?-->\s*", "", TEMPLATE.read_text(), flags=re.S)
    return set(re.findall(r"\{\{([A-Za-z0-9_]+)\}\}", text))


# --------------------------------------------------------------------------
# fast run
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def fast_run(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("repro")
    argv = ["--fast", "--outdir", str(out), "--quiet"]
    for step in QUICK_STEPS:
        argv += ["--only", step]
    rc = main(argv)
    assert rc == 0, "reproduce_paper reported a failed step"
    return out


def test_fast_only_writes_its_outputs(fast_run: Path) -> None:
    assert (fast_run / "results.json").is_file()
    assert (fast_run / "tables" / "T0_labellar_types.csv").is_file()
    assert (fast_run / "tables" / "T0_labellar_types.md").is_file()
    fig = fast_run / "figures" / "F1_architecture.png"
    assert fig.is_file() and fig.stat().st_size > 10_000
    assert (fast_run / "figures" / "F1_architecture.svg").is_file()
    assert (fast_run / "figures" / "captions.md").is_file()


def test_fast_run_writes_nothing_into_the_repo(fast_run: Path) -> None:
    """--outdir must redirect every artefact; nothing leaks into papers/."""
    assert fast_run != PAPERS
    doc = json.loads((fast_run / "results.json").read_text())
    for entry in doc["figures"] + doc["tables"]:
        assert not Path(entry["path"]).is_absolute()
        assert (fast_run / entry["path"]).is_file()


def test_fast_results_validate(fast_run: Path) -> None:
    doc = json.loads((fast_run / "results.json").read_text())
    assert validate_results(doc) == []
    assert doc["fast"] is True
    # provenance is always added, whatever --only asks for
    assert set(doc["steps"]) == set(QUICK_STEPS) | {"provenance"}


def test_fast_run_is_deterministic(tmp_path: Path) -> None:
    """Same seed, same commit, same values (runtimes excepted)."""
    runs = []
    for i in range(2):
        out = tmp_path / f"run{i}"
        assert main(["--fast", "--outdir", str(out), "--quiet", "--only", "census"]) == 0
        doc = json.loads((out / "results.json").read_text())
        runs.append({k: v["value"] for k, v in doc["values"].items()})
    assert runs[0] == runs[1]


def test_library_call_with_no_argv_only_prints_guidance(capsys, tmp_path: Path) -> None:
    """`flylab reproduce-paper` calls main() with no argv.

    It must not inherit the host process's command line and must not spend
    minutes writing into papers/ as a side effect of being imported.
    """
    before = sorted(PAPERS.glob("**/*")) if PAPERS.exists() else []
    assert main() == 0
    out = capsys.readouterr().out
    assert "scripts/reproduce_paper.py" in out
    assert "--fast" in out
    assert sorted(PAPERS.glob("**/*")) == before


def test_unknown_step_is_rejected(tmp_path: Path) -> None:
    assert main(["--outdir", str(tmp_path), "--quiet", "--only", "does_not_exist"]) == 2


# --------------------------------------------------------------------------
# the committed record vs the paper
# --------------------------------------------------------------------------
@pytest.mark.skipif(not COMMITTED_RESULTS.exists(), reason="papers/results.json not built yet")
def test_committed_results_validate() -> None:
    doc = json.loads(COMMITTED_RESULTS.read_text())
    assert validate_results(doc) == []


@pytest.mark.skipif(not COMMITTED_RESULTS.exists(), reason="papers/results.json not built yet")
def test_committed_results_cover_every_key_the_paper_quotes() -> None:
    doc = json.loads(COMMITTED_RESULTS.read_text())
    values = doc["values"]
    missing = sorted(template_keys() - set(values))
    assert not missing, f"paper template references keys absent from results.json: {missing}"


@pytest.mark.skipif(not COMMITTED_RESULTS.exists(), reason="papers/results.json not built yet")
def test_committed_results_cover_the_declared_paper_keys() -> None:
    values = json.loads(COMMITTED_RESULTS.read_text())["values"]
    missing = sorted(set(PAPER_KEYS) - set(values))
    assert not missing, f"PAPER_KEYS absent from results.json: {missing}"


@pytest.mark.skipif(not TEMPLATE.exists(), reason="paper template not present")
def test_declared_paper_keys_are_actually_used() -> None:
    """PAPER_KEYS is the test contract; it must not rot into a wish list."""
    unused = sorted(set(PAPER_KEYS) - template_keys())
    assert not unused, f"PAPER_KEYS entries no longer used by the template: {unused}"


@pytest.mark.skipif(not COMMITTED_RESULTS.exists(), reason="papers/results.json not built yet")
def test_committed_run_was_not_a_fast_run() -> None:
    doc = json.loads(COMMITTED_RESULTS.read_text())
    assert doc["fast"] is False, "papers/results.json came from a --fast run; rebuild it"


@pytest.mark.skipif(not COMMITTED_RESULTS.exists(), reason="papers/results.json not built yet")
def test_rendered_draft_has_no_unresolved_placeholders() -> None:
    draft = PAPERS / "IJRC_FlyLab_draft.md"
    if not draft.exists():
        pytest.skip("draft not rendered yet")
    text = draft.read_text()
    assert "[[MISSING:" not in text
    assert "{{" not in text


# --------------------------------------------------------------------------
# the whole pipeline (minutes)
# --------------------------------------------------------------------------
@pytest.mark.slow
def test_full_fast_pipeline(tmp_path: Path) -> None:
    """Every step, at --fast, into a temp dir. Minutes, not seconds."""
    out = tmp_path / "full"
    assert main(["--fast", "--outdir", str(out), "--quiet"]) == 0
    doc = json.loads((out / "results.json").read_text())
    assert validate_results(doc) == []
    assert [n for n, s in doc["steps"].items() if s["status"] != "ok"] == []
    assert set(doc["steps"]) == set(STEP_NAMES)
    assert len([f for f in doc["figures"] if f["name"].endswith(".png")]) == 10
    missing = sorted(set(PAPER_KEYS) - set(doc["values"]))
    assert not missing, f"a full run did not produce: {missing}"
