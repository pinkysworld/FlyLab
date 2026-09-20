"""Contract tests for the tutorial and the interpretation guide.

Two things are checked, and they are the two that rot first:

1. **Every command the tutorial prints exists.**  `docs/TUTORIAL.md` and
   `docs/INTERPRETATION.md` are full of pasteable `flylab ...` lines, and a
   renamed command would leave a documented command that does not run.  Every
   such line is parsed out of the fenced blocks and checked against the Typer
   app's registered commands, options included.
2. **The lesson runner runs.**  ``flylab tutorial <n>`` executes real commands,
   so a lesson whose step has drifted fails here rather than in front of a
   newcomer.  Two fast lessons are run end to end.

The documents' *numbers* are not asserted: they are a record of a run, and the
tests that pin the science live beside the modules that compute it.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from typer.main import get_command
from typer.testing import CliRunner

from flylab.cli import TUTORIAL_LESSONS, app

runner = CliRunner(env={"COLUMNS": "250", "TERM": "dumb", "NO_COLOR": "1"})

DOCS = Path(__file__).resolve().parents[1] / "docs"
TUTORIAL = DOCS / "TUTORIAL.md"
INTERPRETATION = DOCS / "INTERPRETATION.md"

#: a fenced code block, with its info string
FENCE = re.compile(r"^```(\w*)\n(.*?)^```", re.M | re.S)

#: lines that show a shell prompt as well as the command
PROMPT = re.compile(r"^\s*\$\s+")


def _click_app():
    return get_command(app)


def _registered() -> dict[str, object]:
    """Command name -> click command, including sub-groups."""
    return dict(_click_app().commands)


def _documented_commands(text: str) -> list[list[str]]:
    """Every `flylab ...` invocation inside a fenced block, as argv."""
    found: list[list[str]] = []
    for info, body in FENCE.findall(text):
        if info not in ("", "bash", "sh", "console"):
            continue
        for raw in body.splitlines():
            line = PROMPT.sub("", raw).strip()
            if not line.startswith("flylab "):
                continue
            # drop a trailing comment and any shell redirection: what is being
            # checked is the command and its options, not the shell around it
            line = line.split("#", 1)[0].split(">", 1)[0].split("|", 1)[0]
            if "$" in line:
                # a shell loop variable is not a literal argument; the same
                # command appears unrolled elsewhere in the document
                continue
            found.append(line.split()[1:])
    return found


def _lesson_commands() -> list[list[str]]:
    return [list(step["argv"]) for le in TUTORIAL_LESSONS for step in le["steps"] if "argv" in step]


# --------------------------------------------------------------------------
# the documents exist and cross-link
# --------------------------------------------------------------------------
def test_both_documents_exist():
    assert TUTORIAL.is_file()
    assert INTERPRETATION.is_file()


def test_the_documents_link_to_each_other_and_to_the_guardrail():
    tut = TUTORIAL.read_text()
    interp = INTERPRETATION.read_text()
    assert "INTERPRETATION.md" in tut
    assert "TUTORIAL.md" in interp
    # the interpretation guide must defer to NOVELTY.md rather than restate it
    assert "NOVELTY.md" in interp


def test_the_tutorial_has_a_lesson_anchor_for_every_registered_lesson():
    tut = TUTORIAL.read_text()
    for lesson in TUTORIAL_LESSONS:
        assert f'<a id="lesson-{lesson["id"]}"></a>' in tut, lesson["id"]


# --------------------------------------------------------------------------
# every documented command exists
# --------------------------------------------------------------------------
def test_the_tutorial_actually_shows_commands():
    # a parser that silently found nothing would make the next test vacuous
    assert len(_documented_commands(TUTORIAL.read_text())) >= 10


@pytest.mark.parametrize(
    "argv",
    _documented_commands(TUTORIAL.read_text())
    + _documented_commands(INTERPRETATION.read_text())
    + _lesson_commands(),
    ids=lambda a: " ".join(a) or "<bare>",
)
def test_every_documented_command_exists(argv):
    """The command resolves, and `--help` accepts the options as written."""
    assert argv, "a bare `flylab` in a document is not a lesson"
    commands = _registered()
    assert argv[0] in commands, f"{argv[0]!r} is not a flylab command; have {sorted(commands)}"
    # `--help` type-checks the options without running the analysis.  A few
    # commands take positional arguments that --help does not need, so the
    # options are what is being checked here.
    result = runner.invoke(app, argv + ["--help"])
    assert result.exit_code == 0, result.output


# --------------------------------------------------------------------------
# the lesson list
# --------------------------------------------------------------------------
def test_lesson_ids_and_keys_are_unique_and_sequential():
    ids = [le["id"] for le in TUTORIAL_LESSONS]
    keys = [le["key"] for le in TUTORIAL_LESSONS]
    assert ids == list(range(1, len(ids) + 1))
    assert len(set(keys)) == len(keys)


def test_every_lesson_says_what_it_does_not_tell_you():
    for le in TUTORIAL_LESSONS:
        assert le["teaches"].strip()
        assert le["not_"].strip(), le["key"]
        assert le["steps"], le["key"]


def test_tutorial_list():
    result = runner.invoke(app, ["tutorial", "--list"])
    assert result.exit_code == 0, result.output
    for le in TUTORIAL_LESSONS:
        assert le["title"] in result.output
        assert le["key"] in result.output
    assert "docs/TUTORIAL.md" in result.output
    assert "docs/INTERPRETATION.md" in result.output


def test_tutorial_with_no_argument_lists_the_lessons():
    assert runner.invoke(app, ["tutorial"]).output == runner.invoke(
        app, ["tutorial", "--list"]
    ).output


def test_tutorial_json_index_is_machine_readable():
    import json

    result = runner.invoke(app, ["tutorial", "--json"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert len(payload) == len(TUTORIAL_LESSONS)
    assert payload[0]["id"] == 1
    assert payload[0]["does_not_tell_you"]
    assert all(c.startswith("flylab ") or c == "python" or c.startswith("python ")
               for row in payload for c in row["commands"])


def test_an_unknown_lesson_is_refused_with_the_list():
    result = runner.invoke(app, ["tutorial", "42"])
    assert result.exit_code == 2
    assert "no tutorial lesson" in result.output


# --------------------------------------------------------------------------
# the lesson runner
# --------------------------------------------------------------------------
def test_lesson_7_runs_end_to_end():
    """Typed evidence: a good row is accepted, three bad ones are refused."""
    result = runner.invoke(app, ["tutorial", "7"])
    assert result.exit_code == 0, result.output
    out = result.output
    assert "functional_engagement" in out
    assert "not_modelled" in out
    # the refusals are the point of the lesson
    assert out.count("EvidenceTypeError") >= 5
    assert "ACCEPTED - this is a bug" not in out
    assert TUTORIAL_LESSONS[6]["not_"][:40] in out


def test_lesson_4_runs_the_real_compare():
    """Receptor and circuit selectivity ordering, from the real payload."""
    result = runner.invoke(app, ["tutorial", "4"])
    assert result.exit_code == 0, result.output
    out = result.output
    assert "flylab compare imidacloprid deltamethrin" in out
    assert "Imidacloprid" in out and "Deltamethrin" in out
    assert "what this does NOT tell you" in out


def test_a_lesson_can_be_selected_by_key():
    by_id = runner.invoke(app, ["tutorial", "7"])
    by_key = runner.invoke(app, ["tutorial", "library"])
    assert by_key.exit_code == 0
    assert by_key.output == by_id.output


def test_lesson_json_does_not_execute_anything():
    import json

    result = runner.invoke(app, ["tutorial", "5", "--json"])
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["key"] == "dependence"
    assert payload["commands"][0].startswith("flylab dependence")
