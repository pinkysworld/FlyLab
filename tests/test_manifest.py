"""The repository artifact manifest and the dependency lock.

What is being pinned down: every class of input the committed results depend on
is covered, a change to any of them is detected and named, and a dependency
difference is reported without failing a release -- because a manifest is
verified on machines other than the one that wrote it.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.manifest import (  # noqa: E402
    FATAL_SECTIONS,
    MANIFEST_VERSION,
    SECTIONS,
    build_manifest,
    dependency_closure,
    main,
    repo_root,
    verify_manifest,
    write_lockfile,
    write_manifest,
)

ROOT = repo_root()


@pytest.fixture(scope="module")
def manifest():
    return build_manifest(ROOT)


# --------------------------------------------------------------------------
# coverage
# --------------------------------------------------------------------------
def test_every_artifact_class_is_covered(manifest):
    assert set(manifest["sections"]) == set(SECTIONS) | {"dependencies"}
    assert manifest["flylab_manifest_version"] == MANIFEST_VERSION
    for name, section in manifest["sections"].items():
        assert len(section["sha256"]) == 64, name
        assert section["description"]
        assert section["n_entries"] == len(section["entries"])


def test_the_code_section_covers_the_package_and_the_build_metadata(manifest):
    entries = manifest["sections"]["code"]["entries"]
    for path in ("flylab/cli.py", "flylab/spec.py", "flylab/report/card.py", "pyproject.toml"):
        assert path in entries, path
    assert all(not p.endswith(".pyc") and "__pycache__" not in p for p in entries)
    assert all(len(sha) == 64 for sha in entries.values())


def test_the_data_sections_cover_the_committed_inputs(manifest):
    graphs = manifest["sections"]["graphs"]["entries"]
    assert "data/derived/malecns_named_neighborhood.json" in graphs
    assert all(p.startswith("data/derived/") for p in graphs)

    literature = manifest["sections"]["literature"]["entries"]
    assert len(literature) >= 6
    assert all(p.endswith(".yaml") for p in literature)

    library = manifest["sections"]["library"]["entries"]
    assert list(library) == ["flylab/pharm/library.yaml"]


def test_the_library_hash_is_the_one_notebooks_record(manifest):
    from flylab.pharm.occupancy import library_sha256

    assert manifest["sections"]["library"]["entries"]["flylab/pharm/library.yaml"] == library_sha256()


def test_the_paper_section_covers_the_committed_artifacts(manifest):
    paper = manifest["sections"]["paper"]["entries"]
    assert "papers/results.json" in paper
    assert any(p.startswith("papers/figures/") for p in paper)
    assert any(p.startswith("papers/tables/") for p in paper)


def test_no_feather_is_ever_hashed(manifest):
    for section in manifest["sections"].values():
        assert not any(str(p).endswith(".feather") for p in section["entries"])


def test_which_sections_are_fatal(manifest):
    for name, section in manifest["sections"].items():
        assert section["fatal"] is (name in FATAL_SECTIONS)
    assert "dependencies" not in FATAL_SECTIONS


# --------------------------------------------------------------------------
# verify
# --------------------------------------------------------------------------
def test_a_freshly_written_manifest_verifies(tmp_path, manifest):
    path = write_manifest(manifest, tmp_path / "m.json")
    diff = verify_manifest(path, ROOT)
    assert diff["ok"] is True
    assert all(section["ok"] for section in diff["sections"].values())
    assert all(section["summary"].endswith("entries match") for section in diff["sections"].values())


def test_a_changed_artifact_is_named_and_fails_the_check(tmp_path, manifest):
    tampered = json.loads(json.dumps(manifest))
    entries = tampered["sections"]["paper"]["entries"]
    victim = sorted(entries)[0]
    entries[victim] = "0" * 64
    path = tmp_path / "m.json"
    path.write_text(json.dumps(tampered))

    diff = verify_manifest(path, ROOT)
    assert diff["ok"] is False
    paper = diff["sections"]["paper"]
    assert paper["ok"] is False
    assert [c["path"] for c in paper["changed"]] == [victim]
    assert paper["changed"][0]["committed"] == "0" * 64
    assert "changed" in paper["summary"]
    assert "flylab manifest --write" in diff["how_to_fix"]


def test_added_and_removed_files_are_reported(tmp_path, manifest):
    tampered = json.loads(json.dumps(manifest))
    entries = tampered["sections"]["literature"]["entries"]
    entries["data/literature/does_not_exist.yaml"] = "1" * 64
    removed = sorted(entries)[0]
    del entries[removed]
    path = tmp_path / "m.json"
    path.write_text(json.dumps(tampered))

    diff = verify_manifest(path, ROOT)["sections"]["literature"]
    assert diff["removed"] == ["data/literature/does_not_exist.yaml"]
    assert diff["added"] == [removed]


def test_a_dependency_difference_is_reported_but_not_fatal(tmp_path, manifest):
    tampered = json.loads(json.dumps(manifest))
    packages = tampered["sections"]["dependencies"]["entries"]
    name = sorted(packages)[0]
    packages[name] = "0.0.0-from-another-machine"
    path = tmp_path / "m.json"
    path.write_text(json.dumps(tampered))

    diff = verify_manifest(path, ROOT)
    assert diff["ok"] is True
    assert diff["sections"]["dependencies"]["ok"] is False
    assert any("dependencies" in note for note in diff["notes"])


def test_a_missing_manifest_is_an_explicit_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="flylab manifest --write"):
        verify_manifest(tmp_path / "absent.json", ROOT)


def test_a_manifest_from_another_format_version_is_flagged(tmp_path, manifest):
    other = json.loads(json.dumps(manifest))
    other["flylab_manifest_version"] = "99.0"
    path = tmp_path / "m.json"
    path.write_text(json.dumps(other))
    assert any("format" in note for note in verify_manifest(path, ROOT)["notes"])


# --------------------------------------------------------------------------
# dependency closure and the lock
# --------------------------------------------------------------------------
def test_the_closure_starts_from_the_declared_dependencies():
    closure = dependency_closure(ROOT)
    assert {"numpy", "pyyaml", "typer", "pydantic"} <= set(closure["roots"])
    installed = {name.lower() for name in closure["packages"]}
    assert {"numpy", "pyyaml"} <= installed
    assert all(version for version in closure["packages"].values())


def test_the_closure_is_narrower_than_pip_freeze():
    """A lock filtered to the closure, not a dump of the whole environment."""
    from importlib import metadata

    everything = {d.metadata["Name"].lower() for d in metadata.distributions() if d.metadata["Name"]}
    closure = {name.lower() for name in dependency_closure(ROOT)["packages"]}
    assert closure <= everything
    assert len(closure) < len(everything)


def test_the_lockfile_pins_versions_and_says_what_it_is_not(tmp_path):
    written = write_lockfile(ROOT, tmp_path / "requirements-lock.txt")
    text = Path(written["path"]).read_text()
    assert "WHAT IT IS NOT" in text and "full environment specification" in text
    pins = [line for line in text.splitlines() if line and not line.startswith("#")]
    assert len(pins) == written["n_packages"] >= 5
    assert all("==" in line for line in pins)
    assert any(line.lower().startswith("numpy==") for line in pins)


def test_the_committed_lockfile_exists_and_is_pinned():
    lock = ROOT / "requirements-lock.txt"
    assert lock.exists(), "run `flylab manifest --lock`"
    pins = [line for line in lock.read_text().splitlines() if line and not line.startswith("#")]
    assert pins and all("==" in line for line in pins)


# --------------------------------------------------------------------------
# the script entry point
# --------------------------------------------------------------------------
def test_main_writes_and_then_checks(tmp_path, capsys):
    path = tmp_path / "m.json"
    assert main(["--write", "--path", str(path)]) == 0
    assert path.exists()
    assert main(["--check", "--path", str(path)]) == 0
    out = capsys.readouterr().out
    assert "manifest check" in out


def test_main_exits_non_zero_on_a_mismatch(tmp_path, manifest, capsys):
    tampered = json.loads(json.dumps(manifest))
    tampered["sections"]["graphs"]["entries"]["data/derived/ghost.json"] = "2" * 64
    path = tmp_path / "m.json"
    path.write_text(json.dumps(tampered))
    assert main(["--check", "--path", str(path)]) == 1
    assert "MISMATCH" in capsys.readouterr().out


def test_main_reports_a_missing_manifest_without_a_traceback(tmp_path, capsys):
    assert main(["--check", "--path", str(tmp_path / "absent.json")]) == 2
    assert "does not exist" in capsys.readouterr().err
