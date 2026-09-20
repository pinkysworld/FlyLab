"""Declarative experiment specs: validation, the run directory, determinism.

The contract this file pins down is the one a researcher relies on: a spec file
either runs or fails with a message naming the valid options, a run writes a
directory that describes itself, and the same spec with the same seed produces
the same numbers and the same artifact hashes.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
import yaml

from flylab.spec import (
    ANALYSES,
    ENGINES,
    SPEC_VERSION,
    VOLATILE_KEYS,
    ExperimentSpec,
    SpecError,
    analysis_entry_point,
    estimate_runtime,
    load_spec,
    make_spec,
    run_spec,
    spec_schema,
    spec_sha256,
    strip_volatile,
)

#: the example from the module docstring, shrunk so the suite stays fast
SPEC_YAML = """\
name: fipronil_connectome_dependence
compound: fipronil
concentrations: [1.0e-8, 1.0e-6]
graph: named
readouts: [mean_hz, mn9_hz]
engines: [rate]
"""


def write(tmp_path: Path, text: str, name: str = "spec.yaml") -> Path:
    path = tmp_path / name
    path.write_text(text)
    return path


# --------------------------------------------------------------------------
# loading and validation
# --------------------------------------------------------------------------
def test_load_spec_resolves_defaults(tmp_path):
    spec = load_spec(write(tmp_path, SPEC_YAML))
    assert spec.name == "fipronil_connectome_dependence"
    assert spec.compounds == ["fipronil"]          # `compound:` is folded into `compounds`
    assert spec.concentrations == [1e-8, 1e-6]
    assert spec.engines == ["rate"]
    assert spec.readouts == ["mean_hz", "mn9_hz"]
    # defaults the file never mentioned
    assert spec.replicates == 1 and spec.seed == 0 and spec.include_vehicle is True
    assert spec.flylab_spec_version == SPEC_VERSION
    assert spec.analyses == {}


def test_multiple_compounds_are_supported(tmp_path):
    spec = load_spec(write(tmp_path, "compounds: [fipronil, imidacloprid]\nconcentrations: [1.0e-6]\n"))
    assert spec.compounds == ["fipronil", "imidacloprid"]


def test_compound_and_compounds_merge_without_losing_either(tmp_path):
    spec = load_spec(
        write(tmp_path, "compound: fipronil\ncompounds: [imidacloprid]\nconcentrations: [1.0e-6]\n")
    )
    assert spec.compounds == ["fipronil", "imidacloprid"]


def test_a_v05_design_file_is_a_valid_spec(tmp_path):
    """`flylab experiment run` is an alias, so its design files must still load."""
    spec = load_spec(
        write(tmp_path, "assay: subgraph\ncompounds: [nicotine]\nconcs_M: [1.0e-6]\nreplicates: 2\n")
    )
    assert spec.engines == ["subgraph"] and ENGINES[spec.engines[0]] == "subgraph"
    assert spec.concentrations == [1e-6] and spec.replicates == 2


@pytest.mark.parametrize(
    "bad, must_mention",
    [
        ({"compuonds": ["fipronil"]}, "compounds"),
        ({"compound": "not_a_drug"}, "fipronil"),
        ({"graph": "wholebrain"}, "named"),
        ({"engines": ["quantum"]}, "rate"),
        ({"readouts": ["glow"]}, "mean_hz"),
        ({"analyses": {"dependance": {}}}, "dependence"),
        ({"analyses": {"dependence": {"perms": 10}}}, "permutations"),
        ({"analyses": {"dependence": {"correction": "bonferroni"}}}, "benjamini-hochberg"),
        ({"replicates": 0}, "replicates"),
        ({"concentrations": [-1.0]}, ">= 0"),
        ({"concentrations": []}, "concentration"),
    ],
)
def test_invalid_specs_name_the_valid_options(bad, must_mention):
    with pytest.raises(SpecError) as exc:
        make_spec({"compounds": ["fipronil"], "concentrations": [1e-6], **bad})
    assert must_mention in str(exc.value)


def test_unknown_key_suggests_the_near_miss():
    with pytest.raises(SpecError) as exc:
        make_spec({"concentrasions": [1e-6]})
    assert "concentrations" in str(exc.value)


def test_missing_file_is_a_spec_error(tmp_path):
    with pytest.raises(SpecError, match="not found"):
        load_spec(tmp_path / "nope.yaml")


def test_broken_yaml_is_a_spec_error(tmp_path):
    with pytest.raises(SpecError, match="not valid YAML"):
        load_spec(write(tmp_path, "name: [unclosed\n"))


def test_spec_hash_is_stable_and_sensitive():
    a = make_spec({"compounds": ["fipronil"], "concentrations": [1e-6]})
    b = make_spec({"compound": "fipronil", "concentrations": [1e-6]})
    c = make_spec({"compounds": ["fipronil"], "concentrations": [1e-6], "seed": 1})
    assert spec_sha256(a) == spec_sha256(b)
    assert spec_sha256(a) != spec_sha256(c)


# --------------------------------------------------------------------------
# schema
# --------------------------------------------------------------------------
def test_spec_schema_is_machine_readable_and_complete():
    schema = spec_schema()
    assert schema["type"] == "object" and schema["additionalProperties"] is False
    props = schema["properties"]
    for field in ("compounds", "concentrations", "graph", "engines", "readouts", "analyses", "seed"):
        assert field in props, field
        assert props[field].get("description")
    assert props["engines"]["items"]["enum"] == ["rate", "lif", "taste", "taste_map"]
    assert "fipronil" in props["compounds"]["items"]["enum"]
    assert set(schema["x-analyses"]) == set(ANALYSES)
    for meta in schema["x-analyses"].values():
        assert meta["module"].startswith("flylab.analysis.")
        assert meta["artifact"].endswith(".json")
    assert schema["x-volatile-keys"] == sorted(VOLATILE_KEYS)
    # the example in the schema is itself a valid spec
    assert make_spec(yaml.safe_load(schema["example"])).compounds == ["fipronil"]


def test_schema_documents_the_run_directory():
    layout = spec_schema()["x-run-directory"]
    for name in ("manifest.json", "spec.yaml", "notebooks/", "results.csv", "cards/", "figures/"):
        assert name in layout


# --------------------------------------------------------------------------
# dry run
# --------------------------------------------------------------------------
def test_dry_run_writes_nothing_and_estimates(tmp_path):
    spec = make_spec(
        {
            "compounds": ["fipronil"],
            "concentrations": [1e-8, 1e-6],
            "analyses": {"dependence": {"permutations": 1000}},
        }
    )
    outdir = tmp_path / "run"
    out = run_spec(spec, outdir, dry_run=True)
    assert out["dry_run"] is True
    assert not outdir.exists()
    assert out["runtime_estimate"]["estimate_s"] > 0
    assert "dependence.json" in out["would_write"]
    assert out["spec"]["concentrations"] == [1e-8, 1e-6]


def test_runtime_estimate_scales_with_the_design():
    small = estimate_runtime({"compounds": ["fipronil"], "concentrations": [1e-6]})
    big = estimate_runtime(
        {"compounds": ["fipronil", "imidacloprid"], "concentrations": [1e-8, 1e-6], "replicates": 3}
    )
    assert big["n_assay_runs"] > small["n_assay_runs"]
    assert big["estimate_s"] > small["estimate_s"]


# --------------------------------------------------------------------------
# a real run
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    outdir = tmp_path_factory.mktemp("spec_run") / "run"
    spec = make_spec(yaml.safe_load(SPEC_YAML))
    manifest = run_spec(spec, outdir)
    return outdir, manifest


def test_run_writes_a_self_describing_directory(run_dir):
    outdir, manifest = run_dir
    assert (outdir / "manifest.json").exists()
    assert (outdir / "spec.yaml").exists()
    assert (outdir / "results.csv").exists()
    assert (outdir / "notebooks").is_dir()
    # 2 concentrations + 1 vehicle, one engine
    assert len(list((outdir / "notebooks").glob("*.json"))) == 3
    assert manifest["spec_sha256"] and manifest["determinism"]["digest"]
    assert manifest["code"]["flylab_version"]
    assert manifest["map_id"] == "male-cns:v1.0"
    assert manifest["graphs"]["named"]["n_nodes"] > 0
    assert manifest["seeds"]["global"] == 0
    assert manifest["timings_s"]["total_s"] >= 0


def test_resolved_spec_yaml_carries_the_filled_in_defaults(run_dir):
    outdir, _ = run_dir
    text = (outdir / "spec.yaml").read_text()
    resolved = yaml.safe_load(text)
    assert resolved["replicates"] == 1 and resolved["include_vehicle"] is True
    assert resolved["compounds"] == ["fipronil"]
    assert "spec_sha256:" in text


def test_results_csv_has_one_row_per_run(run_dir):
    outdir, manifest = run_dir
    rows = list(csv.DictReader((outdir / "results.csv").open()))
    assert len(rows) == 3 == manifest["n_rows"]
    assert list(rows[0])[:5] == ["engine", "condition", "compound", "conc_M", "replicate"]
    assert {r["condition"] for r in rows} == {"drug", "vehicle"}
    assert float(rows[0]["mean_hz"]) >= 0.0


def test_only_requested_analyses_are_written(run_dir):
    outdir, manifest = run_dir
    for name, meta in ANALYSES.items():
        assert not (outdir / meta["artifact"]).exists(), name
    assert manifest["analyses_requested"] == []


def test_every_artifact_is_hashed_in_the_manifest(run_dir):
    outdir, manifest = run_dir
    on_disk = {
        p.relative_to(outdir).as_posix()
        for p in outdir.rglob("*")
        if p.is_file() and p.name != "manifest.json"
    }
    assert set(manifest["artifacts"]) == on_disk
    for meta in manifest["artifacts"].values():
        assert len(meta["sha256"]) == 64 and len(meta["content_sha256"]) == 64
        assert meta["bytes"] > 0


def test_cards_are_written_for_the_headline_result(run_dir):
    outdir, manifest = run_dir
    card_path = outdir / "cards" / "fipronil__rate.json"
    assert card_path.exists() and (outdir / "cards" / "fipronil__rate.md").exists()
    card = json.loads(card_path.read_text())
    assert card["subject"]["compound"] == "fipronil"
    # the headline dose is the highest in the spec
    assert card["subject"]["concentration_M"] == pytest.approx(1e-6)
    assert "cards/fipronil__rate.md" in manifest["cards"]


def test_notebooks_keep_their_warnings(run_dir):
    outdir, _ = run_dir
    notebook = json.loads(next((outdir / "notebooks").glob("*fipronil*")).read_text())
    assert notebook["live_lab"] is None
    assert notebook["warnings"]


# --------------------------------------------------------------------------
# determinism
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def twin_runs(tmp_path_factory):
    """The same spec, run twice into two directories. Figures included."""
    root = tmp_path_factory.mktemp("twin")
    spec = make_spec(yaml.safe_load(SPEC_YAML))
    return root, run_spec(spec, root / "a"), run_spec(spec, root / "b")


def test_the_same_spec_and_seed_reproduce_byte_for_byte(twin_runs):
    tmp_path, first, second = twin_runs

    assert (tmp_path / "a" / "results.csv").read_bytes() == (tmp_path / "b" / "results.csv").read_bytes()
    assert (tmp_path / "a" / "spec.yaml").read_bytes() == (tmp_path / "b" / "spec.yaml").read_bytes()
    assert first["determinism"]["digest"] == second["determinism"]["digest"]
    assert set(first["artifacts"]) == set(second["artifacts"])
    for path, meta in first["artifacts"].items():
        assert meta["content_sha256"] == second["artifacts"][path]["content_sha256"], path


def test_a_different_seed_changes_the_spec_hash_but_the_rate_engine_is_deterministic(tmp_path):
    base = make_spec({**yaml.safe_load(SPEC_YAML), "figures": False})
    shifted = make_spec({**yaml.safe_load(SPEC_YAML), "figures": False, "seed": 7})
    assert spec_sha256(base) != spec_sha256(shifted)
    a = run_spec(base, tmp_path / "a")
    b = run_spec(shifted, tmp_path / "b")
    # the rate engine carries no RNG, so the readouts must not move with the seed
    assert (tmp_path / "a" / "results.csv").read_text() == (tmp_path / "b" / "results.csv").read_text()
    assert a["seeds"]["global"] == 0 and b["seeds"]["global"] == 7


def test_volatile_keys_are_what_makes_the_bytes_differ(twin_runs):
    tmp_path, first, _second = twin_runs
    notebooks = [p for p in first["artifacts"] if p.startswith("notebooks/")]
    assert notebooks
    for path in notebooks:
        a = json.loads((tmp_path / "a" / path).read_text())
        b = json.loads((tmp_path / "b" / path).read_text())
        assert strip_volatile(a) == strip_volatile(b)
        assert set(first["artifacts"][path]["volatile_fields_removed"]) <= VOLATILE_KEYS


# --------------------------------------------------------------------------
# analyses and degradation
# --------------------------------------------------------------------------
def test_dependence_is_written_and_carries_its_settings(tmp_path):
    spec = make_spec(
        {
            "name": "dep",
            "compounds": ["fipronil"],
            "concentrations": [1e-6],
            "readouts": ["mean_hz"],
            "figures": False,
            "analyses": {"dependence": {"permutations": 5, "correction": "benjamini-hochberg"}},
        }
    )
    manifest = run_spec(spec, tmp_path / "run")
    payload = json.loads((tmp_path / "run" / "dependence.json").read_text())
    assert payload["settings"]["permutations"] == 5
    assert payload["entry_point"].startswith("flylab.analysis.dependence.")
    assert payload["cells"] and payload["cells"][0]["compound"] == "fipronil"
    assert payload["multiplicity_correction"]["method"] == "benjamini-hochberg"
    assert manifest["analyses_completed"] == ["dependence"]
    # the card picks the analysis up
    card = json.loads((tmp_path / "run" / "cards" / "fipronil__rate.json").read_text())
    assert card["connectome_dependence"]["assessed"] is True


def test_analysis_entry_points_are_discovered_not_hardcoded():
    fn, dotted = analysis_entry_point("dependence")
    assert callable(fn) and dotted == "flylab.analysis.dependence.dependence_profile"
    with pytest.raises(SpecError, match="Valid options"):
        analysis_entry_point("telepathy")


def test_a_renamed_entry_point_produces_an_actionable_error(monkeypatch):
    monkeypatch.setitem(
        ANALYSES,
        "dependence",
        {**ANALYSES["dependence"], "candidates": ("dependence_profile_v2",)},
    )
    with pytest.raises(SpecError) as exc:
        analysis_entry_point("dependence")
    message = str(exc.value)
    assert "dependence_profile_v2" in message and "flylab.analysis.dependence" in message
    assert "currently exports" in message


def test_figures_degrade_gracefully_without_matplotlib(tmp_path, monkeypatch):
    import builtins

    real_import = builtins.__import__

    def no_matplotlib(name, *args, **kw):
        if name.startswith("matplotlib"):
            raise ImportError("no matplotlib in this build")
        return real_import(name, *args, **kw)

    monkeypatch.setattr(builtins, "__import__", no_matplotlib)
    spec = make_spec({"compounds": ["fipronil"], "concentrations": [1e-6], "readouts": ["mean_hz"]})
    manifest = run_spec(spec, tmp_path / "run")
    assert manifest["figures"]["available"] is False
    assert "matplotlib" in manifest["figures"]["reason"]
    assert (tmp_path / "run" / "results.csv").exists()  # everything else still written
    assert any("no figures" in w for w in manifest["warnings"])


def test_figures_are_written_when_matplotlib_is_installed(tmp_path):
    pytest.importorskip("matplotlib")
    spec = make_spec({"compounds": ["fipronil"], "concentrations": [1e-8, 1e-6], "readouts": ["mean_hz"]})
    manifest = run_spec(spec, tmp_path / "run")
    assert manifest["figures"]["available"] is True
    assert manifest["figures"]["files"] == ["figures/rate_mean_hz.png"]
    assert (tmp_path / "run" / "figures" / "rate_mean_hz.png").exists()


def test_the_run_states_that_nothing_was_measured_in_a_fly(run_dir):
    _, manifest = run_dir
    assert manifest["label"] == "model_derived"
    assert "living fly" in manifest["disclaimer"]
    assert any("simulated" in w or "No row is a measurement" in w for w in manifest["warnings"])


def test_spec_model_rejects_unknown_fields():
    with pytest.raises(Exception):
        ExperimentSpec(nonsense=1)


# --------------------------------------------------------------------------
# the pydantic-free twin
# --------------------------------------------------------------------------
def test_the_spec_resolves_identically_without_pydantic():
    """The science core must import and validate where pydantic does not exist.

    Same pattern as ``flylab/assays/experiment.py``: the pydantic model is the
    reference and the dataclass fallback has to match it field for field,
    including the error messages.
    """
    import builtins
    import importlib
    import sys

    from flylab.spec import HAVE_PYDANTIC

    assert HAVE_PYDANTIC, "this test documents the pydantic path"
    raw = {"compound": "fipronil", "concentrations": [1e-6], "engines": ["rate"], "seed": 3}
    reference = make_spec(raw).model_dump()

    blocked = ("pydantic", "fastapi", "starlette", "typer", "uvicorn", "pandas", "pyarrow")
    real_import = builtins.__import__

    def guarded(name, *args, **kw):
        if name.split(".")[0] in blocked:
            raise ImportError(name)
        return real_import(name, *args, **kw)

    saved = {m: sys.modules[m] for m in list(sys.modules) if m.startswith("flylab")}
    try:
        for mod in list(sys.modules):
            if mod.startswith("flylab"):
                del sys.modules[mod]
        builtins.__import__ = guarded
        fallback = importlib.import_module("flylab.spec")
        assert fallback.HAVE_PYDANTIC is False
        assert fallback.make_spec(raw).model_dump() == reference
        with pytest.raises(fallback.SpecError, match="unknown compound"):
            fallback.make_spec({"compound": "not_a_drug"})
        with pytest.raises(fallback.SpecError, match="Valid options"):
            fallback.make_spec({"nonsense": 1})
    finally:
        builtins.__import__ = real_import
        for mod in list(sys.modules):
            if mod.startswith("flylab"):
                del sys.modules[mod]
        sys.modules.update(saved)

