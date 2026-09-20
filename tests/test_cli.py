"""CLI contract tests.

Two things are checked for every command: ``--help`` works (so a typo in an
option annotation cannot ship), and the command actually runs and prints the
numbers it promises.  Commands whose module has not landed yet must exit with a
friendly message, never a traceback.
"""

from __future__ import annotations

import json

import pytest

from flylab import __version__ as flylab_version
from typer.testing import CliRunner

from flylab.cli import app

# Typer renders help through Rich: it wraps at the terminal width and emits ANSI
# codes, so a narrow CI terminal splits long option names like
# "--closure-min-weight" across lines and the plain `in` check fails. Pin a wide,
# colourless terminal for every invocation so help assertions are stable.
runner = CliRunner(env={"COLUMNS": "250", "TERM": "dumb", "NO_COLOR": "1"})

#: every command and sub-command the bench exposes
COMMANDS = [
    ["occupancy"],
    ["compare"],
    ["dashboard"],
    ["dependence"],
    ["ablation"],
    ["stability"],
    ["uncertainty"],
    ["voi"],
    ["claims"],
    ["list-drugs"],
    ["meta"],
    ["assay"],
    ["assay-cns"],
    ["assay-subgraph"],
    ["assay-spiking"],
    ["assay-taste-map"],
    ["ensemble"],
    ["sensitivity"],
    ["ic50"],
    ["exposure"],
    ["null-panel"],
    ["selectivity"],
    ["predictions"],
    ["reproduce-paper"],
    ["graph"],
    ["graph", "info"],
    ["graph", "impact"],
    ["experiment"],
    ["experiment", "run"],
    ["experiment", "example"],
    ["run"],
    ["spec-schema"],
    ["card"],
    ["manifest"],
    ["extract-subgraph"],
    ["download-malecns"],
    ["serve"],
    ["tutorial"],
]


@pytest.mark.parametrize("cmd", COMMANDS, ids=lambda c: " ".join(c))
def test_help_works_for_every_command(cmd):
    result = runner.invoke(app, cmd + ["--help"])
    assert result.exit_code == 0, result.output
    assert "Usage" in result.output


def test_no_args_shows_help():
    result = runner.invoke(app, [])
    assert "Usage" in result.output


def test_tutorial_lists_lessons_and_runs_one():
    """`flylab tutorial` is the executable half of docs/TUTORIAL.md.

    The full contract (every documented command exists, two lessons run end to
    end) lives in tests/test_tutorial.py; this is the CLI-level smoke test.
    """
    listing = runner.invoke(app, ["tutorial", "--list"])
    assert listing.exit_code == 0, listing.output
    assert "docs/TUTORIAL.md" in listing.output

    lesson = runner.invoke(app, ["tutorial", "7"])
    assert lesson.exit_code == 0, lesson.output
    assert "EvidenceTypeError" in lesson.output
    assert "what this does NOT tell you" in lesson.output


# --------------------------------------------------------------------------
# library / occupancy
# --------------------------------------------------------------------------
def test_occupancy_table_and_json():
    result = runner.invoke(app, ["occupancy", "imidacloprid", "--conc", "1e-6"])
    assert result.exit_code == 0
    assert "insect_nAChR" in result.output
    assert "vertebrate_nAChR_a4b2" in result.output

    result = runner.invoke(app, ["occupancy", "imidacloprid", "--conc", "1e-6", "--json"])
    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert payload["key"] == "imidacloprid"
    assert payload["selectivity"]["nAChR"]["insect_occupancy"] > 0.9


def test_occupancy_unknown_compound_is_an_error_not_a_crash():
    result = runner.invoke(app, ["occupancy", "not-a-drug", "--conc", "1e-6"])
    assert result.exit_code != 0


def test_compare_two_compounds():
    result = runner.invoke(
        app, ["compare", "imidacloprid", "fipronil", "--conc", "1e-6", "--no-dependence"]
    )
    assert result.exit_code == 0
    assert "Imidacloprid" in result.output and "Fipronil" in result.output
    # the decision columns, then the per-compound occupancy tables
    assert "recSI" in result.output and "cirSI" in result.output
    assert "insect_RDL" in result.output


def test_compare_json_is_the_compare_endpoint_payload():
    payload = json.loads(
        runner.invoke(
            app,
            [
                "compare",
                "imidacloprid",
                "deltamethrin",
                "--conc",
                "1e-6",
                "--no-dependence",
                "--json",
            ],
        ).output
    )
    assert {r["compound"] for r in payload["rows"]} == {"imidacloprid", "deltamethrin"}
    assert payload["best_receptor_si"] == "imidacloprid"
    assert payload["best_circuit_si"] == "deltamethrin"
    assert payload["runtime_estimate"]["estimate_s"] > 0


# --------------------------------------------------------------------------
# dashboard / decision layer
# --------------------------------------------------------------------------
def test_dashboard_prints_the_decision_blocks():
    result = runner.invoke(app, ["dashboard", "fipronil", "--conc", "1e-6", "--n", "20"])
    assert result.exit_code == 0, result.output
    out = result.output
    assert "Fipronil" in out and "1.00e-06 M" in out
    assert "insect engagement" in out and "vertebrate engagement" in out
    # the vertebrate number the v0.5 correction made unhideable
    assert "0.476" in out
    assert "vertebrate reaches 20% engagement at" in out
    assert "specific wiring evidence: present (topology-dependent)" in out
    assert "what can I trust?" in out
    assert "Live validation" in out and "none" in out
    assert "why this happened" in out


def test_dashboard_json_matches_the_endpoint():
    payload = json.loads(
        runner.invoke(
            app, ["dashboard", "imidacloprid", "--conc", "1e-6", "--no-dependence", "--json"]
        ).output
    )
    assert payload["compound"]["key"] == "imidacloprid"
    assert payload["headline"]["insect_engagement"]["value"] > 0.9
    assert payload["trust"]["blended_confidence"] is None
    assert payload["why"]["inputs"]["potency_is_inert"] is True


def test_dependence_command_names_the_verdict():
    result = runner.invoke(app, ["dependence", "imidacloprid", "--conc", "1e-6", "--n", "20"])
    assert result.exit_code == 0, result.output
    assert "composition-dominated" in result.output
    assert "resolution" in result.output


def test_dependence_landscape_estimates_before_running():
    result = runner.invoke(app, ["dependence", "--landscape"])
    assert result.exit_code == 0, result.output
    assert "estimated" in result.output
    assert "--run" in result.output


def test_ablation_command():
    result = runner.invoke(app, ["ablation", "imidacloprid", "--conc", "1e-6"])
    assert result.exit_code == 0, result.output
    for level in ("A_receptor_only", "B_composition_only", "C_topology_only", "D_full_flylab"):
        assert level in result.output


@pytest.mark.parametrize("cmd", [["stability"], ["uncertainty"], ["voi"]])
def test_expensive_commands_can_print_an_estimate_only(cmd):
    result = runner.invoke(app, cmd + ["--estimate"])
    assert result.exit_code == 0, result.output
    assert "estimated" in result.output


def test_claims_command_prints_the_chain_and_the_three_columns():
    result = runner.invoke(app, ["claims", "imidacloprid", "--conc", "1e-6"])
    assert result.exit_code == 0, result.output
    out = result.output
    for step in ("parameter_source", "mechanism_rule", "malecns_edges", "readout"):
        assert step in out
    assert "FACT" in out and "INFERENCE" in out and "UNKNOWN" in out
    assert "MaleCNS v1.0" in out

    payload = json.loads(
        runner.invoke(app, ["claims", "fipronil", "--conc", "1e-6", "--json"]).output
    )
    assert payload["label_counts"]["OBSERVED"] == 1


def test_list_drugs():
    result = runner.invoke(app, ["list-drugs"])
    assert result.exit_code == 0
    assert "imidacloprid" in result.output

    payload = json.loads(runner.invoke(app, ["list-drugs", "--json"]).output)
    assert any(row["key"] == "fipronil" for row in payload)


def test_meta_prints_graphs_and_library_hash():
    result = runner.invoke(app, ["meta"])
    assert result.exit_code == 0
    # Never assert a literal version: three surfaces had drifted to 0.5.0
    # while the package was 0.6.0, and hardcoded tests hid it.
    assert f"FlyLab {flylab_version}" in result.output
    assert "taste_motor" in result.output
    assert "sha256" in result.output

    payload = json.loads(runner.invoke(app, ["meta", "--json"]).output)
    assert payload["graphs"]["named"]["available"] is True


# --------------------------------------------------------------------------
# assays
# --------------------------------------------------------------------------
def test_assay_taste():
    result = runner.invoke(app, ["assay", "--compound", "imidacloprid", "--conc", "1e-6"])
    assert result.exit_code == 0
    assert "mn9_sugar_hz" in result.output


def test_assay_cns():
    result = runner.invoke(app, ["assay-cns", "--compound", "imidacloprid", "--conc", "1e-6"])
    assert result.exit_code == 0
    assert "cns_excitation_index" in result.output


def test_assay_subgraph_on_both_graphs():
    for graph in ("named", "taste_motor"):
        result = runner.invoke(app, ["assay-subgraph", "--graph", graph, "--conc", "1e-6"])
        assert result.exit_code == 0, result.output
        assert "mean_hz" in result.output


def test_assay_spiking():
    result = runner.invoke(app, ["assay-spiking", "--t-ms", "120", "--conc", "1e-6"])
    assert result.exit_code == 0, result.output
    assert "LIF" in result.output and "mn9_hz" in result.output


def test_assay_taste_map():
    result = runner.invoke(app, ["assay-taste-map", "--conc", "1e-6"])
    assert result.exit_code == 0, result.output
    assert "n_sweet_grn" in result.output
    assert "hypothesis" in result.output


# --------------------------------------------------------------------------
# uncertainty / analysis
# --------------------------------------------------------------------------
def test_ensemble():
    result = runner.invoke(app, ["ensemble", "--n-rep", "3"])
    assert result.exit_code == 0, result.output
    assert "n_rep=3" in result.output


def test_sensitivity():
    result = runner.invoke(app, ["sensitivity", "--readout", "mean_hz", "--conc", "1e-6"])
    assert result.exit_code == 0, result.output
    assert "drive_hz" in result.output


def test_ic50_says_it_is_model_derived():
    result = runner.invoke(app, ["ic50", "--n-boot", "10", "--n-rep", "1", "--readout", "mean_hz"])
    assert result.exit_code == 0, result.output
    assert "model IC50" in result.output
    assert "not an animal IC50" in result.output


def test_exposure():
    result = runner.invoke(app, ["exposure", "--dose", "2", "--route", "bath"])
    assert result.exit_code == 0, result.output
    assert "Cmax" in result.output and "AUC" in result.output


# --------------------------------------------------------------------------
# graph
# --------------------------------------------------------------------------
def test_graph_info():
    result = runner.invoke(app, ["graph", "info", "--graph", "taste_motor"])
    assert result.exit_code == 0, result.output
    assert "LB3b" in result.output
    assert "acetylcholine" in result.output

    payload = json.loads(runner.invoke(app, ["graph", "info", "--json"]).output)
    assert payload["n_nodes"] > 0


def test_graph_impact():
    result = runner.invoke(app, ["graph", "impact", "--top", "5", "--conc", "1e-6"])
    assert result.exit_code == 0, result.output
    assert "gains" in result.output
    assert "delta" in result.output


# --------------------------------------------------------------------------
# experiments
# --------------------------------------------------------------------------
def test_experiment_example_is_a_runnable_design(tmp_path):
    result = runner.invoke(app, ["experiment", "example"])
    assert result.exit_code == 0
    assert "assay:" in result.output and "concs_M:" in result.output

    import yaml

    spec = yaml.safe_load(result.output)
    assert spec["assay"] == "subgraph"
    assert len(spec["compounds"]) >= 2


def test_experiment_run_writes_csv_and_json(tmp_path):
    design = tmp_path / "design.yaml"
    design.write_text(
        "assay: subgraph\ncompounds: [imidacloprid]\nconcs_M: [1.0e-8, 1.0e-6]\n"
        "replicates: 1\nreadouts: [mn9_hz, mean_hz]\ninclude_vehicle: true\n"
    )
    out_csv = tmp_path / "results.csv"
    out_json = tmp_path / "results.json"
    result = runner.invoke(
        app, ["experiment", "run", str(design), "--out", str(out_csv), "--json", str(out_json)]
    )
    assert result.exit_code == 0, result.output
    lines = out_csv.read_text().strip().splitlines()
    assert lines[0].startswith("compound,conc_M,replicate")
    assert len(lines) == 4  # header + 2 drug rows + 1 vehicle row
    payload = json.loads(out_json.read_text())
    assert payload["n_rows"] == 2
    assert payload["label"] == "model_derived"


def test_experiment_run_prints_csv_when_no_output_path(tmp_path):
    design = tmp_path / "design.yaml"
    design.write_text("compounds: [nicotine]\nconcs_M: [1.0e-6]\nreplicates: 1\ninclude_vehicle: false\n")
    result = runner.invoke(app, ["experiment", "run", str(design)])
    assert result.exit_code == 0, result.output
    assert "compound,conc_M,replicate" in result.output


# --------------------------------------------------------------------------
# lazy commands
# --------------------------------------------------------------------------
#: lazy commands, with arguments small enough that a landed module still runs fast
LAZY = [
    ["selectivity", "--conc", "1e-6"],
    ["predictions", "--n-rep", "1"],
    ["null-panel", "--n", "1", "--modes", "weight_permute"],
    ["reproduce-paper"],
]


@pytest.mark.parametrize("cmd", LAZY, ids=lambda c: c[0])
def test_lazy_commands_run_or_say_why_not(cmd):
    """Either the module has landed and the command works, or it exits 2 with a note."""
    result = runner.invoke(app, cmd)
    if result.exit_code == 0:
        return
    assert result.exit_code == 2, result.output
    assert result.exception is None or isinstance(result.exception, SystemExit)
    assert "not" in result.output.lower()


# --------------------------------------------------------------------------
# options that only exist in v0.5
# --------------------------------------------------------------------------
def test_extract_subgraph_has_the_new_options():
    out = runner.invoke(app, ["extract-subgraph", "--help"]).output
    assert "--closure-min-weight" in out
    assert "--out" in out


def test_serve_has_open_flag():
    assert "--open" in runner.invoke(app, ["serve", "--help"]).output


def test_json_flag_is_offered_on_table_commands():
    for cmd in ("occupancy", "list-drugs", "meta", "assay-subgraph", "ensemble", "ic50", "exposure"):
        assert "--json" in runner.invoke(app, [cmd, "--help"]).output, cmd


# --------------------------------------------------------------------------
# specs, cards and the artifact manifest
# --------------------------------------------------------------------------
SPEC_YAML = """\
name: cli_spec
compound: fipronil
concentrations: [1.0e-6]
readouts: [mean_hz]
engines: [rate]
figures: false
"""


def test_spec_schema_prints_the_fields_and_the_analyses():
    out = runner.invoke(app, ["spec-schema"])
    assert out.exit_code == 0, out.output
    for field in ("compounds", "concentrations", "engines", "analyses"):
        assert field in out.output
    for analysis in ("dependence", "robustness", "uncertainty"):
        assert analysis in out.output


def test_spec_schema_json_is_a_schema():
    out = runner.invoke(app, ["spec-schema", "--json"])
    schema = json.loads(out.output)
    assert schema["type"] == "object"
    assert "compounds" in schema["properties"]


def test_spec_schema_example_is_a_runnable_spec(tmp_path):
    out = runner.invoke(app, ["spec-schema", "--example"])
    assert out.exit_code == 0
    path = tmp_path / "spec.yaml"
    path.write_text(out.output)
    dry = runner.invoke(app, ["run", str(path), "--dry-run", "--outdir", str(tmp_path / "run")])
    assert dry.exit_code == 0, dry.output


@pytest.fixture(scope="module")
def cli_run(tmp_path_factory):
    """One `flylab run` for the whole module; the card tests read it back."""
    root = tmp_path_factory.mktemp("cli_run")
    spec = root / "spec.yaml"
    spec.write_text(SPEC_YAML)
    outdir = root / "run"
    result = runner.invoke(app, ["run", str(spec), "--outdir", str(outdir)])
    assert result.exit_code == 0, result.output
    return outdir, result


def test_run_dry_run_validates_and_estimates_without_running(tmp_path):
    spec = tmp_path / "spec.yaml"
    spec.write_text(SPEC_YAML)
    outdir = tmp_path / "run"
    result = runner.invoke(app, ["run", str(spec), "--outdir", str(outdir), "--dry-run"])
    assert result.exit_code == 0, result.output
    assert "estimate:" in result.output
    assert "spec_sha256" in result.output
    assert not outdir.exists()


def test_run_writes_a_run_directory(cli_run):
    outdir, result = cli_run
    assert (outdir / "manifest.json").exists()
    assert (outdir / "results.csv").exists()
    assert (outdir / "cards").is_dir()
    assert "digest" in result.output


def test_run_reports_an_invalid_spec_without_a_traceback(tmp_path):
    spec = tmp_path / "spec.yaml"
    spec.write_text("compound: not_a_real_drug\nconcentrations: [1.0e-6]\n")
    result = runner.invoke(app, ["run", str(spec), "--outdir", str(tmp_path / "run")])
    assert result.exit_code == 2
    assert result.exception is None or isinstance(result.exception, SystemExit)
    assert "unknown compound" in result.output


def test_experiment_run_is_an_alias_of_run(tmp_path):
    """The v0.5 entry point keeps working and gains the new flags."""
    design = tmp_path / "design.yaml"
    design.write_text("assay: subgraph\ncompounds: [fipronil]\nconcs_M: [1.0e-6]\n")
    outdir = tmp_path / "run"
    result = runner.invoke(app, ["experiment", "run", str(design), "--outdir", str(outdir), "--dry-run"])
    assert result.exit_code == 0, result.output
    assert "estimate:" in result.output
    assert not outdir.exists()


def test_card_prints_every_section():
    """The standalone path: no run directory, the assay is run on the spot."""
    result = runner.invoke(app, ["card", "fipronil", "--conc", "1e-6"])
    assert result.exit_code == 0, result.output
    from flylab.report.card import SECTIONS

    for section in SECTIONS:
        assert section in result.output
    assert "not assessed in this run" in result.output


def test_card_reads_a_run_directory_as_json(cli_run):
    outdir, _ = cli_run
    result = runner.invoke(app, ["card", "--run", str(outdir), "--json"])
    assert result.exit_code == 0, result.output
    card = json.loads(result.output)
    assert card["subject"]["compound"] == "fipronil"
    assert card["independent_biological_validation"]["status"] == "none"


def test_card_renders_markdown(cli_run):
    outdir, _ = cli_run
    result = runner.invoke(app, ["card", "--run", str(outdir), "--markdown"])
    assert result.exit_code == 0, result.output
    assert result.output.startswith("# Claim card")
    assert "## Independent biological validation" in result.output


def test_card_on_a_missing_run_directory_is_an_error_not_a_crash(tmp_path):
    result = runner.invoke(app, ["card", "--run", str(tmp_path / "nothing")])
    assert result.exit_code == 2
    assert result.exception is None or isinstance(result.exception, SystemExit)


def test_manifest_check_reports_rather_than_crashes():
    """The committed manifest may legitimately be stale, or absent in a slim tree.

    The contract is only that --check reports -- never a traceback -- and exits
    0 when it matches, 2 when it does not or when there is nothing to check.
    """
    result = runner.invoke(app, ["manifest", "--check"])
    assert result.exit_code in (0, 2), result.output
    assert result.exception is None or isinstance(result.exception, SystemExit)
    if "does not exist" in result.output:
        assert "--write" in result.output
        return
    assert "manifest check" in result.output
    for section in ("code", "library", "graphs", "literature", "paper"):
        assert section in result.output


def test_manifest_write_and_check_agree(tmp_path):
    path = tmp_path / "m.json"
    written = runner.invoke(app, ["manifest", "--write", "--path", str(path)])
    assert written.exit_code == 0, written.output
    assert path.exists()
    checked = runner.invoke(app, ["manifest", "--check", "--path", str(path)])
    assert checked.exit_code == 0, checked.output
    assert "MISMATCH" not in checked.output


def test_manifest_check_without_a_manifest_says_how_to_make_one(tmp_path):
    result = runner.invoke(app, ["manifest", "--check", "--path", str(tmp_path / "absent.json")])
    assert result.exit_code == 2
    assert "--write" in result.output


