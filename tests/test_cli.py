"""CLI contract tests.

Two things are checked for every command: ``--help`` works (so a typo in an
option annotation cannot ship), and the command actually runs and prints the
numbers it promises.  Commands whose module has not landed yet must exit with a
friendly message, never a traceback.
"""

from __future__ import annotations

import json

import pytest
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
    ["extract-subgraph"],
    ["download-malecns"],
    ["serve"],
]


@pytest.mark.parametrize("cmd", COMMANDS, ids=lambda c: " ".join(c))
def test_help_works_for_every_command(cmd):
    result = runner.invoke(app, cmd + ["--help"])
    assert result.exit_code == 0, result.output
    assert "Usage" in result.output


def test_no_args_shows_help():
    result = runner.invoke(app, [])
    assert "Usage" in result.output


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
    result = runner.invoke(app, ["compare", "imidacloprid", "fipronil", "--conc", "1e-6"])
    assert result.exit_code == 0
    assert "Imidacloprid" in result.output and "Fipronil" in result.output


def test_list_drugs():
    result = runner.invoke(app, ["list-drugs"])
    assert result.exit_code == 0
    assert "imidacloprid" in result.output

    payload = json.loads(runner.invoke(app, ["list-drugs", "--json"]).output)
    assert any(row["key"] == "fipronil" for row in payload)


def test_meta_prints_graphs_and_library_hash():
    result = runner.invoke(app, ["meta"])
    assert result.exit_code == 0
    assert "FlyLab 0.5.0" in result.output
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
