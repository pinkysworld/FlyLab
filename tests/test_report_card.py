"""Claim cards: nine sections, always all nine, never a quiet omission.

The point of a card is that the sections a claim would rather not have are the
ones it cannot drop.  These tests fail if a section disappears, if an
unassessed analysis is silently omitted instead of being marked
"not assessed in this run", or if the independent-validation section ever
reports anything other than "none".
"""

from __future__ import annotations

import json

import pytest

from flylab.report import card as card_mod
from flylab.report.card import (
    CARD_VERSION,
    SECTIONS,
    SECTION_TITLES,
    cards_for_run,
    claim_card,
    to_markdown,
)
from flylab.spec import make_spec, run_spec

COMPOUND = "fipronil"
CONC = 1e-6


@pytest.fixture(scope="module")
def notebook():
    from flylab.assays.ensemble import run_assay

    return run_assay("subgraph", COMPOUND, CONC, seed=0)


@pytest.fixture(scope="module")
def vehicle():
    from flylab.assays.ensemble import run_assay

    return run_assay("subgraph", None, 0.0, seed=0)


@pytest.fixture(scope="module")
def card(notebook, vehicle):
    return claim_card(notebook, compound=COMPOUND, conc_M=CONC, engine="rate", vehicle=vehicle)


# --------------------------------------------------------------------------
# shape
# --------------------------------------------------------------------------
def test_a_card_has_exactly_the_nine_sections(card):
    assert card["sections"] == list(SECTIONS)
    assert len(SECTIONS) == 9
    for key in SECTIONS:
        assert key in card, key
        assert isinstance(card[key], dict)
    assert card["flylab_card_version"] == CARD_VERSION
    assert set(SECTION_TITLES) == set(SECTIONS)


def test_the_claim_is_one_sentence_naming_compound_dose_and_readout(card):
    sentence = card["claim"]["sentence"]
    assert sentence.count(". ") <= 1 and sentence.endswith(".")
    assert COMPOUND in sentence and "MaleCNS" in sentence
    assert card["claim"]["readout"] in sentence
    assert card["claim"]["value"] is not None
    assert card["claim"]["vehicle"] is not None        # the vehicle was supplied
    assert card["claim"]["relative_change"] is not None


def test_status_is_model_derived_and_says_so(card):
    assert card["status"]["status"] == "model-derived"
    assert card["status"]["label"] == "model_derived"
    assert "living fly" in card["status"]["statement"]


def test_evidence_inputs_carry_parameter_type_and_evidence_distance(card):
    section = card["evidence_inputs"]
    assert section["n_rows"] >= 1
    assert section["n_modelled"] + section["n_not_modelled"] == section["n_rows"]
    modelled = [r for r in section["rows"] if r["modelled"]]
    assert modelled, "fipronil must have at least one typed row"
    for row in modelled:
        assert row["parameter_type"] in ("Kd", "Ki", "EC50", "IC50", "Kb")
        assert row["evidence_distance"]                 # the source relation
        assert row["source"]
        assert row["engagement"] is not None
    for row in section["rows"]:
        if not row["modelled"]:
            assert row["engagement"] is None            # N/A, never a small number


def test_measured_not_assumed_is_the_malecns_edges(card):
    section = card["measured_not_assumed"]
    assert section["what"] == "the MaleCNS v1.0 edges"
    assert section["n_nodes"] and section["n_edges"]
    assert "MaleCNS" in (section["citation"] or "")
    assert "one link" in section["statement"]


def test_the_three_named_assumptions_are_always_present(card):
    ids = {item["id"] for item in card["assumptions_introduced"]["items"]}
    assert {"transmitter_sign", "gain_mapping", "free_concentration"} <= ids


def test_independent_biological_validation_is_none(card):
    section = card["independent_biological_validation"]
    assert section["status"] == "none"
    assert section["statement"].startswith("None.")
    assert section["what_would_change_it"]


def test_named_unknowns_are_named_not_parameterised(card):
    items = card["named_unknowns"]["items"]
    ids = {item["id"] for item in items}
    assert {"cns_free_concentration", "mn9_receptor_expression", "biological_variance"} <= ids
    for item in items:
        assert item["would_resolve"]


# --------------------------------------------------------------------------
# the sections that depend on an analysis
# --------------------------------------------------------------------------
def test_unassessed_sections_say_so_rather_than_disappearing(card):
    for key in ("connectome_dependence", "specification_stability"):
        section = card[key]
        assert section["assessed"] is False
        assert "not assessed in this run" in section["statement"]
        assert section["how_to_assess"]


def test_a_requested_but_uncovered_analysis_reads_differently(notebook):
    payload = {"analysis": "dependence", "cells": [{"compound": "imidacloprid", "conc_M": 1e-6, "modes": []}]}
    section = claim_card(
        notebook,
        compound=COMPOUND,
        conc_M=CONC,
        dependence=payload,
        requested_analyses=["dependence"],
    )["connectome_dependence"]
    assert section["assessed"] is False
    assert "requested" in section["statement"]


def test_dependence_lists_the_nulls_it_is_distinguishable_from(notebook):
    payload = {
        "analysis": "dependence",
        "cells": [
            {
                "compound": COMPOUND,
                "conc_M": CONC,
                "class": "topology-dependent",
                "n": 200,
                "readout": "mean_hz",
                "necessary_information_level": {"level": "wiring_without_transmitter_identity"},
                "modes": [
                    {"mode": "weight_permute", "p_two_sided": 0.006, "beats_null": True, "n": 200,
                     "information_kept": "the real edge list and the transmitters"},
                    {"mode": "sign_permute", "p_two_sided": 0.42, "beats_null": False, "n": 200,
                     "information_kept": "the real edges"},
                ],
            }
        ],
    }
    section = claim_card(
        notebook, compound=COMPOUND, conc_M=CONC, dependence=payload, requested_analyses=["dependence"]
    )["connectome_dependence"]
    assert section["assessed"] is True
    assert [m["null_model"] for m in section["distinguishable_from"]] == ["weight_permute"]
    assert [m["null_model"] for m in section["not_distinguishable_from"]] == ["sign_permute"]
    assert section["distinguishable_from"][0]["p_two_sided"] == pytest.approx(0.006)
    assert section["class"] == "topology-dependent"
    # the level block is flattened to its name, and kept whole beside it
    assert section["necessary_information_level"] == "wiring_without_transmitter_identity"
    assert isinstance(section["necessary_information_level_detail"], dict)
    assert "weight_permute" in section["statement"] and "0.006" in section["statement"]


def test_dependence_with_no_beaten_null_says_the_connectome_was_not_needed(notebook):
    payload = {
        "cells": [
            {
                "compound": COMPOUND,
                "conc_M": CONC,
                "class": "composition-dominated",
                "modes": [
                    {"mode": "weight_permute", "p_two_sided": 0.58, "beats_null": False},
                    {"mode": "erdos_renyi", "p_two_sided": 0.47, "beats_null": False},
                ],
            }
        ]
    }
    section = claim_card(
        notebook, compound=COMPOUND, conc_M=CONC, dependence=payload, requested_analyses=["dependence"]
    )["connectome_dependence"]
    assert section["distinguishable_from"] == []
    assert "does not need the measured wiring" in section["statement"]


def test_specification_stability_reports_retention_out_of_the_family(notebook):
    payload = {
        "analysis": "robustness",
        "result": {
            "family_size": 25,
            "default_spec": "flylab_biphasic@1.00x",
            "rows": [
                {
                    "conclusion": "C3_rdl_disinhibition",
                    "claim": "fipronil (1 uM) disinhibits the circuit",
                    "n_specs": 25,
                    "n_retained": 25,
                    "fraction_retained": 1.0,
                    "fragile": False,
                },
                {
                    "conclusion": "C1_nicotinic_suppression",
                    "claim": "a nicotinic agonist suppresses network activity",
                    "n_specs": 25,
                    "n_retained": 15,
                    "fraction_retained": 0.6,
                    "fragile": True,
                    "failing_specs": ["linear@1.00x"],
                },
            ],
        },
    }
    section = claim_card(
        notebook, compound=COMPOUND, conc_M=CONC, robustness=payload, requested_analyses=["robustness"]
    )["specification_stability"]
    assert section["assessed"] is True
    assert section["family_size"] == 25
    retained = {c["conclusion"]: c["retained"] for c in section["conclusions"]}
    assert retained == {"C3_rdl_disinhibition": 25, "C1_nicotinic_suppression": 15}
    assert section["conclusions_mentioning_compound"] == ["C3_rdl_disinhibition"]
    assert "25/25" in section["statement"]


# --------------------------------------------------------------------------
# markdown
# --------------------------------------------------------------------------
def test_markdown_renders_every_section_in_order(card):
    text = to_markdown(card)
    positions = [text.index(f"## {SECTION_TITLES[key]}") for key in SECTIONS if key != "claim"]
    assert positions == sorted(positions)
    assert text.startswith("# Claim card -- fipronil")
    assert "> " in text                      # the claim is a blockquote
    assert "| receptor | parameter type |" in text
    assert "not assessed in this run" in text
    assert text.endswith("\n")


def test_markdown_shortens_long_sources_but_json_keeps_them(card):
    text = to_markdown(card)
    long_sources = [r["source"] for r in card["evidence_inputs"]["rows"] if len(str(r["source"])) > 200]
    assert long_sources, "fipronil carries at least one long provenance note"
    assert long_sources[0] not in text
    assert "…" in text
    for line in text.splitlines():
        if line.startswith("| insect_RDL"):
            assert line.count("|") == 7      # the table stays a table


def test_markdown_survives_a_card_with_no_evidence_rows():
    minimal = claim_card({"compound": None, "concentration_M": None})
    text = to_markdown(minimal)
    for key in SECTIONS:
        assert f"## {SECTION_TITLES[key]}" in text


# --------------------------------------------------------------------------
# cards for a run directory
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    outdir = tmp_path_factory.mktemp("cards_run") / "run"
    spec = make_spec(
        {
            "name": "cards",
            "compounds": ["fipronil", "imidacloprid"],
            "concentrations": [1e-8, 1e-6],
            "readouts": ["mean_hz"],
            "figures": False,
        }
    )
    run_spec(spec, outdir)
    return outdir


def test_cards_for_run_reads_the_written_cards(run_dir):
    cards = cards_for_run(run_dir)
    assert len(cards) == 2
    assert sorted(c["subject"]["compound"] for c in cards) == ["fipronil", "imidacloprid"]
    for c in cards:
        assert c["sections"] == list(SECTIONS)
        assert c["subject"]["concentration_M"] == pytest.approx(1e-6)


def test_cards_for_run_rebuilds_from_notebooks_when_no_cards_were_written(run_dir, tmp_path):
    import shutil

    copy = tmp_path / "copy"
    shutil.copytree(run_dir, copy)
    shutil.rmtree(copy / "cards")
    cards = cards_for_run(copy)
    assert len(cards) == 2
    assert all(c["claim"]["sentence"] for c in cards)
    # the vehicle notebook in the directory is found, so the claim is a change
    assert all(c["claim"]["vehicle"] is not None for c in cards)


def test_cards_for_run_rejects_a_missing_directory(tmp_path):
    with pytest.raises(FileNotFoundError):
        cards_for_run(tmp_path / "nothing_here")


def test_card_json_round_trips(card):
    assert json.loads(json.dumps(card, default=str))["sections"] == list(SECTIONS)


def test_the_module_is_importable_without_the_server():
    """The card is a library function; no web framework may be required."""
    import importlib

    assert callable(getattr(importlib.import_module("flylab.report"), "claim_card"))
    assert card_mod.__name__ == "flylab.report.card"
