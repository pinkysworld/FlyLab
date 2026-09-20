"""The claim-provenance layer.

``flylab/analysis/claims.py`` is what stops a FlyLab number being quoted
without its chain: every result is walked back through the parameter source,
the engagement transformation, the mechanism rule, the expression assumption,
the MaleCNS edges, the transmitter predictions, the engine and the readout,
and every link says how strong it is.  These tests pin the vocabulary, the
order, and the three things the summary must never lose sight of -- that the
connectome is a measurement, that the gain rule is not, and that the CNS free
concentration is unknown.
"""

from __future__ import annotations

import pytest

from flylab.analysis.claims import (
    CHAIN_STEPS,
    CLASSIFICATIONS,
    LABELS,
    claim_audit,
    fact_inference_unknown,
    label_counts,
    to_markdown,
)
from flylab.assays.subgraph import run_subgraph_assay


@pytest.fixture(scope="module")
def notebook():
    return run_subgraph_assay("imidacloprid", 1e-6)


@pytest.fixture(scope="module")
def audit(notebook):
    return claim_audit(notebook)


# --------------------------------------------------------------------------
# vocabulary
# --------------------------------------------------------------------------
def test_labels_and_classifications_are_the_documented_five():
    assert LABELS == (
        "OBSERVED",
        "LITERATURE-DERIVED",
        "MODEL-ASSUMPTION",
        "COMPUTED",
        "UNKNOWN",
    )
    assert CLASSIFICATIONS == (
        "LITERATURE",
        "MODEL-DERIVED",
        "MODEL-ASSUMPTION",
        "PREDICTION",
        "NOT MODELLED",
    )


def test_chain_covers_every_documented_step():
    assert CHAIN_STEPS == (
        "result",
        "parameter_source",
        "engagement_transformation",
        "mechanism_rule",
        "expression_assumption",
        "malecns_edges",
        "transmitter_predictions",
        "engine",
        "readout",
    )


# --------------------------------------------------------------------------
# the chain
# --------------------------------------------------------------------------
def test_audit_emits_the_whole_chain_in_order(audit):
    steps = [link["step"] for link in audit["chain"]]
    assert steps == list(CHAIN_STEPS)
    assert [link["order"] for link in audit["chain"]] == list(range(len(CHAIN_STEPS)))


def test_every_link_is_labelled_and_chipped(audit):
    for link in audit["chain"]:
        assert link["label"] in LABELS, link
        assert link["classification"] in CLASSIFICATIONS, link
        assert link["statement"].strip()


def test_only_the_connectome_is_an_observation_of_this_system(audit):
    observed = [link["step"] for link in audit["chain"] if link["label"] == "OBSERVED"]
    assert observed == ["malecns_edges"]
    edges = next(link for link in audit["chain"] if link["step"] == "malecns_edges")
    assert edges["classification"] == "LITERATURE"
    assert edges["detail"]["n_nodes"] > 0 and edges["detail"]["n_edges"] > 0
    assert "MaleCNS" in (edges["detail"].get("citation") or "")


def test_the_gain_rule_is_never_labelled_as_evidence(audit):
    rule = next(link for link in audit["chain"] if link["step"] == "mechanism_rule")
    assert rule["label"] == "MODEL-ASSUMPTION"
    assert rule["classification"] == "MODEL-ASSUMPTION"
    assert any("fitted" in a for a in rule["assumptions"])


def test_transmitters_are_flagged_as_predictions(audit):
    nt = next(link for link in audit["chain"] if link["step"] == "transmitter_predictions")
    assert nt["classification"] == "PREDICTION"
    assert nt["unknowns"]


def test_expression_link_carries_the_coverage_and_the_mn9_gap(audit):
    expr = next(link for link in audit["chain"] if link["step"] == "expression_assumption")
    fraction = expr["detail"]["fraction_cells_mapped"]
    assert fraction is None or 0.0 < fraction < 1.0
    assert any("MN9" in u for u in expr["unknowns"])


def test_parameter_sources_carry_their_type_and_citation(audit):
    sources = next(link for link in audit["chain"] if link["step"] == "parameter_source")["sources"]
    assert sources
    modelled = [s for s in sources if s["classification"] != "NOT MODELLED"]
    assert modelled, "imidacloprid has sourced rows"
    for row in modelled:
        assert row["param_type"] in ("EC50", "IC50", "Kd", "Ki", "Kb")
        assert row["param_value_M"] > 0
        assert row["source"]
        assert row["species"]
    not_modelled = [s for s in sources if s["classification"] == "NOT MODELLED"]
    for row in not_modelled:
        assert row["param_value_M"] is None, "a placeholder must not report a number"


def test_label_counts_add_up(audit):
    counts = label_counts(audit["chain"])
    assert sum(counts.values()) == len(audit["chain"])
    assert counts["OBSERVED"] == 1


# --------------------------------------------------------------------------
# fact / inference / unknown
# --------------------------------------------------------------------------
def test_fact_inference_unknown_is_attached_to_every_audit(audit):
    fiu = audit["fact_inference_unknown"]
    assert set(fiu) >= {"facts", "model_inference", "unknown", "counts"}


def test_facts_are_the_topology_and_the_sourced_parameters(audit):
    facts = audit["fact_inference_unknown"]["facts"]
    kinds = [f["kind"] for f in facts]
    assert "connectome" in kinds
    assert "parameter" in kinds
    connectome = next(f for f in facts if f["kind"] == "connectome")
    assert "MaleCNS" in connectome["statement"]
    for row in (f for f in facts if f["kind"] == "parameter"):
        assert row["label"] == "LITERATURE-DERIVED"
        assert row["source"]


def test_inference_is_the_predicted_rate_change(audit):
    inference = audit["fact_inference_unknown"]["model_inference"]
    rates = [i for i in inference if i["kind"] == "rate_change"]
    assert rates, "a circuit notebook must yield predicted rate changes"
    assert {"mean_hz", "mn9_hz", "dnp01_hz"} <= {i["readout"] for i in rates}
    for row in rates:
        assert row["classification"] == "PREDICTION"
    mean = next(i for i in rates if i["readout"] == "mean_hz")
    assert mean["vehicle"] > mean["treated"], "imidacloprid suppresses the mean rate"


def test_circuit_indices_are_inference_not_fact():
    out = fact_inference_unknown(
        {"compound": "imidacloprid", "concentration_M": 1e-6},
        indices={"receptor_si_log10": 2.7, "circuit_si_log10": 1.79},
    )
    indices = [i for i in out["model_inference"] if i["kind"] == "index"]
    assert len(indices) == 2
    assert all(i["classification"] == "MODEL-DERIVED" for i in indices)
    assert not any(f["kind"] == "index" for f in out["facts"])


def test_unknown_names_the_three_gaps_the_model_cannot_fill(audit):
    ids = {u["id"] for u in audit["fact_inference_unknown"]["unknown"]}
    assert {"cns_free_concentration", "mn9_receptor_expression", "biological_variance"} <= ids
    for row in audit["fact_inference_unknown"]["unknown"]:
        assert row["label"] == "UNKNOWN"
        assert row["would_resolve"]


def test_placeholder_receptors_become_unknowns(audit):
    ids = {u["id"] for u in audit["fact_inference_unknown"]["unknown"]}
    assert any(i.startswith("no_value_") for i in ids)


# --------------------------------------------------------------------------
# inputs other than a notebook
# --------------------------------------------------------------------------
def test_audit_works_from_a_compound_and_a_dose_alone():
    audit = claim_audit(compound="fipronil", conc_M=1e-6)
    assert audit["subject"]["compound"] == "fipronil"
    assert audit["subject"]["concentration_M"] == 1e-6
    assert [link["step"] for link in audit["chain"]] == list(CHAIN_STEPS)
    sources = next(link for link in audit["chain"] if link["step"] == "parameter_source")["sources"]
    gaba = next(s for s in sources if s["receptor"] == "vertebrate_GABA_A")
    assert gaba["param_type"] == "IC50"
    assert gaba["param_value_M"] == pytest.approx(1.1e-6)


def test_a_compound_with_no_sourced_row_is_unknown_not_zero():
    audit = claim_audit(compound="fipronil", conc_M=1e-6)
    sources = next(link for link in audit["chain"] if link["step"] == "parameter_source")["sources"]
    nach = next(s for s in sources if s["receptor"] == "insect_nAChR")
    assert nach["classification"] == "NOT MODELLED"
    assert nach["param_value_M"] is None


def test_audit_never_claims_a_live_measurement(audit):
    text = " ".join(w for w in audit["warnings"])
    assert "live_lab" in text
    assert audit["disclaimer"]
    for link in audit["chain"]:
        assert link["label"] != "OBSERVED" or link["step"] == "malecns_edges"


def test_markdown_render_mentions_every_step(audit):
    md = to_markdown(audit)
    for step in CHAIN_STEPS:
        assert step in md
    assert "## Fact" in md and "## Unknown" in md


def test_empty_input_does_not_raise():
    audit = claim_audit({})
    assert audit["subject"]["compound"] is None
    assert [link["step"] for link in audit["chain"]] == list(CHAIN_STEPS)
    param = next(link for link in audit["chain"] if link["step"] == "parameter_source")
    assert param["label"] == "UNKNOWN"
