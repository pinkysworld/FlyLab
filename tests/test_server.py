"""HTTP contract tests for the bench API.

Every endpoint in ``docs/DESIGN_v0.5.md`` is exercised once.  The
forward-compatible endpoints (null models, selectivity, predictions,
genotypes, mixtures, expression, validation) accept either 200 or 501: the
modules behind them land separately, and the bench must keep serving either
way.  Run times are kept small (short LIF windows, few bootstraps) so the whole
file stays well inside the suite budget.
"""

from __future__ import annotations

import json

import pytest
from fastapi.testclient import TestClient

from flylab.server import MAX_EXPERIMENT_ROWS, VERSION, app

#: endpoints whose backing module may not exist yet
FORWARD_OK = (200, 501)


@pytest.fixture(scope="module")
def client() -> TestClient:
    with TestClient(app) as c:
        yield c


# --------------------------------------------------------------------------
# page + meta
# --------------------------------------------------------------------------
def test_index_serves_the_bench(client):
    r = client.get("/")
    assert r.status_code == 200
    assert "FLY" in r.text and "app.js" in r.text


def test_static_assets_are_mounted(client):
    for path in ("/static/index.html", "/static/app.js", "/static/styles.css"):
        assert client.get(path).status_code == 200, path


def test_health(client):
    body = client.get("/api/health").json()
    assert body == {"ok": True, "version": VERSION}


def test_meta_has_everything_the_ui_needs(client):
    r = client.get("/api/meta")
    assert r.status_code == 200
    m = r.json()
    assert m["version"] == "0.5.0"
    assert m["map"]["id"] and m["map"]["citation"]

    # library identity
    assert len(m["library"]["sha256"]) == 64
    assert m["library"]["n_compounds"] >= 5

    # compounds carry class and tiered insect targets
    keys = {c["key"] for c in m["compounds"]}
    assert {"imidacloprid", "nicotine", "fipronil"} <= keys
    imi = next(c for c in m["compounds"] if c["key"] == "imidacloprid")
    assert imi["class"] == "neonicotinoid"
    assert any(t["receptor"] == "insect_nAChR" for t in imi["insect_targets"])
    assert all("evidence_tier" in t for t in imi["insect_targets"])

    # receptors, mechanism rules and gain keys
    assert m["receptors"]["insect_nAChR"]["organism"] == "insect"
    assert any(row["gain"] == "g_ach" for row in m["mechanisms"])
    assert set(m["gain_keys"]) == {"g_ach", "g_gaba", "g_glu", "g_oct", "g_nav", "ach_tone"}

    # both graphs, with node/edge/seed counts
    for name in ("named", "taste_motor"):
        g = m["graphs"][name]
        assert g["available"] is True
        assert g["n_nodes"] > 0 and g["n_edges"] > 0
        assert g["map"] and g["citation"]
        assert g["seeds"]["MN9"] >= 1
    assert m["graphs"]["taste_motor"]["seeds"]["LB3b"] >= 1

    # gustatory seed hypothesis is carried, with its caveat
    assert "hypothesis" in json.dumps(m["gustatory_hypothesis"]).lower()
    assert m["limits"]["max_experiment_rows"] == MAX_EXPERIMENT_ROWS


def test_census(client):
    body = client.get("/api/census").json()
    assert body["n_traced"] > 1000
    assert body["neurotransmitter_counts"]


# --------------------------------------------------------------------------
# library / occupancy
# --------------------------------------------------------------------------
def test_drugs_list_and_detail(client):
    keys = client.get("/api/drugs").json()["compounds"]
    assert "imidacloprid" in keys
    entry = client.get("/api/drugs/imidacloprid").json()
    assert entry["key"] == "imidacloprid"
    assert entry["receptors"]["insect_nAChR"]["ec50_M"] > 0
    assert client.get("/api/drugs/not-a-drug").status_code == 404


def test_occupancy_and_selectivity(client):
    body = client.get("/api/occupancy", params={"compound": "imidacloprid", "conc_M": 1e-6}).json()
    assert body["compound"] == "Imidacloprid"
    by = {r["receptor"]: r for r in body["receptors"]}
    # the headline claim of the project: high insect, low vertebrate
    assert by["insect_nAChR"]["occupancy"] > 0.9
    assert by["vertebrate_nAChR_a4b2"]["occupancy"] < 0.5
    assert body["selectivity"]["nAChR"]["log10_ec50_ratio_vert_over_insect"] > 0


def test_occupancy_rejects_unknown_compound_and_negative_dose(client):
    assert client.get("/api/occupancy", params={"compound": "nope", "conc_M": 1e-6}).status_code == 404
    assert client.get("/api/occupancy", params={"compound": "nicotine", "conc_M": -1}).status_code == 400


def test_occupancy_curve(client):
    body = client.get("/api/occupancy/curve", params={"compound": "nicotine"}).json()
    points = body["points"]
    assert len(points) > 8
    series = [p["receptors"]["insect_nAChR"] for p in points]
    assert series == sorted(series)  # monotone in concentration
    assert body["receptors"]["insect_nAChR"]["family"]


def test_occupancy_ci(client):
    body = client.get(
        "/api/occupancy/ci", params={"compound": "imidacloprid", "conc_M": 1e-8, "n": 40}
    ).json()
    row = body["receptors"]["insect_nAChR"]
    assert row["p2_5"] <= row["point"] <= row["p97_5"]
    assert client.get("/api/occupancy/ci", params={"compound": "nicotine", "n": 1}).status_code == 400


# --------------------------------------------------------------------------
# assays
# --------------------------------------------------------------------------
def test_assay_taste(client):
    nb = client.post("/api/assay/taste", json={"compound": "imidacloprid", "conc_M": 1e-6}).json()
    assert nb["assay"] == "taste_mn9"
    assert nb["readouts"]["mn9_sugar_hz"] >= 0
    assert nb["warnings"]
    assert nb["live_lab"] is None


def test_assay_taste_dose_response(client):
    body = client.get("/api/assay/taste/dose-response", params={"compound": "imidacloprid"}).json()
    assert len(body["points"]) == 6


def test_assay_cns(client):
    nb = client.post("/api/assay/cns", json={"compound": "imidacloprid", "conc_M": 1e-6}).json()
    assert nb["readouts"]["n_traced"] > 1000
    assert "cns_excitation_index" in nb["readouts"]["treated"]


def test_assay_subgraph_on_both_graphs(client):
    for graph in ("named", "taste_motor"):
        nb = client.post(
            "/api/assay/subgraph", json={"compound": "imidacloprid", "conc_M": 1e-6, "graph": graph}
        ).json()
        assert nb["readouts"]["n_nodes"] > 0
        assert nb["readouts"]["runtime_ms"] >= 0
        assert nb["gains"]["g_ach"] < 1.0  # imidacloprid desensitises at this dose


def test_assay_subgraph_dose_response(client):
    body = client.get("/api/assay/subgraph/dose-response", params={"compound": "imidacloprid"}).json()
    assert len(body["points"]) == 6


def test_assay_spiking(client):
    nb = client.post(
        "/api/assay/spiking", json={"compound": "imidacloprid", "conc_M": 1e-6, "t_ms": 150, "seed": 0}
    ).json()
    r = nb["readouts"]
    assert r["engine"] == "lif"
    assert r["t_ms"] == 150
    assert isinstance(r["raster"], list)
    assert r["psth"]["bin_ms"] == 10.0
    assert "MN9" in r["psth"]["cells"]


def test_assay_spiking_bounds(client):
    assert client.post("/api/assay/spiking", json={"t_ms": 99999}).status_code == 422
    assert client.post("/api/assay/spiking", json={"conc_M": -1}).status_code == 422


def test_assay_taste_map(client):
    nb = client.post(
        "/api/assay/taste-map",
        json={"compound": "imidacloprid", "conc_M": 1e-6, "sugar_hz": 150, "bitter_hz": 0,
              "engine": "rate", "seed": 0},
    ).json()
    r = nb["readouts"]
    assert r["n_sweet_grn"] > 0 and r["n_bitter_grn"] > 0
    assert r["mn9_sugar_hz"] is not None
    assert any("hypothesis" in w for w in nb["warnings"])


def test_assay_ensemble(client):
    nb = client.post(
        "/api/assay/ensemble",
        json={"assay": "subgraph", "compound": "imidacloprid", "conc_M": 1e-6, "n_rep": 3},
    ).json()
    u = nb["uncertainty"]
    assert u["n_rep"] == 3
    assert u["ci"]["mean_hz"][0] <= u["mean"]["mean_hz"] <= u["ci"]["mean_hz"][1]


def test_ensemble_rejects_too_many_replicates(client):
    assert client.post("/api/assay/ensemble", json={"n_rep": 999}).status_code == 422


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------
def test_analysis_sensitivity(client):
    body = client.post(
        "/api/analysis/sensitivity",
        json={"assay": "subgraph", "compound": "imidacloprid", "conc_M": 1e-6, "readout": "mean_hz"},
    ).json()
    assert body["rows"]
    assert {"param", "low", "high", "base", "span"} <= set(body["rows"][0])
    assert body["label"] == "model_derived"


def test_analysis_ic50_is_labelled_model_derived(client):
    body = client.post(
        "/api/analysis/ic50",
        json={
            "assay": "subgraph",
            "compound": "imidacloprid",
            "readout": "mean_hz",
            "concs_M": [1e-9, 1e-7, 1e-5],
            "n_boot": 10,
            "n_rep": 1,
        },
    ).json()
    assert body["label"] == "model_derived"
    assert body["fit"]["ic50"] > 0
    assert len(body["points"]) == 3
    assert any("not an animal IC50" in w for w in body["warnings"])


def test_analysis_ic50_rejects_empty_ladder(client):
    assert client.post("/api/analysis/ic50", json={"concs_M": []}).status_code == 400


# --------------------------------------------------------------------------
# graph viewer
# --------------------------------------------------------------------------
def test_graph_viewer_payload(client):
    body = client.get(
        "/api/graph", params={"graph": "named", "min_weight": 5, "max_edges": 400,
                              "compound": "imidacloprid", "conc_M": 1e-6}
    ).json()
    assert body["graph"] == "named"
    assert body["n_nodes"] > 0 and body["n_edges"] > 0
    node = body["nodes"][0]
    assert {"id", "x", "y", "nt", "rate", "rate_vehicle", "is_seed"} <= set(node)
    edge = body["edges"][0]
    assert {"source", "target", "weight", "eff_weight", "eff_weight_vehicle", "eff_weight_delta"} <= set(edge)
    assert any(n["is_seed"] for n in body["nodes"])
    # the drug must actually move at least one cholinergic edge
    assert any(abs(e["eff_weight_delta"]) > 0 for e in body["edges"])


def test_graph_viewer_rejects_bad_bounds(client):
    assert client.get("/api/graph", params={"max_edges": 0}).status_code == 400
    assert client.get("/api/graph", params={"min_weight": -1}).status_code == 400
    assert client.get("/api/graph", params={"graph": "nope"}).status_code == 422


def test_graph_impact_named(client):
    body = client.post(
        "/api/graph/impact", json={"compound": "imidacloprid", "conc_M": 1e-6, "graph": "named"}
    ).json()
    assert body["top_nodes"] and body["top_edges"]
    assert body["gains"]["g_ach"] < 1.0
    assert "grn_to_mn9_paths" not in body


def test_graph_impact_taste_motor_adds_paths(client):
    body = client.post(
        "/api/graph/impact",
        json={"compound": "imidacloprid", "conc_M": 1e-6, "graph": "taste_motor", "top_paths": 5},
    ).json()
    paths = body["grn_to_mn9_paths"]
    assert "sweet" in paths and "bitter" in paths
    assert "hypothesis" in paths["note"]


# --------------------------------------------------------------------------
# exposure
# --------------------------------------------------------------------------
def test_exposure(client):
    body = client.post(
        "/api/exposure", json={"compound": "imidacloprid", "dose": 1.0, "route": "feeding", "t_h": 12}
    ).json()
    assert body["cmax"] > 0 and body["tmax"] > 0 and body["auc"] > 0
    assert len(body["t_h"]) == len(body["conc_M"])
    assert any("placeholder" in w for w in body["warnings"])


def test_exposure_rejects_unknown_route_and_bad_step(client):
    assert client.post("/api/exposure", json={"route": "injection"}).status_code == 422
    assert client.post("/api/exposure", json={"t_h": 1, "dt_h": 5}).status_code == 400


# --------------------------------------------------------------------------
# experiment
# --------------------------------------------------------------------------
def test_experiment_runs_a_design(client):
    body = client.post(
        "/api/experiment",
        json={
            "assay": "subgraph",
            "compounds": ["imidacloprid", "fipronil"],
            "concs_M": [1e-8, 1e-6],
            "replicates": 1,
            "readouts": ["mn9_hz", "mean_hz", "g_ach"],
            "include_vehicle": False,
        },
    ).json()
    assert body["n_rows"] == 2 * 2 * 1
    assert len(body["rows"]) == 4
    assert "notebooks" not in body  # omitted unless keep_notebooks
    assert body["summary"][0]["n"] >= 1
    assert body["label"] == "model_derived"


def test_experiment_keeps_notebooks_when_asked(client):
    body = client.post(
        "/api/experiment",
        json={"compounds": ["imidacloprid"], "concs_M": [1e-6], "replicates": 1, "keep_notebooks": True},
    ).json()
    assert body["notebooks"]


def test_experiment_row_cap(client):
    r = client.post(
        "/api/experiment",
        json={"compounds": ["imidacloprid"] * 20, "concs_M": [1e-6] * 20, "replicates": 10},
    )
    assert r.status_code == 400
    assert str(MAX_EXPERIMENT_ROWS) in r.json()["detail"]


def test_experiment_csv(client):
    r = client.post(
        "/api/experiment",
        json={"compounds": ["imidacloprid"], "concs_M": [1e-8, 1e-6], "replicates": 1},
        headers={"Accept": "text/csv"},
    )
    assert r.status_code == 200
    r = client.post(
        "/api/experiment/csv",
        json={"compounds": ["imidacloprid"], "concs_M": [1e-8, 1e-6], "replicates": 1},
    )
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/csv")
    assert "attachment" in r.headers["content-disposition"]
    lines = r.text.strip().splitlines()
    assert lines[0].startswith("compound,conc_M,replicate")
    assert len(lines) >= 3


# --------------------------------------------------------------------------
# notebook: live_lab import
# --------------------------------------------------------------------------
def test_live_lab_import_attaches_a_parsed_table(client):
    nb = client.post("/api/assay/taste", json={"compound": "nicotine", "conc_M": 1e-6}).json()
    assert nb["live_lab"] is None
    out = client.post(
        "/api/notebook/live-lab",
        json={
            "notebook": nb,
            "csv_text": "fly_id,conc_M,per_score\n1,1e-6,0.4\n2,1e-6,0.65\n3,1e-6,\n",
            "assay_name": "PER",
            "source": "bench notebook p.12",
        },
    ).json()
    live = out["live_lab"]
    assert live["assay_name"] == "PER"
    assert live["source"] == "bench notebook p.12"
    assert live["n_rows"] == 3
    assert live["columns"] == ["fly_id", "conc_M", "per_score"]
    # numeric where possible, None where the cell was blank
    assert live["table"][0] == {"fly_id": 1, "conc_M": 1e-6, "per_score": 0.4}
    assert live["table"][2]["per_score"] is None
    assert "imported_utc" in live
    assert any("never generates" in w for w in out["warnings"])
    # the endpoint invents nothing: the readouts are untouched
    assert out["readouts"] == nb["readouts"]


def test_live_lab_rejects_empty_and_headerless_input(client):
    nb = {"warnings": []}
    assert client.post("/api/notebook/live-lab", json={"notebook": nb, "csv_text": "  "}).status_code == 400
    assert (
        client.post("/api/notebook/live-lab", json={"notebook": nb, "csv_text": "only,a,header\n"}).status_code
        == 400
    )


# --------------------------------------------------------------------------
# forward-compatible endpoints
# --------------------------------------------------------------------------
def test_analysis_null(client):
    r = client.post(
        "/api/analysis/null",
        json={"compound": "imidacloprid", "conc_M": 1e-6, "mode": "weight_permute", "n": 2},
    )
    assert r.status_code in FORWARD_OK
    if r.status_code == 501:
        assert r.json()["detail"].startswith("module")
    else:
        body = r.json()
        assert "real_effect" in body and "null_effects" in body


def test_analysis_null_panel(client):
    r = client.post(
        "/api/analysis/null-panel",
        json={"compound": "imidacloprid", "n": 2, "modes": ["weight_permute"],
              "include_taste_map": False},
    )
    assert r.status_code in FORWARD_OK
    if r.status_code == 200:
        assert r.json()["rows"]


def test_analysis_selectivity_table(client):
    r = client.get("/api/analysis/selectivity", params={"conc_M": 1e-6})
    assert r.status_code in FORWARD_OK
    assert client.get("/api/analysis/selectivity", params={"conc_M": -1}).status_code == 400


def test_analysis_selectivity_landscape(client):
    r = client.post(
        "/api/analysis/selectivity-landscape",
        json={"assay": "subgraph", "graphs": ["named"], "compounds": ["imidacloprid"]},
    )
    assert r.status_code in FORWARD_OK
    if r.status_code == 200:
        assert r.json()["rows"]


def test_analysis_circuit_si(client):
    r = client.post("/api/analysis/circuit-si", json={"compound": "imidacloprid", "graph": "named"})
    assert r.status_code in FORWARD_OK
    if r.status_code == 200:
        assert "circuit_si_log10" in r.json()


def test_predictions(client):
    r = client.get("/api/predictions", params={"n_rep": 1})
    assert r.status_code in FORWARD_OK
    assert client.get("/api/predictions", params={"n_rep": 0}).status_code == 400


def test_genotypes(client):
    r = client.get("/api/genotypes")
    assert r.status_code in FORWARD_OK


def test_genotype_kwarg_is_forwarded_or_501(client):
    r = client.post(
        "/api/assay/subgraph", json={"compound": "imidacloprid", "conc_M": 1e-6, "genotype": "Rdl-MDRR"}
    )
    assert r.status_code in FORWARD_OK
    if r.status_code == 501:
        assert r.json()["detail"] == "genotype not supported yet"
    # no genotype means no keyword, so the plain path always works
    assert client.post("/api/assay/subgraph", json={"compound": "imidacloprid", "conc_M": 1e-6}).status_code == 200


def test_mixture(client):
    r = client.post(
        "/api/mixture",
        json={"components": [{"compound": "imidacloprid", "conc_M": 1e-7},
                             {"compound": "fipronil", "conc_M": 1e-7}], "model": "bliss"},
    )
    assert r.status_code in FORWARD_OK
    assert client.post("/api/mixture", json={"components": []}).status_code == 400


def test_expression(client):
    assert client.get("/api/expression").status_code in FORWARD_OK


def test_validation(client):
    assert client.get("/api/validation").status_code in FORWARD_OK


def test_at_least_one_forward_endpoint_is_documented_either_way(client):
    """The contract is 200-or-501; a 500 from any of them is a bug."""
    paths = [
        ("GET", "/api/genotypes", None),
        ("GET", "/api/expression", None),
        ("GET", "/api/validation", None),
        ("GET", "/api/analysis/selectivity", None),
    ]
    for method, path, body in paths:
        r = client.get(path) if method == "GET" else client.post(path, json=body)
        assert r.status_code in FORWARD_OK, (path, r.status_code, r.text[:200])


# --------------------------------------------------------------------------
# cross-cutting
# --------------------------------------------------------------------------
def test_cors_allows_a_localhost_origin(client):
    r = client.get("/api/health", headers={"Origin": "http://localhost:5173"})
    assert r.headers.get("access-control-allow-origin") == "http://localhost:5173"


def test_missing_graph_file_is_a_404_not_a_500(client, monkeypatch):
    import flylab.server as server

    def boom(*_a, **_kw):
        raise FileNotFoundError("neighborhood JSON missing")

    monkeypatch.setattr(server, "rate_network", boom)
    r = client.get("/api/graph", params={"graph": "named"})
    assert r.status_code == 404
    assert "missing" in r.json()["detail"]
