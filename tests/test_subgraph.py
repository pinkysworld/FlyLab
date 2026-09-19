import json
from pathlib import Path
from flylab.assays.subgraph import run_subgraph_assay

GRAPH = Path("data/derived/malecns_named_neighborhood.json")

def test_graph_has_seeds():
    assert GRAPH.exists()
    g = json.loads(GRAPH.read_text())
    assert g["n_nodes"] > 10
    assert g["seeds"]["MN9"]
    assert g["seeds"]["DNp01"]
    assert g["n_edges"] > 10

def test_subgraph_assay_runs():
    nb = run_subgraph_assay("imidacloprid", 1e-6, GRAPH)
    assert nb["readouts"]["named"]["MN9"][0]["hz"] is not None
    assert nb["map"]["name"] == "male-cns:v1.0"
