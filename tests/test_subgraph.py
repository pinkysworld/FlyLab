import json
from pathlib import Path
import pytest
GRAPH = Path("data/derived/malecns_named_neighborhood.json")
@pytest.mark.skipif(not GRAPH.exists(), reason="neighborhood JSON not built")
def test_graph_has_seeds():
    g = json.loads(GRAPH.read_text())
    assert g["n_nodes"] > 10
    assert g["seeds"]["MN9"]
    assert g["seeds"]["DNp01"]
