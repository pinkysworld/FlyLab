"""The browser bridge must answer exactly what the HTTP API answers.

The static (Pyodide) build runs ``flylab.browser.bridge`` instead of
``flylab.server``; if the two ever disagree, a notebook exported from the
GitHub Pages bench stops being the same object as one exported from the served
bench, and the reproducibility claim in ``docs/PAGES.md`` is void.  So every
representative route is run through both transports in the same process and
compared value by value.

Only transport-level noise is ignored: wall-clock fields (``runtime_ms``) and
the creation timestamps (``created_utc`` / ``imported_utc``), which differ
between two runs of the *same* backend as well.
"""
from __future__ import annotations

import builtins
import importlib
import sys

import pytest
from fastapi.testclient import TestClient

from flylab.browser import bridge
from flylab.server import app

# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
#: keys that legitimately differ between two runs of anything
#: (``runtime_s`` is the analysis layer's own wall-clock field)
VOLATILE = {"created_utc", "imported_utc", "runtime_ms", "runtime_s"}

TOL = 1e-9


def scrub(obj):
    if isinstance(obj, dict):
        return {k: scrub(v) for k, v in obj.items() if k not in VOLATILE}
    if isinstance(obj, list):
        return [scrub(v) for v in obj]
    return obj


def diff(a, b, path=""):
    """Every place ``a`` (bridge) and ``b`` (HTTP) disagree, floats to ``TOL``."""
    out: list[str] = []
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b)):
            if key not in a:
                out.append(f"{path}.{key}: missing from the bridge answer")
            elif key not in b:
                out.append(f"{path}.{key}: missing from the HTTP answer")
            else:
                out += diff(a[key], b[key], f"{path}.{key}")
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: length {len(a)} vs {len(b)}")
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                out += diff(x, y, f"{path}[{i}]")
    elif isinstance(a, bool) or isinstance(b, bool):
        if a is not b:
            out.append(f"{path}: {a!r} vs {b!r}")
    elif isinstance(a, float) or isinstance(b, float):
        if a is None or b is None:
            if a is not b:
                out.append(f"{path}: {a!r} vs {b!r}")
        elif abs(float(a) - float(b)) > TOL * max(1.0, abs(float(b))):
            out.append(f"{path}: {a!r} vs {b!r}")
    elif a != b:
        out.append(f"{path}: {a!r} vs {b!r}")
    return out


@pytest.fixture(scope="module")
def client() -> TestClient:
    with TestClient(app) as c:
        yield c


def http_json(client: TestClient, method: str, route: str, payload):
    if method == "GET":
        r = client.get(route, params=payload or {})
    else:
        r = client.post(route, json=payload or {})
    return r


#: (id, method, route, payload) -- one representative call per contract area.
#: Everything here is sized to keep the whole file inside the suite budget.
CASES = [
    ("meta", "GET", "/api/meta", None),
    ("health", "GET", "/api/health", None),
    ("census", "GET", "/api/census", None),
    ("drugs", "GET", "/api/drugs", None),
    ("drug", "GET", "/api/drugs/imidacloprid", None),
    ("occupancy", "GET", "/api/occupancy", {"compound": "imidacloprid", "conc_M": 1e-6}),
    ("occupancy-curve", "GET", "/api/occupancy/curve", {"compound": "nicotine"}),
    (
        "occupancy-ci",
        "GET",
        "/api/occupancy/ci",
        {"compound": "imidacloprid", "conc_M": 1e-6, "n": 50, "seed": 3},
    ),
    ("assay-taste", "POST", "/api/assay/taste", {"compound": "imidacloprid", "conc_M": 1e-6}),
    ("assay-cns", "POST", "/api/assay/cns", {"compound": "imidacloprid", "conc_M": 1e-6}),
    (
        "assay-subgraph",
        "POST",
        "/api/assay/subgraph",
        {"compound": "imidacloprid", "conc_M": 1e-6, "graph": "named"},
    ),
    (
        "assay-spiking",
        "POST",
        "/api/assay/spiking",
        {"compound": "fipronil", "conc_M": 1e-6, "t_ms": 120, "seed": 1},
    ),
    (
        "assay-taste-map",
        "POST",
        "/api/assay/taste-map",
        {"compound": "imidacloprid", "conc_M": 1e-6, "engine": "rate"},
    ),
    (
        "ensemble",
        "POST",
        "/api/assay/ensemble",
        {"assay": "subgraph", "compound": "nicotine", "conc_M": 1e-6, "n_rep": 3, "seed": 2},
    ),
    (
        "sensitivity",
        "POST",
        "/api/analysis/sensitivity",
        {"assay": "subgraph", "compound": "imidacloprid", "conc_M": 1e-6, "readout": "mean_hz"},
    ),
    (
        "ic50",
        "POST",
        "/api/analysis/ic50",
        {
            "assay": "subgraph",
            "compound": "imidacloprid",
            "readout": "mean_hz",
            "concs_M": [1e-9, 1e-7, 1e-6, 1e-5],
            "n_boot": 0,
            "n_rep": 1,
            "seed": 0,
        },
    ),
    (
        "graph",
        "GET",
        "/api/graph",
        {
            "graph": "named",
            "min_weight": 5,
            "max_edges": 300,
            "compound": "imidacloprid",
            "conc_M": 1e-6,
        },
    ),
    (
        "graph-impact",
        "POST",
        "/api/graph/impact",
        {
            "compound": "imidacloprid",
            "conc_M": 1e-6,
            "graph": "named",
            "top_nodes": 5,
            "top_edges": 5,
            "top_paths": 3,
        },
    ),
    (
        "exposure",
        "POST",
        "/api/exposure",
        {"compound": "imidacloprid", "dose": 1.0, "route": "feeding", "t_h": 6, "dt_h": 0.5},
    ),
    (
        "experiment",
        "POST",
        "/api/experiment",
        {"compounds": ["imidacloprid"], "concs_M": [1e-8, 1e-6], "replicates": 1, "seed": 0},
    ),
    ("genotypes", "GET", "/api/genotypes", None),
    (
        "mixture",
        "POST",
        "/api/mixture",
        {
            "components": [
                {"compound": "imidacloprid", "conc_M": 1e-7},
                {"compound": "fipronil", "conc_M": 1e-7},
            ],
            "assay": "subgraph",
        },
    ),
    ("selectivity", "GET", "/api/analysis/selectivity", {"conc_M": 1e-6}),
    ("predictions", "GET", "/api/predictions", {"n_rep": 1, "seed": 0}),
    ("expression", "GET", "/api/expression", None),
    ("validation", "GET", "/api/validation", None),
    (
        "null",
        "POST",
        "/api/analysis/null",
        {
            "assay": "subgraph",
            "compound": "imidacloprid",
            "conc_M": 1e-6,
            "mode": "sign_permute",
            "n": 5,
            "seed": 0,
        },
    ),
    (
        "live-lab",
        "POST",
        "/api/notebook/live-lab",
        {"notebook": {"warnings": []}, "csv_text": "fly,per\n1,0.4\n2,0.6\n"},
    ),
    # -- the v0.6 dashboard / decision layer.  Everything here is pinned to the
    # cheapest settings the route allows: the point is that both transports run
    # the same function, not that the analysis is well resolved.
    (
        "dashboard",
        "GET",
        "/api/dashboard",
        {"compound": "fipronil", "conc_M": 1e-6, "n": 2, "graph": "named"},
    ),
    (
        "dashboard-estimate",
        "GET",
        "/api/dashboard",
        {"compound": "imidacloprid", "estimate_only": True},
    ),
    (
        "compare",
        "POST",
        "/api/compare",
        {
            "compounds": ["imidacloprid", "fipronil"],
            "conc_M": 1e-6,
            "include_dependence": False,
        },
    ),
    ("claims", "POST", "/api/claims", {"compound": "imidacloprid", "conc_M": 1e-6}),
    (
        "dependence",
        "POST",
        "/api/dependence",
        {"compound": "imidacloprid", "conc_M": 1e-6, "n": 2, "seed": 0},
    ),
    ("dependence-landscape", "POST", "/api/dependence/landscape", {"n": 2}),
    ("ablation", "POST", "/api/ablation", {"compound": "imidacloprid", "conc_M": 1e-6}),
    ("stability-estimate", "POST", "/api/robustness/stability", {"estimate_only": True}),
    ("thresholds-estimate", "POST", "/api/robustness/thresholds", {"estimate_only": True}),
    ("uncertainty-estimate", "POST", "/api/uncertainty/global", {"estimate_only": True}),
    ("voi-estimate", "POST", "/api/voi", {"estimate_only": True}),
    (
        "genotype-panel",
        "POST",
        "/api/genotype/panel",
        {"compound": "deltamethrin", "conc_M": 1e-6},
    ),
    (
        "isobologram",
        "POST",
        "/api/mixture/isobologram",
        {"compound_a": "imidacloprid", "compound_b": "fipronil", "n": 3},
    ),
]


# --------------------------------------------------------------------------
# parity
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name,method,route,payload", CASES, ids=[c[0] for c in CASES])
def test_bridge_matches_http(client, name, method, route, payload):
    response = http_json(client, method, route, payload)
    got = bridge.call(route, payload)
    if response.status_code == 501:
        assert got.get("error", {}).get("status") == 501, f"{name}: bridge should also be 501"
        return
    assert response.status_code == 200, f"{name}: HTTP said {response.status_code}"
    assert "error" not in got or not isinstance(got.get("error"), dict), f"{name}: {got}"
    problems = diff(scrub(got), scrub(response.json()), route)
    assert not problems, f"{name}: {len(problems)} differences\n" + "\n".join(problems[:20])


def test_experiment_csv_text_is_the_same_table(client):
    payload = {"compounds": ["imidacloprid"], "concs_M": [1e-6], "replicates": 1, "seed": 0}
    http = client.post("/api/experiment/csv", json=payload)
    assert http.status_code == 200
    out = bridge.call("/api/experiment/csv", payload)
    assert out["csv"] == http.text
    assert out["content_type"].startswith("text/csv")
    assert out["filename"].startswith("flylab_experiment_")


def test_null_panel_matches_http(client):
    """The slowest route in the contract: one mode, two shuffles."""
    payload = {
        "compound": "imidacloprid",
        "conc_M": 1e-6,
        "n": 2,
        "seed": 0,
        "modes": ["sign_permute"],
        "include_taste_map": False,
    }
    http = client.post("/api/analysis/null-panel", json=payload)
    got = bridge.call("/api/analysis/null-panel", payload)
    if http.status_code == 501:
        assert got["error"]["status"] == 501
        return
    assert http.status_code == 200
    problems = diff(scrub(got), scrub(http.json()), "/api/analysis/null-panel")
    assert not problems, "\n".join(problems[:20])


# --------------------------------------------------------------------------
# contract
# --------------------------------------------------------------------------
def test_every_server_route_is_mirrored():
    """Nothing may be added to server.py without being added to the bridge."""
    served = set()
    for route in app.routes:
        path = getattr(route, "path", "")
        if path.startswith("/api/"):
            served.add(path)
    mirrored = set(bridge.routes())
    assert served <= mirrored, f"not mirrored by the bridge: {sorted(served - mirrored)}"
    assert mirrored <= served, f"bridge invents routes: {sorted(mirrored - served)}"
    assert len(mirrored) == 46


def test_version_and_limits_match_the_server():
    from flylab import server

    assert bridge.VERSION == server.VERSION
    assert bridge.MAX_N_REP == server.MAX_N_REP
    assert bridge.MAX_T_MS == server.MAX_T_MS
    assert bridge.MAX_EXPERIMENT_ROWS == server.MAX_EXPERIMENT_ROWS
    assert bridge.MAX_CONC_M == server.MAX_CONC_M
    v = bridge.version()
    assert v["version"] == server.VERSION and v["n_routes"] == 46


def test_errors_are_returned_not_raised():
    assert bridge.call("/api/nope") == {
        "error": {"status": 404, "detail": "no such route '/api/nope'"}
    }
    assert bridge.call("/api/drugs/unobtainium")["error"]["status"] == 404
    assert bridge.call("/api/occupancy", {"compound": "unobtainium"})["error"]["status"] == 404
    assert bridge.call("/api/occupancy", {"compound": "imidacloprid", "conc_M": -1})["error"][
        "status"
    ] == 400
    assert bridge.call("/api/assay/subgraph", {"steps": 10**9})["error"]["status"] == 400
    bad = bridge.call("/api/experiment", {"compounds": ["a"] * 50, "concs_M": [1e-6] * 50, "replicates": 5})
    assert bad["error"]["status"] == 400 and "limit" in bad["error"]["detail"]
    assert bridge.call("/api/notebook/live-lab", {"notebook": {}, "csv_text": " "})["error"][
        "status"
    ] == 400


def test_query_string_and_payload_are_the_same_request():
    a = bridge.call("/api/occupancy?compound=nicotine&conc_M=1e-6")
    b = bridge.call("/api/occupancy", {"compound": "nicotine", "conc_M": 1e-6})
    assert scrub(a) == scrub(b)


def test_genotype_requests_report_501_like_the_server(client):
    payload = {"compound": "imidacloprid", "conc_M": 1e-6, "genotype": "Rdl-A301S"}
    http = client.post("/api/assay/subgraph", json=payload)
    got = bridge.call("/api/assay/subgraph", payload)
    assert http.status_code == got.get("error", {}).get("status", 200)


# --------------------------------------------------------------------------
# the browser LIF kernel
# --------------------------------------------------------------------------
def test_browser_lif_constants_keep_the_step_and_shorten_the_window():
    from flylab.circuit import lif

    assert lif.BROWSER_T_MS == 200.0
    assert lif.BROWSER_DT_MS == lif.DT_MS == 0.1, "the browser must not integrate more coarsely"
    assert lif.BROWSER_SPARSE is True
    # the native default is still the dense kernel
    assert lif.LIFNetwork.__init__.__defaults__[-1] is False


@pytest.mark.parametrize("graph_name", ["named", "taste_motor"])
def test_sparse_kernel_gives_identical_spike_trains(graph_name):
    """The browser kernel must be the same simulation, not a similar one."""
    import numpy as np

    from flylab.circuit.lif import LIFNetwork
    from flylab.circuit.rate import compute_gains, load_graph

    graph = load_graph(graph_name)
    net = LIFNetwork(graph, seed=0)
    drive = {n["bodyId"]: 26.0 for n in graph["nodes"]}
    for ids in graph["seeds"].values():
        for body_id in ids:
            drive[body_id] = 40.0

    assert net.nnz < len(net) ** 2 / 10, "the graph is meant to be sparse"
    indptr, indices, data = net.csr()
    assert indptr.size == len(net) + 1 and indices.size == data.size == net.nnz
    assert int(indptr[-1]) == net.nnz

    for gains in (None, compute_gains("imidacloprid", 1e-6)[0], compute_gains("fipronil", 1e-6)[0]):
        for seed in (0, 7):
            net.seed = seed
            dense = net.run(drive, gains, t_ms=200.0, sparse=False)
            net.seed = seed
            sparse = net.run(drive, gains, t_ms=200.0, sparse=True)
            assert dense.spikes == sparse.spikes
            assert np.array_equal(dense.rates_hz, sparse.rates_hz)
            assert dense.window_ms == sparse.window_ms


def test_lif_network_helper_honours_the_sparse_flag():
    from flylab.circuit.lif import lif_network
    from flylab.circuit.rate import load_graph

    graph = load_graph("named")
    assert lif_network(graph, seed=1, sparse=True).sparse is True
    assert lif_network(graph, seed=1).sparse is False


# --------------------------------------------------------------------------
# no pydantic, no fastapi
# --------------------------------------------------------------------------
BLOCKED = ("pydantic", "fastapi", "starlette", "typer", "uvicorn", "pandas", "pyarrow")


@pytest.fixture
def without_web_stack(monkeypatch):
    """Make ``import pydantic`` (etc.) fail, as it does inside Pyodide."""
    real_import = builtins.__import__

    def guarded(name, *args, **kw):
        if name.split(".")[0] in BLOCKED:
            raise ImportError(f"{name} is not available in this environment")
        return real_import(name, *args, **kw)

    for mod in list(sys.modules):
        if mod.split(".")[0] in BLOCKED or mod.startswith("flylab"):
            monkeypatch.delitem(sys.modules, mod, raising=False)
    monkeypatch.setattr(builtins, "__import__", guarded)
    yield
    for mod in list(sys.modules):
        if mod.startswith("flylab"):
            sys.modules.pop(mod, None)


def test_experiment_design_falls_back_to_a_dataclass(without_web_stack):
    experiment = importlib.import_module("flylab.assays.experiment")
    assert experiment.HAVE_PYDANTIC is False
    design = experiment.ExperimentDesign(
        compounds=["imidacloprid"], concs_M=[1e-6], replicates=2, seed=1
    )
    dumped = design.model_dump()
    assert dumped["compounds"] == ["imidacloprid"]
    assert dumped["concs_M"] == [1e-6]
    assert dumped["replicates"] == 2 and dumped["seed"] == 1
    assert dumped["readouts"] == list(experiment.READOUT_KEYS)
    assert dumped["assay"] == "subgraph" and dumped["options"] == {}
    # unknown keys are ignored, exactly as pydantic's extra="ignore" does
    assert "genotype" not in experiment.ExperimentDesign(genotype="Rdl").model_dump()

    for kwargs, message in (
        ({"compounds": []}, "design needs at least one compound"),
        ({"concs_M": []}, "design needs at least one concentration"),
        ({"concs_M": [-1e-6]}, "concentrations must be >= 0"),
        ({"replicates": 0}, "replicates must be >= 1"),
        ({"assay": "nope"}, "assay must be one of"),
    ):
        with pytest.raises(ValueError) as exc:
            experiment.ExperimentDesign(**kwargs)
        assert message in str(exc.value)


def test_the_same_design_dumps_the_same_dict_with_and_without_pydantic():
    """The pydantic path is the reference; the fallback must match it exactly."""
    from flylab.assays.experiment import HAVE_PYDANTIC, ExperimentDesign

    assert HAVE_PYDANTIC, "this test documents the pydantic path"
    kwargs = dict(assay="spiking", compounds=["nicotine"], concs_M=[1e-7], replicates=3, seed=5)
    reference = ExperimentDesign(**kwargs).model_dump()

    real_import = builtins.__import__

    def guarded(name, *args, **kw):
        if name.split(".")[0] in BLOCKED:
            raise ImportError(name)
        return real_import(name, *args, **kw)

    saved = {m: sys.modules[m] for m in list(sys.modules) if m.startswith("flylab")}
    try:
        for mod in list(sys.modules):
            if mod.startswith("flylab"):
                del sys.modules[mod]
        builtins.__import__ = guarded
        fallback = importlib.import_module("flylab.assays.experiment")
        assert fallback.HAVE_PYDANTIC is False
        assert fallback.ExperimentDesign(**kwargs).model_dump() == reference
    finally:
        builtins.__import__ = real_import
        for mod in list(sys.modules):
            if mod.startswith("flylab"):
                del sys.modules[mod]
        sys.modules.update(saved)


def test_bridge_runs_the_whole_core_without_pydantic_or_fastapi(without_web_stack):
    """``python -c "import flylab.browser.bridge"`` with the web stack hidden."""
    module = importlib.import_module("flylab.browser.bridge")
    assert "pydantic" not in sys.modules and "fastapi" not in sys.modules

    nb = module.call("/api/assay/subgraph", {"compound": "imidacloprid", "conc_M": 1e-6})
    assert nb["readouts"]["mean_hz"] > 0
    assert nb["provenance"]["library_sha256"]

    out = module.call(
        "/api/experiment", {"compounds": ["imidacloprid"], "concs_M": [1e-6], "replicates": 1}
    )
    assert out["n_rows"] == 1 and out["design"]["assay"] == "subgraph"
    assert module.call("/api/census")["map"].startswith("male-cns")
    assert "pydantic" not in sys.modules


def test_bridge_import_is_cheap():
    """A cold browser tab pays for this import before anything else happens."""
    import subprocess

    code = (
        "import sys, flylab.browser.bridge as b;"
        "assert 'numpy' not in sys.modules, 'numpy must not be imported at module level';"
        "assert len(b.routes()) == 46;"
        "print(len([m for m in sys.modules if m.startswith('flylab')]))"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    assert int(out.stdout.strip()) <= 4, "the bridge pulled in more of flylab than it needs"


def test_browser_defaults_to_the_fast_dependence_setting():
    """The static build must not silently run the paper's shuffle count."""
    from flylab.analysis.dependence import FAST_N

    assert bridge.FAST_DEPENDENCE_N == FAST_N
    out = bridge.call("/api/dashboard", {"compound": "imidacloprid", "estimate_only": True})
    assert out["runtime_estimate"]["n"] == FAST_N
    assert out["runtime_estimate"]["estimate_s"] > 0


def test_long_routes_answer_an_estimate_without_running():
    """Every expensive route states a runtime before it spends one."""
    for route in (
        "/api/robustness/stability",
        "/api/robustness/thresholds",
        "/api/uncertainty/global",
        "/api/voi",
        "/api/ablation",
    ):
        out = bridge.call(route, {"estimate_only": True})
        assert out["estimate_only"] is True, route
        assert out["runtime_estimate"]["estimate_s"] > 0, route
    # the landscape is minutes of compute, so estimate-only is its *default*
    out = bridge.call("/api/dependence/landscape", {})
    assert out["estimate_only"] is True
    assert out["runtime_estimate"]["estimate_s"] > 0


def test_the_dashboard_builder_is_shared_with_the_server():
    """Not merely equal answers: literally the same function object."""
    from flylab import server

    assert server.bridge is bridge
    assert callable(bridge.build_dashboard)
    assert callable(bridge.build_compare)
    assert callable(bridge.build_claims)
