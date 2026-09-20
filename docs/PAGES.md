# The static bench (GitHub Pages)

FlyLab runs in two places and is **one program**. The served bench is
`flylab/server.py` (FastAPI) on your machine. The static bench is the same
`flylab` wheel running inside your browser tab on Pyodide, answering the same
routes through `flylab/browser/bridge.py`. There is no second implementation of
the pharmacology or the circuit models: `flylab-boot.js` is a loader, `app.js`
is the same UI file in both builds, and every number on the page comes out of
Python.

```
                     ┌───────────────────────────── same package ─────────┐
browser tab          │ flylab.pharm / circuit / assays / analysis / maps   │
  app.js  ──►  window.flylabCall  ──►  flylab.browser.bridge.call(route, payload)
                                                      ▲
localhost      app.js  ──►  fetch("/api/...")  ──►  flylab.server (FastAPI)
```

## What runs where

| | served bench | static bench |
|---|---|---|
| transport | FastAPI + pydantic over HTTP | `bridge.call(route, payload)` in Pyodide |
| Python | your interpreter | Pyodide 0.28.3 (CPython 3.13, WebAssembly) |
| packages | full dependency list | `numpy`, `pyyaml` only (`pip install flylab[browser]`) |
| data files | repo `data/` | fetched at boot into `<site-packages>/data/` |
| compute | your CPU, threads available | one WebAssembly thread in the tab |
| your data | stays local | stays local — nothing is uploaded, there is no server |

`fastapi`, `pydantic`, `typer`, `uvicorn`, `pandas` and `pyarrow` are **not**
installed in the browser. `pyarrow` has no WebAssembly build at all, which is
why the wheel is installed with `deps: false` and why the science core had to
stop importing pydantic:

* `flylab/assays/experiment.py` was the only core module that imported
  pydantic. It now imports it inside a `try`, and falls back to a dataclass
  `ExperimentDesign` with the same field names, defaults, coercion, error
  messages and `.model_dump()`. `tests/test_browser_bridge.py` asserts the two
  paths dump identical dicts.
* `flylab/maps/malecns.py` imports `pandas` only *inside* `load_census`, and
  only when a local MaleCNS atlas is present. In the browser there is none, so
  the committed `data/derived/malecns_census_v1.json` fallback is used — the
  same fallback CI uses.
* `flylab/notebook/schema.py` shells out to `git rev-parse` for provenance.
  That call fails inside Pyodide and is already wrapped in `try/except`, so a
  browser notebook carries `git_sha: null` and every other provenance field
  (flylab version, `library.yaml` SHA-256, platform, RNG seed) intact.

Nothing else in `flylab/` (outside `server.py` and `cli.py`) imports a web
framework. `python -c "import flylab.browser.bridge"` works with pydantic,
fastapi, typer, pandas and pyarrow all hidden.

## The bridge

```python
from flylab.browser import bridge

bridge.call("/api/meta")
bridge.call("/api/assay/subgraph", {"compound": "imidacloprid", "conc_M": 1e-6})
bridge.call("/api/occupancy?compound=nicotine&conc_M=1e-6")   # query string also works
bridge.routes()    # the 34 server paths
bridge.version()   # version, limits, data root
```

It mirrors all 34 routes of `flylab/server.py` — same handler bodies, same
defaults, same bounds (`MAX_N_REP`, `MAX_T_MS`, `MAX_EXPERIMENT_ROWS`,
`MAX_CONC_M`), same forward-compatible 501s. Two deliberate transport-level
differences:

* errors are returned, not raised: `{"error": {"status": 400|404|501, "detail": ...}}`.
  `KeyError` / `FileNotFoundError` → 404, `ValueError` and out-of-range fields
  → 400 (FastAPI answers 422 for a malformed body; the bridge has no request
  objects, so it calls that 400 too).
* `POST /api/experiment/csv` returns `{"csv": ..., "filename": ..., "content_type": ...}`
  instead of a file download.

**Data files.** Every module except `flylab/maps/malecns.py` looks for its data
relative to the CWD as well as to the package, so a wheel layout works as long
as the files sit where `Path(__file__).parents[2] / "data"` points — i.e.
`<site-packages>/data/…`, which is exactly where `flylab-boot.js` writes them.
`malecns.py` computes `CENSUS_FALLBACK` / `GUSTATORY_SEEDS` as absolute paths at
import time with no fallback; if those paths are wrong for the layout in use,
`bridge._prepare()` repoints those two module attributes (and, if needed,
chdirs) rather than editing the module. `FLYLAB_DATA_DIR` overrides the search
and names the directory *containing* `data/`.

## LIF in WebAssembly

The browser must not silently run a *different* model, so the integration step
stays at `DT_MS = 0.1 ms`. What the browser build shortens is the window:
`flylab.circuit.lif.BROWSER_T_MS = 200.0` ms (the Spikes panel defaults to it
and says so in one line under the heading), and the rate engine stays the
default everywhere else.

`lif.py` also grew a CSR-style sparse kernel, built once from the same signed
matrix the dense kernel uses (`LIFNetwork.csr()`, `run(..., sparse=True)`,
`BROWSER_SPARSE`). Both kernels reduce each column in ascending row order, so
**the spike trains are identical, not merely similar** — `tests/test_browser_bridge.py`
asserts equality of `spikes` and `rates_hz` for both graphs, two seeds, vehicle
and two drugs, and the existing `tests/test_circuit_lif.py` regressions pin the
dense default.
Native default behaviour is unchanged (`sparse=False`).

Measured here (Python 3.11, single core, 500 ms window, `dt = 0.1 ms`):

| graph | nodes | non-zeros | fill | dense | sparse | speed-up |
|---|---|---|---|---|---|---|
| `named` | 1126 | 1230 | 0.10 % | 384 ms | 352 ms | 1.09× |
| `named`, fipronil | 1126 | 1230 | 0.10 % | 412 ms | 359 ms | 1.15× |
| `taste_motor` | 1841 | 18 345 | 0.54 % | 630 ms | 492 ms | 1.28× |
| `taste_motor`, fipronil | 1841 | 18 345 | 0.54 % | 953 ms | 582 ms | 1.64× |

The win grows with the spike count, because the sparse kernel only ever touches
the entries of the columns that fired; most of the remaining cost is the 5000
Python-level time steps, which neither kernel avoids.

## Building the site

```bash
python scripts/build_pages.py            # -> dist/pages (gitignored)
python scripts/build_pages.py --serve    # build, then preview on :8000
python scripts/build_pages.py --vendor-pyodide   # no CDN at run time
```

`dist/pages/` contains `index.html` (the repo's own, with `/static/...` asset
paths made relative and `flylab-boot.js` injected before `app.js`), `app.js` and
`styles.css` copied byte for byte, the wheel built by `pip wheel . --no-deps`,
the data files, and `manifest.json` with the git SHA, the wheel and
`library.yaml` SHA-256s, every file size and the pinned Pyodide version. The
build prints a size report and **fails if the total exceeds 25 MB**.

Size of this build:

| item | size |
|---|---|
| `data/derived/malecns_taste_motor_neighborhood.json` | 1.03 MB |
| `wheels/flylab-0.5.0-py3-none-any.whl` | 191 KB |
| `data/derived/malecns_named_neighborhood.json` | 169 KB |
| `app.js` | 90 KB |
| `index.html` | 27 KB |
| `data/derived/malecns_census_v1.json` | 26 KB |
| six `data/literature/*.yaml` | 93 KB |
| `styles.css` | 15 KB |
| `flylab-boot.js` | 9 KB |
| **total (16 files)** | **1.64 MB** |

Pyodide itself is loaded from the pinned CDN
`https://cdn.jsdelivr.net/pyodide/v0.28.3/full/` and cached by the browser, so
it is not part of that 1.64 MB: `pyodide.asm.wasm` 8.25 MB, `python_stdlib.zip`
2.30 MB, `pyodide.asm.js` 1.02 MB, the numpy wheel 2.97 MB, pyyaml 117 KB,
micropip + packaging 182 KB.

`--vendor-pyodide` copies exactly those files into the build — the loader, the
interpreter and only the wheels this bench imports, not the other ~340 packages
of a Pyodide release — for networks that block the CDN. That build is **16.61 MB
in 26 files**, still inside the 25 MB budget, and needs no third-party host at
run time. Every network step in `flylab-boot.js` is retried three times and each
package is verified by importing it, so a dropped transfer costs seconds rather
than the page; if it fails anyway, the page says plainly that the runtime could
not start and points the visitor at the served bench.

## Measured in a real browser

Headless Chromium, `dist/pages/` over `python -m http.server`, one core, cold
cache, all CDN traffic through a TLS-intercepting proxy — treat these as an
upper bound, not a benchmark:

| step | browser (Pyodide 0.28.3) | same call natively | note |
|---|---|---|---|
| first load, runtime vendored in the build | **8.6 s** (3.5 s in `flylab-boot.js`) | — | no third-party host; zero console errors, zero failed requests |
| first load from the CDN, cold cache | **37.6 s** (31.8 s in `flylab-boot.js`) | — | includes one automatic retry of `python_stdlib.zip`; an uninterrupted CDN boot measured **22.9 s** |
| `/api/assay/subgraph` (imidacloprid 1e-6, `named`) | **278 ms** | 36 ms | ~8x, and the returned `mean_hz` is `0.45846211979249263` in both, to the last digit |
| `/api/assay/spiking`, 200 ms window | **736 ms** | 295 ms | ~2.5x; `mn9_hz = 12.5 Hz`, 826 vehicle / 563 treated spikes in both |
| Spikes panel (two 200 ms LIF runs + plots) | **1.3 s** | — | |
| Scorecard tab | 1.0 s | | |
| Curves tab (occupancy curve + circuit IC50 ladder) | 15.3 s | | the heaviest routine panel |
| Circuit tab (viewer + impact) | 2.6 s | | |
| Taste tab (reduced control + map path) | 3.0 s | | |
| Notebook tab | 0.9 s | | |

Zero console errors, zero page errors and zero failed requests in the vendored
run; Plotly and cytoscape load from their own CDNs exactly as in the served
bench. The panel timings above are the same either way — once Python is up, the
runtime's origin makes no difference.

The environment these were measured in routes all HTTPS through an intercepting
proxy that intermittently drops a large transfer, which is why `flylab-boot.js`
retries every network step three times and proves each package by importing it
before continuing: a dropped `python_stdlib.zip`, `micropip` or `pyyaml` wheel
then costs a few seconds instead of the whole page. On a network that drops
those repeatedly, ship the vendored build.


## Publishing

`.github/workflows/pages.yml` builds on push to `main` and on
`workflow_dispatch` (**not** on pull requests — a PR must never publish), then
deploys with `actions/deploy-pages` using `pages: write` + `id-token: write`.

**One-time repository setting:** the owner must go to *Settings → Pages →
Build and deployment → Source* and choose **GitHub Actions**. Until that is
done the workflow builds and uploads the artifact but the deploy step fails.

## The reproducibility claim (for the paper)

> The bench in the browser and the bench on a workstation are the same program.
> Both import the same `flylab` wheel; the browser reaches it through
> `flylab.browser.bridge`, the workstation through `flylab.server`, and
> `tests/test_browser_bridge.py` asserts that the two transports return equal
> JSON — every readout to 1e-9, every warning, every provenance field — for a
> representative call to each route family, including all five assays, the
> ensemble, sensitivity, circuit IC50, the graph viewer, mixtures, genotypes,
> selectivity, the pre-registered predictions and a null-model panel. A
> notebook exported from the GitHub Pages page therefore carries the same
> `library.yaml` SHA-256 and the same numbers as one exported locally, and a
> reader can re-run any figure in this paper without installing anything.

Two caveats that belong next to that claim: `git_sha` is `null` in browser
notebooks (there is no git in WebAssembly — the build's SHA is in
`manifest.json` instead), and the browser's Spikes panel defaults to a 200 ms
window rather than 500 ms. The step size, the parameters, the map, the library
and the patch rules are identical.
