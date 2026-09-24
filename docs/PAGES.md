# Browser bench and GitHub Pages

FlyLab has two interactive transports around the same Python science core:

1. the local FastAPI bench,
2. the static browser bench running the FlyLab wheel inside Pyodide.

The browser build is not a JavaScript reimplementation of the pharmacology or circuit models.

## Architecture

~~~text
                         shared FlyLab Python package
                 pharm | circuit | assays | analysis | maps
                                /                 \
                               /                   \
                  local server                    browser
                  FastAPI HTTP               Pyodide bridge
                       |                           |
                       +----------- app.js --------+
~~~

The local application reaches scientific functions through flylab/server.py.

The static application reaches the same core through flylab/browser/bridge.py.

The user interface is shared.

The bench stores the primary run and experiment-design settings in the current
browser. **Export setup** writes a versioned JSON preset that can be imported
in another browser build; it contains settings, not assay results. At phone
widths, the FlyLab and Bench tab rows become a single grouped native panel
picker, synchronized with the active panel.

## What "same core" means

The browser and server builds share:

- pharmacology library parsing,
- typed receptor engagement,
- mechanism rules,
- derived MaleCNS circuit data,
- rate-model assays,
- supported LIF assays,
- analysis functions exposed by the bridge,
- notebook and provenance structures.

There is no second set of receptor rules or circuit equations in JavaScript.

## What it does not mean

The browser build should not be described as reproducing the entire publication pipeline.

The native paper pipeline also uses publication tooling such as Matplotlib and writes figures, tables, results.json and rendered manuscript artifacts.

The accurate reproducibility claim is:

> The browser executes the same FlyLab scientific core for supported interactive endpoints, while the native reproduction pipeline regenerates the complete publication artifact set.

## Browser dependencies

The browser build intentionally keeps its Python dependency surface small.

The science core used by Pyodide must not require server-only packages at module import time.

Server and CLI dependencies such as FastAPI, pydantic, typer, pandas and pyarrow are therefore isolated from the browser-safe core wherever possible.

The browser package subset is declared in pyproject.toml.

## Bridge

The bridge accepts route-style calls and dispatches them to Python handlers.

Example:

~~~python
from flylab.browser import bridge

bridge.call("/api/meta")
bridge.call(
    "/api/assay/subgraph",
    {"compound": "imidacloprid", "conc_M": 1e-6},
)
~~~

The exact route set is defined by the current implementation. Do not copy a route count into durable documentation because the count changes as the API grows.

## Transport parity

tests/test_browser_bridge.py compares representative browser-bridge and server calls, and separately asserts that the bridge mirrors the server's route set **exactly in both directions**: every /api/ route the server serves is mirrored by the bridge, and the bridge invents none the server does not serve. The per-call comparison is representative; the route-set contract is not.

The intended contract is:

- same scientific inputs,
- same scientific implementation,
- equal or numerically equivalent supported outputs,
- equivalent warnings and provenance semantics.

This is transport parity, not evidence that every possible route, parameter combination or publication artifact has been exhaustively cross-validated.

## Provenance differences

A browser has no Git repository checkout in the usual sense.

As a result:

- notebook git_sha may be null in Pyodide,
- the static-site build records the build revision in its manifest,
- library and data hashes remain available where supported.

These transport differences should be explicit rather than hidden.

## LIF execution

The browser may use browser-oriented execution defaults to keep interactive latency acceptable. Concretely, the static build defaults the spiking window to 200 ms against 500 ms natively, with the same 0.1 ms integration step: a shorter run of the same model, not a coarser one.

Any such difference must be visible in the UI and documentation. A shorter simulation window is not the same experiment and must never be silently presented as identical to a longer native run.

The underlying model equations and parameters should remain shared unless a difference is explicitly documented.

## Building the static site

~~~bash
python scripts/build_pages.py
~~~

Preview locally:

~~~bash
python scripts/build_pages.py --serve
~~~

Build with a vendored Pyodide runtime:

~~~bash
python scripts/build_pages.py --vendor-pyodide
~~~

The generated site is written under dist/pages and is not committed as source.

## GitHub Pages deployment

The Pages workflow is defined in:

~~~text
.github/workflows/pages.yml
~~~

It builds and publishes from main and can also be started manually.

Pull requests do not publish directly to the live site.

The deployment workflow intentionally separates:

- build,
- artifact upload,
- Pages deployment.

## Local browser bench

For local development, the simplest supported path is:

~~~bash
python -m pip install -e ".[dev,viz]"
flylab serve
~~~

This runs the FastAPI transport and the shared UI at http://127.0.0.1:8765. It does not open a browser: pass --open if you want one launched, or use --host and --port to change the address.

## Zero-install bench

The public static bench is intended to be available at:

https://pinkysworld.github.io/FlyLab/

The GitHub Pages workflow is the source of truth for deployment status.

## Testing

Fast test suite:

~~~bash
python -m pytest -q
~~~

Longer end-to-end browser and static-build checks (the browser smoke test
requires Playwright and Chromium):

~~~bash
python -m pip install -e ".[dev,viz,browser-test]"
python -m playwright install chromium
pytest -q -m slow
~~~

The scheduled and manually dispatched slow CI job installs Chromium and runs
the static browser smoke test. It builds the Pages site, boots Pyodide in a
headless browser, runs the dashboard through the bridge, switches panels, and
blocks the chart scripts to check that their data tables remain available.
Pyodide and its packages still load from the configured runtime URLs in this
test; it does not verify a fully offline browser session.

## Performance

Browser performance depends strongly on:

- Pyodide startup,
- browser engine,
- cache state,
- CPU,
- selected assay,
- simulation length,
- graph size,
- whether dependencies are vendored.

Do not treat old timing numbers in issue threads or historical commits as stable performance specifications.

If reproducible performance numbers are needed for a paper or release, record:

- FlyLab commit,
- browser version,
- operating system,
- CPU,
- cold or warm cache,
- graph,
- assay parameters,
- number of repetitions.

## Offline and restricted-network use

The chart libraries are vendored by default; --no-vendor-js opts out. The Pages workflow builds with --vendor-pyodide, so the scientific runtime does not depend on a third-party Python CDN at execution time either. The generated manifest.json records both decisions as js_vendored and pyodide_vendored, so a published build states whether it is self-contained rather than leaving it to be inferred.

The published site is therefore self-contained as built. That self-containment is currently verified by inspecting the build and its manifest, not by an automated browser test, so a stronger claim than this one needs such a test first: a check that a loaded page issues no required off-origin request.

## Reproducibility wording for papers

Recommended wording:

> FlyLab exposes the same Python scientific core through a local server and a Pyodide browser bridge. Transport-parity tests compare representative route families and verify equivalent scientific outputs. The complete publication artifact set is regenerated by the native reproduction pipeline from committed inputs.

Avoid:

> The complete paper regenerates in the browser.

Avoid:

> Every browser and native route is proven identical.

Those statements are stronger than the current automated test contract.

## Related documentation

- [Documentation index](README.md)
- [Architecture](ARCHITECTURE.md)
- [Software stack](STACK.md)
- [Research positioning](NOVELTY.md)
- [Main README](../README.md)
