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

tests/test_browser_bridge.py compares representative browser-bridge and server calls.

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

The browser may use browser-oriented execution defaults to keep interactive latency acceptable.

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

This runs the FastAPI transport and shared UI in the default browser.

## Zero-install bench

The public static bench is intended to be available at:

https://pinkysworld.github.io/FlyLab/

The GitHub Pages workflow is the source of truth for deployment status.

## Testing

Fast test suite:

~~~bash
python -m pytest -q
~~~

Longer end-to-end browser and static-build checks:

~~~bash
pytest -m slow
~~~

The ordinary tests workflow currently runs the fast suite. Slow checks should be run before tagged research releases and can also be placed in a scheduled or release-specific CI job.

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

The build system can vendor Pyodide assets so the scientific runtime does not depend on a third-party Python CDN at execution time.

Front-end asset policy should be checked against the generated build before making a "fully offline" claim. A release should only use that wording after an automated browser test confirms there are no required off-origin requests.

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
