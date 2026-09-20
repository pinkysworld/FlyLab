# Software stack

Goal: one lab that feels like software on Windows and macOS (Linux comes free), not a notebook you have to babysit.

## Decision

| Layer | Choice | Why |
|---|---|---|
| Science core | Python 3.11+ | Connectome loaders, typed engagement, permutation analyses, pytest; every fly sim already lives here |
| Numerics | numpy only (no scipy) | Spearman, Kendall, Nelder–Mead and the bootstrap are ~200 lines each; scipy has no WebAssembly guarantee and the browser build has to import the same package |
| Circuit kernel | numpy rate model + numpy LIF | the rate model is deterministic and ~50 ms per run, which is what makes a 1000-shuffle permutation profile and an 84-cell FDR-corrected landscape on two cuts affordable; Brian2 stays optional and unused. The engine row-normalises its weight matrix by default (`normalise="row_abs"`, with `none` and `degree` selectable) and that choice is load-bearing for the composition-versus-topology question, so it is swept rather than assumed |
| API | FastAPI + pydantic | local server, OpenAPI, one place to hang the UI and the browser bridge; the browser bridge mirrors every route |
| Browser | Pyodide + `flylab/browser/bridge.py` | the *same wheel* answers the same routes in a tab; no second implementation (`docs/PAGES.md`) |
| UI | plain HTML + `app.js`, CDN libs, no build step | one file serves both the local server and the static site |
| Plots | Plotly.js + cytoscape.js (UI) / matplotlib (paper figures) | the paper needs 300 dpi PNG+SVG from the same numbers the UI shows |
| Packaging | `flylab serve` opens the browser; GitHub Pages for the zero-install bench | no Electron, no Tauri yet |
| Data | user-downloaded Feather; only small derived JSON in git | CC-BY maps stay off git; the 1.1 GB matrix is never needed to reproduce the paper |
| Reproduction | `scripts/reproduce_paper.py` (argparse, named steps, `--fast`, `--only`, `--outdir`) | one command regenerates every figure, table and number; `papers/results.json` is the machine-readable record, and a test asserts the committed record was produced at the shipped statistical effort rather than merely being self-consistent |
| Declarative runs | `flylab run` over a YAML spec (`flylab/spec.py`) | a whole run as one file: resolved spec, notebooks, analyses, figures, claim cards and a manifest with two hashes per artifact; the same spec and seed reproduce it byte for byte (`docs/WORKFLOW.md`) |
| CI | GitHub Actions: `tests.yml`, `pages.yml`, `malecns-subgraph.yml` | tests on push, static bench on `main`, the 1.1 GB cut on demand |

Ship Gate 3 as **localhost in the default browser**. That is already cross-platform, and since v0.5 there is a second, zero-install route: the same wheel on Pyodide, published from GitHub Pages. Wrap in [Tauri](https://tauri.app/) only when someone asks for a Dock / Start-menu icon.

### Two constraints the browser build imposed on the core

Worth knowing before adding a dependency: `pyarrow` has no WebAssembly build and `pydantic`, `fastapi`, `typer` and `pandas` are not installed in the tab. The science core therefore may not import any of them at module scope. `flylab/assays/experiment.py` imports pydantic inside a `try` and falls back to a dataclass with identical behaviour; `flylab/maps/malecns.py` imports pandas only inside `load_census`. Keep it that way: `python -c "import flylab.browser.bridge"` must work with all five hidden.

## What researchers actually click

```
+--------------------------------------------------------------+
| FlyLab          compound [ imidacloprid v ]  conc [ 1e-6 ]   |
|                 [ Run assay ]  [ Export notebook ]           |
+------------------------------+-------------------------------+
| FLY CIRCUIT                  | VERTEBRATE PANEL              |
| named cells + rates          | engagement bars + Hill curves |
| MN9 / GRNs / DNp01           | nAChR a4b2 / GABA-A           |
+------------------------------+-------------------------------+
| notebook preview (JSON) + warnings                           |
+--------------------------------------------------------------+
```

If the window needs a tutorial, the lab is not ready.

## Do not use for the product UI

| Tool | Why not as the lab face |
|---|---|
| Streamlit / Gradio | Fine for a weekend prototype. Rerun-the-script model feels cheap for a dual-pane bench and fights custom layout |
| Jupyter as the only UI | Great for us. Terrible as the thing you hand a visiting student |
| Electron | Two Chromiums plus a 166k-neuron runtime |
| Unity / Unreal / PharmVR-style 3D | Different product (procedure training), months before the first EC50 |
| MATLAB App Designer | Not a shippable Win/Mac lab without licences |
| Qt / PySide only | Native and solid, but you then maintain two front ends the moment you want remote or a paper supplement |
| Pure JS spike engine | You will re-implement Brian2 badly |

Streamlit is allowed as a *throwaway* Gate-3 spike if FastAPI+React slips. It must not become the architecture.

## Local run

```bash
python -m pip install -e ".[dev,viz]"
python -m pytest -q                  # fast suite; `-m slow` adds the end-to-end checks
flylab serve                         # http://127.0.0.1:8765
python scripts/reproduce_paper.py    # ~35-45 min: 16 figures, T0-T25, results.json
```

Same command on Windows PowerShell and macOS Terminal. Python from python.org or conda-forge, not the Windows Store stub.

## Extras

`pip install -e ".[dev]"` for pytest and httpx, `".[viz]"` for matplotlib (needed by `scripts/reproduce_paper.py`), `".[browser]"` for the Pyodide-safe subset.

BLAS threading: the matrices here are small (1126 and 1841 nodes) and an unpinned thread pool costs roughly 10× on a contended machine. `scripts/reproduce_paper.py` pins the thread-count environment variables before importing numpy and re-execs once if numpy is already loaded. If you write another long-running driver, do the same.

## Optional later

- Tauri 2: `FlyLab.app` / `FlyLab.exe` wrapping localhost
- neuPrint only as a browser link for cell lookup, not as a runtime dependency
- WebGL instancing of somata
- Linux: already works if the two supported platforms work
