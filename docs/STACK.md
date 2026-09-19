# Software stack

Goal: one lab that feels like software on Windows and macOS (Linux comes free), not a notebook you have to babysit.

## Decision

| Layer | Choice | Why |
|---|---|---|
| Science core | Python 3.11+ | Connectome loaders, occupancy, pytest, every fly sim already lives here |
| Circuit kernel | numpy LIF first, Brian2 optional | flypoke / Shiu models are Python; do not rewrite spikes in JS |
| API | FastAPI | local server, OpenAPI, easy to hang a UI or a script on |
| UI | Vite + React + TypeScript | comfortable research chrome: two panes, sliders, tables, export |
| Plots | Plotly.js (UI) / matplotlib (tests) | dose–response and occupancy traces |
| Circuit view v1 | SVG / 2D canvas of named cells | 166k WebGL dots are a later vanity feature |
| Packaging | `flylab serve` opens the browser; optional Tauri shell later | no Electron |
| Data | user-downloaded Feather in `data/raw/` | CC-BY maps stay off git |

Ship Gate 3 as **localhost in the default browser**. That is already cross-platform. Wrap in [Tauri](https://tauri.app/) only when someone asks for a Dock / Start-menu icon.

## What researchers actually click

```
+--------------------------------------------------------------+
| FlyLab          compound [ imidacloprid v ]  conc [ 1e-6 ]   |
|                 [ Run assay ]  [ Export notebook ]           |
+------------------------------+-------------------------------+
| FLY CIRCUIT                  | VERTEBRATE PANEL              |
| named cells + rates          | occupancy bars + Hill curves  |
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

## Local run (target)

```bash
python -m pip install -e ".[dev,ui]"
flylab serve          # http://127.0.0.1:8765
```

Same command on Windows PowerShell and macOS Terminal. Python from python.org or conda-forge, not the Windows Store stub.

## Optional later

- Tauri 2: `FlyLab.app` / `FlyLab.exe` wrapping localhost
- neuPrint only as a browser link for cell lookup, not as a runtime dependency
- WebGL instancing of somata after Gate 4
- Linux: already works if the two supported platforms work
