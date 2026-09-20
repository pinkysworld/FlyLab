#!/usr/bin/env python
"""Assemble the static (GitHub Pages) build of the FlyLab bench.

The static build is the *same* bench: the same ``flylab`` package, built as a
wheel and installed into Pyodide, answering the same routes through
``flylab.browser.bridge`` instead of FastAPI.  No science is reimplemented in
JavaScript -- ``dist/pages/`` is a loader plus the UI plus the wheel plus the
data files.

    python scripts/build_pages.py                 # build into dist/pages
    python scripts/build_pages.py --serve         # build, then preview on :8000
    python scripts/build_pages.py --skip-wheel    # reuse the wheel already there

Layout produced::

    dist/pages/
      index.html          flylab/static/index.html, asset paths made relative
      app.js styles.css   copied verbatim
      flylab-boot.js      Pyodide loader, exposes window.flylabCall
      wheels/flylab-*.whl the package under test
      data/derived/*.json data/literature/*.yaml
      manifest.json       git sha, library sha256, sizes, Pyodide version

See ``docs/PAGES.md``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
STATIC = REPO / "flylab" / "static"
DEFAULT_OUT = REPO / "dist" / "pages"

#: chart libraries the UI loads from a CDN when served by FastAPI.  The
#: published site vendors them instead: a visitor behind a proxy that blocks or
#: throttles third-party CDNs would otherwise get a bench whose Python works
#: but whose plots and circuit graph never render.
JS_VENDOR = {
    "plotly.min.js": "https://cdn.plot.ly/plotly-2.35.2.min.js",
    "cytoscape.min.js": "https://cdn.jsdelivr.net/npm/cytoscape@3.30.2/dist/cytoscape.min.js",
}

#: pinned Pyodide release (jsDelivr mirrors the official pyodide CDN paths)
PYODIDE_VERSION = "0.28.3"
PYODIDE_CDN = "https://cdn.jsdelivr.net/pyodide/v{version}/full/"

#: packages the wheel needs at runtime, resolved by micropip from Pyodide's
#: own lock file (flylab's other dependencies -- fastapi, pydantic, typer,
#: pandas, pyarrow -- are server/CLI-only and are never installed here)
RUNTIME_PACKAGES = ["numpy", "pyyaml"]

#: distribution name -> module the boot script imports to prove the install
#: really landed. micropip logs a failed wheel fetch and carries on, so
#: without this check a flaky CDN would only surface on the first assay.
RUNTIME_IMPORTS = {"numpy": "numpy", "pyyaml": "yaml"}

#: data files the bridge needs; copied to ``<site-packages>/<path>`` at boot,
#: which is exactly where every module's own resolver already looks
DATA_FILES = [
    "data/derived/malecns_named_neighborhood.json",
    "data/derived/malecns_taste_motor_neighborhood.json",
    "data/derived/malecns_census_v1.json",
    "data/derived/malecns_gustatory_seeds.json",
    "data/literature/behavioral_assays.yaml",
    "data/literature/published_rank_orders.yaml",
    "data/literature/receptor_expression_by_class.yaml",
    "data/literature/resistance_alleles.yaml",
    "data/literature/mixtures.yaml",
    "data/literature/fly_pharmacokinetics.yaml",
]

MAX_TOTAL_MB = 25.0


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_sha() -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO,
            capture_output=True,
            text=True,
            timeout=10,
            check=True,
        )
        return out.stdout.strip() or None
    except Exception:
        return None


def _git_dirty() -> bool | None:
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=REPO,
            capture_output=True,
            text=True,
            timeout=10,
            check=True,
        )
        return bool(out.stdout.strip())
    except Exception:
        return None


def _human(n: int) -> str:
    return f"{n / 1024:.0f} KB" if n < 1024 * 1024 else f"{n / 1024 / 1024:.2f} MB"


# --------------------------------------------------------------------------
# steps
# --------------------------------------------------------------------------
def vendor_js(out: Path) -> list[str]:
    """Download the chart libraries into ``vendor/`` so the site needs no CDN."""
    import urllib.request

    dest = out / "vendor"
    dest.mkdir(parents=True, exist_ok=True)
    written = []
    for name, url in JS_VENDOR.items():
        with urllib.request.urlopen(url, timeout=120) as r:
            (dest / name).write_bytes(r.read())
        written.append(f"vendor/{name}")
    return written


def copy_ui(out: Path, vendor_js_assets: bool = True) -> list[str]:
    """Copy the three static files, making the two asset paths relative.

    GitHub Pages serves a project site from ``/<repo>/``, so the served page
    cannot keep ``/static/app.js``: the rewrite is purely about the URL prefix
    and changes no markup, no styling and no behaviour.  ``flylab-boot.js`` is
    injected before ``app.js`` so that ``window.FLYLAB_STATIC`` and
    ``window.flylabReady`` exist by the time the app script is parsed.
    """
    written = []
    for name in ("app.js", "styles.css"):
        shutil.copy2(STATIC / name, out / name)
        written.append(name)

    html = (STATIC / "index.html").read_text()
    html = html.replace('src="/static/app.js"', 'src="app.js"')
    html = html.replace('href="/static/styles.css"', 'href="styles.css"')
    html = html.replace("/static/styles.css", "styles.css")
    if 'src="app.js"' not in html:
        raise SystemExit("index.html no longer references /static/app.js; update build_pages.py")
    if vendor_js_assets:
        # Point the two chart libraries at the local copies.  app.js keeps its
        # own fallback loader, so a missing vendored file still degrades to the
        # table views rather than breaking the page.
        for name, url in JS_VENDOR.items():
            html = html.replace(url, f"vendor/{name}")
    html = html.replace('<script src="app.js">', '<script src="flylab-boot.js"></script>\n<script src="app.js">')
    (out / "index.html").write_text(html)
    written.append("index.html")
    return written


def build_wheel(out: Path, skip: bool = False, prebuilt: Path | None = None) -> Path:
    """``pip wheel . --no-deps`` into ``dist/pages/wheels``.

    ``prebuilt`` copies an existing wheel instead of building one, which is how
    the cheap test in ``tests/test_build_pages.py`` avoids a ``pip wheel`` run.
    """
    wheels = out / "wheels"
    wheels.mkdir(parents=True, exist_ok=True)
    existing = sorted(wheels.glob("flylab-*.whl"))
    if prebuilt is not None:
        for old_wheel in existing:
            old_wheel.unlink()
        dest = wheels / Path(prebuilt).name
        shutil.copy2(prebuilt, dest)
        print(f"  (--wheel) using {dest.name}")
        return dest
    if skip and existing:
        print(f"  (--skip-wheel) reusing {existing[-1].name}")
        return existing[-1]
    for old in existing:
        old.unlink()
    cmd = [sys.executable, "-m", "pip", "wheel", str(REPO), "--no-deps", "-w", str(wheels), "-q"]
    subprocess.run(cmd, check=True, cwd=str(REPO))
    built = sorted(wheels.glob("flylab-*.whl"))
    if not built:
        raise SystemExit("pip wheel produced no flylab wheel")
    return built[-1]


def wheel_contains(wheel: Path, member_suffix: str) -> bool:
    with zipfile.ZipFile(wheel) as z:
        return any(name.endswith(member_suffix) for name in z.namelist())


def copy_data(out: Path, wheel: Path) -> list[str]:
    copied = []
    for rel in DATA_FILES:
        src = REPO / rel
        if not src.exists():
            print(f"  ! missing {rel} (the bench will answer 404 for it)")
            continue
        dest = out / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        copied.append(rel)
    # library.yaml is package data, so normally it rides inside the wheel
    if not wheel_contains(wheel, "pharm/library.yaml"):
        dest = out / "data" / "library.yaml"
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO / "flylab" / "pharm" / "library.yaml", dest)
        copied.append("data/library.yaml")
        print("  ! library.yaml is not inside the wheel; shipped separately")
    return copied


BOOT_JS = r"""/* FlyLab static build: boot Pyodide, install the flylab wheel, expose
 * window.flylabCall(route, payload) -> Promise<parsed JSON>.
 *
 * Generated by scripts/build_pages.py -- edit that, not this file.
 *
 * This file contains no pharmacology and no circuit maths: every answer comes
 * from flylab.browser.bridge.call, the same function the FastAPI server's
 * routes call.  app.js sees one difference only: window.FLYLAB_STATIC.
 */
(function () {
  "use strict";

  var CFG = __CONFIG__;

  window.FLYLAB_STATIC = true;
  window.FLYLAB_BUILD = CFG;

  /* ---------------------------------------------------------- overlay */
  var overlay = null;
  var bar = null;
  var label = null;
  var detail = null;

  function makeOverlay() {
    overlay = document.createElement("div");
    overlay.id = "flylab-boot";
    overlay.innerHTML =
      '<div class="fb-card" role="status" aria-live="polite">' +
      '<div class="fb-title">FlyLab</div>' +
      '<div class="fb-sub">virtual Drosophila pharmacology bench &mdash; static build</div>' +
      '<div class="fb-track"><div class="fb-bar"></div></div>' +
      '<div class="fb-label">starting Python in your browser&hellip;</div>' +
      '<div class="fb-detail">Pyodide ' + CFG.pyodide_version +
      ' &middot; flylab ' + CFG.flylab_version + " &middot; nothing is uploaded: every number is computed on this machine</div>" +
      "</div>";
    var css = document.createElement("style");
    css.textContent =
      "#flylab-boot{position:fixed;inset:0;z-index:9999;display:flex;align-items:center;" +
      "justify-content:center;background:var(--bg,#0f1115);color:var(--text-1,#e8eaed);" +
      "font:14px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;padding:16px}" +
      "#flylab-boot .fb-card{max-width:520px;width:100%;text-align:left}" +
      "#flylab-boot .fb-title{font-size:28px;font-weight:700;letter-spacing:.08em}" +
      "#flylab-boot .fb-sub{opacity:.75;margin:.25rem 0 1.2rem}" +
      "#flylab-boot .fb-track{height:6px;border-radius:99px;background:rgba(127,127,127,.25);overflow:hidden}" +
      "#flylab-boot .fb-bar{height:100%;width:4%;border-radius:99px;background:currentColor;transition:width .25s ease}" +
      "#flylab-boot .fb-label{margin-top:.75rem;font-variant-numeric:tabular-nums}" +
      "#flylab-boot .fb-detail{margin-top:.5rem;opacity:.6;font-size:12px}" +
      "#flylab-boot.fb-error .fb-bar{background:#d9534f;width:100%}" +
      "@media (prefers-color-scheme: light){#flylab-boot{background:#fbfbfd;color:#1b1c1f}}";
    document.head.appendChild(css);
    document.body.appendChild(overlay);
    bar = overlay.querySelector(".fb-bar");
    label = overlay.querySelector(".fb-label");
    detail = overlay.querySelector(".fb-detail");
  }

  function step(pct, text) {
    if (!overlay) return;
    if (pct !== null && pct !== undefined) bar.style.width = Math.max(0, Math.min(100, pct)) + "%";
    if (text) label.textContent = text;
  }

  /* A CDN having a bad minute must not cost a visitor the whole bench: every
     network step is retried before the page gives up on it. */
  async function withRetry(what, fn, attempts) {
    var last = null;
    var n = attempts || 3;
    for (var i = 1; i <= n; i++) {
      try {
        return await fn();
      } catch (err) {
        last = err;
        console.warn("flylab: " + what + " failed (attempt " + i + "/" + n + ")", err);
        if (i < n) {
          step(null, what + " failed, retrying (" + i + "/" + n + ")\u2026");
          await new Promise(function (r) { setTimeout(r, 900 * i); });
        }
      }
    }
    throw new Error(what + " failed after " + n + " attempts: " + String((last && last.message) || last));
  }

  function fail(message, hint) {
    if (!overlay) return;
    overlay.classList.add("fb-error");
    label.textContent = message;
    detail.innerHTML = hint || "";
  }

  function done() {
    if (!overlay) return;
    overlay.remove();
    overlay = null;
  }

  /* ---------------------------------------------------------- loading */
  function loadScript(url) {
    return new Promise(function (resolve, reject) {
      var tag = document.createElement("script");
      tag.src = url;
      tag.onload = function () { resolve(true); };
      tag.onerror = function () { reject(new Error("could not load " + url)); };
      document.head.appendChild(tag);
    });
  }

  var bridgeCall = null;

  async function start() {
    makeOverlay();
    var t0 = performance.now();
    step(5, "downloading the Python runtime…");
    await withRetry("loading pyodide.js", function () {
      return loadScript(CFG.pyodide_url + "pyodide.js");
    });

    step(18, "starting Python…");
    var pyodide = await withRetry("starting Pyodide", function () {
      return window.loadPyodide({
        indexURL: CFG.pyodide_url,
        stdout: function (s) { console.log("[py]", s); },
        stderr: function (s) { console.warn("[py]", s); },
      });
    });
    window.pyodide = pyodide;

    step(38, "installing numpy and pyyaml…");
    var micropip = await withRetry("loading micropip", async function () {
      await pyodide.loadPackage("micropip", {
        messageCallback: function (m) { step(42, String(m)); },
      });
      return pyodide.pyimport("micropip");
    });
    /* one string per call: micropip.install is a Python function, and a JS
       array would arrive as a JsProxy rather than a list */
    for (var i = 0; i < CFG.runtime_packages.length; i++) {
      await (function (name) {
        var module = CFG.runtime_imports[name] || name;
        /* install AND import inside one retry: micropip reports a failed
           wheel fetch on the console and resolves anyway, so the import is
           what actually proves the package is there. numpy is not optional --
           the rate model, the LIF kernel and every readout are numpy. */
        return withRetry("installing " + name, async function () {
          step(null, "installing " + name + "…");
          await micropip.install(name);
          await pyodide.runPythonAsync("import " + module);
        });
      })(CFG.runtime_packages[i]);
    }

    step(62, "installing the flylab wheel…");
    /* deps:false -- the wheel's metadata also lists the server and CLI
       dependencies (fastapi, pydantic, typer, pandas, pyarrow); the science
       core imports none of them and pyarrow has no WebAssembly build.
       callKwargs, because a trailing JS object is not a Python kwargs dict. */
    await withRetry("installing the flylab wheel", function () {
      return micropip.install.callKwargs(new URL(CFG.wheel, document.baseURI).href, { deps: false });
    });

    step(74, "copying the MaleCNS neighborhood and literature tables…");
    pyodide.globals.set("_flylab_data_files", CFG.data_files);
    pyodide.globals.set("_flylab_base_url", new URL(".", document.baseURI).href);
    await pyodide.runPythonAsync(
      "import pathlib, flylab\n" +
      "from pyodide.http import pyfetch\n" +
      "root = pathlib.Path(flylab.__file__).resolve().parents[1]\n" +
      "for rel in list(_flylab_data_files):\n" +
      "    dest = root / rel\n" +
      "    dest.parent.mkdir(parents=True, exist_ok=True)\n" +
      "    resp = await pyfetch(_flylab_base_url + rel)\n" +
      "    dest.write_bytes(await resp.bytes())\n"
    );

    step(90, "warming the bridge…");
    bridgeCall = pyodide.runPython(
      "from flylab.browser.bridge import call_json, version\n" +
      "import json as _json\n" +
      "print('flylab bridge', _json.dumps(version()))\n" +
      "call_json\n"
    );

    step(100, "ready");
    window.FLYLAB_BOOT_MS = performance.now() - t0;
    console.log("flylab: static bench ready in " + Math.round(window.FLYLAB_BOOT_MS) + " ms");
    done();
  }

  /* window.flylabCall is the whole transport: one JSON in, one JSON out. */
  window.flylabCall = function (route, payload) {
    return window.flylabReady.then(function () {
      var text;
      try {
        text = bridgeCall(route, payload === undefined || payload === null ? null : JSON.stringify(payload));
      } catch (err) {
        return { error: { status: 500, detail: String((err && err.message) || err) } };
      }
      try {
        return JSON.parse(text);
      } catch (err) {
        return { error: { status: 500, detail: "the bridge returned malformed JSON" } };
      }
    });
  };

  window.flylabReady = new Promise(function (resolve, reject) {
    function go() {
      start().then(resolve, function (err) {
        console.error(err);
        fail(
          "could not start the in-browser Python runtime",
          "This build loads Pyodide " + CFG.pyodide_version + " from <code>" + CFG.pyodide_url +
          "</code>. If that host is blocked on this network, the bench cannot run here: " +
          "use the served build instead (<code>pip install flylab &amp;&amp; flylab serve</code>). " +
          "<br>Details: " + String((err && err.message) || err)
        );
        reject(err);
      });
    }
    if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", go);
    else go();
  });
})();
"""


def write_boot(out: Path, wheel: Path, data_files: list[str], pyodide_version: str,
               flylab_version: str, pyodide_url: str) -> None:
    cfg = {
        "pyodide_version": pyodide_version,
        "pyodide_url": pyodide_url,
        "wheel": f"wheels/{wheel.name}",
        "runtime_packages": RUNTIME_PACKAGES,
        "runtime_imports": RUNTIME_IMPORTS,
        "data_files": data_files,
        "flylab_version": flylab_version,
    }
    (out / "flylab-boot.js").write_text(BOOT_JS.replace("__CONFIG__", json.dumps(cfg, indent=2)))


def vendor_pyodide(out: Path, version: str) -> str:
    """Copy the Pyodide runtime into the build (for networks that block the CDN).

    Only the loader, the interpreter and the wheels this bench imports are
    vendored; the other ~340 packages in a Pyodide release are not.
    """
    import urllib.request

    base = PYODIDE_CDN.format(version=version)
    dest = out / "pyodide"
    dest.mkdir(parents=True, exist_ok=True)
    lock_path = dest / "pyodide-lock.json"
    core = [
        "pyodide.js",
        "pyodide.asm.js",
        "pyodide.asm.wasm",
        "pyodide.mjs",
        "python_stdlib.zip",
        "pyodide-lock.json",
    ]
    for name in core:
        with urllib.request.urlopen(base + name, timeout=120) as r:
            (dest / name).write_bytes(r.read())
    lock = json.loads(lock_path.read_text())
    wanted = set(RUNTIME_PACKAGES) | {"micropip", "packaging"}
    queue = list(wanted)
    seen: set[str] = set()
    while queue:
        key = queue.pop().lower()
        if key in seen:
            continue
        seen.add(key)
        entry = lock["packages"].get(key)
        if not entry:
            continue
        queue.extend(entry.get("depends") or [])
        with urllib.request.urlopen(base + entry["file_name"], timeout=300) as r:
            (dest / entry["file_name"]).write_bytes(r.read())
    return "pyodide/"


def size_report(out: Path) -> tuple[list[tuple[str, int]], int]:
    rows: list[tuple[str, int]] = []
    total = 0
    for path in sorted(out.rglob("*")):
        if path.is_file():
            n = path.stat().st_size
            rows.append((str(path.relative_to(out)), n))
            total += n
    return rows, total


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def build(
    out: Path = DEFAULT_OUT,
    skip_wheel: bool = False,
    pyodide_version: str = PYODIDE_VERSION,
    vendor: bool = False,
    max_mb: float = MAX_TOTAL_MB,
    clean: bool = True,
    wheel_path: Path | None = None,
    vendor_js_assets: bool = True,
) -> dict:
    """Build the static site and return its manifest."""
    if clean and out.exists() and not skip_wheel:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    print(f"FlyLab static build -> {out}")
    ui = copy_ui(out, vendor_js_assets=vendor_js_assets)
    print(f"  ui: {', '.join(ui)}")
    if vendor_js_assets:
        ui += vendor_js(out)
        print(f"  vendored chart libraries: {', '.join(JS_VENDOR)}")
    wheel = build_wheel(out, skip=skip_wheel, prebuilt=wheel_path)
    print(f"  wheel: {wheel.name} ({_human(wheel.stat().st_size)})")
    data = copy_data(out, wheel)
    print(f"  data: {len(data)} files")

    pyodide_url = PYODIDE_CDN.format(version=pyodide_version)
    if vendor:
        print("  vendoring the Pyodide runtime (this downloads ~15 MB)")
        pyodide_url = vendor_pyodide(out, pyodide_version)

    flylab_version = wheel.name.split("-")[1]
    write_boot(out, wheel, data, pyodide_version, flylab_version, pyodide_url)

    lib = REPO / "flylab" / "pharm" / "library.yaml"
    rows, total = size_report(out)
    manifest = {
        "name": "flylab-static",
        "flylab_version": flylab_version,
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_sha": _git_sha(),
        "git_dirty": _git_dirty(),
        "pyodide_version": pyodide_version,
        "pyodide_url": pyodide_url,
        "js_vendored": vendor_js_assets,
        "pyodide_vendored": bool(vendor),
        "wheel": f"wheels/{wheel.name}",
        "wheel_sha256": _sha256(wheel),
        "library_sha256": _sha256(lib),
        "runtime_packages": RUNTIME_PACKAGES,
        "runtime_imports": RUNTIME_IMPORTS,
        "data_files": data,
        "files": [{"path": p, "bytes": n} for p, n in rows],
        "total_bytes": total,
        "max_total_bytes": int(max_mb * 1024 * 1024),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    rows, total = size_report(out)
    manifest["files"] = [{"path": p, "bytes": n} for p, n in rows]
    manifest["total_bytes"] = total
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print("\n  size report")
    for p, n in sorted(rows, key=lambda r: -r[1])[:14]:
        print(f"    {_human(n):>10}  {p}")
    print(f"    {'-' * 10}")
    print(f"    {_human(total):>10}  TOTAL ({len(rows)} files, budget {max_mb:.0f} MB)")
    if total > max_mb * 1024 * 1024:
        raise SystemExit(
            f"static build is {_human(total)}, over the {max_mb:.0f} MB budget. "
            "Drop a data file, or vendor fewer Pyodide packages."
        )
    return manifest


def serve(out: Path, port: int = 8000) -> None:
    print(f"\nserving {out} at http://127.0.0.1:{port}/  (ctrl-c to stop)")
    subprocess.run(
        [sys.executable, "-m", "http.server", str(port), "--bind", "127.0.0.1"],
        cwd=str(out),
        check=False,
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output directory (default dist/pages)")
    ap.add_argument("--skip-wheel", action="store_true", help="reuse the wheel already in the output dir")
    ap.add_argument("--wheel", type=Path, default=None,
                    help="copy this prebuilt flylab wheel instead of running pip wheel")
    ap.add_argument("--pyodide-version", default=PYODIDE_VERSION)
    ap.add_argument("--vendor-pyodide", action="store_true",
                    help="download the Pyodide runtime into the build instead of using the CDN")
    ap.add_argument("--no-vendor-js", action="store_true",
                    help="load Plotly and cytoscape from their CDNs instead of vendoring them")
    ap.add_argument("--max-mb", type=float, default=MAX_TOTAL_MB)
    ap.add_argument("--serve", action="store_true", help="preview the build over python -m http.server")
    ap.add_argument("--port", type=int, default=8000)
    args = ap.parse_args(argv)

    manifest = build(
        out=args.out,
        skip_wheel=args.skip_wheel,
        pyodide_version=args.pyodide_version,
        vendor=args.vendor_pyodide,
        max_mb=args.max_mb,
        wheel_path=args.wheel,
        vendor_js_assets=not args.no_vendor_js,
    )
    print(f"\n  git sha      {manifest['git_sha']}")
    print(f"  library      {manifest['library_sha256'][:16]}…")
    print(f"  pyodide      {manifest['pyodide_version']} ({'vendored' if manifest['pyodide_vendored'] else 'CDN'})")
    if args.serve:
        serve(args.out, args.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
