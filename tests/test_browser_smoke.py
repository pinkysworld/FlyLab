"""A real-browser smoke test for the static Pyodide bench.

This is intentionally opt-in: it builds a site, starts a local HTTP server,
downloads the pinned Pyodide runtime in Chromium, and runs one dashboard call.
"""
from __future__ import annotations

import functools
import importlib.util
import json
import sys
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def _load_build_pages():
    spec = importlib.util.spec_from_file_location("build_pages", REPO / "scripts" / "build_pages.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["build_pages"] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class _QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, _format, *args):
        pass


@pytest.mark.slow
def test_static_bench_starts_runs_dashboard_and_falls_back_without_charts(tmp_path):
    playwright_api = pytest.importorskip(
        "playwright.sync_api",
        reason="install the browser-test extra and Chromium to run the static-browser smoke test",
    )

    site = tmp_path / "pages"
    _load_build_pages().build(out=site)

    handler = functools.partial(_QuietHandler, directory=str(site))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    server_thread = threading.Thread(target=server.serve_forever, daemon=True)
    server_thread.start()
    url = f"http://127.0.0.1:{server.server_port}/"

    try:
        with playwright_api.sync_playwright() as playwright:
            browser = playwright.chromium.launch(headless=True)
            try:
                page = browser.new_page()
                page.set_default_timeout(30_000)

                def block_chart_assets(route):
                    request = route.request
                    path = request.url.lower()
                    is_chart_script = request.resource_type == "script" and any(
                        name in path for name in ("plotly", "cytoscape")
                    )
                    if is_chart_script:
                        route.abort()
                    else:
                        route.continue_()

                page.route("**/*", block_chart_assets)
                page.goto(url, wait_until="domcontentloaded", timeout=60_000)

                # The overlay remains visible on boot failure. Wait until the
                # loader removes it, or surfaces its explicit error state.
                page.wait_for_function(
                    "() => !document.querySelector('#flylab-boot') || "
                    "document.querySelector('#flylab-boot.fb-error') !== null",
                    timeout=10 * 60_000,
                )
                boot_error = page.locator("#flylab-boot.fb-error .fb-label").count()
                assert boot_error == 0, page.locator("#flylab-boot").inner_text()
                assert page.evaluate("window.FLYLAB_STATIC") is True

                # The first dashboard activation is automatic. Wait for its
                # data and verify the chart-less path leaves the equivalent table.
                page.wait_for_function(
                    "() => document.querySelectorAll('#dash-tiles .tile').length >= 4 && "
                    "document.querySelector('#status-text')?.textContent === 'dashboard ready'",
                    timeout=180_000,
                )
                assert page.locator("#tbl-dash-ladder tbody tr").count() >= 1
                assert "Plotly could not be loaded" in page.locator("#plot-dash-ladder").inner_text()

                # Exercise a deliberate dashboard run and prove it traversed
                # the Pyodide bridge, then switch panels and return.
                page.evaluate(
                    """() => {
                      const original = window.flylabCall;
                      window.__smokeCalls = {};
                      window.flylabCall = function (route, payload) {
                        const key = route.split('?')[0];
                        window.__smokeCalls[key] = (window.__smokeCalls[key] || 0) + 1;
                        return original.call(this, route, payload);
                      };
                    }"""
                )
                page.locator("button[data-run='dashboard']").click()
                page.wait_for_function(
                    "() => (window.__smokeCalls['/api/dashboard'] || 0) === 1 && "
                    "document.querySelector('#status-text')?.textContent === 'dashboard ready'",
                    timeout=180_000,
                )

                # Portable setup files round-trip the run controls and the
                # experiment form without carrying stale calculated output.
                page.locator("#f-compound").select_option("fipronil")
                page.locator("#f-conc").fill("3e-7")
                page.locator("#f-conc").dispatch_event("change")
                page.locator("#f-exp-compounds").select_option(["fipronil"])
                page.locator("#f-exp-start").fill("1e-8")
                page.locator("#f-exp-start").dispatch_event("change")
                page.locator("#f-exp-stop").fill("1e-5")
                page.locator("#f-exp-stop").dispatch_event("change")
                page.locator("#f-exp-ppd").fill("2")
                page.locator("#f-exp-ppd").dispatch_event("change")
                page.locator("#f-exp-reps").fill("4")
                page.locator("#f-exp-reps").dispatch_event("change")
                page.locator("#f-exp-assay").select_option("taste_map")
                page.locator("#f-exp-readouts").select_option(["mn9_hz", "mean_hz"])
                page.locator("#f-exp-randomize").check()
                page.locator("#f-exp-blind").check()

                with page.expect_download() as download_info:
                    page.locator("#btn-preset-export").click()
                download = download_info.value
                preset_path = tmp_path / download.suggested_filename
                download.save_as(preset_path)
                preset = json.loads(preset_path.read_text())
                assert preset["format"] == "flylab.design-preset"
                assert preset["version"] == 1
                assert preset["design"]["compound"] == "fipronil"
                assert preset["experiment"]["assay"] == "taste_map"
                assert preset["experiment"]["replicates"] == 4

                page.locator("#f-compound").select_option("imidacloprid")
                page.locator("#f-exp-assay").select_option("subgraph")
                page.locator("#f-exp-reps").fill("2")
                page.locator("#f-exp-reps").dispatch_event("change")
                page.locator("#f-preset-file").set_input_files(str(preset_path))
                page.locator("#preset-status").wait_for(state="visible")
                page.wait_for_function("() => document.querySelector('#preset-status').textContent.startsWith('Imported preset.')")
                assert page.locator("#f-compound").input_value() == "fipronil"
                assert float(page.locator("#f-conc").input_value()) == pytest.approx(3e-7)
                assert page.locator("#f-exp-assay").input_value() == "taste_map"
                assert page.locator("#f-exp-reps").input_value() == "4"
                assert page.locator("#f-exp-randomize").is_checked()
                assert page.locator("#f-exp-blind").is_checked()

                # Invalid files must not partially change the active design.
                invalid_path = tmp_path / "invalid-preset.json"
                invalid_path.write_text(json.dumps({"format": "wrong"}))
                page.locator("#f-preset-file").set_input_files(str(invalid_path))
                page.wait_for_function("() => document.querySelector('#preset-status').textContent.startsWith('Import failed:')")
                assert page.locator("#f-compound").input_value() == "fipronil"

                page.locator("#tab-notebook").click()
                assert page.locator("#tab-notebook").get_attribute("aria-selected") == "true"
                assert "active" in page.locator("#panel-notebook").get_attribute("class")
                page.locator("#tab-dashboard").click()
                assert page.locator("#tab-dashboard").get_attribute("aria-selected") == "true"
                assert "active" in page.locator("#panel-dashboard").get_attribute("class")

                # At phone width, the grouped native picker replaces the two
                # horizontal tab rows and stays synced with the active panel.
                page.set_viewport_size({"width": 390, "height": 844})
                assert page.locator("#mobile-panel-picker").is_visible()
                assert page.locator("#main > .tabs").first.is_hidden()
                page.locator("#mobile-panel-picker").select_option("controls")
                assert "active" in page.locator("#panel-controls").get_attribute("class")
                assert page.locator("#tab-controls").get_attribute("aria-selected") == "true"
            finally:
                browser.close()
    finally:
        server.shutdown()
        server.server_close()
        server_thread.join(timeout=5)
