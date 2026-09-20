"""The static build must assemble a complete, self-describing site.

The cheap tests here never run ``pip wheel``: they hand ``build_pages`` a
prebuilt (stub) wheel and check the layout, the manifest and the boot script.
The real end-to-end build is marked ``slow`` and is excluded from ``pytest -q``
by the marker filter in ``pyproject.toml`` (run it with ``pytest -m slow``).
"""
from __future__ import annotations

import importlib.util
import json
import sys
import zipfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

spec = importlib.util.spec_from_file_location("build_pages", REPO / "scripts" / "build_pages.py")
build_pages = importlib.util.module_from_spec(spec)
sys.modules["build_pages"] = build_pages
assert spec.loader is not None
spec.loader.exec_module(build_pages)


@pytest.fixture
def stub_wheel(tmp_path: Path) -> Path:
    """A wheel-shaped zip: enough for the layout checks, free to make."""
    path = tmp_path / "flylab-0.5.0-py3-none-any.whl"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("flylab/__init__.py", "__version__ = '0.5.0'\n")
        z.writestr("flylab/pharm/library.yaml", "schema_version: 2\n")
        z.writestr("flylab-0.5.0.dist-info/METADATA", "Name: flylab\nVersion: 0.5.0\n")
    return path


@pytest.fixture
def site(tmp_path: Path, stub_wheel: Path) -> tuple[Path, dict]:
    out = tmp_path / "pages"
    manifest = build_pages.build(out=out, wheel_path=stub_wheel)
    return out, manifest


def test_site_has_every_file_the_bench_needs(site):
    out, _ = site
    for name in ("index.html", "app.js", "styles.css", "flylab-boot.js", "manifest.json"):
        assert (out / name).is_file(), name
    assert (out / "wheels" / "flylab-0.5.0-py3-none-any.whl").is_file()
    for rel in build_pages.DATA_FILES:
        assert (out / rel).is_file(), rel
    assert (out / "data" / "derived" / "malecns_named_neighborhood.json").stat().st_size > 1000


def test_manifest_records_the_build(site):
    out, manifest = site
    on_disk = json.loads((out / "manifest.json").read_text())
    assert on_disk == manifest
    assert manifest["flylab_version"] == "0.5.0"
    assert len(manifest["library_sha256"]) == 64
    assert len(manifest["wheel_sha256"]) == 64
    assert manifest["pyodide_version"] == build_pages.PYODIDE_VERSION
    assert manifest["pyodide_url"].startswith("https://")
    assert manifest["git_sha"] is None or len(manifest["git_sha"]) == 40
    assert manifest["total_bytes"] == sum(f["bytes"] for f in manifest["files"])
    assert manifest["total_bytes"] < manifest["max_total_bytes"]
    paths = {f["path"] for f in manifest["files"]}
    assert "index.html" in paths and "flylab-boot.js" in paths
    assert manifest["runtime_packages"] == ["numpy", "pyyaml"]


def test_page_loads_the_boot_script_before_the_app(site):
    out, _ = site
    html = (out / "index.html").read_text()
    assert "/static/app.js" not in html and "/static/styles.css" not in html
    assert html.index("flylab-boot.js") < html.index('src="app.js"')
    # the UI itself is copied byte for byte: the static build changes transport,
    # not the bench
    assert (out / "app.js").read_bytes() == (REPO / "flylab" / "static" / "app.js").read_bytes()
    assert (out / "styles.css").read_bytes() == (
        REPO / "flylab" / "static" / "styles.css"
    ).read_bytes()


def test_boot_script_is_configured_and_carries_no_science(site):
    out, manifest = site
    boot = (out / "flylab-boot.js").read_text()
    assert "__CONFIG__" not in boot
    assert manifest["wheel"] in boot
    assert build_pages.PYODIDE_VERSION in boot
    assert "window.flylabCall" in boot and "window.FLYLAB_STATIC" in boot
    assert "flylab.browser.bridge" in boot
    for forbidden in ("ec50", "occupancy", "g_ach", "hill"):
        assert forbidden not in boot.lower(), f"{forbidden!r} must not be reimplemented in JS"


def test_build_refuses_to_go_over_budget(tmp_path, stub_wheel):
    with pytest.raises(SystemExit) as exc:
        build_pages.build(out=tmp_path / "tiny", wheel_path=stub_wheel, max_mb=0.001)
    assert "budget" in str(exc.value)


def test_data_files_listed_for_the_boot_script_exist_in_the_repo():
    for rel in build_pages.DATA_FILES:
        assert (REPO / rel).is_file(), rel


@pytest.mark.slow
def test_end_to_end_build_with_a_real_wheel(tmp_path):
    """The whole thing: pip wheel, data copy, manifest, size budget."""
    out = tmp_path / "pages"
    manifest = build_pages.build(out=out)
    wheel = out / manifest["wheel"]
    assert wheel.is_file()
    with zipfile.ZipFile(wheel) as z:
        names = z.namelist()
    assert "flylab/browser/bridge.py" in names
    assert any(n.endswith("pharm/library.yaml") for n in names)
    assert manifest["total_bytes"] < 25 * 1024 * 1024
    assert manifest["flylab_version"] == "0.5.0"
