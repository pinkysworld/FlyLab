"""FlyLab notebook schema.

0.3 adds ``created_utc``, a ``provenance`` block, a ``gains`` block and an
``uncertainty`` slot. Every key present in 0.2 is still there, so code written
against the older schema keeps working.
"""

from __future__ import annotations

import copy
import platform as _platform
import subprocess
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any

NOTEBOOK_VERSION = "0.3"

_REPO_ROOT = Path(__file__).resolve().parents[2]


def _flylab_version() -> str | None:
    try:
        from importlib.metadata import PackageNotFoundError, version

        try:
            return version("flylab")
        except PackageNotFoundError:
            pass
    except Exception:  # pragma: no cover - importlib.metadata is stdlib
        pass
    try:
        import flylab

        return getattr(flylab, "__version__", None)
    except Exception:  # pragma: no cover
        return None


@lru_cache(maxsize=1)
def git_sha() -> str | None:
    """Current commit, or None outside a git checkout."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=_REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=5,
            check=True,
        )
        return out.stdout.strip() or None
    except Exception:
        return None


@lru_cache(maxsize=1)
def _library_hash() -> str | None:
    try:
        from flylab.pharm.occupancy import library_sha256

        return library_sha256()
    except Exception:
        return None


def provenance(rng_seed: int | None = None) -> dict[str, Any]:
    """Provenance block: versions, commit, library hash, platform, seed."""
    return {
        "flylab_version": _flylab_version(),
        "git_sha": git_sha(),
        "library_sha256": _library_hash(),
        "map": {"name": None, "version": None, "citation": None},
        "rng_seed": rng_seed,
        "platform": f"{_platform.platform()} python-{_platform.python_version()}",
    }


def empty_notebook(assay: str, seed: int | None = None) -> dict[str, Any]:
    """A blank notebook for one assay run."""
    return {
        "flylab_notebook_version": NOTEBOOK_VERSION,
        "assay": assay,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "provenance": provenance(seed),
        "map": {"name": None, "version": None, "citation": None},
        "compound": None,
        "concentration_M": None,
        "occupancy": [],
        "gains": {},
        "readouts": {},
        "uncertainty": None,
        "live_lab": None,
        "warnings": [],
    }


def set_map(notebook: dict[str, Any], name: str | None, version: str | None, citation: str | None) -> dict[str, Any]:
    """Set the map block in both the legacy slot and the provenance block."""
    block = {"name": name, "version": version, "citation": citation}
    notebook["map"] = block
    notebook.setdefault("provenance", provenance())["map"] = copy.deepcopy(block)
    return notebook
