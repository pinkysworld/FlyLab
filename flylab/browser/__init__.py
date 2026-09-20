"""Backend-agnostic entry point for the FlyLab science core.

``flylab.server`` (FastAPI) and ``flylab.browser.bridge`` (Pyodide /
WebAssembly) are two transports over *the same* library functions: the routes,
bounds and error contract are mirrored, no pharmacology or circuit maths lives
in either one, and neither is allowed to be a second implementation of the
science.  See ``docs/PAGES.md``.
"""
from __future__ import annotations

__all__ = ["bridge"]
