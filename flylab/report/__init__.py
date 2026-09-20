"""Reporting artifacts: the human-readable end of a FlyLab run.

:mod:`flylab.report.card` turns one result into a **claim card** -- a fixed
nine-section artifact that says what is claimed, what backs it, what was
assumed, whether the claim needed the connectome at all, and what is still
unknown.  The card is emitted as JSON and as Markdown; nothing in this package
computes pharmacology or circuit maths.

Imports are function-local so ``import flylab.report`` stays cheap and works in
the browser build, where neither pydantic nor a web framework exists.
"""

from __future__ import annotations

__all__ = ["claim_card", "to_markdown", "cards_for_run", "SECTIONS", "CARD_VERSION"]


def __getattr__(name: str):  # pragma: no cover - thin lazy re-export
    if name in __all__:
        from flylab.report import card as _card

        return getattr(_card, name)
    raise AttributeError(name)
