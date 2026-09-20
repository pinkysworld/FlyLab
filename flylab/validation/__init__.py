"""Retrospective validation of FlyLab against published results.

Nothing in this package is fitted: it compares the model's own output with
orderings and classifications taken from the literature and reports the
mismatches as findings.
"""

from flylab.validation.rank import (
    KNOWN_DISCREPANCIES,
    kendall_tau,
    load_rank_orders,
    spearman_rho,
    validate_all,
    validate_entry,
    validation_report_markdown,
)

__all__ = [
    "load_rank_orders",
    "validate_entry",
    "validate_all",
    "validation_report_markdown",
    "spearman_rho",
    "kendall_tau",
    "KNOWN_DISCREPANCIES",
]
