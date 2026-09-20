"""Hill fits of the *model's own* dose-response curves.

Nothing in here is fitted to animal data.  Every number these functions return
is a descriptive summary of a FlyLab simulation and must be labelled
``model_derived`` in any notebook or figure (see ``docs/DESIGN_v0.5.md``).

scipy is deliberately **not** a dependency: the 4-parameter Hill fit uses a
coarse grid over ``log10(IC50)`` and slope followed by a small pure-numpy
Nelder-Mead refinement.
"""
from __future__ import annotations

from typing import Any, Callable, Sequence

import numpy as np

__all__ = ["hill4", "hill4_fit", "bootstrap_fit", "nelder_mead"]


def hill4(x: np.ndarray, top: float, bottom: float, ic50: float, slope: float) -> np.ndarray:
    """4-parameter Hill / logistic in concentration space.

    ``y = bottom + (top - bottom) / (1 + (x / ic50) ** slope)``

    With ``slope > 0`` the curve falls from ``top`` (low concentration) to
    ``bottom`` (high concentration), i.e. an inhibition curve.
    """
    x = np.asarray(x, dtype=float)
    ic50 = max(float(ic50), 1e-30)
    with np.errstate(over="ignore", invalid="ignore"):
        ratio = np.power(np.maximum(x, 1e-30) / ic50, float(slope))
    out = bottom + (top - bottom) / (1.0 + ratio)
    return np.nan_to_num(out, nan=float(bottom), posinf=float(top), neginf=float(bottom))


def _hill4_log(u: np.ndarray, top: float, bottom: float, log_ic50: float, slope: float) -> np.ndarray:
    """Same curve parameterised by ``u = log10(x)`` (numerically nicer)."""
    z = np.clip(float(slope) * (np.asarray(u, dtype=float) - float(log_ic50)), -60.0, 60.0)
    return bottom + (top - bottom) / (1.0 + np.power(10.0, z))


def nelder_mead(
    fun: Callable[[np.ndarray], float],
    x0: Sequence[float],
    step: Sequence[float] | None = None,
    max_iter: int = 800,
    tol: float = 1e-10,
) -> tuple[np.ndarray, float, int]:
    """Minimal Nelder-Mead simplex minimiser (pure numpy)."""
    x0 = np.asarray(x0, dtype=float)
    n = x0.size
    step = np.asarray(step if step is not None else np.where(np.abs(x0) > 1e-8, 0.1 * np.abs(x0), 0.1), dtype=float)
    sim = np.vstack([x0] + [x0 + np.eye(n)[i] * step[i] for i in range(n)])
    f = np.array([fun(p) for p in sim])
    it = 0
    for it in range(1, int(max_iter) + 1):
        order = np.argsort(f)
        sim, f = sim[order], f[order]
        if np.max(np.abs(f[1:] - f[0])) <= tol and np.max(np.abs(sim[1:] - sim[0])) <= tol:
            break
        centroid = sim[:-1].mean(axis=0)
        xr = centroid + 1.0 * (centroid - sim[-1])
        fr = fun(xr)
        if fr < f[0]:
            xe = centroid + 2.0 * (centroid - sim[-1])
            fe = fun(xe)
            sim[-1], f[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < f[-2]:
            sim[-1], f[-1] = xr, fr
        else:
            xc = centroid + 0.5 * (sim[-1] - centroid)
            fc = fun(xc)
            if fc < f[-1]:
                sim[-1], f[-1] = xc, fc
            else:
                sim[1:] = sim[0] + 0.5 * (sim[1:] - sim[0])
                f[1:] = np.array([fun(p) for p in sim[1:]])
    order = np.argsort(f)
    return sim[order][0], float(f[order][0]), it


def _r2(y: np.ndarray, yhat: np.ndarray) -> float:
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    if ss_tot <= 1e-12:
        return 1.0 if ss_res <= 1e-12 else 0.0
    return 1.0 - ss_res / ss_tot


def hill4_fit(
    x: Sequence[float],
    y: Sequence[float],
    slope_bounds: tuple[float, float] = (0.15, 8.0),
    max_iter: int = 400,
) -> dict[str, Any]:
    """Fit ``top``, ``bottom``, ``ic50``, ``slope`` by least squares in log-x space.

    Returns ``{top, bottom, ic50, log10_ic50, slope, r2, n, converged,
    sse, in_range}``.  ``in_range`` says whether the fitted IC50 lies inside the
    span of the tested concentrations; if not, treat the value as an
    extrapolation.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x, y = x[keep], y[keep]
    n = int(x.size)
    if n < 4:
        raise ValueError("hill4_fit needs at least 4 finite points with x > 0")
    u = np.log10(x)
    order = np.argsort(u)
    u, y, x = u[order], y[order], x[order]

    y_lo, y_hi = float(y.min()), float(y.max())
    if abs(y_hi - y_lo) < 1e-12:  # flat curve: no IC50 to find
        return {
            "top": y_hi,
            "bottom": y_lo,
            "ic50": float("nan"),
            "log10_ic50": float("nan"),
            "slope": float("nan"),
            "r2": 0.0,
            "sse": 0.0,
            "n": n,
            "converged": False,
            "in_range": False,
            "note": "response is flat across the tested range; no IC50 is defined",
        }

    top0 = float(np.mean(y[: max(1, n // 4)]))
    bot0 = float(np.mean(y[-max(1, n // 4) :]))

    def sse(p: np.ndarray) -> float:
        slope = float(np.clip(p[3], *slope_bounds))
        return float(np.sum((y - _hill4_log(u, p[0], p[1], p[2], slope)) ** 2))

    # coarse vectorised grid over log10(IC50) x slope, top/bottom from the ends
    grid_p = np.linspace(u.min() - 1.0, u.max() + 1.0, 41)
    grid_s = np.array([0.3, 0.6, 1.0, 1.5, 2.5, 4.0])
    zz = np.clip(grid_s[None, :, None] * (u[None, None, :] - grid_p[:, None, None]), -60.0, 60.0)
    pred = bot0 + (top0 - bot0) / (1.0 + np.power(10.0, zz))
    grid_sse = np.sum((pred - y[None, None, :]) ** 2, axis=2)
    gi, gj = np.unravel_index(int(np.argmin(grid_sse)), grid_sse.shape)
    p0 = np.array([top0, bot0, float(grid_p[gi]), float(grid_s[gj])])
    p, f, it = nelder_mead(
        sse,
        p0,
        step=[max(abs(top0), 1.0) * 0.2, max(abs(bot0), 1.0) * 0.2, 0.5, 0.4],
        max_iter=max_iter,
    )
    slope = float(np.clip(p[3], *slope_bounds))
    log_ic50 = float(p[2])
    yhat = _hill4_log(u, p[0], p[1], log_ic50, slope)
    return {
        "top": float(p[0]),
        "bottom": float(p[1]),
        "ic50": float(10.0**log_ic50),
        "log10_ic50": log_ic50,
        "slope": slope,
        "r2": _r2(y, yhat),
        "sse": float(f),
        "n": n,
        "converged": bool(it < max_iter),
        "in_range": bool(u.min() <= log_ic50 <= u.max()),
    }


def bootstrap_fit(
    x: Sequence[float],
    y_reps: Sequence[Sequence[float]],
    n_boot: int = 200,
    seed: int = 0,
    ci: tuple[float, float] = (2.5, 97.5),
) -> dict[str, Any]:
    """Bootstrap a Hill fit by resampling replicates within each concentration.

    ``y_reps`` is ``n_points x n_replicates``.  Returns the fit on the replicate
    means plus percentile CIs for every parameter.
    """
    x = np.asarray(x, dtype=float)
    Y = np.asarray([np.asarray(r, dtype=float) for r in y_reps], dtype=float)
    if Y.ndim == 1:
        Y = Y[:, None]
    base = hill4_fit(x, Y.mean(axis=1))
    rng = np.random.default_rng(seed)
    n_pts, n_rep = Y.shape
    keys = ("ic50", "log10_ic50", "slope", "top", "bottom", "r2")
    draws: dict[str, list[float]] = {k: [] for k in keys}
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n_rep, size=(n_pts, n_rep))
        ys = np.take_along_axis(Y, idx, axis=1).mean(axis=1)
        try:
            fit = hill4_fit(x, ys)
        except ValueError:
            continue
        if not np.isfinite(fit["ic50"]):
            continue
        for k in keys:
            draws[k].append(float(fit[k]))
    out_ci: dict[str, list[float]] = {}
    for k, vals in draws.items():
        if len(vals) >= 5:
            lo, hi = np.percentile(vals, ci)
            out_ci[k] = [float(lo), float(hi)]
    return {
        "fit": base,
        "ci": out_ci,
        "n_boot": int(n_boot),
        "n_ok": len(draws["ic50"]),
        "method": "nonparametric bootstrap over model replicates, percentile CI",
    }
