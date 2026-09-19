from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt

def save_wholens_figures(nb: dict, outdir: Path) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    n = nb["readouts"]["neurotransmitter_counts"]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    keys = list(n.keys())
    ax.bar(keys, [n[k] for k in keys], color="#1f4e5f")
    ax.set_ylabel("traced neurons (MaleCNS v1.0)")
    ax.set_title("Predicted transmitter census")
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    p = outdir / "fig_nt_census.png"
    fig.savefig(p, dpi=160); plt.close(fig); paths.append(p)
    occ = nb.get("occupancy") or []
    if occ:
        fig, ax = plt.subplots(figsize=(7.2, 3.6))
        labels = [r["receptor"] for r in occ]
        vals = [r["occupancy"] for r in occ]
        colors = ["#d9772c" if l.startswith("insect") else "#3d5a5b" for l in labels]
        ax.bar(labels, vals, color=colors)
        ax.set_ylim(0, 1); ax.set_ylabel("occupancy")
        ax.set_title(f"{nb.get('compound')} @ {nb.get('concentration_M'):.1e} M")
        ax.tick_params(axis="x", rotation=25)
        fig.tight_layout()
        p = outdir / "fig_occupancy.png"
        fig.savefig(p, dpi=160); plt.close(fig); paths.append(p)
    return paths
