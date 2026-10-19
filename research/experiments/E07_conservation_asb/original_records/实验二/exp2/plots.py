from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_metaprofile(x: np.ndarray, y: np.ndarray, out: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x, y, color="#3C5488", linewidth=1.8)
    ax.axhline(0.0, color="#999999", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Relative position (motif-centered)")
    ax.set_ylabel("Mean PhyloP")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_core_vs_flank_violin(inst: pd.DataFrame, out: Path) -> None:
    if inst.empty:
        return
    df = inst[["mean_phyloP_core", "mean_phyloP_flank"]].melt(var_name="region", value_name="PhyloP")
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.violinplot(data=df, x="region", y="PhyloP", ax=ax, hue="region", palette=["#4DBBD5", "#E64B35"], legend=False)
    ax.set_title("Core vs flank mean PhyloP (instances)")
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_motif_position_heatmap(pos: pd.DataFrame, out: Path, max_motifs: int = 40) -> None:
    if pos.empty:
        return
    if "task_id" in pos.columns:
        pos = pos.copy()
        pos["motif_key"] = pos["task_id"].astype(str) + "|" + pos["motif_id"].astype(str)
        keycol = "motif_key"
    else:
        pos = pos.copy()
        pos["motif_key"] = pos["motif_id"].astype(str)
        keycol = "motif_key"
    motifs = pos[keycol].unique().tolist()
    if len(motifs) > max_motifs:
        vol = pos.groupby(keycol)["mean_phyloP"].count().sort_values(ascending=False)
        motifs = vol.head(max_motifs).index.tolist()
        sub = pos[pos[keycol].isin(motifs)]
    else:
        sub = pos
    pivot = sub.pivot_table(index=keycol, columns="rel_pos", values="mean_phyloP", aggfunc="mean")
    fig, ax = plt.subplots(figsize=(min(14, 0.35 * pivot.shape[1] + 4), max(4, 0.25 * pivot.shape[0] + 2)))
    sns.heatmap(pivot, cmap="viridis", ax=ax, cbar_kws={"label": "mean PhyloP"})
    ax.set_title("Motif × relative position")
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_ic_cwm_vs_phylop_scatter(pos: pd.DataFrame, out: Path) -> None:
    if pos.empty:
        return
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(pos["abs_cwm"], pos["mean_phyloP"], s=12, alpha=0.5, c="#3C5488", label="|CWM|")
    ax.set_xlabel("|CWM| (pattern position)")
    ax.set_ylabel("mean PhyloP")
    ax.legend()
    ax.set_title("Position-wise |CWM| vs mean PhyloP")
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
