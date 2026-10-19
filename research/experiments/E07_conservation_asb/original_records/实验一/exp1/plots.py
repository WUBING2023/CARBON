from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_peak_vs_matched(df: Path | pd.DataFrame, out_path: Path, *, max_points: int = 50000) -> None:
    if isinstance(df, Path):
        peak_df = pd.read_csv(df)
    else:
        peak_df = df
    peak_df = peak_df[np.isfinite(peak_df["mean_phyloP_core"]) & np.isfinite(peak_df["mean_phyloP_matched_bg"])]
    if peak_df.empty:
        return
    if len(peak_df) > max_points:
        peak_df = peak_df.sample(max_points, random_state=0)
    plt.figure(figsize=(8, 5))
    data = pd.DataFrame(
        {
            "PhyloP": np.concatenate(
                [peak_df["mean_phyloP_core"].values, peak_df["mean_phyloP_matched_bg"].values]
            ),
            "Region": ["peak core"] * len(peak_df) + ["matched non-peak"] * len(peak_df),
        }
    )
    sns.violinplot(data=data, x="Region", y="PhyloP", inner="box")
    plt.title("Mean PhyloP: peak core vs matched non-peak (Tier A)")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150)
    plt.close()


def plot_metaplot(x: np.ndarray, y: np.ndarray, out_path: Path, *, title: str) -> None:
    plt.figure(figsize=(8, 4))
    plt.plot(x, y, color="C0", linewidth=2)
    plt.axvline(0, color="gray", linestyle="--", linewidth=1)
    plt.xlabel("Distance to summit (bp)")
    plt.ylabel("Mean PhyloP")
    plt.title(title)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150)
    plt.close()


def plot_task_delta_heatmap(task_summary: pd.DataFrame, out_path: Path) -> None:
    if task_summary.empty or "delta_median" not in task_summary.columns:
        return
    sub = task_summary[["task_id", "cell_line", "delta_median"]].dropna()
    if sub.empty:
        return
    sub = sub.sort_values("delta_median", ascending=False)
    mat = sub["delta_median"].to_numpy(dtype=np.float64).reshape(-1, 1)
    plt.figure(figsize=(4, max(6, len(sub) * 0.08)))
    sns.heatmap(mat, cmap="RdBu_r", center=0, cbar_kws={"label": "ΔPhyloP median"})
    plt.yticks(np.arange(len(sub)) + 0.5, sub["task_id"].tolist(), fontsize=4)
    plt.xticks([0.5], ["Δmedian"])
    plt.title("Per-task ΔPhyloP (core − matched)")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
