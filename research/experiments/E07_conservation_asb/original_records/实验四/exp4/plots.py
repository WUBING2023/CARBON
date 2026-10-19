from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_delta_pcc_by_bin(bin_summary: pd.DataFrame, out_path: Path, *, strat: str = "window_phyloP") -> None:
    """跨 task：每个 bin 上 median_delta_pcc 的中位数及四分位距示意。"""
    sub = bin_summary[bin_summary["stratification"] == strat].copy()
    if sub.empty:
        return
    bins = sorted(sub["bin_id"].unique())
    medians = []
    lows = []
    highs = []
    for b in bins:
        chunk = sub.loc[sub["bin_id"] == b, "median_delta_pcc"].astype(float)
        chunk = chunk[np.isfinite(chunk)]
        if chunk.size == 0:
            medians.append(float("nan"))
            lows.append(float("nan"))
            highs.append(float("nan"))
        else:
            medians.append(float(np.median(chunk)))
            lows.append(float(np.percentile(chunk, 25)))
            highs.append(float(np.percentile(chunk, 75)))

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(bins), dtype=np.float64)
    ax.axhline(0.0, color="0.5", linewidth=0.8, linestyle="--")
    ax.plot(x, medians, marker="o", color="#E64B35", label="median of task medians")
    yerr = np.vstack([np.array(medians) - np.array(lows), np.array(highs) - np.array(medians)])
    m = np.isfinite(medians)
    ax.errorbar(x[m], np.array(medians)[m], yerr=yerr[:, m], fmt="none", ecolor="#E64B35", alpha=0.5, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(b)) for b in bins])
    ax.set_xlabel("conservation bin (0=low)")
    ax.set_ylabel("Δ PCC (MARCP − no_conservation)")
    ax.set_title(f"Stratification: {strat}")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_delta_spearman_by_bin(bin_summary: pd.DataFrame, out_path: Path, *, strat: str = "window_phyloP") -> None:
    sub = bin_summary[bin_summary["stratification"] == strat].copy()
    if sub.empty:
        return
    bins = sorted(sub["bin_id"].unique())
    medians = []
    lows = []
    highs = []
    for b in bins:
        chunk = sub.loc[sub["bin_id"] == b, "median_delta_spearman"].astype(float)
        chunk = chunk[np.isfinite(chunk)]
        if chunk.size == 0:
            medians.append(float("nan"))
            lows.append(float("nan"))
            highs.append(float("nan"))
        else:
            medians.append(float(np.median(chunk)))
            lows.append(float(np.percentile(chunk, 25)))
            highs.append(float(np.percentile(chunk, 75)))

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(bins), dtype=np.float64)
    ax.axhline(0.0, color="0.5", linewidth=0.8, linestyle="--")
    ax.plot(x, medians, marker="o", color="#4DBBD5", label="median of task medians")
    yerr = np.vstack([np.array(medians) - np.array(lows), np.array(highs) - np.array(medians)])
    m = np.isfinite(medians)
    ax.errorbar(x[m], np.array(medians)[m], yerr=yerr[:, m], fmt="none", ecolor="#4DBBD5", alpha=0.5, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels([str(int(b)) for b in bins])
    ax.set_xlabel("conservation bin (0=low)")
    ax.set_ylabel("Δ Spearman (MARCP − no_conservation)")
    ax.set_title(f"Stratification: {strat}")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
