#!/usr/bin/env python3
"""
Scan test_metrics_summary.txt under 各个具体实验的指标/MARCP, build a tidy table,
compute deltas vs baseline, export CSV, and write figures to figures/.
"""
from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
plt.rcParams.update(
    {
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
    }
)

ABLATION_ORDER = [
    "baseline",
    "no_se",
    "no_multiscale",
    "no_gru",
    "no_combined_loss",
    "no_conservative",
]

METRIC_KEYS = [
    "Mean_PCC",
    "Spearman",
    "Signal_AUC",
    "Signal_AP",
    "Peak_AUC",
    "Peak_AP",
    "Training_time_sec",
    "Testing_time_sec",
]


def _parse_summary_file(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Metric"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            key, val = parts[0].strip(), parts[1].strip()
            try:
                out[key] = float(val)
            except ValueError:
                continue
    return out


def _task_cell_line(task: str) -> str:
    parts = task.split("_")
    if len(parts) >= 2 and parts[1] in ("HepG2", "K562"):
        return parts[1]
    return "other"


def load_runs(metric_base: Path) -> pd.DataFrame:
    rows = []
    pattern = metric_base / "**" / "test_metrics_summary.txt"
    for path in sorted(metric_base.glob("**/test_metrics_summary.txt")):
        rel = path.relative_to(metric_base)
        parts = rel.parts
        if len(parts) < 3:
            continue
        task, ablation = parts[0], parts[1]
        m = _parse_summary_file(path)
        row = {
            "task": task,
            "ablation": ablation,
            "cell_line": _task_cell_line(task),
            "path": str(path),
        }
        for k in METRIC_KEYS:
            row[k] = m.get(k, np.nan)
        rows.append(row)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    cat = pd.Categorical(df["ablation"], categories=ABLATION_ORDER, ordered=True)
    df["ablation"] = cat
    return df.sort_values(["task", "ablation"]).reset_index(drop=True)


def paired_with_baseline(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (task, non-baseline ablation) with baseline metrics and deltas."""
    base = df[df["ablation"].astype(str) == "baseline"][
        ["task"] + METRIC_KEYS
    ].rename(columns={k: f"baseline_{k}" for k in METRIC_KEYS})
    non = df[df["ablation"].astype(str) != "baseline"].copy()
    merged = non.merge(base, on="task", how="inner")
    for k in METRIC_KEYS:
        merged[f"delta_{k}"] = merged[k] - merged[f"baseline_{k}"]
    return merged


def bootstrap_mean_ci(
    x: np.ndarray, n_boot: int = 2000, seed: int = 42, alpha: float = 0.05
) -> tuple[float, float, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    means = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        means[i] = np.mean(rng.choice(x, size=n, replace=True))
    lo, hi = np.percentile(means, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(np.mean(x)), float(lo), float(hi)


def plot_distribution(df: pd.DataFrame, out_dir: Path, metric: str) -> None:
    sub = df[df["ablation"].astype(str).isin(ABLATION_ORDER)].copy()
    sub = sub.dropna(subset=[metric])
    if sub.empty:
        return
    data = [sub[sub["ablation"].astype(str) == a][metric].values for a in ABLATION_ORDER]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    positions = range(1, len(ABLATION_ORDER) + 1)
    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("#7eb0d5")
        patch.set_alpha(0.85)
    ax.set_xticks(positions)
    ax.set_xticklabels(ABLATION_ORDER, rotation=25, ha="right")
    ax.set_ylabel(metric)
    ax.set_title(f"Distribution of {metric} by ablation (N runs per box)")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"01_box_{metric}.png")
    plt.close(fig)


def plot_mean_delta_bars(paired: pd.DataFrame, out_dir: Path, metric_key: str = "Mean_PCC") -> None:
    col = f"delta_{metric_key}"
    abs_order = [a for a in ABLATION_ORDER if a != "baseline"]
    means, los, his = [], [], []
    ns = []
    for a in abs_order:
        sel = paired[paired["ablation"].astype(str) == a][col].dropna()
        ns.append(len(sel))
        m, lo, hi = bootstrap_mean_ci(sel.values)
        means.append(m)
        los.append(m - lo)
        his.append(hi - m)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(abs_order))
    ax.bar(x, means, yerr=[los, his], capsize=4, color="#41ae76", edgecolor="#238b45", alpha=0.9)
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(abs_order, rotation=25, ha="right")
    ax.set_ylabel(f"Δ {metric_key} (ablation − baseline)")
    ax.set_title(
        f"Mean paired Δ vs baseline (95% bootstrap CI); N per bar: {dict(zip(abs_order, ns))}"
    )
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"02_mean_delta_{metric_key}.png")
    plt.close(fig)


def plot_win_rate(paired: pd.DataFrame, out_dir: Path, metric_key: str = "Mean_PCC") -> None:
    col = f"delta_{metric_key}"
    abs_order = [a for a in ABLATION_ORDER if a != "baseline"]
    rates = []
    ns = []
    for a in abs_order:
        sel = paired[paired["ablation"].astype(str) == a][col].dropna()
        n = len(sel)
        ns.append(n)
        if n == 0:
            rates.append(0.0)
        else:
            rates.append(float(np.mean(sel.values > 0)))
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(abs_order))
    ax.bar(x, rates, color="#8856a7", edgecolor="#542788", alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(abs_order, rotation=25, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of tasks with Δ > 0")
    ax.set_title(f"Win rate vs baseline on {metric_key} (paired tasks; N={dict(zip(abs_order, ns))})")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"03_win_rate_{metric_key}.png")
    plt.close(fig)


def plot_baseline_vs_delta(
    paired: pd.DataFrame, out_dir: Path, metric_key: str = "Mean_PCC"
) -> None:
    col = f"delta_{metric_key}"
    bcol = f"baseline_{metric_key}"
    abs_order = [a for a in ABLATION_ORDER if a != "baseline"]
    ncols, nrows = 3, 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 7))
    axes = axes.flatten()
    for i, a in enumerate(abs_order):
        ax = axes[i]
        sub = paired[paired["ablation"].astype(str) == a]
        sub = sub.dropna(subset=[bcol, col])
        ax.scatter(sub[bcol], sub[col], s=12, alpha=0.45, c="#08519c")
        ax.axhline(0, color="k", linewidth=0.6)
        ax.set_title(f"{a} (N={len(sub)})")
        ax.set_xlabel(f"baseline {metric_key}")
        ax.set_ylabel(f"Δ {metric_key}")
        ax.grid(True, alpha=0.25)
    for j in range(len(abs_order), len(axes)):
        fig.delaxes(axes[j])
    fig.suptitle(f"Baseline strength vs paired Δ ({metric_key})", y=1.02)
    fig.tight_layout()
    fig.savefig(out_dir / f"04_baseline_vs_delta_{metric_key}.png")
    plt.close(fig)


def plot_cell_line_delta(
    paired: pd.DataFrame, out_dir: Path, metric_key: str = "Mean_PCC"
) -> None:
    col = f"delta_{metric_key}"
    abs_order = [a for a in ABLATION_ORDER if a != "baseline"]
    lines = ["HepG2", "K562"]
    fig, ax = plt.subplots(figsize=(10, 5))
    width = 0.35
    x = np.arange(len(abs_order))
    for li, line in enumerate(lines):
        means = []
        for a in abs_order:
            sub = paired[
                (paired["ablation"].astype(str) == a) & (paired["cell_line"] == line)
            ][col].dropna()
            means.append(sub.mean() if len(sub) else np.nan)
        offset = (li - 0.5) * width
        ax.bar(x + offset, means, width, label=line, alpha=0.88)
    ax.axhline(0, color="k", linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(abs_order, rotation=25, ha="right")
    ax.set_ylabel(f"Mean Δ {metric_key}")
    ax.set_title(f"Paired mean Δ by cell line (HepG2 vs K562 from task name)")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / f"05_cell_line_delta_{metric_key}.png")
    plt.close(fig)


def plot_metric_correlation_heatmap(df: pd.DataFrame, out_dir: Path) -> None:
    cols = [c for c in METRIC_KEYS if c not in ("Training_time_sec", "Testing_time_sec")]
    sub = df[cols].dropna()
    if sub.shape[0] < 10:
        return
    cm = sub.corr(method="pearson")
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(cm.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_yticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=35, ha="right")
    ax.set_yticklabels(cols)
    for i in range(len(cols)):
        for j in range(len(cols)):
            ax.text(j, i, f"{cm.values[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title("Pearson correlation across runs (all ablations pooled)")
    fig.tight_layout()
    fig.savefig(out_dir / "06_metric_correlation.png")
    plt.close(fig)


def plot_training_time(df: pd.DataFrame, out_dir: Path) -> None:
    col = "Training_time_sec"
    sub = df.dropna(subset=[col])
    if len(sub) < 50 or sub[col].max() <= 0:
        return
    data = [sub[sub["ablation"].astype(str) == a][col].values / 60.0 for a in ABLATION_ORDER]
    if all(len(d) == 0 for d in data):
        return
    fig, ax = plt.subplots(figsize=(10, 4.5))
    positions = range(1, len(ABLATION_ORDER) + 1)
    ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
    ax.set_xticks(positions)
    ax.set_xticklabels(ABLATION_ORDER, rotation=25, ha="right")
    ax.set_ylabel("Training time (minutes)")
    ax.set_title("Training time distribution by ablation (from test_metrics_summary)")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "07_training_time_box.png")
    plt.close(fig)


def print_stats(df: pd.DataFrame, paired: pd.DataFrame) -> str:
    lines = []
    lines.append(f"summary_files: {len(df)}")
    lines.append(f"unique_tasks: {df['task'].nunique()}")
    lines.append("runs_per_ablation:")
    for a in ABLATION_ORDER:
        n = len(df[df['ablation'].astype(str) == a])
        lines.append(f"  {a}: {n}")
    lines.append(f"paired_rows (non-baseline with baseline): {len(paired)}")
    abs_order = [a for a in ABLATION_ORDER if a != "baseline"]
    for a in abs_order:
        sub = paired[paired["ablation"].astype(str) == a]
        lines.append(f"  paired_N[{a}]: {len(sub)}")
    return "\n".join(lines)


def main() -> None:
    root = Path(__file__).resolve().parent
    metric_base = root.parent / "各个具体实验的指标" / "MARCP"
    data_dir = root / "data"
    fig_dir = root / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = load_runs(metric_base)
    if df.empty:
        raise SystemExit(f"No data under {metric_base}")

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    df["generated_at_utc"] = ts

    df.to_csv(data_dir / "ablation_runs.csv", index=False)

    paired = paired_with_baseline(df)
    paired.to_csv(data_dir / "ablation_paired_vs_baseline.csv", index=False)

    stats_text = print_stats(df, paired)
    (data_dir / "parse_stats.txt").write_text(stats_text + "\n", encoding="utf-8")
    print(stats_text)

    # Figures
    for m in ("Mean_PCC", "Peak_AUC", "Spearman"):
        if df[m].notna().any():
            plot_distribution(df, fig_dir, m)

    # Primary paper-style comparisons: Spearman; Mean_PCC retained for reference
    plot_mean_delta_bars(paired, fig_dir, "Spearman")
    plot_win_rate(paired, fig_dir, "Spearman")
    plot_baseline_vs_delta(paired, fig_dir, "Spearman")
    plot_cell_line_delta(paired, fig_dir, "Spearman")
    plot_mean_delta_bars(paired, fig_dir, "Mean_PCC")
    plot_win_rate(paired, fig_dir, "Mean_PCC")
    plot_baseline_vs_delta(paired, fig_dir, "Mean_PCC")
    plot_cell_line_delta(paired, fig_dir, "Mean_PCC")
    plot_metric_correlation_heatmap(df, fig_dir)
    plot_training_time(df, fig_dir)

    meta = {
        "generated_at_utc": ts,
        "metric_base": str(metric_base),
        "n_summary_files": len(df),
        "n_paired_rows": len(paired),
    }
    (data_dir / "run_meta.txt").write_text(
        "\n".join(f"{k}\t{v}" for k, v in meta.items()) + "\n", encoding="utf-8"
    )
    print(f"Wrote CSV and figures under {root}")


if __name__ == "__main__":
    main()
