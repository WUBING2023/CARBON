#!/usr/bin/env python3
"""
Analyze paired ablation vs baseline CSV for outliers, metric conflicts,
cross-ablation patterns, leverage points, and training-time quirks.
Reads data/ablation_paired_vs_baseline.csv; writes data/flagged_tasks.csv,
消融异常与特殊模式报告.md, and figures/08_delta_pcc_vs_peak_scatter.png.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update(
    {
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "font.size": 9,
    }
)

NON_BASELINE = ["no_se", "no_multiscale", "no_gru", "no_combined_loss", "no_conservative"]

# Sign conflict: ignore tiny deltas
MIN_ABS_DELTA = 0.01

# Bipolar pattern: strong combined-loss gain + strong GRU removal loss
BIPOLAR_NCL_MIN = 0.05
BIPOLAR_NO_GRU_MAX = -0.05

# Time anomaly: large time shift but tiny PCC change
TIME_ABS_LARGE = 300.0  # seconds
PCC_SMALL_FOR_TIME = 0.01
TIME_MAD_Z = 4.0


def _mad_z(x: np.ndarray) -> np.ndarray:
    """Robust z-score using median and MAD; 0.6745 makes MAD comparable to std for Normal."""
    x = np.asarray(x, dtype=float)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad <= 0:
        return np.zeros_like(x)
    return 0.6745 * (x - med) / mad


def _iqr_bounds(x: np.ndarray) -> tuple[float, float, float, float, float]:
    q1, q3 = np.percentile(x, [25, 75])
    iqr = q3 - q1
    lo = q1 - 1.5 * iqr
    hi = q3 + 1.5 * iqr
    return float(q1), float(q3), float(iqr), float(lo), float(hi)


def collect_iqr_outliers(paired: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    rows = []
    summary = []
    for ab in NON_BASELINE:
        sub = paired[paired["ablation"].astype(str) == ab]["delta_Mean_PCC"].dropna()
        arr = sub.values
        if len(arr) < 5:
            continue
        q1, q3, iqr, lo, hi = _iqr_bounds(arr)
        summary.append(
            {
                "ablation": ab,
                "n": len(arr),
                "Q1": q1,
                "Q3": q3,
                "IQR": iqr,
                "fence_low": lo,
                "fence_high": hi,
            }
        )
        for _, r in paired[paired["ablation"].astype(str) == ab].iterrows():
            d = r["delta_Mean_PCC"]
            if not np.isfinite(d):
                continue
            if d < lo or d > hi:
                rows.append(
                    {
                        "task": r["task"],
                        "ablation": ab,
                        "flag": "IQR_outlier_delta_Mean_PCC",
                        "detail": f"delta={d:.6f} fence=[{lo:.4f},{hi:.4f}]",
                        "delta_Mean_PCC": d,
                    }
                )
    return pd.DataFrame(rows), summary


def collect_sign_conflicts(paired: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in paired.iterrows():
        dp = r["delta_Mean_PCC"]
        ds = r["delta_Spearman"]
        dk = r["delta_Peak_AUC"]
        if not (np.isfinite(dp) and np.isfinite(ds) and np.isfinite(dk)):
            continue
        ab = str(r["ablation"])
        if abs(dp) >= MIN_ABS_DELTA and abs(ds) >= MIN_ABS_DELTA and np.sign(dp) != np.sign(ds):
            rows.append(
                {
                    "task": r["task"],
                    "ablation": ab,
                    "flag": "sign_conflict_PCC_vs_Spearman",
                    "detail": f"dPCC={dp:.5f} dSp={ds:.5f}",
                    "delta_Mean_PCC": dp,
                }
            )
        if abs(dp) >= MIN_ABS_DELTA and abs(dk) >= MIN_ABS_DELTA and np.sign(dp) != np.sign(dk):
            rows.append(
                {
                    "task": r["task"],
                    "ablation": ab,
                    "flag": "sign_conflict_PCC_vs_PeakAUC",
                    "detail": f"dPCC={dp:.5f} dPeakAUC={dk:.5f}",
                    "delta_Mean_PCC": dp,
                }
            )
    return pd.DataFrame(rows)


def pivot_deltas(paired: pd.DataFrame) -> pd.DataFrame:
    sub = paired[paired["ablation"].astype(str).isin(NON_BASELINE)].copy()
    wide = sub.pivot_table(
        index="task",
        columns="ablation",
        values="delta_Mean_PCC",
        aggfunc="first",
    )
    return wide


def cross_ablation_patterns(wide: pd.DataFrame) -> dict[str, Any]:
    cols = NON_BASELINE
    missing = [c for c in cols if c not in wide.columns]
    if missing:
        return {"error": f"missing columns: {missing}"}
    m = wide.reindex(columns=cols).dropna(how="any")
    all_neg = m.apply(lambda row: (row[cols] < 0).all(), axis=1)
    all_pos = m.apply(lambda row: (row[cols] > 0).all(), axis=1)
    ncl = "no_combined_loss"
    gru = "no_gru"
    bipolar = m.apply(
        lambda row: (
            np.isfinite(row.get(ncl, np.nan))
            and np.isfinite(row.get(gru, np.nan))
            and row[ncl] >= BIPOLAR_NCL_MIN
            and row[gru] <= BIPOLAR_NO_GRU_MAX
        ),
        axis=1,
    )
    return {
        "all_five_negative_count": int(all_neg.sum()),
        "all_five_negative_tasks": all_neg[all_neg].index.tolist()[:15],
        "all_five_positive_count": int(all_pos.sum()),
        "all_five_positive_tasks": all_pos[all_pos].index.tolist()[:15],
        "bipolar_ncl_vs_nogru_count": int(bipolar.sum()),
        "bipolar_ncl_vs_nogru_tasks": bipolar[bipolar].index.tolist()[:15],
    }


def no_combined_loss_leverage(paired: pd.DataFrame, topk: int = 15) -> dict[str, Any]:
    sub = paired[paired["ablation"].astype(str) == "no_combined_loss"][
        ["task", "delta_Mean_PCC"]
    ].dropna()
    sub = sub.sort_values("delta_Mean_PCC", ascending=False)
    top = sub.head(topk)
    bot = sub.tail(topk)
    vals = sub["delta_Mean_PCC"].values
    mean_all = float(np.mean(vals))
    trim_k = 5
    trimmed = vals.copy()
    if len(trimmed) > 2 * trim_k:
        trimmed = np.sort(trimmed)[trim_k:-trim_k]
    mean_trim = float(np.mean(trimmed)) if len(trimmed) else mean_all
    mean_wo_top5 = float(np.mean(np.sort(vals)[:-5])) if len(vals) > 5 else mean_all
    return {
        "mean_all": mean_all,
        "median": float(np.median(vals)),
        "mean_trim_inner": mean_trim,
        "mean_after_drop_top5_largest": mean_wo_top5,
        "top_positive_tasks": list(zip(top["task"].tolist(), top["delta_Mean_PCC"].tolist())),
        "bottom_tasks": list(zip(bot["task"].tolist(), bot["delta_Mean_PCC"].tolist())),
    }


def extreme_by_ablation(paired: pd.DataFrame, k: int = 10) -> dict[str, dict[str, Any]]:
    out = {}
    for ab in NON_BASELINE:
        sub = paired[paired["ablation"].astype(str) == ab][["task", "delta_Mean_PCC"]].dropna()
        if sub.empty:
            continue
        sub = sub.sort_values("delta_Mean_PCC")
        out[ab] = {
            "worst_k": list(zip(sub.head(k)["task"].tolist(), sub.head(k)["delta_Mean_PCC"].tolist())),
            "best_k": list(zip(sub.tail(k)["task"].tolist(), sub.tail(k)["delta_Mean_PCC"].tolist())),
        }
    return out


def training_time_flags(paired: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ab in NON_BASELINE:
        sub = paired[paired["ablation"].astype(str) == ab]
        dt = sub["delta_Training_time_sec"].dropna()
        dp = sub["delta_Mean_PCC"].dropna()
        if len(dt) < 5:
            continue
        z = _mad_z(sub["delta_Training_time_sec"].values)
        for i, (_, r) in enumerate(sub.iterrows()):
            dtr = r["delta_Training_time_sec"]
            dpc = r["delta_Mean_PCC"]
            if not np.isfinite(dtr):
                continue
            zv = z[i] if len(z) == len(sub) else np.nan
            large_time = abs(dtr) >= TIME_ABS_LARGE or (
                np.isfinite(zv) and abs(zv) >= TIME_MAD_Z
            )
            small_pcc = np.isfinite(dpc) and abs(dpc) < PCC_SMALL_FOR_TIME
            if large_time and small_pcc:
                rows.append(
                    {
                        "task": r["task"],
                        "ablation": ab,
                        "flag": "time_shift_large_PCC_small",
                        "detail": f"dTrain_sec={dtr:.1f} z_mad~{zv:.2f} dPCC={dpc}",
                        "delta_Mean_PCC": dpc,
                    }
                )
    return pd.DataFrame(rows)


def plot_scatter_with_outliers(
    paired: pd.DataFrame, iqr_flags: pd.DataFrame, out_path: Path
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    axes = axes.flatten()
    out_set = set()
    if not iqr_flags.empty and "task" in iqr_flags.columns and "ablation" in iqr_flags.columns:
        for _, r in iqr_flags.iterrows():
            out_set.add((r["task"], str(r["ablation"])))
    for i, ab in enumerate(NON_BASELINE):
        ax = axes[i]
        sub = paired[paired["ablation"].astype(str) == ab].dropna(
            subset=["delta_Mean_PCC", "delta_Peak_AUC"]
        )
        ax.scatter(
            sub["delta_Peak_AUC"],
            sub["delta_Mean_PCC"],
            s=14,
            alpha=0.45,
            c="#3182bd",
        )
        ax.axhline(0, color="k", linewidth=0.5)
        ax.axvline(0, color="k", linewidth=0.5)
        for _, r in sub.iterrows():
            if (r["task"], ab) in out_set:
                ax.scatter(
                    [r["delta_Peak_AUC"]],
                    [r["delta_Mean_PCC"]],
                    s=40,
                    facecolors="none",
                    edgecolors="#cb181d",
                    linewidths=1.2,
                )
        ax.set_xlabel("Δ Peak_AUC")
        ax.set_ylabel("Δ Mean_PCC")
        ax.set_title(ab)
        ax.grid(True, alpha=0.25)
    fig.suptitle("Δ Mean_PCC vs Δ Peak_AUC (red ring = IQR outlier on Δ Mean_PCC)", y=1.01)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def merge_flags(*dfs: pd.DataFrame) -> pd.DataFrame:
    parts = [d for d in dfs if d is not None and not d.empty]
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True).drop_duplicates(
        subset=["task", "ablation", "flag"], keep="first"
    )


def write_report(
    path: Path,
    paired_display: str,
    iqr_summary: list[dict],
    iqr_outliers: pd.DataFrame,
    conflicts: pd.DataFrame,
    cross: dict[str, Any],
    leverage: dict[str, Any],
    extremes: dict[str, dict],
    time_flags: pd.DataFrame,
    all_flags: pd.DataFrame,
    generated: str,
) -> None:
    lines = [
        "# 全量消融 — 异常与特殊模式分析",
        "",
        f"**生成时间（UTC）**：{generated}",
        "",
        "**数据**：配对表来自",
        f"`{paired_display}`。",
        "",
        "## 1. 方法摘要",
        "",
        "- **离群（主指标）**：对每个非 baseline 消融，在 `delta_Mean_PCC` 上使用 **Tukey 栅栏**：",
        "  下栅栏 = Q1 − 1.5×IQR，上栅栏 = Q3 + 1.5×IQR；落在外侧记为离群。",
        "- **符号冲突**：同一 `(task, ablation)` 上，若 |Δ|≥0.01 且 `delta_Mean_PCC` 与 `delta_Spearman`（或与 `delta_Peak_AUC`）**符号相反**，记为冲突。",
        "- **跨消融模式**：将 5 种消融的 `delta_Mean_PCC` 按 task 对齐；统计「五者全负」「五者全正」；",
        f"  「两极」：`no_combined_loss` Δ≥{BIPOLAR_NCL_MIN} 且同 task 上 `no_gru` Δ≤{BIPOLAR_NO_GRU_MAX}。",
        "- **no_combined_loss 杠杆点**：列出 Δ 最大/最小的任务；并报告去掉最大 5 个 Δ 后的均值。",
        "- **训练时间异常**：|ΔTraining_sec|≥300 或 |MAD 稳健 z|≥4，且 |ΔMean_PCC|<0.01。",
        "",
        "> 离群为统计意义上的稀有，不等于跑错；需核对请查 `logs/MARCP_<task>_<ablation>.log`。",
        "",
        "## 2. IQR 栅栏与离群数（delta_Mean_PCC）",
        "",
        "| ablation | N | Q1 | Q3 | IQR | 下栅栏 | 上栅栏 |",
        "|----------|---|----|----|-----|--------|--------|",
    ]
    for s in iqr_summary:
        lines.append(
            f"| {s['ablation']} | {s['n']} | {s['Q1']:.5f} | {s['Q3']:.5f} | "
            f"{s['IQR']:.5f} | {s['fence_low']:.5f} | {s['fence_high']:.5f} |"
        )
    lines.extend(
        [
            "",
            f"标为离群的行数：**{len(iqr_outliers)}**（见 `data/flagged_tasks.csv` 中 `IQR_outlier_*`）。",
            "",
            "## 3. 指标符号冲突（PCC vs Spearman / PCC vs Peak_AUC）",
            "",
            f"- PCC vs Spearman 冲突行数：**{len(conflicts[conflicts['flag'].str.contains('Spearman', na=False)]) if not conflicts.empty else 0}**",
            f"- PCC vs Peak_AUC 冲突行数：**{len(conflicts[conflicts['flag'].str.contains('PeakAUC', na=False)]) if not conflicts.empty else 0}**",
            "",
        ]
    )
    if not conflicts.empty:
        lines.append("示例（最多 20 行）：")
        lines.append("")
        lines.append("| task | ablation | flag | detail |")
        lines.append("|------|----------|------|--------|")
        for _, r in conflicts.head(20).iterrows():
            lines.append(f"| {r['task']} | {r['ablation']} | {r['flag']} | {r['detail']} |")
        lines.append("")

    lines.extend(
        [
            "## 4. 跨消融模式（同一 task 上 5 种消融的 Δ Mean_PCC）",
            "",
        ]
    )
    if cross:
        lines.append(f"- **五消融全部劣于 baseline（全负）**：{cross.get('all_five_negative_count', 0)} 个任务")
        ex = cross.get("all_five_negative_tasks") or []
        if ex:
            lines.append(f"  - 示例：`{', '.join(ex[:5])}` …")
        lines.append(f"- **五消融全部优于 baseline（全正）**：{cross.get('all_five_positive_count', 0)} 个任务")
        ex2 = cross.get("all_five_positive_tasks") or []
        if ex2:
            lines.append(f"  - 示例：`{', '.join(ex2[:5])}` …")
        lines.append(
            f"- **两极（no_combined_loss 强正且 no_gru 强负）**：{cross.get('bipolar_ncl_vs_nogru_count', 0)} 个任务"
        )
        ex3 = cross.get("bipolar_ncl_vs_nogru_tasks") or []
        if ex3:
            lines.append(f"  - 示例：`{', '.join(ex3[:5])}` …")
        lines.append("")

    lines.extend(
        [
            "## 5. no_combined_loss：均值与「杠杆」任务",
            "",
            f"- 全体均值 Δ Mean_PCC：**{leverage['mean_all']:.6f}**",
            f"- 中位数：**{leverage['median']:.6f}**",
            f"- 去掉最大 5 个 Δ 后的均值：**{leverage['mean_after_drop_top5_largest']:.6f}**",
            f"- 去掉内外各 5 个后的截尾均值（若可算）：**{leverage['mean_trim_inner']:.6f}**",
            "",
            "**Δ 最大的若干任务（推高均值）**：",
            "",
            "| task | Δ Mean_PCC |",
            "|------|------------|",
        ]
    )
    for t, v in leverage.get("top_positive_tasks", [])[:15]:
        lines.append(f"| {t} | {v:.6f} |")
    lines.extend(["", "**Δ 最小的若干任务**：", "", "| task | Δ Mean_PCC |", "|------|------------|"])
    for t, v in leverage.get("bottom_tasks", [])[:15]:
        lines.append(f"| {t} | {v:.6f} |")
    lines.append("")

    lines.extend(["## 6. 各消融下 Δ Mean_PCC 最差 / 最好任务（各 10 个）", ""])
    for ab, d in extremes.items():
        lines.append(f"### {ab}")
        lines.append("")
        lines.append("最差（最负）：")
        for t, v in d.get("worst_k", [])[:10]:
            lines.append(f"- `{t}` — {v:.6f}")
        lines.append("")
        lines.append("最好（最正）：")
        for t, v in d.get("best_k", [])[:10]:
            lines.append(f"- `{t}` — {v:.6f}")
        lines.append("")

    lines.extend(
        [
            "## 7. 训练时间波动大且 PCC 变化小",
            "",
            f"标记行数：**{len(time_flags)}**（详见 CSV `time_shift_large_PCC_small`）。",
            "",
            "## 8. 汇总：flagged_tasks.csv",
            "",
            f"合并去重后共 **{len(all_flags)}** 条记录；完整列表见 `data/flagged_tasks.csv`。",
            "",
            "## 9. 图",
            "",
            "![ΔPCC vs ΔPeakAUC](figures/08_delta_pcc_vs_peak_scatter.png)",
            "",
            "红色空心圈：该 `(task, ablation)` 在 `delta_Mean_PCC` 上为 IQR 离群。",
            "",
            "---",
            "",
            "**复现**：`../venv/bin/python analyze_ablation_anomalies.py`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    root = Path(__file__).resolve().parent
    data_dir = root / "data"
    fig_dir = root / "figures"
    paired_path = data_dir / "ablation_paired_vs_baseline.csv"
    if not paired_path.is_file():
        raise SystemExit(f"Missing {paired_path}; run plot_ablation_summary.py first.")

    paired = pd.read_csv(paired_path)
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    iqr_df, iqr_summary = collect_iqr_outliers(paired)
    conflicts = collect_sign_conflicts(paired)
    wide = pivot_deltas(paired)
    cross = cross_ablation_patterns(wide)
    leverage = no_combined_loss_leverage(paired)
    extremes = extreme_by_ablation(paired)
    time_df = training_time_flags(paired)

    for df in (iqr_df, conflicts, time_df):
        if not df.empty and "flag" in df.columns:
            pass

    all_flags = merge_flags(iqr_df, conflicts, time_df)

    out_csv = data_dir / "flagged_tasks.csv"
    all_flags.to_csv(out_csv, index=False)

    plot_scatter_with_outliers(paired, iqr_df, fig_dir / "08_delta_pcc_vs_peak_scatter.png")

    write_report(
        root / "消融异常与特殊模式报告.md",
        "data/ablation_paired_vs_baseline.csv",
        iqr_summary,
        iqr_df,
        conflicts,
        cross,
        leverage,
        extremes,
        time_df,
        all_flags,
        generated,
    )

    print(f"Wrote {out_csv}, 消融异常与特殊模式报告.md, figures/08_delta_pcc_vs_peak_scatter.png")
    print(f"Total flags: {len(all_flags)}, IQR outliers: {len(iqr_df)}")


if __name__ == "__main__":
    main()
