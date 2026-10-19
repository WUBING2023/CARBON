from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats


def safe_wilcoxon(x: np.ndarray) -> tuple[float, float]:
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size < 6:
        return float("nan"), float("nan")
    try:
        w = stats.wilcoxon(x, zero_method="wilcox", alternative="two-sided", method="auto")
        return float(w.statistic), float(w.pvalue)
    except ValueError:
        return float("nan"), float("nan")


def spearman_trend(bin_index: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """bin_index 单调整数 0..K-1；y 为每 bin 的标量统计。"""
    m = np.isfinite(y) & np.isfinite(bin_index.astype(np.float64))
    if int(m.sum()) < 3:
        return {"spearman_r": float("nan"), "spearman_p": float("nan"), "kendall_tau": float("nan"), "kendall_p": float("nan")}
    bx = bin_index[m].astype(np.float64)
    yy = y[m].astype(np.float64)
    sr = stats.spearmanr(bx, yy, nan_policy="omit")
    kt = stats.kendalltau(bx, yy, nan_policy="omit")
    s_r = getattr(sr, "correlation", None)
    s_p = getattr(sr, "pvalue", None)
    if s_r is None:
        s_r, s_p = sr[0], sr[1]
    k_t = getattr(kt, "statistic", None)
    k_p = getattr(kt, "pvalue", None)
    if k_t is None:
        k_t, k_p = kt[0], kt[1]
    return {
        "spearman_r": float(s_r) if s_r is not None and np.isfinite(s_r) else float("nan"),
        "spearman_p": float(s_p) if s_p is not None and np.isfinite(s_p) else float("nan"),
        "kendall_tau": float(k_t) if k_t is not None and np.isfinite(k_t) else float("nan"),
        "kendall_p": float(k_p) if k_p is not None and np.isfinite(k_p) else float("nan"),
    }


def cross_task_enrichment_vs_gain(task_table: pd.DataFrame) -> dict[str, Any]:
    """
    task_table 列需含 conservation_enrichment, ablation_gain（全为有限值的任务参与）。
    """
    if task_table.empty or not {"conservation_enrichment", "ablation_gain"}.issubset(task_table.columns):
        return {"n_tasks": 0, "spearman_r": float("nan"), "spearman_p": float("nan")}
    t = task_table.dropna(subset=["conservation_enrichment", "ablation_gain"])
    if len(t) < 3:
        return {"n_tasks": int(len(t)), "spearman_r": float("nan"), "spearman_p": float("nan")}
    sr = stats.spearmanr(
        t["conservation_enrichment"].to_numpy(dtype=np.float64),
        t["ablation_gain"].to_numpy(dtype=np.float64),
        nan_policy="omit",
    )
    s_r = getattr(sr, "correlation", None)
    s_p = getattr(sr, "pvalue", None)
    if s_r is None:
        s_r, s_p = sr[0], sr[1]
    return {
        "n_tasks": int(len(t)),
        "spearman_r": float(s_r) if s_r is not None and np.isfinite(s_r) else float("nan"),
        "spearman_p": float(s_p) if s_p is not None and np.isfinite(s_p) else float("nan"),
    }
