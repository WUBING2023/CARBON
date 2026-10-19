from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats


def instance_core_vs_flank_wilcoxon(inst: pd.DataFrame) -> dict[str, float]:
    if inst.empty or "mean_phyloP_core" not in inst.columns:
        return {"statistic": float("nan"), "pvalue": float("nan"), "n": 0.0}
    a = inst["mean_phyloP_core"].to_numpy(dtype=np.float64)
    b = inst["mean_phyloP_flank"].to_numpy(dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if a.size < 3:
        return {"statistic": float("nan"), "pvalue": float("nan"), "n": float(a.size)}
    try:
        w = stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided", method="auto")
        return {"statistic": float(w.statistic), "pvalue": float(w.pvalue), "n": float(a.size)}
    except ValueError:
        return {"statistic": float("nan"), "pvalue": float("nan"), "n": float(a.size)}


def motif_ci_summary(inst: pd.DataFrame) -> pd.DataFrame:
    if inst.empty:
        return pd.DataFrame()
    inst = inst.copy()
    inst["CI"] = inst["mean_phyloP_core"] - inst["mean_phyloP_flank"]
    g = inst.groupby("motif_id", as_index=False).agg(
        n_instances=("instance_id", "count"),
        mean_CI=("CI", "mean"),
        median_CI=("CI", "median"),
    )
    return g


def position_spearman(position_df: pd.DataFrame) -> dict[str, Any]:
    if position_df.empty or len(position_df) < 4:
        return {"rho_cwm": float("nan"), "p_cwm": float("nan"), "rho_ic": float("nan"), "p_ic": float("nan")}
    x = position_df["abs_cwm"].to_numpy(dtype=np.float64)
    y = position_df["mean_phyloP"].to_numpy(dtype=np.float64)
    ic = position_df["ppm_ic"].to_numpy(dtype=np.float64)
    r1, p1 = stats.spearmanr(x, y, nan_policy="omit")
    r2, p2 = stats.spearmanr(ic, y, nan_policy="omit")
    return {
        "rho_cwm": float(r1) if r1 == r1 else float("nan"),
        "p_cwm": float(p1) if p1 == p1 else float("nan"),
        "rho_ic": float(r2) if r2 == r2 else float("nan"),
        "p_ic": float(p2) if p2 == p2 else float("nan"),
    }


def permutation_spearman_labels(
    position_df: pd.DataFrame,
    *,
    n_perm: int,
    seed: int,
) -> dict[str, float]:
    """打乱 rel_pos 与 mean_phyloP 的配对（保持 marginal），检验 Spearman(abs_cwm, phyloP)。"""
    if position_df.empty or n_perm <= 0:
        return {"perm_p_cwm": float("nan")}
    df = position_df.copy()
    rng = np.random.default_rng(seed)
    x = df["abs_cwm"].to_numpy(dtype=np.float64)
    y = df["mean_phyloP"].to_numpy(dtype=np.float64)
    r0, _ = stats.spearmanr(x, y, nan_policy="omit")
    r0 = float(r0) if r0 == r0 else 0.0
    cnt = 0
    for _ in range(n_perm):
        yp = rng.permutation(y)
        r, _ = stats.spearmanr(x, yp, nan_policy="omit")
        r = float(r) if r == r else 0.0
        if abs(r) >= abs(r0):
            cnt += 1
    return {"perm_p_cwm": float((cnt + 1) / (n_perm + 1))}


def cross_task_ci_stats(ci_values: np.ndarray) -> dict[str, float]:
    ci_values = ci_values[np.isfinite(ci_values)]
    if ci_values.size == 0:
        return {}
    try:
        w = stats.wilcoxon(ci_values, zero_method="wilcox", alternative="two-sided", method="auto")
        wp, ws = float(w.pvalue), float(w.statistic)
    except ValueError:
        wp, ws = float("nan"), float("nan")
    return {
        "n": float(ci_values.size),
        "median_CI": float(np.median(ci_values)),
        "mean_CI": float(np.mean(ci_values)),
        "wilcoxon_vs0_pvalue": wp,
        "wilcoxon_vs0_statistic": ws,
    }
