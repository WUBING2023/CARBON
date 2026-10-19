from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn import metrics

from exp4.paths import marcp_root


def _ensure_marcp_utils():
    root = str(marcp_root().resolve())
    if root not in sys.path:
        sys.path.insert(0, root)
    from MARCP.utils.functions import find_peaks as marcp_find_peaks

    return marcp_find_peaks


def row_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """a,b: (N, L) -> (N,)"""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    am = a.mean(axis=1, keepdims=True)
    bm = b.mean(axis=1, keepdims=True)
    ac = a - am
    bc = b - bm
    num = (ac * bc).sum(axis=1)
    den = np.sqrt(np.maximum((ac * ac).sum(axis=1) * (bc * bc).sum(axis=1), 1e-18))
    r = num / den
    r[~np.isfinite(r)] = np.nan
    return r


def row_spearman(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """逐行 Spearman = 对秩做 Pearson。"""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ra = np.apply_along_axis(stats.rankdata, 1, a)
    rb = np.apply_along_axis(stats.rankdata, 1, b)
    return row_pearson(ra, rb)


def assign_quantile_bins(values: np.ndarray, n_bins: int) -> tuple[np.ndarray, np.ndarray]:
    """
    返回 (bin_id 0..B-1, n_effective_bins)。
    对非有限值赋 -1；qcut duplicates 时实际 bin 数可能 < n_bins。
    """
    v = np.asarray(values, dtype=np.float64)
    n = v.size
    bin_id = np.full(n, -1, dtype=np.int32)
    m = np.isfinite(v)
    if int(m.sum()) < max(10, n_bins * 2):
        return bin_id, 0
    try:
        cats = pd.qcut(v[m], q=n_bins, labels=False, duplicates="drop")
    except ValueError:
        return bin_id, 0
    cats = np.asarray(cats, dtype=np.int32)
    bin_id[m] = cats
    n_eff = int(np.max(bin_id)) + 1 if np.any(bin_id >= 0) else 0
    return bin_id, n_eff


def flattened_signal_peak_ap(
    y_flat: np.ndarray,
    pred_flat: np.ndarray,
) -> tuple[float, float]:
    """
    与 MARCP calculate_batch_metrics_whole 一致：Signal AP（y>median），Peak AP（find_peaks 于 y_flat，阈值 2/2）。
    """
    marcp_find_peaks = _ensure_marcp_utils()
    y_flat = np.asarray(y_flat, dtype=np.float64).ravel()
    pred_flat = np.asarray(pred_flat, dtype=np.float64).ravel()
    if y_flat.size == 0:
        return float("nan"), float("nan")

    sp_whole = stats.spearmanr(pred_flat, y_flat, nan_policy="omit").correlation
    _ = sp_whole

    peak_mark, _ = marcp_find_peaks(y_flat, signal_max_threshold=2, signal_sum_threshold=2)
    peak_mark = np.asarray(peak_mark, dtype=np.int32)
    n_pos = int(peak_mark.sum())
    n_neg = int(peak_mark.size - n_pos)
    if n_pos > 0 and n_neg > 0:
        peak_ap = float(metrics.average_precision_score(y_true=peak_mark, y_score=pred_flat))
    else:
        peak_ap = float("nan")

    med = float(np.median(y_flat))
    signal_mark = (y_flat > med).astype(np.int32)
    n_sp = int(signal_mark.sum())
    n_sn = int(signal_mark.size - n_sp)
    if n_sp > 0 and n_sn > 0:
        signal_ap = float(metrics.average_precision_score(y_true=signal_mark, y_score=pred_flat))
    else:
        signal_ap = float("nan")

    return signal_ap, peak_ap
