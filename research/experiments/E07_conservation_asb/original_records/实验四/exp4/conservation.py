from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from exp4.data_loader import norm_chrom
from exp4.peaks import find_peaks, peak_segments


def window_mean_phylop(x5: np.ndarray) -> float:
    """x5: (101, 5) 最后一维为 PhyloP。"""
    p = np.asarray(x5[:, 4], dtype=np.float64)
    return float(np.mean(p))


def peak_core_mean_phylop(
    y: np.ndarray,
    phylop: np.ndarray,
    *,
    signal_max_threshold: float = 2.0,
    signal_sum_threshold: float = 4.0,
    signal_length: int = 1,
) -> float:
    """
    无峰时返回 nan；该窗不参与 peak_core 分层（调用方应 mask）。
    若有多个 peak 段，合并所有段上的 PhyloP 均值（段并集上的 mean）。
    """
    yv = np.asarray(y, dtype=np.float64).ravel()
    phy = np.asarray(phylop, dtype=np.float64).ravel()
    mask, pc = find_peaks(
        yv,
        signal_max_threshold=signal_max_threshold,
        signal_sum_threshold=signal_sum_threshold,
        signal_length=signal_length,
    )
    if pc == 0:
        return float("nan")
    segs = peak_segments(mask)
    vals: list[float] = []
    for s, e in segs:
        vals.append(float(np.mean(phy[s : e + 1])))
    return float(np.mean(vals)) if vals else float("nan")


def _clamp_interval(start: int, end: int, lo: int = 0, hi: int = 100) -> tuple[int, int]:
    return max(lo, start), min(hi, end)


def motif_core_mean_phylop(
    x5: np.ndarray,
    core_start: int,
    core_end: int,
) -> float:
    """[core_start, core_end] 闭区间，与窗内 0..100 索引一致。"""
    phy = np.asarray(x5[:, 4], dtype=np.float64)
    s, e = _clamp_interval(int(core_start), int(core_end))
    if s > e:
        return float("nan")
    return float(np.mean(phy[s : e + 1]))


def load_motif_core_table(path: Path) -> pd.DataFrame:
    """
    CSV/TSV：至少含 task_id, chr, chromStart, chromEnd, strand, motif_core_start, motif_core_end
    （后两者为窗内 0-based 闭区间）。
    """
    path = Path(path)
    if path.suffix.lower() in (".parquet", ".pq"):
        try:
            return pd.read_parquet(path)
        except ImportError:
            raise RuntimeError("读取 parquet 需要安装 pyarrow") from None
    sep = "\t" if path.suffix.lower() in (".tsv", ".bed") else ","
    return pd.read_csv(path, sep=sep)


def attach_motif_core_means(
    meta: pd.DataFrame,
    motif_df: pd.DataFrame,
    X: np.ndarray,
) -> np.ndarray:
    """
    meta 与 X 行对齐，X: (N, L, C)；返回 motif_core_mean_phyloP（无匹配为 nan）。
    meta 需含 task_id, chr, chromStart, chromEnd, strand。
    """
    n = len(meta)
    out = np.full(n, np.nan, dtype=np.float64)
    if motif_df.empty or n == 0:
        return out

    m = motif_df.copy()
    m["_k"] = list(
        zip(
            m["task_id"].astype(str),
            m["chr"].map(norm_chrom),
            m["chromStart"].astype(int),
            m["chromEnd"].astype(int),
            m["strand"].astype(str),
        )
    )
    m = m.drop_duplicates(subset=["_k"], keep="first")
    lut = m.set_index("_k")
    keys = list(
        zip(
            meta["task_id"].astype(str),
            meta["chr"].map(norm_chrom),
            meta["chromStart"].astype(int),
            meta["chromEnd"].astype(int),
            meta["strand"].astype(str),
        )
    )
    for i, k in enumerate(keys):
        if k not in lut.index:
            continue
        row = lut.loc[k]
        cs = int(row["motif_core_start"])
        ce = int(row["motif_core_end"])
        out[i] = motif_core_mean_phylop(X[i], cs, ce)
    return out
