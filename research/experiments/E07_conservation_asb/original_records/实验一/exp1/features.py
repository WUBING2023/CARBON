from __future__ import annotations

import numpy as np


def onehot_to_seq(x4: np.ndarray) -> list[str]:
    """x4 shape (L,4) A,C,G,U -> 碱基字符。"""
    m = np.asarray(x4, dtype=np.float32)
    idx = np.argmax(m, axis=1)
    mp = np.array(["A", "C", "G", "U"])
    return [str(mp[i]) for i in idx]


def gc_content(x4: np.ndarray, start: int, end: int) -> float:
    """[start, end] 闭区间。"""
    sl = onehot_to_seq(x4[start : end + 1])
    if not sl:
        return float("nan")
    gc = sum(1 for b in sl if b in ("G", "C"))
    return gc / len(sl)


def mean_phylop(phylop: np.ndarray, start: int, end: int) -> float:
    if start > end:
        return float("nan")
    return float(np.mean(phylop[start : end + 1]))


def clamp_interval(start: int, end: int, lo: int = 0, hi: int = 100) -> tuple[int, int]:
    return max(lo, start), min(hi, end)


def flank_mean_phylop(phylop: np.ndarray, seg_start: int, seg_end: int, flank: int = 20) -> float:
    """peak 段两侧各 flank nt（窗内截断），合并两侧均值；缺失侧跳过。"""
    left_s, left_e = clamp_interval(seg_start - flank, seg_start - 1)
    right_s, right_e = clamp_interval(seg_end + 1, seg_end + flank)
    vals: list[float] = []
    if left_s <= left_e:
        vals.append(mean_phylop(phylop, left_s, left_e))
    if right_s <= right_e:
        vals.append(mean_phylop(phylop, right_s, right_e))
    if not vals:
        return float("nan")
    return float(np.mean(vals))


def summit_window_mean(phylop: np.ndarray, summit: int, half: int) -> float:
    s, e = clamp_interval(summit - half, summit + half)
    return mean_phylop(phylop, s, e)
