from __future__ import annotations

import numpy as np

from exp1.features import gc_content


def pick_matched_controls(
    y: np.ndarray,
    yctl: np.ndarray,
    x: np.ndarray,
    peak_mask: np.ndarray,
    seg_start: int,
    seg_end: int,
    *,
    n_controls: int = 3,
    rng: np.random.Generator | None = None,
) -> list[tuple[int, int, float]]:
    """
    在 peak_mask==0 的位置上滑动等长窗口，按 GC / sum(Y) / sum(Y_ctl) 距离选对照。
    返回 [(start, end, score), ...] 最多 n_controls 个；不足则返回更少。
    """
    if rng is None:
        rng = np.random.default_rng(0)

    peak_len = seg_end - seg_start + 1
    L = len(y)
    phylop = x[:, 4]
    x4 = x[:, :4]

    def seg_stats(s: int, e: int) -> tuple[float, float, float, float]:
        gcv = gc_content(x4, s, e)
        ys = float(np.sum(y[s : e + 1]))
        cs = float(np.sum(yctl[s : e + 1]))
        phy = float(np.mean(phylop[s : e + 1]))
        return gcv, ys, cs, phy

    g0, ys0, cs0, _ = seg_stats(seg_start, seg_end)

    candidates: list[tuple[int, float]] = []
    for s in range(0, L - peak_len + 1):
        e = s + peak_len - 1
        if s == seg_start:
            continue
        if np.any(peak_mask[s : e + 1] != 0):
            continue
        g1, ys1, cs1, _ = seg_stats(s, e)
        score = abs(g1 - g0) + abs(ys1 - ys0) / (abs(ys0) + 1.0) + abs(cs1 - cs0) / (abs(cs0) + 1.0)
        candidates.append((s, score))

    if not candidates:
        return []

    candidates.sort(key=lambda t: t[1])
    top = candidates[: max(n_controls * 5, 20)]
    rng.shuffle(top)
    top.sort(key=lambda t: t[1])
    out: list[tuple[int, int, float]] = []
    for s, sc in top[:n_controls]:
        out.append((s, s + peak_len - 1, sc))
    return out
