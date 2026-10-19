"""与 MARCP/MARCP/utils/functions.py find_peaks 语义一致的窗内 peak 标记。"""
from __future__ import annotations

import numpy as np


def find_peaks(
    signal_list: np.ndarray,
    *,
    signal_max_threshold: float = 2.0,
    signal_sum_threshold: float = 4.0,
    signal_length: int = 1,
) -> tuple[np.ndarray, int]:
    """
    连续非零段；段内 max>=signal_max_threshold 或 sum>=signal_sum_threshold 则整段标为 1。
    返回 (text_label 长度与输入相同, peak 段计数)。
    """
    y = np.asarray(signal_list, dtype=np.float64).ravel()
    text_label = np.zeros(len(y), dtype=np.int32)
    peak_count = 0

    sub_signal_start: int | None = None
    sub_signal_max = 0.0
    sub_signal_sum = 0.0

    for base_pos, signal_value in enumerate(y):
        if signal_value != 0:
            if sub_signal_start is None:
                sub_signal_start = base_pos
            sub_signal_max = max(sub_signal_max, float(signal_value))
            sub_signal_sum += float(signal_value)
        elif sub_signal_start is not None:
            sub_signal_end = base_pos - 1
            seg_len = sub_signal_end - sub_signal_start + 1
            if seg_len >= signal_length and (
                sub_signal_max >= signal_max_threshold or sub_signal_sum >= signal_sum_threshold
            ):
                text_label[sub_signal_start : sub_signal_end + 1] = 1
                peak_count += 1
            sub_signal_start = None
            sub_signal_max = 0.0
            sub_signal_sum = 0.0

    if sub_signal_start is not None:
        sub_signal_end = len(y) - 1
        seg_len = sub_signal_end - sub_signal_start + 1
        if seg_len >= signal_length and (
            sub_signal_max >= signal_max_threshold or sub_signal_sum >= signal_sum_threshold
        ):
            text_label[sub_signal_start : sub_signal_end + 1] = 1
            peak_count += 1

    return text_label, peak_count


def peak_segments(mask: np.ndarray) -> list[tuple[int, int]]:
    """mask 为 0/1 时返回 [(start, end), ...] 闭区间。"""
    m = np.asarray(mask, dtype=np.int32).ravel()
    out: list[tuple[int, int]] = []
    i = 0
    n = len(m)
    while i < n:
        if m[i] != 1:
            i += 1
            continue
        j = i
        while j < n and m[j] == 1:
            j += 1
        out.append((i, j - 1))
        i = j
    return out
