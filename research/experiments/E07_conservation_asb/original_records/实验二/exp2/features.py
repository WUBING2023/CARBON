from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def mean_phylop_slice(phylop: np.ndarray, a: int, b: int) -> float:
    """[a, b) 半开区间。"""
    a = int(max(0, a))
    b = int(min(int(phylop.shape[0]), b))
    if b <= a:
        return float("nan")
    return float(np.mean(phylop[a:b]))


def flank_mean_outside_core(
    phylop: np.ndarray,
    core_start: int,
    core_end: int,
    flank_bp: int,
) -> float:
    """core [cs, ce) 两侧各 flank_bp（截断到窗缘），合并两侧均值。"""
    L = int(phylop.shape[0])
    left_a = max(0, core_start - flank_bp)
    left_b = core_start
    right_a = core_end
    right_b = min(L, core_end + flank_bp)
    vals: list[float] = []
    if left_a < left_b:
        vals.append(mean_phylop_slice(phylop, left_a, left_b))
    if right_a < right_b:
        vals.append(mean_phylop_slice(phylop, right_a, right_b))
    if not vals:
        return float("nan")
    return float(np.mean(vals))


def random_non_core_segment_mean_phylop(
    phylop: np.ndarray,
    core_start: int,
    core_end: int,
    rng: np.random.Generator,
    max_tries: int = 40,
) -> float:
    """
    在窗内取与 core 等长的区间，尽量不与 core 重叠，作为对照 mean PhyloP。
    若窗太短无法避免重叠则退化为任意等长区间。
    """
    L = int(phylop.shape[0])
    cs, ce = int(core_start), int(core_end)
    ln = ce - cs
    if ln <= 0 or ln > L:
        return float("nan")
    for _ in range(max_tries):
        s = int(rng.integers(0, max(1, L - ln + 1)))
        e = s + ln
        if e <= cs or s >= ce:
            return mean_phylop_slice(phylop, s, e)
    s = int(rng.integers(0, max(1, L - ln + 1)))
    return mean_phylop_slice(phylop, s, s + ln)


def ppm_information_content_row(ppm_row: np.ndarray, eps: float = 1e-9) -> float:
    """每位点 IC（bits），RNA 4 字母背景 0.25。"""
    p = np.asarray(ppm_row, dtype=np.float64).ravel()
    p = np.clip(p, eps, 1.0)
    p = p / p.sum()
    return float(np.sum(p * np.log2(p + eps) - p * np.log2(0.25)))


@dataclass
class MetaprofileAccumulator:
    """以 motif 中心对齐的窗内相对坐标 metaprofile。"""
    half: int
    sum_w: np.ndarray = field(init=False)
    cnt_w: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        w = 2 * int(self.half) + 1
        self.sum_w = np.zeros(w, dtype=np.float64)
        self.cnt_w = np.zeros(w, dtype=np.float64)

    def add_instance(self, phylop: np.ndarray, core_start: int, core_end: int) -> None:
        if core_end <= core_start:
            return
        center = 0.5 * (core_start + core_end - 1)
        for j in range(int(core_start), int(core_end)):
            rel = int(round(j - center))
            if abs(rel) > self.half:
                continue
            idx = rel + self.half
            self.sum_w[idx] += float(phylop[j])
            self.cnt_w[idx] += 1.0

    def curve(self) -> tuple[np.ndarray, np.ndarray]:
        with np.errstate(invalid="ignore"):
            y = np.where(self.cnt_w > 0, self.sum_w / np.maximum(self.cnt_w, 1e-9), np.nan)
        x = np.arange(-self.half, self.half + 1, dtype=np.int32)
        return x, y


def mean_abs_cwm_on_core(cwm_L4: np.ndarray, core_start: int, core_end: int, pattern_len: int) -> float:
    """将 pattern 长度 L 与窗内 core 对齐：假设 core 长度 == pattern_len，取逐位 |CWM| 平均。"""
    if cwm_L4 is None or cwm_L4.ndim != 2 or cwm_L4.shape[1] != 4:
        return float("nan")
    Lp = int(cwm_L4.shape[0])
    cs, ce = int(core_start), int(core_end)
    if ce - cs != Lp:
        return float("nan")
    v = np.sum(np.abs(cwm_L4), axis=-1)
    return float(np.mean(v))


def mean_ic_on_core(ppm_L4: np.ndarray, core_start: int, core_end: int, pattern_len: int) -> float:
    if ppm_L4 is None or ppm_L4.ndim != 2 or ppm_L4.shape[1] != 4:
        return float("nan")
    Lp = int(ppm_L4.shape[0])
    cs, ce = int(core_start), int(core_end)
    if ce - cs != Lp:
        return float("nan")
    s = ppm_L4.sum(axis=-1, keepdims=True)
    s = np.where(s < 1e-12, 1.0, s)
    ppm = ppm_L4 / s
    ics = [ppm_information_content_row(ppm[i]) for i in range(Lp)]
    return float(np.mean(ics))
