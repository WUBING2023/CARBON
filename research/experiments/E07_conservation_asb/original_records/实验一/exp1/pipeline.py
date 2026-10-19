from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from exp1.data_loader import filter_df_indices, load_chr_arrays, read_test_bed_keys
from exp1.features import flank_mean_phylop, gc_content, mean_phylop, summit_window_mean
from exp1.matching import pick_matched_controls
from exp1.paths import TEST_CHROMOSOMES, default_bed_path
from exp1.peaks import find_peaks, peak_segments
from exp1.task_ids import parse_rbp_cell_line


@dataclass
class MetaplotAccumulator:
    half: int = 25
    sum_w: np.ndarray = field(init=False)
    cnt_w: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        w = 2 * self.half + 1
        self.sum_w = np.zeros(w, dtype=np.float64)
        self.cnt_w = np.zeros(w, dtype=np.float64)

    def add_row(self, phylop: np.ndarray, summit: int) -> None:
        if summit < self.half or summit > 100 - self.half:
            return
        for d in range(-self.half, self.half + 1):
            j = d + self.half
            self.sum_w[j] += float(phylop[summit + d])
            self.cnt_w[j] += 1.0

    def curve(self) -> tuple[np.ndarray, np.ndarray]:
        with np.errstate(invalid="ignore"):
            y = np.where(self.cnt_w > 0, self.sum_w / np.maximum(self.cnt_w, 1e-9), np.nan)
        x = np.arange(-self.half, self.half + 1, dtype=np.int32)
        return x, y


def process_task_chrom(
    task_id: str,
    chrom: str,
    h5data: Path,
    *,
    use_bed: bool,
    bed_path: Path | None,
    signal_sum_threshold: float,
    n_controls: int,
    rng: np.random.Generator,
    meta: MetaplotAccumulator,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    """
    返回 (peak 行字典列表, QC 计数)。
    """
    loaded = load_chr_arrays(task_id, chrom, h5data)
    qc: dict[str, float] = {
        "n_rows": 0.0,
        "n_no_peak": 0.0,
        "n_peak_segments": 0.0,
    }
    if loaded is None:
        return [], qc
    df, x_all, y_all, yctl_all = loaded
    test_keys = read_test_bed_keys(bed_path) if (use_bed and bed_path) else set()
    idxs = filter_df_indices(df, test_keys)
    rows_out: list[dict[str, Any]] = []
    rbp, cell_line = parse_rbp_cell_line(task_id)

    for ii in idxs:
        i = int(ii)
        row = df.iloc[i]
        y = y_all[i]
        x = x_all[i]
        yctl = yctl_all[i]
        qc["n_rows"] += 1.0

        peak_mask, pc = find_peaks(
            y,
            signal_max_threshold=2.0,
            signal_sum_threshold=signal_sum_threshold,
            signal_length=1,
        )
        if pc == 0:
            qc["n_no_peak"] += 1.0
            continue

        summit = int(np.argmax(y))
        phylop = np.asarray(x[:, 4], dtype=np.float64)
        meta.add_row(phylop, summit)

        ch = str(row["chr"])
        cstart = int(row["chromStart"])
        cend = int(row["chromEnd"])
        strand = str(row["strand"])
        peak_name = row["peak"] if "peak" in row.index and not pd.isna(row.get("peak")) else ""

        for si, (seg_start, seg_end) in enumerate(peak_segments(peak_mask)):
            qc["n_peak_segments"] += 1.0
            peak_len = seg_end - seg_start + 1
            seg_y = y[seg_start : seg_end + 1]
            peak_max_signal = float(np.max(seg_y))
            peak_sum_signal = float(np.sum(seg_y))
            mean_core = mean_phylop(phylop, seg_start, seg_end)
            mean_fl = flank_mean_phylop(phylop, seg_start, seg_end, flank=20)
            mean_s5 = summit_window_mean(phylop, summit, 5)
            mean_s10 = summit_window_mean(phylop, summit, 10)

            ctrls = pick_matched_controls(
                y, yctl, x, peak_mask, seg_start, seg_end, n_controls=n_controls, rng=rng
            )
            if not ctrls:
                mean_bg = float("nan")
            else:
                bg_means = []
                for cs, ce, _ in ctrls:
                    bg_means.append(mean_phylop(phylop, cs, ce))
                mean_bg = float(np.mean(bg_means))

            gc_peak = gc_content(x[:, :4], seg_start, seg_end)

            peak_id = f"{task_id}_{norm_chrom_key(ch)}_{cstart}_{cend}_{strand}_seg{si}"
            rows_out.append(
                {
                    "task_id": task_id,
                    "rbp": rbp,
                    "cell_line": cell_line,
                    "chrom": chrom,
                    "chr": ch,
                    "chromStart": cstart,
                    "chromEnd": cend,
                    "strand": strand,
                    "peak_row_peak": str(peak_name),
                    "peak_id": peak_id,
                    "row_index": i,
                    "summit_idx": summit,
                    "seg_start": seg_start,
                    "seg_end": seg_end,
                    "peak_len": peak_len,
                    "peak_max_signal": peak_max_signal,
                    "peak_sum_signal": peak_sum_signal,
                    "mean_phyloP_core": mean_core,
                    "mean_phyloP_flank": mean_fl,
                    "mean_phyloP_summit_pm5": mean_s5,
                    "mean_phyloP_summit_pm10": mean_s10,
                    "mean_phyloP_matched_bg": mean_bg,
                    "gc_content_peak": gc_peak,
                    "region_type": "NA",
                    "signal_sum_threshold": signal_sum_threshold,
                }
            )

    return rows_out, qc


def norm_chrom_key(ch: str) -> str:
    return str(ch).replace("chr", "")


def run_task(
    task_id: str,
    h5data: Path,
    *,
    use_bed: bool,
    signal_sum_threshold: float,
    n_controls: int,
    rng: np.random.Generator,
    meta: MetaplotAccumulator,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    bed_path = default_bed_path(task_id) if use_bed else None
    all_rows: list[dict[str, Any]] = []
    qc_tot = {"n_rows": 0.0, "n_no_peak": 0.0, "n_peak_segments": 0.0}
    for chrom in TEST_CHROMOSOMES:
        rows, qc = process_task_chrom(
            task_id,
            chrom,
            h5data,
            use_bed=use_bed,
            bed_path=bed_path,
            signal_sum_threshold=signal_sum_threshold,
            n_controls=n_controls,
            rng=rng,
            meta=meta,
        )
        all_rows.extend(rows)
        for k in qc_tot:
            qc_tot[k] += qc[k]
    return all_rows, qc_tot


def summarize_task_peaks(peak_df: pd.DataFrame) -> dict[str, Any] | None:
    """单 task 层级的 Δ 与 Wilcoxon（core vs matched_bg）。"""
    if peak_df.empty:
        return None
    a = peak_df["mean_phyloP_core"].to_numpy(dtype=np.float64)
    b = peak_df["mean_phyloP_matched_bg"].to_numpy(dtype=np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    if not np.any(m):
        return None
    a, b = a[m], b[m]
    diff = a - b
    out: dict[str, Any] = {
        "n_peaks": int(len(diff)),
        "delta_mean": float(np.mean(diff)),
        "delta_median": float(np.median(diff)),
        "frac_delta_gt0": float(np.mean(diff > 0)),
    }
    try:
        w = stats.wilcoxon(diff, zero_method="wilcox", alternative="two-sided", method="auto")
        out["wilcoxon_statistic"] = float(w.statistic)
        out["wilcoxon_pvalue"] = float(w.pvalue)
    except ValueError:
        out["wilcoxon_statistic"] = float("nan")
        out["wilcoxon_pvalue"] = float("nan")
    return out
