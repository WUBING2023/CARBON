"""
单 task 工作单元，供 ProcessPoolExecutor 调用（须为模块顶层函数以便 pickle）。
"""
from __future__ import annotations

import os
import zlib
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from exp1.pipeline import MetaplotAccumulator, run_task, summarize_task_peaks


def init_thread_env(threads_per_worker: int) -> None:
    """每个子进程内限制 BLAS/OpenMP 线程数，避免 70×多核过订阅。"""
    s = str(max(1, int(threads_per_worker)))
    for k in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[k] = s


def _task_rng_seed(global_seed: int, task_id: str) -> int:
    h = zlib.adler32(task_id.encode("utf-8"))
    return int((global_seed + h) & 0xFFFFFFFF)


def run_one_task(
    packed: tuple[str, str, bool, float, int, int, int],
) -> dict[str, Any]:
    """
    packed = (task_id, h5data_dir, use_bed, signal_sum_threshold, n_controls, seed, meta_half)
    返回可 pickle 字典：rows, qc, summary, meta_sum, meta_cnt
    """
    task_id, h5_str, use_bed, signal_sum_threshold, n_controls, seed, meta_half = packed
    h5 = Path(h5_str)
    rng = np.random.default_rng(_task_rng_seed(seed, task_id))
    meta = MetaplotAccumulator(half=int(meta_half))
    rows, qc = run_task(
        task_id,
        h5,
        use_bed=use_bed,
        signal_sum_threshold=float(signal_sum_threshold),
        n_controls=int(n_controls),
        rng=rng,
        meta=meta,
    )
    qc["task_id"] = task_id
    pdf = pd.DataFrame(rows)
    summ = summarize_task_peaks(pdf)
    if summ is None:
        summary = {
            "task_id": task_id,
            "n_peaks": 0,
            "delta_mean": float("nan"),
            "delta_median": float("nan"),
            "frac_delta_gt0": float("nan"),
            "wilcoxon_statistic": float("nan"),
            "wilcoxon_pvalue": float("nan"),
        }
    else:
        summary = dict(summ)
        summary["task_id"] = task_id

    return {
        "task_id": task_id,
        "rows": rows,
        "qc": qc,
        "summary": summary,
        "meta_sum": meta.sum_w.copy(),
        "meta_cnt": meta.cnt_w.copy(),
        "meta_half": int(meta_half),
    }
