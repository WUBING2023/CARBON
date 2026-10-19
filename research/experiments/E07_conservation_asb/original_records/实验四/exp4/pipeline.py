from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from exp4.data_loader import filter_df_indices, load_chr_arrays, read_test_bed_keys
from exp4.paths import TEST_CHROMOSOMES, default_bed_path


def collect_task_windows(
    task_id: str,
    h5data: Path,
    *,
    use_bed: bool,
    chroms: tuple[str, ...] = TEST_CHROMOSOMES,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray] | None:
    """
    合并 test 染色体上的 test 窗（可选 BED 子集）。
    返回 (meta, X, Y)；meta 含原 df 列 + task_id + chrom。
    """
    bed_path = default_bed_path(task_id) if use_bed else None
    test_keys = read_test_bed_keys(bed_path) if (use_bed and bed_path) else set()

    parts_meta: list[pd.DataFrame] = []
    parts_x: list[np.ndarray] = []
    parts_y: list[np.ndarray] = []

    for chrom in chroms:
        loaded = load_chr_arrays(task_id, chrom, h5data)
        if loaded is None:
            continue
        df, x_all, y_all, _yctl = loaded
        idxs = filter_df_indices(df, test_keys)
        if idxs.size == 0:
            continue
        dfp = df.iloc[idxs].copy()
        dfp["task_id"] = task_id
        dfp["chrom"] = str(chrom)
        parts_meta.append(dfp.reset_index(drop=True))
        parts_x.append(np.asarray(x_all[idxs], dtype=np.float32))
        parts_y.append(np.asarray(y_all[idxs], dtype=np.float32))

    if not parts_meta:
        return None
    meta = pd.concat(parts_meta, ignore_index=True)
    X = np.concatenate(parts_x, axis=0)
    Y = np.concatenate(parts_y, axis=0)
    return meta, X, Y
