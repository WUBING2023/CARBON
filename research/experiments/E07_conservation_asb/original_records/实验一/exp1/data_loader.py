from __future__ import annotations

from pathlib import Path
import h5py
import numpy as np
import pandas as pd

from exp1.paths import default_bed_path


def norm_chrom(c: str) -> str:
    return str(c).replace("chr", "")


def read_test_bed_keys(bed_path: Path) -> set[tuple[str, int, int, str]]:
    keys: set[tuple[str, int, int, str]] = set()
    if not bed_path.is_file():
        return keys
    for line in bed_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        chrom = norm_chrom(parts[0])
        start, end = int(parts[1]), int(parts[2])
        strand = parts[5]
        keys.add((chrom, start, end, strand))
    return keys


def load_chr_arrays(
    task: str,
    chrom: str,
    h5data: Path,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray] | None:
    """
    读取单条染色体的元数据 df 与 _5ch 中 X,Y,Y_ctl。
    若文件缺失返回 None。
    """
    c5 = h5data / task / "_5ch" / f"{task}_chr{chrom}_5ch.h5"
    meta_path = h5data / task / f"{task}_chr{chrom}.h5"
    if not c5.is_file() or not meta_path.is_file():
        return None
    df = pd.read_hdf(meta_path, key="df")
    with h5py.File(c5, "r") as hf:
        x = np.asarray(hf["X"], dtype=np.float32)
        y = np.asarray(hf["Y"], dtype=np.float32)
        yctl = np.asarray(hf["Y_ctl"], dtype=np.float32)
    if len(df) != len(y):
        raise ValueError(f"{task} chr{chrom}: df rows {len(df)} != Y rows {len(y)}")
    return df, x, y, yctl


def row_key_series(df: pd.DataFrame) -> pd.Series:
    return pd.Series(
        list(
            zip(
                df["chr"].map(norm_chrom),
                df["chromStart"].astype(int),
                df["chromEnd"].astype(int),
                df["strand"].astype(str),
            )
        ),
        index=df.index,
    )


def filter_df_indices(df: pd.DataFrame, test_keys: set[tuple[str, int, int, str]]) -> np.ndarray:
    if not test_keys:
        return np.arange(len(df), dtype=np.int64)
    ck = row_key_series(df)
    mask = ck.apply(lambda t: t in test_keys)
    return np.flatnonzero(mask.to_numpy())
