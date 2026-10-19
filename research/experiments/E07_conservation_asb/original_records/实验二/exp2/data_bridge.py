from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from exp2.genome_coords import norm_chrom


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


def row_key(chrom: str, cstart: int, cend: int, strand: str) -> tuple[str, int, int, str]:
    return (norm_chrom(chrom), int(cstart), int(cend), str(strand))


def resolve_5ch_paths(h5data: Path, task_id: str, h5_mode: str, single_h5: str | None) -> list[Path]:
    h5data = Path(h5data)
    if h5_mode == "all_chr":
        d = h5data / task_id / "_5ch"
        if not d.is_dir():
            return []
        return sorted(d.glob(f"{task_id}_chr*_5ch.h5"))
    if single_h5:
        parts = [Path(p.strip()) for p in str(single_h5).split(",") if p.strip()]
        return [p for p in parts if p.is_file()]
    return []


def h5_concat_cumsum(paths: list[Path]) -> tuple[np.ndarray, list[Path]]:
    sizes: list[int] = []
    for p in paths:
        with h5py.File(p, "r") as hf:
            sizes.append(int(hf["X"].shape[0]))
    cum = np.zeros(len(sizes) + 1, dtype=np.int64)
    if sizes:
        cum[1:] = np.cumsum(np.array(sizes, dtype=np.int64))
    return cum, paths


def global_index_to_file_local(cumsum: np.ndarray, global_idx: int) -> tuple[int, int]:
    import bisect

    fi = bisect.bisect_right(cumsum, int(global_idx)) - 1
    if fi < 0:
        fi = 0
    local = int(global_idx - cumsum[fi])
    return fi, local


def parse_chr_from_5ch_path(path: Path) -> str:
    import re

    m = re.search(r"_chr(\d+)_5ch\.h5$", path.name, re.I)
    return m.group(1) if m else ""


class MetaTableCache:
    """按 5ch 文件缓存 meta df，避免重复读 HDF。"""

    def __init__(self, h5data: Path, task_id: str) -> None:
        self.h5data = Path(h5data)
        self.task_id = task_id
        self._dfs: dict[str, pd.DataFrame] = {}

    def row_for_5ch_row(self, chr5_path: Path, local_row: int) -> pd.Series | None:
        p = Path(chr5_path)
        chr_id = parse_chr_from_5ch_path(p)
        meta_path = self.h5data / self.task_id / f"{self.task_id}_chr{chr_id}.h5"
        key = str(meta_path.resolve())
        if key not in self._dfs:
            if not Path(meta_path).is_file():
                return None
            self._dfs[key] = pd.read_hdf(meta_path, key="df")
        df = self._dfs[key]
        if local_row < 0 or local_row >= len(df):
            return None
        return df.iloc[int(local_row)]


def load_x_phylop_row(chr5_path: Path, local_row: int) -> np.ndarray | None:
    p = Path(chr5_path)
    if not p.is_file():
        return None
    with h5py.File(p, "r") as hf:
        x = hf["X"][int(local_row)]
    arr = np.asarray(x, dtype=np.float32)
    if arr.ndim != 2 or arr.shape[1] < 5:
        return None
    return np.asarray(arr[:, 4], dtype=np.float64)


def load_ohe_row(chr5_path: Path, local_row: int) -> np.ndarray | None:
    p = Path(chr5_path)
    if not p.is_file():
        return None
    with h5py.File(p, "r") as hf:
        x = hf["X"][int(local_row)]
    arr = np.asarray(x, dtype=np.float32)
    if arr.ndim != 2 or arr.shape[1] < 4:
        return None
    return np.asarray(arr[:, :4], dtype=np.float64)
