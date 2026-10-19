from __future__ import annotations

import os
import zlib
from pathlib import Path
from typing import Any

import numpy as np

from exp2.pipeline import run_task_exp2


def init_thread_env(threads_per_worker: int) -> None:
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
    packed: tuple[str, str, str, bool, int, int, int, int],
) -> dict[str, Any]:
    """
    packed = (task_id, h5data_dir, motif_task_dir, use_bed, flank_bp, meta_half, min_inst, seed)
    """
    task_id, h5_str, motif_dir, use_bed, flank_bp, meta_half, min_inst, seed = packed
    h5 = Path(h5_str)
    mdir = Path(motif_dir)
    rng = np.random.default_rng(_task_rng_seed(int(seed), task_id))
    return run_task_exp2(
        task_id,
        mdir,
        h5,
        use_bed=bool(use_bed),
        flank_bp=int(flank_bp),
        meta_half=int(meta_half),
        min_instances_per_motif=int(min_inst),
        rng=rng,
    )


def task_has_modisco(motif_task_dir: Path) -> bool:
    p = Path(motif_task_dir) / "tfmodisco" / "modisco_results.h5"
    q = Path(motif_task_dir) / "predictions" / "sample_indices.npy"
    r = Path(motif_task_dir) / "attributions" / "ig_hyp_4ch.npz"
    return p.is_file() and q.is_file() and r.is_file()
