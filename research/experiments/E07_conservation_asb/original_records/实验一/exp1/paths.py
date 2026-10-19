from __future__ import annotations

import os
from pathlib import Path


def rbp_root() -> Path:
    return Path(os.environ.get("EXP1_RBP_ROOT", "/data/CARBON_PROJECT")).resolve()


def h5data_root() -> Path:
    return Path(os.environ.get("EXP1_H5DATA", str(rbp_root() / "数据" / "h5data"))).resolve()


def bed_root() -> Path:
    return Path(os.environ.get("EXP1_BED_ROOT", str(rbp_root() / "数据" / "new_beddata"))).resolve()


def default_bed_path(task: str) -> Path:
    return bed_root() / task / f"{task}.101.peak_test.seqBw.treat_ctl.bed"


TEST_CHROMOSOMES = ("1", "8", "15")
