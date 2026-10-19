from __future__ import annotations

import os
from pathlib import Path


def rbp_root() -> Path:
    return Path(os.environ.get("EXP4_RBP_ROOT", "/data/CARBON_PROJECT")).resolve()


def h5data_root() -> Path:
    return Path(os.environ.get("EXP4_H5DATA", str(rbp_root() / "数据" / "h5data"))).resolve()


def bed_root() -> Path:
    return Path(os.environ.get("EXP4_BED_ROOT", str(rbp_root() / "数据" / "new_beddata"))).resolve()


def ablation_weights_root() -> Path:
    return Path(
        os.environ.get(
            "EXP4_MARCP_WEIGHTS_ROOT",
            str(rbp_root() / "消融实验" / "全量消融运行" / "模型权重" / "MARCP"),
        )
    ).resolve()


def marcp_root() -> Path:
    return Path(os.environ.get("EXP4_MARCP_CODE", str(rbp_root() / "MARCP"))).resolve()


def default_bed_path(task: str) -> Path:
    return bed_root() / task / f"{task}.101.peak_test.seqBw.treat_ctl.bed"


def weights_baseline(task: str) -> Path:
    return ablation_weights_root() / task / "baseline" / "model_params.pth"


def weights_no_conservative(task: str) -> Path:
    return ablation_weights_root() / task / "no_conservative" / "model_params.pth"


TEST_CHROMOSOMES = ("1", "8", "15")
