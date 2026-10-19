from __future__ import annotations

import os
from pathlib import Path


def rbp_root() -> Path:
    return Path(os.environ.get("EXP2_RBP_ROOT", os.environ.get("EXP1_RBP_ROOT", "/data/CARBON_PROJECT"))).resolve()


def h5data_root() -> Path:
    return Path(os.environ.get("EXP2_H5DATA", os.environ.get("EXP1_H5DATA", str(rbp_root() / "数据" / "h5data")))).resolve()


def bed_root() -> Path:
    return Path(os.environ.get("EXP2_BED_ROOT", os.environ.get("EXP1_BED_ROOT", str(rbp_root() / "数据" / "new_beddata")))).resolve()


def motif_outputs_root() -> Path:
    return Path(
        os.environ.get(
            "EXP2_MOTIF_OUTPUTS",
            str(rbp_root() / "可解释性分析" / "Motif" / "outputs"),
        )
    ).resolve()


def motif_marcp_root() -> Path:
    """用于 sys.path 导入 motif_analysis。"""
    return Path(os.environ.get("EXP2_MOTIF_MARCP", str(rbp_root() / "可解释性分析" / "Motif"))).resolve()


def default_bed_path(task: str) -> Path:
    return bed_root() / task / f"{task}.101.peak_test.seqBw.treat_ctl.bed"


def default_phylo_bw() -> Path | None:
    p = Path(os.environ.get("EXP2_PHYLOP_BW", str(rbp_root() / "数据" / "phyloP" / "hg38.phyloP100way.bw")))
    return p if p.is_file() else None


def motif_task_dir(task_id: str) -> Path:
    return motif_outputs_root() / task_id


TEST_CHROMOSOMES = ("1", "8", "15")
