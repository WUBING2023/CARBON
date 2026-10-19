from __future__ import annotations

from pathlib import Path


def parse_rbp_cell_line(task_id: str) -> tuple[str, str]:
    """
    从 task 目录名解析 RBP 与 cell line。
    约定：HepG2 或 K562 为细胞系标记，其前一段（可含下划线）为 RBP 名。
    例：AARS_K562_ENCSR825SVO -> (AARS, K562)
        AGGF1_HepG2_ENCSR543TPH -> (AGGF1, HepG2)
    """
    parts = str(task_id).split("_")
    for i, p in enumerate(parts):
        if p in ("HepG2", "K562"):
            rbp = "_".join(parts[:i]) if i > 0 else parts[0]
            return rbp, p
    return task_id, "NA"


def list_tasks(h5data: Path) -> list[str]:
    h5data = Path(h5data)
    if not h5data.is_dir():
        return []
    names: list[str] = []
    for p in sorted(h5data.iterdir()):
        if p.is_dir() and (p / "_5ch").is_dir():
            names.append(p.name)
    return names
