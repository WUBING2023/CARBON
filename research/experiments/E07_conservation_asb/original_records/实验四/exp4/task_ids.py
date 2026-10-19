from __future__ import annotations

from pathlib import Path


def parse_rbp_cell_line(task_id: str) -> tuple[str, str]:
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
