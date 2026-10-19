from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np


@dataclass
class SeqletRecord:
    pattern_id: str
    h5_group_path: str
    instance_index: int
    example_idx: int
    start: int
    end: int
    is_revcomp: bool
    seqlet_score: float


def load_pattern_ranking(task_dir: Path) -> list[dict[str, Any]] | None:
    p = task_dir / "tfmodisco_export" / "pattern_ranking.json"
    if not p.is_file():
        return None
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None
    rows = data.get("patterns")
    if not isinstance(rows, list):
        return None
    return [r for r in rows if isinstance(r, dict) and r.get("h5_group_path")]


def iter_seqlets_from_group(
    hf: h5py.File,
    pattern_id: str,
    group_path: str,
) -> list[SeqletRecord]:
    gpath = group_path.rstrip("/")
    if f"{gpath}/seqlets" not in hf:
        return []
    sl = hf[f"{gpath}/seqlets"]
    if not isinstance(sl, h5py.Group):
        return []
    req = ("start", "end", "example_idx")
    if not all(k in sl for k in req):
        return []
    starts = np.asarray(sl["start"][:], dtype=np.int64)
    ends = np.asarray(sl["end"][:], dtype=np.int64)
    examples = np.asarray(sl["example_idx"][:], dtype=np.int64)
    n = int(starts.shape[0])
    if ends.shape[0] != n or examples.shape[0] != n:
        return []
    rev = np.asarray(sl["is_revcomp"][:], dtype=bool) if "is_revcomp" in sl else np.zeros(n, dtype=bool)
    contrib = None
    if "contrib_scores" in sl:
        try:
            contrib = np.asarray(sl["contrib_scores"][:], dtype=np.float64)
        except Exception:
            contrib = None
    out: list[SeqletRecord] = []
    for i in range(n):
        st, en, ex = int(starts[i]), int(ends[i]), int(examples[i])
        rc = bool(rev[i]) if i < rev.size else False
        score = float("nan")
        if contrib is not None and contrib.ndim == 3 and contrib.shape[0] > i:
            score = float(np.sum(np.abs(contrib[i])))
        out.append(
            SeqletRecord(
                pattern_id=pattern_id,
                h5_group_path=gpath,
                instance_index=i,
                example_idx=ex,
                start=st,
                end=en,
                is_revcomp=rc,
                seqlet_score=score,
            )
        )
    return out


def load_all_seqlets(modisco_h5: Path, patterns: list[dict[str, Any]]) -> list[SeqletRecord]:
    modisco_h5 = Path(modisco_h5)
    all_s: list[SeqletRecord] = []
    with h5py.File(modisco_h5, "r") as hf:
        for row in patterns:
            pid = str(row.get("pattern_id", ""))
            gp = str(row.get("h5_group_path", ""))
            if not gp or gp not in hf:
                continue
            all_s.extend(iter_seqlets_from_group(hf, pid, gp))
    return all_s


def load_pattern_contrib_and_sequence(
    modisco_h5: Path,
    group_path: str,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """返回 (contrib_scores Lx4, sequence Lx4) 原始 h5 浮点。"""
    modisco_h5 = Path(modisco_h5)
    with h5py.File(modisco_h5, "r") as hf:
        if group_path not in hf:
            return None, None
        g = hf[group_path]
        if not isinstance(g, h5py.Group):
            return None, None
        c = g["contrib_scores"][()] if "contrib_scores" in g else None
        s = g["sequence"][()] if "sequence" in g else None
        if c is not None:
            c = np.asarray(c, dtype=np.float64)
        if s is not None:
            s = np.asarray(s, dtype=np.float64)
        return c, s
