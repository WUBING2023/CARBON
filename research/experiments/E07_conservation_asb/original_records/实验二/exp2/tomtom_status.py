from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


def best_tomtom_q_by_query(task_dir: Path) -> dict[str, float]:
    """
    解析 Tomtom 输出，对每个 Query_ID 取最小 q-value。
    优先 tomtom_with_rbp.tsv，其次 tomtom.tsv。
    """
    out: dict[str, float] = {}
    for rel in (
        "tomtom_cisbp/tomtom_out/tomtom_with_rbp.tsv",
        "tomtom_cisbp/tomtom_out/tomtom.tsv",
    ):
        p = task_dir / rel
        if not p.is_file():
            continue
        try:
            df = pd.read_csv(p, sep="\t")
        except Exception:
            continue
        if "Query_ID" not in df.columns or "q-value" not in df.columns:
            continue
        g = df.groupby("Query_ID")["q-value"].min()
        for qid, qv in g.items():
            try:
                out[str(qid)] = float(qv)
            except (TypeError, ValueError):
                continue
        if out:
            break
    return out


def motif_match_status(q: float | None, thr: float = 0.05) -> str:
    if q is None or not (q == q):
        return "no_tomtom"
    if q < thr:
        return f"q_lt_{thr:g}"
    return "ns"


def annotate_patterns(
    patterns: list[dict[str, Any]],
    qmap: dict[str, float],
) -> list[dict[str, Any]]:
    rows = []
    for row in patterns:
        pid = str(row.get("pattern_id", ""))
        q = qmap.get(pid)
        rows.append({**row, "tomtom_min_q": q, "motif_match_status": motif_match_status(q)})
    return rows
