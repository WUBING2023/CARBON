#!/usr/bin/env python3
"""
扫描 outputs/<TASK>/ 下 Motif 主流程与 CISBP Tomtom 完成状态（仅检查文件是否存在）。

pipeline：以 tfmodisco/modisco_results.h5 是否存在为准。

  python3 scripts/query_motif_batch_status.py
  python3 scripts/query_motif_batch_status.py --tsv --summary
  python3 scripts/query_motif_batch_status.py --missing cisbp
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SKIP_DIR_NAMES = frozenset(
    {
        "batch_motif_logs",
        "fixtures",
        "demo",
        "workspace",  # 阶段 A 池的临时目录，非 <TASK> 名
        ".git",
    }
)


def _has(p: Path) -> bool:
    return p.is_file()


def scan(outputs_root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not outputs_root.is_dir():
        return rows
    for child in sorted(outputs_root.iterdir()):
        if not child.is_dir():
            continue
        name = child.name
        if name.startswith(".") or name in SKIP_DIR_NAMES:
            continue
        modisco_h5 = child / "tfmodisco" / "modisco_results.h5"
        qm = child / "tfmodisco_export" / "queries.meme"
        cis = child / "tomtom_cisbp" / "tomtom_out" / "tomtom.tsv"
        rows.append(
            {
                "task": name,
                "pipeline_ok": "1" if _has(modisco_h5) else "0",
                "queries_ok": "1" if _has(qm) else "0",
                "cisbp_ok": "1" if _has(cis) else "0",
            }
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="Query per-task Motif + CISBP Tomtom status")
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument("--tsv", action="store_true")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument(
        "--missing",
        type=str,
        choices=("none", "pipeline", "queries", "cisbp"),
        default="none",
    )
    args = ap.parse_args()
    out_root = args.outputs_root.resolve()
    rows = scan(out_root)

    key_map = {
        "pipeline": "pipeline_ok",
        "queries": "queries_ok",
        "cisbp": "cisbp_ok",
    }

    if args.missing != "none":
        k = key_map[args.missing]
        for r in rows:
            if r[k] == "0":
                print(r["task"])
        return

    if args.summary:
        n = len(rows)
        np = sum(1 for r in rows if r["pipeline_ok"] == "1")
        nq = sum(1 for r in rows if r["queries_ok"] == "1")
        nc = sum(1 for r in rows if r["cisbp_ok"] == "1")
        print(f"tasks={n} pipeline={np}/{n} queries={nq}/{n} cisbp={nc}/{n}")
        return

    fields = ["task", "pipeline_ok", "queries_ok", "cisbp_ok"]
    if args.tsv:
        w = csv.DictWriter(sys.stdout, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
        return

    headers = ["task", "pipe", "qry", "cis"]
    widths = [max(len(h), 4) for h in headers]
    for r in rows:
        widths[0] = max(widths[0], len(r["task"]))
    line = f"{'task':<{widths[0]}}  {'pipe':>4}  {'qry':>3}  {'cis':>3}"
    print(line)
    print("-" * len(line))
    for r in rows:
        print(
            f"{r['task']:<{widths[0]}}  {r['pipeline_ok']:>4}  {r['queries_ok']:>3}  {r['cisbp_ok']:>3}"
        )


if __name__ == "__main__":
    main()
