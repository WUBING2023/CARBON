#!/usr/bin/env python3
"""
扫描 outputs/<TASK>/tomtom_cisbp/tomtom_out/tomtom_with_rbp.tsv，
统计 cisbp_priority != match_task_gene（Tomtom 命中 CISBP 但 RBP 名称与任务基因不一致）的汇总数据。
不绘图。

  python3 scripts/stats_non_match_task_gene.py
  python3 scripts/stats_non_match_task_gene.py --out-txt reference_data/motif_p_value/stats_non_match_task_gene.txt
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SKIP_DIRS = frozenset({"batch_motif_logs", "fixtures", "demo", "workspace", ".git"})
MATCH_PRIORITY = "match_task_gene"


def _parse_float(s: str) -> float:
    s = (s or "").strip()
    if not s or s.lower() in ("nan", "na", ""):
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def _p_from_row(row: dict[str, str]) -> float:
    for k in ("_p_num", "p-value", "P-value", "p_value", "_p"):
        if k in row and row[k] is not None:
            v = _parse_float(str(row[k]))
            if v == v:
                return v
    return float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument(
        "--out-txt",
        type=Path,
        default=ROOT / "reference_data" / "motif_p_value" / "stats_non_match_task_gene.txt",
    )
    ap.add_argument(
        "--out-tasks-zero",
        type=Path,
        default=ROOT / "reference_data" / "motif_p_value" / "tasks_zero_match_task_gene.txt",
        help="列出整份 Tomtom 中没有任何 match_task_gene 行的任务名（每行一个）",
    )
    ap.add_argument("--sig-p", type=float, default=0.05)
    args = ap.parse_args()

    out_root = args.outputs_root.resolve()
    sig = args.sig_p

    total_rows = 0
    match_rows = 0
    non_rows = 0
    non_p_lt = 0

    tasks_with_tsv = 0
    tasks_zero_match: list[str] = []
    tasks_with_match: list[str] = []

    per_task_match: dict[str, int] = {}
    per_task_non: dict[str, int] = {}

    pattern_all: set[tuple[str, str]] = set()
    pattern_has_match: set[tuple[str, str]] = set()
    non_pattern_p_lt: set[tuple[str, str]] = set()

    non_rbp_names: set[str] = set()

    for child in sorted(out_root.iterdir()):
        if not child.is_dir() or child.name in SKIP_DIRS or child.name.startswith("."):
            continue
        tsv = child / "tomtom_cisbp" / "tomtom_out" / "tomtom_with_rbp.tsv"
        if not tsv.is_file():
            continue
        tasks_with_tsv += 1
        task = child.name
        tm = tn = 0
        with open(tsv, encoding="utf-8", errors="replace", newline="") as f:
            rd = csv.DictReader(f, delimiter="\t")
            fieldnames = rd.fieldnames or []
            pri_col = next((h for h in fieldnames if h and str(h).strip() == "cisbp_priority"), None)
            if pri_col is None:
                print(f"[warn] skip {tsv}: no cisbp_priority", file=sys.stderr)
                continue
            for row in rd:
                total_rows += 1
                qid = (row.get("Query_ID") or "").strip()
                pattern_all.add((task, qid))
                pri = (row.get(pri_col) or "").strip()
                pv = _p_from_row(row)
                if pri == MATCH_PRIORITY:
                    match_rows += 1
                    tm += 1
                    pattern_has_match.add((task, qid))
                else:
                    non_rows += 1
                    tn += 1
                    name = (row.get("cisbp_RBP_Name") or "").strip()
                    if name:
                        non_rbp_names.add(name)
                    if pv == pv and pv < sig:
                        non_p_lt += 1
                        non_pattern_p_lt.add((task, qid))
        per_task_match[task] = tm
        per_task_non[task] = tn
        if tm == 0:
            tasks_zero_match.append(task)
        else:
            tasks_with_match.append(task)

    tasks_zero_match.sort()
    tasks_with_match.sort()
    patterns_no_match = len(pattern_all - pattern_has_match)

    lines = [
        "Tomtom + CISBP 联表统计：未按「任务基因」匹配（cisbp_priority != match_task_gene）",
        "",
        "说明：match_task_gene 表示 cisbp_RBP_Name 与任务期望基因（run_pipeline 传入的 expected_gene）一致；",
        "       其余行记为「其他」RBP 的 CISBP motif 命中，即未按本任务的蛋白质名称匹配。",
        "",
        f"扫描目录：{out_root}",
        f"含 tomtom_with_rbp.tsv 的任务数：{tasks_with_tsv}",
        "",
        "【全表行数】",
        f"  Tomtom 联表总行数（所有命中）：{total_rows}",
        f"  match_task_gene 行数：{match_rows}",
        f'  非任务基因匹配（"other"）行数：{non_rows}',
        "",
        f"【显著性 p < {sig}，且为非任务基因匹配】",
        f"  行数：{non_p_lt}",
        f"  不重复 (task, Query_ID) 数：{len(non_pattern_p_lt)}",
        "",
        "【按任务】",
        f"  至少存在 1 条 match_task_gene 的任务数：{len(tasks_with_match)}",
        f"  整份 Tomtom 中没有任何 match_task_gene 的任务数：{len(tasks_zero_match)}",
        "",
        "【按 pattern（task + Query_ID）】",
        f"  不重复 pattern 总数：{len(pattern_all)}",
        f"  至少有一条 match_task_gene 的 pattern 数：{len(pattern_has_match)}",
        f"  完全没有 match_task_gene 的 pattern 数：{patterns_no_match}",
        "",
        "【非任务基因匹配行中的 CISBP RBP 名称】",
        f"  不重复 cisbp_RBP_Name 数：{len(non_rbp_names)}",
        "",
    ]

    args.out_txt.parent.mkdir(parents=True, exist_ok=True)
    args.out_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")

    args.out_tasks_zero.parent.mkdir(parents=True, exist_ok=True)
    args.out_tasks_zero.write_text("\n".join(tasks_zero_match) + "\n", encoding="utf-8")

    print(f"Wrote {args.out_txt}")
    print(f"Wrote {args.out_tasks_zero} ({len(tasks_zero_match)} tasks)")


if __name__ == "__main__":
    main()
