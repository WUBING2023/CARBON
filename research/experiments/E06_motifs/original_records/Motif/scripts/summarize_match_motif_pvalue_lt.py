#!/usr/bin/env python3
"""
读取 reference_data/motif_p_value/match_motif_p_sorted.csv，
统计 p_value < 阈值时的行数、唯一 (task, Query_ID) 数、唯一 cisbp_RBP_Name 数，
并写出 summary_p_lt_0.05.txt（默认阈值 0.05，即 Tomtom 常规显著性）。

  python3 scripts/summarize_match_motif_pvalue_lt.py
  python3 scripts/summarize_match_motif_pvalue_lt.py --threshold 0.05
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CSV = ROOT / "reference_data/motif_p_value/match_motif_p_sorted.csv"
DEFAULT_OUT = ROOT / "reference_data/motif_p_value/summary_p_lt_0.05.txt"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--threshold", type=float, default=0.05)
    args = ap.parse_args()
    t = args.threshold

    rows: list[dict[str, str]] = []
    with args.csv.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)

    def pval(r: dict[str, str]) -> float:
        return float(r["p_value"])

    sel = [r for r in rows if pval(r) < t]
    motifs = sorted({(r["task"], r["Query_ID"]) for r in sel})
    prots = sorted({r["cisbp_RBP_Name"] for r in sel})

    lines = [
        "match_motif_p_sorted.csv 汇总（cisbp_priority == match_task_gene 的 Tomtom 行）",
        f"阈值：p_value < {t}",
        "",
        f"匹配行总数（数据行）：{len(rows)}",
        f"p < {t} 的行数：{len(sel)}",
        "",
        f"p < {t} 的不重复 motif 数（task + Query_ID 组合）：{len(motifs)}",
        f"p < {t} 的不重复蛋白质数（cisbp_RBP_Name）：{len(prots)}",
        "",
        "口径说明",
        "- motif：每条分析任务下的 TF-Modisco 查询模式，用 (task, Query_ID) 去重。",
        "- 蛋白质：CISBP 侧 RBP 名称（cisbp_RBP_Name），与任务目标基因一致。",
        "",
        f"p < {t} 涉及的蛋白质（按字母序）：",
        ", ".join(prots),
        "",
        f"p < {t} 的不重复 motif（task, Query_ID）：",
    ]
    for task, qid in motifs:
        lines.append(f"- {task}, {qid}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
