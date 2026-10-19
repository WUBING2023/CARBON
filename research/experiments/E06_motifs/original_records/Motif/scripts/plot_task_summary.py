#!/usr/bin/env python3
"""
生成任务汇总 PNG（不写 HTML）。

  cd /data/CARBON_PROJECT/可解释性分析/Motif
  python3 scripts/plot_task_summary.py --task-dir outputs/AGGF1_HepG2_ENCSR543TPH

输出：<task-dir>/figures/*.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from motif_analysis.task_summary_figures import save_all_summary_figures


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=None, help="默认 <task-dir>/figures")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    paths = save_all_summary_figures(args.task_dir, out_dir=args.out_dir, dpi=args.dpi)
    for p in paths:
        print(p)


if __name__ == "__main__":
    main()
