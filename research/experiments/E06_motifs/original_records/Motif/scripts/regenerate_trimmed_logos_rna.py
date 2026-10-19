#!/usr/bin/env python3
"""
将各任务 tfmodisco/report_simple/trimmed_logos 下 CWM logo 按 RNA 重绘（第 4 碱基 **U** 而非 T）。

依据 modisco_results.h5 中 contrib_scores，trim 逻辑与 MEME modiscolite.report.create_modisco_logos 一致。

  cd 可解释性分析/Motif
  python3 scripts/regenerate_trimmed_logos_rna.py --summary
  python3 scripts/regenerate_trimmed_logos_rna.py
  python3 scripts/regenerate_trimmed_logos_rna.py --task-dir outputs/RBFOX2_K562_ENCSR756CKJ
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SKIP_NAMES = frozenset({"batch_motif_logs", "fixtures", "demo", "workspace", ".git"})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument("--task-dir", type=Path, default=None, help="仅处理该任务目录（含 tfmodisco/modisco_results.h5）")
    ap.add_argument("--trim-threshold", type=float, default=0.3, help="与 modisco report-simple 默认 trim 一致")
    ap.add_argument("--summary", action="store_true", help="只统计将处理的 h5 数量")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    from motif_analysis.rna_cwm_trimmed_logos import write_trimmed_logos_from_h5

    h5_list: list[Path] = []
    if args.task_dir is not None:
        td = args.task_dir.resolve()
        h5 = td / "tfmodisco" / "modisco_results.h5"
        if h5.is_file():
            h5_list.append(h5)
        else:
            print(f"[error] missing {h5}", file=sys.stderr)
            raise SystemExit(1)
    else:
        out_root = args.outputs_root.resolve()
        if not out_root.is_dir():
            print(f"[error] not a directory: {out_root}", file=sys.stderr)
            raise SystemExit(1)
        for child in sorted(out_root.iterdir()):
            if not child.is_dir() or child.name in SKIP_NAMES or child.name.startswith("."):
                continue
            h5 = child / "tfmodisco" / "modisco_results.h5"
            if h5.is_file():
                h5_list.append(h5)

    print(f"[rna_trimmed_logos] tasks_with_h5={len(h5_list)}", flush=True)
    if args.summary:
        return

    n_ok = 0
    n_fail = 0
    for h5 in h5_list:
        task_root = h5.parent.parent
        logo_dir = task_root / "tfmodisco" / "report_simple" / "trimmed_logos"
        try:
            if args.dry_run:
                print(f"  would write: {logo_dir}", flush=True)
            else:
                tags = write_trimmed_logos_from_h5(
                    h5,
                    logo_dir,
                    trim_threshold=args.trim_threshold,
                )
                print(f"  OK {task_root.name}  logos={len(tags)}", flush=True)
            n_ok += 1
        except Exception as e:
            n_fail += 1
            print(f"  FAIL {task_root.name}: {e}", flush=True)

    print(f"[rna_trimmed_logos] done ok={n_ok} fail={n_fail}", flush=True)
    raise SystemExit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
