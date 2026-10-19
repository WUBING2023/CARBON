#!/usr/bin/env python3
"""
将 tfmodisco_export/queries.meme 与 CISBP-RNA Homo sapiens 子库（合并为单个 MEME）做 Tomtom，
并对 tomtom.tsv 联表 RBP_Information_all_motifs.txt，写出 tomtom_with_rbp.tsv 与 comparison_summary_cisbp.txt。

依赖：MEME Suite 的 `tomtom`（PATH 或 `--tomtom-binary`）。输出目录使用 `-oc` 可覆盖重跑。
RNA 单链比对须关闭互补链搜索：实际调用带 **`-no-ssc`**（以及 `-dist pearson -min-overlap 4`，与手动示例一致）。

  python3 scripts/tomtom_modisco_vs_cisbp.py \\
    --task-dir outputs/AGGF1_HepG2_ENCSR543TPH \\
    --expected-gene AGGF1

--max-targets 0 表示使用 pwms 目录下全部条目。
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from motif_analysis.cisbp_compare import load_all_cisbp_pwms
from motif_analysis.cisbp_tomtom_annotate import annotate_and_export_cisbp
from motif_analysis.meme_write import write_meme_file

_DEFAULT_HS_PWM = ROOT / "reference_data" / "CISBP-RNA" / "homo_sapiens" / "pwms"
_DEFAULT_HS_RBP = (
    ROOT / "reference_data" / "CISBP-RNA" / "homo_sapiens" / "RBP_Information_all_motifs.txt"
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-dir", type=Path, required=True)
    ap.add_argument(
        "--cisbp-pwms",
        type=Path,
        default=_DEFAULT_HS_PWM,
        help=f"默认: {_DEFAULT_HS_PWM}",
    )
    ap.add_argument(
        "--cisbp-rbp-table",
        type=Path,
        default=_DEFAULT_HS_RBP,
        help="CISBP RBP–Motif 展开表（默认 homo_sapiens）",
    )
    ap.add_argument(
        "--max-targets",
        type=int,
        default=0,
        help="靶库最多条数；0=全部",
    )
    ap.add_argument(
        "--queries-meme",
        type=Path,
        default=None,
        help="默认 <task-dir>/tfmodisco_export/queries.meme",
    )
    ap.add_argument(
        "--tomtom-binary",
        type=Path,
        default=None,
        help="tomtom 可执行文件（默认从 PATH 查找）",
    )
    ap.add_argument(
        "--expected-gene",
        type=str,
        default=None,
        help="任务基因符号：RBP_Name 匹配时在该 Query 下置顶（如 RBFOX2）",
    )
    ap.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="comparison_summary_cisbp.txt 中每个 query 列出的命中条数上限",
    )
    ap.add_argument(
        "--skip-annotate",
        action="store_true",
        help="仅运行 tomtom，不写联表与 summary",
    )
    ap.add_argument(
        "--annotate-only",
        type=Path,
        default=None,
        help="跳过 Tomtom；对已有 tomtom.tsv 做 CISBP 联表（路径通常为 .../tomtom_out/tomtom.tsv）",
    )
    args = ap.parse_args()
    task_dir = args.task_dir.resolve()
    q_meme = args.queries_meme or (task_dir / "tfmodisco_export" / "queries.meme")
    q_meme = q_meme.resolve()
    out_base = task_dir / "tomtom_cisbp"
    out_base.mkdir(parents=True, exist_ok=True)
    t_out = out_base / "tomtom_out"
    tsv_path = t_out / "tomtom.tsv"

    if args.annotate_only:
        tsv = args.annotate_only.resolve()
        if not tsv.is_file():
            print(f"[error] {tsv} not found", file=sys.stderr)
            raise SystemExit(1)
        annotate_and_export_cisbp(
            tsv,
            args.cisbp_rbp_table.resolve(),
            t_out,
            expected_gene=args.expected_gene,
            top_n=args.top_n,
            note_lines=[
                "NOTE: --annotate-only（未重新运行 tomtom）。",
            ],
        )
        print(f"[ok] CISBP 注释: {t_out / 'tomtom_with_rbp.tsv'}", flush=True)
        return

    if not q_meme.is_file():
        print(
            f"[error] 缺少 {q_meme}，请先运行 scripts/modisco_pattern_rank.py（推荐）或 scripts/export_modisco_motifs.py",
            file=sys.stderr,
        )
        raise SystemExit(1)

    cisbp = load_all_cisbp_pwms(args.cisbp_pwms.resolve())
    if not cisbp:
        print(f"[error] CISBP 目录无有效 .txt: {args.cisbp_pwms}", file=sys.stderr)
        raise SystemExit(1)
    if args.max_targets and args.max_targets > 0:
        cisbp = cisbp[: args.max_targets]

    target_meme = out_base / "cisbp_targets.meme"
    write_meme_file(target_meme, cisbp)
    print(f"[ok] 靶库 MEME（{len(cisbp)} 条）: {target_meme}", flush=True)

    if args.tomtom_binary:
        exe = str(args.tomtom_binary.resolve())
        if not Path(exe).is_file():
            print(f"[error] not a file: {exe}", file=sys.stderr)
            raise SystemExit(2)
    else:
        exe = shutil.which("tomtom")
    if not exe:
        print(
            "[error] 未找到 tomtom（MEME Suite）。已写出 queries 与 cisbp_targets.meme，可手动执行：\n"
            f"  tomtom -no-ssc -dist pearson -min-overlap 4 \\\n"
            f"    {q_meme} {target_meme} -oc {out_base / 'tomtom_out'}\n"
            "安装说明: https://meme-suite.org/meme/doc/download.html",
            file=sys.stderr,
        )
        raise SystemExit(2)
    t_out.mkdir(parents=True, exist_ok=True)

    cmd = [
        exe,
        "-no-ssc",
        "-dist",
        "pearson",
        "-min-overlap",
        "4",
        str(q_meme),
        str(target_meme),
        "-oc",
        str(t_out),
    ]
    print("[run]", " ".join(cmd), flush=True)
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise SystemExit(r.returncode)
    print(f"[ok] Tomtom 输出: {tsv_path}", flush=True)

    if args.skip_annotate:
        return

    if not tsv_path.is_file():
        print(f"[error] 缺少 {tsv_path}", file=sys.stderr)
        raise SystemExit(3)

    annotate_and_export_cisbp(
        tsv_path,
        args.cisbp_rbp_table.resolve(),
        t_out,
        expected_gene=args.expected_gene,
        top_n=args.top_n,
    )
    print(f"[ok] CISBP 注释: {t_out / 'tomtom_with_rbp.tsv'}", flush=True)


if __name__ == "__main__":
    main()
