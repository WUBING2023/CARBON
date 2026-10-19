#!/usr/bin/env python3
"""
按 seqlet 支持度对 TF-MoDISco 的 pattern 排序，并导出 JSON/CSV 与供 Tomtom/CISBP 比对的 MEME。

**top_patterns_for_benchmark.meme**：h5 中各 pattern 的 ``sequence``（seqlet 平均 one-hot），
加 ``--meme-pseudocount`` 后行归一化；全长，**不使用** ``contrib_scores``（CWM）冒充 PPM。

**queries.meme**（默认）：在已加 ``--meme-pseudocount`` 并归一化的 PPM 上计算 IC、做 **IC Span 截断**（允许 bipartite 间隔处 IC 凹陷）；
若未回退，**直接截取** ``PPM[start:end]`` 子矩阵写入 MEME（不再二次 pseudocount）；若未找到高 IC 位点或跨度小于 ``--ic-trim-min-positions``，则**回退为全长**
（与 top 中该 motif 一致）。``--no-ic-trim-queries`` 时 queries 与 top 文件内容相同。

trimmed_logos/ 由 `modisco report-simple` 生成，通常对应 CWM 可视化；与 MEME 中 PPM 来源不同。

  python3 scripts/modisco_pattern_rank.py --task-dir outputs/AGGF1_HepG2_ENCSR543TPH
  python3 scripts/modisco_pattern_rank.py --h5 path/to/modisco_results.h5 --include-subpatterns
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from motif_analysis.meme_write import write_meme_file
from motif_analysis.modisco_h5_patterns import (
    PatternRecord,
    collect_pattern_records,
    load_sequence_ppm_for_meme,
    sort_records_by_seqlets,
)
from motif_analysis.ppm_ic_trim import trim_ppm_ic_span


def _records_to_jsonable(records: list[PatternRecord]) -> list[dict]:
    out = []
    for r in records:
        out.append(
            {
                "pattern_id": r.pattern_id,
                "h5_group_path": r.h5_group_path,
                "track": r.track,
                "n_seqlets": r.n_seqlets,
                "contrib_shape": list(r.contrib_shape) if r.contrib_shape else None,
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
说明（trimmed_logos）:
  report_simple/trimmed_logos 下的 PNG 由 MEME modisco 的 report-simple 子命令生成，
  一般每个顶层 pos_patterns/pattern_k 会有 fwd/rev 的 CWM logo；子模式 subpattern_* 是否另有
  独立 trimmed 图取决于 modisco 版本。排序与 seqlet 数以 modisco_results.h5 为准。

IC Span（queries.meme 默认）:
  在归一化 PPM 上逐列算 IC(bits)；取所有 IC>=阈值的行下标，跨度为 [min, max]（含中间低 IC 间隔）。
""",
    )
    ap.add_argument("--task-dir", type=Path, default=None, help="任务目录（默认读 tfmodisco/modisco_results.h5）")
    ap.add_argument("--h5", type=Path, default=None, help="直接指定 modisco_results.h5")
    ap.add_argument("--include-neg", action="store_true", help="包含 neg_patterns")
    ap.add_argument("--include-subpatterns", action="store_true", help="列出 pattern_k/subpattern_j")
    ap.add_argument("--out-dir", type=Path, default=None, help="默认 <task-dir>/tfmodisco_export 或 h5 同目录")
    ap.add_argument("--top-n", type=int, default=10, help="写入 top_patterns_for_benchmark.meme 的前 N 条（0=不写）")
    ap.add_argument("--no-meme", action="store_true", help="不写出 MEME")
    ap.add_argument(
        "--meme-pseudocount",
        type=float,
        default=0.001,
        help="写入 MEME 前对 sequence PPM 每行每碱基加的伪计数，再归一化（默认 0.001；0=仅归一化原始平均 one-hot）",
    )
    ap.add_argument(
        "--no-ic-trim-queries",
        action="store_true",
        help="queries.meme 与 top_patterns_for_benchmark.meme 相同（不做 IC Span 截断）",
    )
    ap.add_argument(
        "--ic-core-min-bits",
        type=float,
        default=0.4,
        help="Span 截断：列 IC(bits) 阈值，默认 0.4（常用 0.3–0.4）",
    )
    ap.add_argument(
        "--ic-trim-min-positions",
        type=int,
        default=6,
        help="Span 宽度小于该值则回退全长，默认 6",
    )
    args = ap.parse_args()

    if args.h5:
        h5_path = args.h5.resolve()
    elif args.task_dir:
        h5_path = (args.task_dir / "tfmodisco" / "modisco_results.h5").resolve()
    else:
        ap.error("需要 --task-dir 或 --h5")

    if not h5_path.is_file():
        print(f"[error] missing {h5_path}", file=sys.stderr)
        raise SystemExit(1)

    out_dir = args.out_dir
    if out_dir is None:
        out_dir = h5_path.parent.parent / "tfmodisco_export" if args.task_dir else h5_path.parent
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    records = collect_pattern_records(
        h5_path,
        include_neg=args.include_neg,
        include_subpatterns=args.include_subpatterns,
    )
    sorted_recs = sort_records_by_seqlets(records, descending=True)

    js = {
        "modisco_results_h5": str(h5_path),
        "n_nodes": len(sorted_recs),
        "sorted_by": "n_seqlets_desc",
        "patterns": _records_to_jsonable(sorted_recs),
        "ic_trim_queries": {
            "enabled": not args.no_ic_trim_queries,
            "ic_core_min_bits": args.ic_core_min_bits,
            "min_positions": args.ic_trim_min_positions,
            "meme_pseudocount_applied": args.meme_pseudocount,
        },
    }
    json_path = out_dir / "pattern_ranking.json"

    csv_path = out_dir / "pattern_ranking.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["pattern_id", "track", "n_seqlets", "h5_group_path", "contrib_shape"],
        )
        w.writeheader()
        for r in sorted_recs:
            w.writerow(
                {
                    "pattern_id": r.pattern_id,
                    "track": r.track,
                    "n_seqlets": r.n_seqlets,
                    "h5_group_path": r.h5_group_path,
                    "contrib_shape": str(r.contrib_shape) if r.contrib_shape else "",
                }
            )

    print(f"[ok] {csv_path}")

    if args.no_meme or args.top_n <= 0:
        json_path.write_text(json.dumps(js, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"[ok] {json_path}")
        return

    motifs_top: list[tuple[str, object]] = []
    motifs_query: list[tuple[str, object]] = []
    trim_by_id: dict[str, dict] = {}

    for r in sorted_recs[: args.top_n]:
        m_full = load_sequence_ppm_for_meme(
            h5_path, r.h5_group_path, pseudocount=args.meme_pseudocount
        )
        if m_full is None:
            print(
                f"[warn] skip {r.pattern_id}: no sequence (L,4) PPM in h5 group",
                file=sys.stderr,
            )
            continue
        motifs_top.append((r.pattern_id, m_full))

        if args.no_ic_trim_queries:
            motifs_query.append((r.pattern_id, m_full.copy()))
            trim_by_id[r.pattern_id] = {
                "ic_trim_start": 0,
                "ic_trim_end": int(m_full.shape[0]),
                "ic_trim_fallback": True,
                "ic_trim_note": "no_ic_trim_queries",
            }
            continue

        m_span, (ts, te), fb = trim_ppm_ic_span(
            m_full,
            ic_core_min=args.ic_core_min_bits,
            min_len=args.ic_trim_min_positions,
        )
        if fb:
            m_q = m_full.copy()
        else:
            m_q = m_span.copy()
        motifs_query.append((r.pattern_id, m_q))
        trim_by_id[r.pattern_id] = {
            "ic_trim_start": ts,
            "ic_trim_end": te,
            "ic_trim_fallback": fb,
        }

    if motifs_top:
        meme_path = out_dir / "top_patterns_for_benchmark.meme"
        write_meme_file(meme_path, motifs_top)
        print(f"[ok] {meme_path} ({len(motifs_top)} motifs)")
        q_path = out_dir / "queries.meme"
        if args.no_ic_trim_queries:
            shutil.copyfile(meme_path, q_path)
            print(f"[ok] {q_path}（与 top 相同，--no-ic-trim-queries）")
        else:
            write_meme_file(q_path, motifs_query)
            print(f"[ok] {q_path} ({len(motifs_query)} motifs，IC Span Tomtom query）")

    for p in js["patterns"]:
        pid = p["pattern_id"]
        if pid in trim_by_id:
            p.update(trim_by_id[pid])

    json_path.write_text(json.dumps(js, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[ok] {json_path}")


if __name__ == "__main__":
    main()
