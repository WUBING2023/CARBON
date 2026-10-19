#!/usr/bin/env python3
"""
按「消融实验」目录结构一键跑 Motif 流水线：

  权重: /data/CARBON_PROJECT/消融实验/全量消融运行/模型权重/MARCP/<任务名>/baseline/model_params.pth
  数据: /data/CARBON_PROJECT/数据/h5data/<任务名>/_5ch/<任务名>_chr<chr>_5ch.h5

示例:
  python scripts/run_from_ablation_task.py --task AGGF1_HepG2_ENCSR543TPH --chr 22 --max_samples 256
  python scripts/run_from_ablation_task.py --task AGGF1_HepG2_ENCSR543TPH --all-chromosomes --select_top 10000
  python scripts/run_from_ablation_task.py --task FUS_K562_ENCSR069EVH --h5 /path/to/custom.h5
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from motif_analysis.user_paths import baseline_checkpoint, default_chr_h5


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", type=str, required=True, help="如 AGGF1_HepG2_ENCSR543TPH")
    ap.add_argument(
        "--all-chromosomes",
        action="store_true",
        help="合并该任务 _5ch 下全部 chr*_5ch.h5（等价 run_pipeline --h5_mode all_chr）",
    )
    ap.add_argument("--chr", type=str, default="22", help="单染色体模式：使用 {task}_chr{chr}_5ch.h5")
    ap.add_argument("--h5", type=str, default=None, help="覆盖：任意 5ch h5 路径（可逗号分隔多文件）")
    ap.add_argument(
        "--out",
        type=str,
        default=None,
        help="输出目录，默认 outputs/<task>/（实验级目录）",
    )
    ap.add_argument("--max_samples", type=int, default=None, help="合并后最多样本数；默认不截断")
    ap.add_argument("--select_top", type=int, default=10000, help="归因前 Top-K；<=0 为全样本")
    ap.add_argument(
        "--rank_by",
        type=str,
        default="pred_mean",
        choices=["pred_mean", "softmax_mean", "spearman_y"],
        help="Top-K 筛选指标",
    )
    ap.add_argument("--batch_size", type=int, default=16)
    ap.add_argument("--ig_steps", type=int, default=32)
    ap.add_argument("--skip_conv", action="store_true", help="跳过卷积核 Motif")
    ap.add_argument("--skip_meme", action="store_true")
    ap.add_argument("--skip_ism", action="store_true")
    ap.add_argument("--skip_tfmodisco", action="store_true", help="仅写 attributions，不调用 modisco CLI")
    ap.add_argument("--legacy_kmeans", action="store_true", help="额外跑 KMeans seqlet 到 seqlet_motifs_legacy/")
    ap.add_argument(
        "--forward-chunk-rows",
        type=int,
        default=None,
        help="run_pipeline --forward_chunk_rows（all_chr 分块 Top-K）",
    )
    ap.add_argument(
        "--no-chunked-rank",
        action="store_true",
        help="run_pipeline --no_chunked_rank",
    )
    ap.add_argument("--max-seqlets", type=int, default=None, help="run_pipeline --max_seqlets")
    ap.add_argument(
        "--attr_target",
        type=str,
        default="binding_topk",
        choices=[
            "binding_topk",
            "binding_sum",
            "binding_max",
            "logits_topk",
            "logits_sum",
            "logits_max",
        ],
        help="IG 归因标量，默认 binding_topk（与训练解码一致）",
    )
    args = ap.parse_args()

    ckpt = baseline_checkpoint(args.task)

    if args.h5:
        h5_mode = "single"
        h5 = args.h5
    elif args.all_chromosomes:
        h5_mode = "all_chr"
        h5 = ""
    else:
        h5_mode = "single"
        h5 = str(default_chr_h5(args.task, chrom=args.chr))

    out = args.out or str(ROOT / "outputs" / args.task)

    env = {
        **os.environ,
        "PYTHONPATH": "/data/CARBON_PROJECT/MARCP:" + os.environ.get("PYTHONPATH", ""),
        "PYTHONUNBUFFERED": "1",
    }
    cmd = [
        sys.executable,
        "-u",
        str(ROOT / "scripts" / "run_pipeline.py"),
        "--checkpoint",
        str(ckpt),
        "--task",
        args.task,
        "--h5_mode",
        h5_mode,
        "--out",
        out,
        "--batch_size",
        str(args.batch_size),
        "--ig_steps",
        str(args.ig_steps),
        "--select_top",
        str(args.select_top),
        "--rank_by",
        args.rank_by,
    ]
    if h5_mode == "single":
        cmd.extend(["--h5", h5])
    if args.max_samples is not None:
        cmd.extend(["--max_samples", str(args.max_samples)])
    if args.skip_conv:
        cmd.append("--skip_conv")
    if args.skip_meme:
        cmd.append("--skip_meme")
    if args.skip_ism:
        cmd.append("--skip_ism")
    if args.skip_tfmodisco:
        cmd.append("--skip_tfmodisco")
    if args.legacy_kmeans:
        cmd.append("--legacy_kmeans")
    if args.forward_chunk_rows is not None:
        cmd.extend(["--forward_chunk_rows", str(args.forward_chunk_rows)])
    if args.no_chunked_rank:
        cmd.append("--no_chunked_rank")
    if args.max_seqlets is not None:
        cmd.extend(["--max_seqlets", str(args.max_seqlets)])
    cmd.extend(["--attr_target", args.attr_target])

    print("checkpoint:", ckpt, flush=True)
    print("h5_mode:", h5_mode, "h5:", h5 or "(all_chr)", flush=True)
    print("out:", out, flush=True)
    print("cmd:", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=str(ROOT), env=env)
    raise SystemExit(r.returncode)


if __name__ == "__main__":
    main()
