#!/usr/bin/env python3
"""
对「模型权重/MARCP/*/baseline/model_params.pth」下所有存在 baseline 的任务，
依次运行 Motif 流水线。结果写在 **outputs/<TASK_NAME>/**（按实验任务，不按染色体命名）。

**默认（论文级全量）**：合并每任务全部 `chr*_5ch.h5`，`--select-top 5000`，不截断 `max_samples`，
且**不**跳过 TF-MoDISco（除非显式 `--skip-tfmodisco`）。单染色体快速试跑请加 `--single-chromosome`。

  python scripts/run_all_tasks.py
  python scripts/run_all_tasks.py --resume   # 跳过已有 tfmodisco/modisco_results.h5
  python scripts/run_all_tasks.py --single-chromosome --chr 22 --max-samples 128
  python scripts/run_all_tasks.py --dry-run
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from motif_analysis.user_paths import H5DATA_ROOT, MARCP_WEIGHTS_ROOT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chr", type=str, default="22")
    ap.add_argument(
        "--single-chromosome",
        action="store_true",
        help="仅使用各任务 chr{--chr} 单文件（快速试跑）；默认合并全部染色体",
    )
    ap.add_argument("--max-samples", type=int, default=None, help="合并后最多样本数；默认不截断")
    ap.add_argument(
        "--select-top",
        type=int,
        default=10000,
        help="归因前 Top-K（与 run_from_ablation_task 一致）",
    )
    ap.add_argument(
        "--rank-by",
        type=str,
        default="pred_mean",
        choices=["pred_mean", "softmax_mean", "spearman_y"],
        help="Top-K 筛选指标",
    )
    ap.add_argument(
        "--forward-chunk-rows",
        type=int,
        default=512,
        help="all_chr 分块前向时每块行数（run_pipeline --forward_chunk_rows）",
    )
    ap.add_argument(
        "--no-chunked-rank",
        action="store_true",
        help="禁用分块 Top-K，强制整表 concat（易 OOM；仅调试用）",
    )
    ap.add_argument(
        "--max-seqlets",
        type=int,
        default=None,
        help="传给 run_pipeline --max_seqlets（可选）",
    )
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--ig-steps", type=int, default=24)
    ap.add_argument(
        "--resume",
        action="store_true",
        help="若 outputs/<TASK>/tfmodisco/modisco_results.h5 已存在则跳过（视为 pipeline 已跑过）",
    )
    ap.add_argument("--dry-run", action="store_true", help="只列出将运行的任务")
    ap.add_argument("--limit", type=int, default=0, help="仅处理前 N 个任务（调试用，0=不限制）")
    ap.add_argument("--skip-conv", action="store_true", help="跳过卷积核 Motif（透传 run_pipeline --skip_conv）")
    ap.add_argument("--skip-meme", action="store_true")
    ap.add_argument("--skip-ism", action="store_true")
    ap.add_argument(
        "--skip-tfmodisco",
        action="store_true",
        help="传给 run_from_ablation_task：不调用 modisco CLI（批量时更快）",
    )
    args = ap.parse_args()
    all_chromosomes = not args.single_chromosome

    tasks = sorted(p.name for p in MARCP_WEIGHTS_ROOT.iterdir() if p.is_dir())
    env = {**os.environ, "PYTHONPATH": "/data/CARBON_PROJECT/MARCP:" + os.environ.get("PYTHONPATH", "")}

    todo = []
    for task in tasks:
        ckpt = MARCP_WEIGHTS_ROOT / task / "baseline" / "model_params.pth"
        if not ckpt.is_file():
            continue
        if all_chromosomes:
            from motif_analysis.user_paths import list_chr5ch

            try:
                chr_paths = list_chr5ch(task)
            except FileNotFoundError:
                continue
            if not chr_paths:
                continue
        else:
            h5path = H5DATA_ROOT / task / "_5ch" / f"{task}_chr{args.chr}_5ch.h5"
            if not h5path.is_file():
                continue
        out = ROOT / "outputs" / task
        modisco_h5 = out / "tfmodisco" / "modisco_results.h5"
        if args.resume and modisco_h5.is_file():
            continue
        todo.append(task)

    if args.limit and args.limit > 0:
        todo = todo[: args.limit]

    log_dir = ROOT / "outputs" / "batch_motif_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"batch_{ts}.log"

    skipped_no_ckpt = sum(
        1
        for t in tasks
        if not (MARCP_WEIGHTS_ROOT / t / "baseline" / "model_params.pth").is_file()
    )
    skipped_no_h5 = 0
    for task in tasks:
        if not (MARCP_WEIGHTS_ROOT / task / "baseline" / "model_params.pth").is_file():
            continue
        if all_chromosomes:
            from motif_analysis.user_paths import list_chr5ch

            try:
                if not list_chr5ch(task):
                    skipped_no_h5 += 1
            except FileNotFoundError:
                skipped_no_h5 += 1
        elif not (H5DATA_ROOT / task / "_5ch" / f"{task}_chr{args.chr}_5ch.h5").is_file():
            skipped_no_h5 += 1

    summary = {
        "started_at": ts,
        "total_planned": len(todo),
        "skipped_no_baseline_ckpt": skipped_no_ckpt,
        "skipped_no_chr_h5": skipped_no_h5,
        "chr": args.chr,
        "all_chromosomes": all_chromosomes,
        "max_samples": args.max_samples,
        "select_top": args.select_top,
        "rank_by": args.rank_by,
        "forward_chunk_rows": args.forward_chunk_rows,
        "no_chunked_rank": args.no_chunked_rank,
        "ok": [],
        "fail": [],
    }

    def log(msg: str) -> None:
        line = f"[{datetime.now().isoformat(timespec='seconds')}] {msg}"
        print(line, flush=True)
        with open(log_file, "a", encoding="utf-8") as f:
            f.write(line + "\n")

    planned_msg = f"Planned tasks: {len(todo)} (log: {log_file})"
    if args.dry_run:
        # stdout 仅任务名一行一个，便于管道 / run_all_tasks_pool 解析；摘要打到 stderr
        line = f"[{datetime.now().isoformat(timespec='seconds')}] {planned_msg}"
        print(line, file=sys.stderr, flush=True)
        with open(log_file, "a", encoding="utf-8") as f:
            f.write(line + "\n")
        for t in todo:
            print(t)
        return

    log(planned_msg)

    for i, task in enumerate(todo, 1):
        log(f"[{i}/{len(todo)}] START {task}")
        t0 = time.perf_counter()
        cmd = [
            sys.executable,
            str(ROOT / "scripts" / "run_from_ablation_task.py"),
            "--task",
            task,
            "--batch_size",
            str(args.batch_size),
            "--ig_steps",
            str(args.ig_steps),
            "--select_top",
            str(args.select_top),
            "--rank_by",
            args.rank_by,
            "--forward-chunk-rows",
            str(args.forward_chunk_rows),
            "--out",
            str(ROOT / "outputs" / task),
        ]
        if args.max_samples is not None:
            cmd.extend(["--max_samples", str(args.max_samples)])
        if all_chromosomes:
            cmd.append("--all-chromosomes")
        else:
            cmd.extend(["--chr", args.chr])
        if args.no_chunked_rank:
            cmd.append("--no-chunked-rank")
        if args.max_seqlets is not None:
            cmd.extend(["--max-seqlets", str(args.max_seqlets)])
        if args.skip_conv:
            cmd.append("--skip_conv")
        if args.skip_meme:
            cmd.append("--skip_meme")
        if args.skip_ism:
            cmd.append("--skip_ism")
        if args.skip_tfmodisco:
            cmd.append("--skip_tfmodisco")
        try:
            r = subprocess.run(cmd, cwd=str(ROOT), env=env, timeout=None)
            dt = time.perf_counter() - t0
            if r.returncode == 0:
                summary["ok"].append(task)
                log(f"[{i}/{len(todo)}] OK {task} ({dt:.1f}s)")
            else:
                summary["fail"].append({"task": task, "code": r.returncode})
                log(f"[{i}/{len(todo)}] FAIL {task} exit={r.returncode} ({dt:.1f}s)")
        except subprocess.TimeoutExpired:
            summary["fail"].append({"task": task, "code": "timeout"})
            log(f"[{i}/{len(todo)}] TIMEOUT {task}")
        except Exception as e:
            summary["fail"].append({"task": task, "error": str(e)})
            log(f"[{i}/{len(todo)}] ERROR {task}: {e}")

    summary["finished_at"] = datetime.now().isoformat(timespec="seconds")
    summary["n_ok"] = len(summary["ok"])
    summary["n_fail"] = len(summary["fail"])
    sum_path = log_dir / f"batch_summary_{ts}.json"
    sum_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    log(f"Done. summary: {sum_path}")
    raise SystemExit(0 if not summary["fail"] else 1)


if __name__ == "__main__":
    main()
