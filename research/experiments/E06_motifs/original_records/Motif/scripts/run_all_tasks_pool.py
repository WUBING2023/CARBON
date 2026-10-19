#!/usr/bin/env python3
"""
不依赖 GNU parallel 的多 GPU 批量 Motif：按 GPU 启动多个工作线程，
从同一队列取任务，子进程环境变量 CUDA_VISIBLE_DEVICES 仅暴露该卡。

默认每 GPU 1 个线程；可用 --slots 指定每卡并行进程数（多进程共享同一块 GPU，显存与算力会争用）。

阶段 A 输出策略（与 PIPELINE / RUN_PIPELINE 一致）：
  - 运行中写入 outputs/workspace/<TASK>/；成功且存在 tfmodisco/modisco_results.h5 时，
    整目录移至 outputs/<TASK>/（若目标已存在则先删除）。
  - 启动某一任务前若 outputs/workspace/<TASK>/ 已存在（上次中断残留），先删除再跑。
  - 失败或非零退出时保留 workspace 便于排查。

默认任务筛选与 run_all_tasks.py --resume 一致（跳过已有 outputs/<TASK>/tfmodisco/modisco_results.h5）；
加 --no-resume 则与不带 --resume 的 dry-run 一致（凡有权重的任务都入队，用于强制全量重跑）。

每槽位 CPU 线程数由 --threads-per-slot（或环境变量 POOL_THREADS_PER_SLOT，默认 3）写入
OMP_NUM_THREADS / MKL_NUM_THREADS / OPENBLAS_NUM_THREADS / NUMEXPR_NUM_THREADS。

同阶段仅允许一个池实例：默认占用 outputs/batch_motif_logs/locks/stage_a_pipeline_pool.lock
（DISABLE_STAGE_LOCK=1 可跳过，不推荐）。

示例:
  python3 scripts/run_all_tasks_pool.py --dry-run
  python3 scripts/run_all_tasks_pool.py --gpus 0,1,2,3 --slots 20,9,9,9 --threads-per-slot 3
  nohup python3 scripts/run_all_tasks_pool.py --gpus 0,1,2,3 --slots 20,9,9,9 > outputs/batch_motif_logs/pool_nohup.log 2>&1 &
"""
from __future__ import annotations

import argparse
import csv
import os
import queue
import shutil
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from motif_analysis.stage_lock import acquire_lock, lock_disabled_from_env, release_lock

MARCP_DEFAULT = Path("/data/CARBON_PROJECT/MARCP")

WORKSPACE_ROOT = ROOT / "outputs" / "workspace"
MODISCO_H5_REL = Path("tfmodisco") / "modisco_results.h5"


def _ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _todo_tasks(*, resume: bool) -> list[str]:
    cmd = [sys.executable, str(ROOT / "scripts" / "run_all_tasks.py"), "--dry-run"]
    if resume:
        cmd.append("--resume")
    r = subprocess.run(
        cmd,
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        print(r.stderr or r.stdout, file=sys.stderr)
        raise SystemExit(r.returncode)
    out: list[str] = []
    for ln in r.stdout.splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        if "Planned tasks:" in s:
            continue
        if s.startswith("[") and len(s) > 10 and s[1:5].isdigit() and "T" in s[:20]:
            continue
        out.append(s)
    return out


def _default_gpus() -> list[str]:
    try:
        r = subprocess.run(
            ["nvidia-smi", "-L"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        if r.returncode != 0:
            return ["0"]
        n = sum(1 for ln in r.stdout.splitlines() if ln.strip().startswith("GPU "))
        if n < 1:
            n = len([ln for ln in r.stdout.splitlines() if "UUID" in ln or "GPU" in ln[:4]])
        if n < 1:
            return ["0"]
        return [str(i) for i in range(n)]
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return ["0"]


def _promote_workspace_to_outputs(task: str) -> None:
    """Move outputs/workspace/<task> -> outputs/<task> if modisco h5 exists."""
    ws = WORKSPACE_ROOT / task
    modisco = ws / MODISCO_H5_REL
    if not modisco.is_file():
        raise FileNotFoundError(str(modisco))
    final = ROOT / "outputs" / task
    if final.exists():
        shutil.rmtree(final)
    shutil.move(str(ws), str(final))


def _worker(
    gpu: str,
    slot: int,
    q: queue.Queue,
    marcp: Path,
    job_rows: list,
    lock: threading.Lock,
    threads_per_slot: int,
) -> None:
    tp = str(max(1, int(threads_per_slot)))
    env_base = {
        **os.environ,
        "PYTHONPATH": f"{marcp}:{os.environ.get('PYTHONPATH', '')}",
        "CUDA_VISIBLE_DEVICES": gpu,
        "PYTHONUNBUFFERED": "1",
        "OMP_NUM_THREADS": tp,
        "MKL_NUM_THREADS": tp,
        "OPENBLAS_NUM_THREADS": tp,
        "NUMEXPR_NUM_THREADS": tp,
    }
    while True:
        item = q.get()
        try:
            if item is None:
                return
            task = item
            t0 = time.perf_counter()
            ws = WORKSPACE_ROOT / task
            if ws.exists():
                shutil.rmtree(ws)
            WORKSPACE_ROOT.mkdir(parents=True, exist_ok=True)

            cmd = [
                sys.executable,
                "-u",
                str(ROOT / "scripts" / "run_from_ablation_task.py"),
                "--task",
                task,
                "--all-chromosomes",
                "--select_top",
                "10000",
                "--rank_by",
                "pred_mean",
                "--forward-chunk-rows",
                "512",
                "--batch_size",
                "16",
                "--ig_steps",
                "24",
                "--out",
                str(ws),
            ]
            print(f"{_ts()} [gpu {gpu} slot {slot}] START {task} -> workspace", flush=True)
            r = subprocess.run(cmd, cwd=str(ROOT), env=env_base)
            dt = time.perf_counter() - t0
            row = {
                "gpu": gpu,
                "slot": str(slot),
                "task": task,
                "exit": str(r.returncode),
                "seconds": f"{dt:.1f}",
                "time": datetime.now().isoformat(timespec="seconds"),
            }
            promote_note = ""
            if r.returncode == 0:
                modisco = ws / MODISCO_H5_REL
                if modisco.is_file():
                    try:
                        _promote_workspace_to_outputs(task)
                        promote_note = " promoted->outputs"
                    except OSError as e:
                        promote_note = f" PROMOTE_ERR {e}"
                        row["exit"] = "2"
                else:
                    promote_note = " WARN:no_modisco_h5_in_workspace"
                    print(
                        f"{_ts()} [gpu {gpu} slot {slot}] WARN exit=0 but missing {modisco}",
                        flush=True,
                    )
            else:
                print(
                    f"{_ts()} [gpu {gpu} slot {slot}] FAIL {task} workspace kept at {ws}",
                    flush=True,
                )
            with lock:
                job_rows.append(row)
            tag = "OK" if row["exit"] == "0" else "FAIL"
            print(
                f"{_ts()} [gpu {gpu} slot {slot}] {tag} {task} ({dt:.1f}s) exit={row['exit']}{promote_note}",
                flush=True,
            )
        finally:
            q.task_done()


def main() -> None:
    ap = argparse.ArgumentParser(description="Multi-GPU Motif batch without GNU parallel.")
    ap.add_argument("--dry-run", action="store_true", help="仅列出待跑任务并退出")
    ap.add_argument(
        "--gpus",
        type=str,
        default="",
        help="逗号分隔物理 GPU 编号，如 0,1,2,3；默认按 nvidia-smi -L 数量使用 0..N-1",
    )
    ap.add_argument(
        "--slots",
        type=str,
        default="",
        help="与 --gpus 同序的每卡工作线程数，逗号分隔，如 20,9,9,9；空则每 GPU 1 线程",
    )
    ap.add_argument(
        "--threads-per-slot",
        type=int,
        default=int(os.environ.get("POOL_THREADS_PER_SLOT", "3")),
        help="每槽位 OpenMP/BLAS 线程上限（默认 3，或环境变量 POOL_THREADS_PER_SLOT）",
    )
    ap.add_argument(
        "--marcp-root",
        type=Path,
        default=Path(os.environ.get("MARCP_ROOT", str(MARCP_DEFAULT))),
    )
    ap.add_argument(
        "--no-resume",
        action="store_true",
        help="不按 modisco_results.h5 跳过；凡 run_all_tasks 会跑的任务都入队（强制重跑）",
    )
    args = ap.parse_args()

    tasks = _todo_tasks(resume=not args.no_resume)
    if args.dry_run:
        for t in tasks:
            print(t)
        print(f"# count={len(tasks)}", file=sys.stderr)
        return

    if not tasks:
        print(f"{_ts()} No tasks to run.", flush=True)
        return

    lock_fd = None
    exit_code = 0
    try:
        lock_fd = acquire_lock(ROOT, "stage_a_pipeline_pool", disable=lock_disabled_from_env())

        if args.gpus.strip():
            gpu_ids = [x.strip() for x in args.gpus.split(",") if x.strip()]
        else:
            gpu_ids = _default_gpus()

        if not gpu_ids:
            gpu_ids = ["0"]

        if args.slots.strip():
            slot_parts = [x.strip() for x in args.slots.split(",") if x.strip()]
            try:
                slot_counts = [max(1, int(x)) for x in slot_parts]
            except ValueError:
                print("[error] --slots must be comma-separated positive integers", file=sys.stderr)
                raise SystemExit(2)
            if len(slot_counts) != len(gpu_ids):
                print(
                    f"[error] --slots count ({len(slot_counts)}) must match --gpus count ({len(gpu_ids)})",
                    file=sys.stderr,
                )
                raise SystemExit(2)
        else:
            slot_counts = [1] * len(gpu_ids)

        total_workers = sum(slot_counts)
        worker_specs: list[tuple[str, int]] = []
        for gpu, n in zip(gpu_ids, slot_counts):
            for s in range(n):
                worker_specs.append((gpu, s))

        log_dir = ROOT / "outputs" / "batch_motif_logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        joblog = log_dir / f"pool_joblog_{stamp}.csv"

        print(
            f"{_ts()} Tasks={len(tasks)}, gpus={gpu_ids}, slots={slot_counts}, total_workers={total_workers}, "
            f"threads_per_slot={args.threads_per_slot}, workspace={WORKSPACE_ROOT}",
            flush=True,
        )
        print(f"{_ts()} Joblog: {joblog}", flush=True)

        q: queue.Queue[str | None] = queue.Queue()
        for t in tasks:
            q.put(t)
        for _ in range(total_workers):
            q.put(None)

        job_rows: list[dict[str, str]] = []
        lock = threading.Lock()
        threads = [
            threading.Thread(
                target=_worker,
                args=(gpu, slot, q, args.marcp_root.resolve(), job_rows, lock, args.threads_per_slot),
                daemon=False,
                name=f"worker-gpu{gpu}-{slot}",
            )
            for gpu, slot in worker_specs
        ]
        t_wall0 = time.perf_counter()
        for th in threads:
            th.start()
        for th in threads:
            th.join()
        wall = time.perf_counter() - t_wall0

        with open(joblog, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["gpu", "slot", "task", "exit", "seconds", "time"])
            w.writeheader()
            w.writerows(job_rows)

        n_ok = sum(1 for r in job_rows if r["exit"] == "0")
        n_fail = len(job_rows) - n_ok
        print(f"{_ts()} Done wall={wall:.1f}s ok={n_ok} fail={n_fail} -> {joblog}", flush=True)
        exit_code = 0 if n_fail == 0 else 1
    finally:
        release_lock(lock_fd)
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
