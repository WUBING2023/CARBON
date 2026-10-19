#!/usr/bin/env python3
"""
对多个 outputs/<TASK>/ 并行运行 modisco_pattern_rank.py，生成 tfmodisco_export/queries.meme
（需已有 tfmodisco/modisco_results.h5）。

  cd 可解释性分析/Motif
  python3 scripts/batch_modisco_queries_export.py --jobs 40 --summary
  nohup python3 scripts/batch_modisco_queries_export.py --jobs 40 \\
    >> outputs/batch_motif_logs/queries_export.log 2>&1 &
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from motif_analysis.stage_lock import acquire_lock, lock_disabled_from_env, release_lock

SKIP_DIR_NAMES = frozenset(
    {
        "batch_motif_logs",
        "fixtures",
        "demo",
        ".git",
    }
)


def _discover_task_dirs(outputs_root: Path) -> list[Path]:
    out: list[Path] = []
    for child in sorted(outputs_root.iterdir()):
        if not child.is_dir() or child.name in SKIP_DIR_NAMES or child.name.startswith("."):
            continue
        out.append(child)
    return out


def _needs_run(task_dir: Path, *, skip_existing: bool) -> tuple[bool, str]:
    h5 = task_dir / "tfmodisco" / "modisco_results.h5"
    qm = task_dir / "tfmodisco_export" / "queries.meme"
    if not h5.is_file():
        return False, "no_h5"
    if skip_existing and qm.is_file() and qm.stat().st_size > 0:
        return False, "has_queries"
    return True, "queue"


def _run_one(payload: dict) -> dict:
    task_dir = Path(payload["task_dir"])
    name = task_dir.name
    py = payload["python"]
    rank_script = payload["rank_script"]
    extra = payload.get("extra_args", [])
    t0 = datetime.now().isoformat(timespec="seconds")
    cmd = [py, str(rank_script), "--task-dir", str(task_dir), *extra]
    try:
        r = subprocess.run(
            cmd,
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            timeout=payload.get("timeout_sec", 7200),
        )
    except subprocess.TimeoutExpired:
        return {
            "task": name,
            "exit": "timeout",
            "stdout": "",
            "stderr": "timeout",
            "time": t0,
        }
    err = (r.stderr or "") + (r.stdout or "")
    return {
        "task": name,
        "exit": str(r.returncode),
        "stdout": (r.stdout or "")[-4000:],
        "stderr": (r.stderr or "")[-4000:],
        "time": t0,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument("--jobs", type=int, default=40, help="并行子进程数（默认 40）")
    ap.add_argument(
        "--no-skip-existing",
        action="store_true",
        help="即使已有 queries.meme 也重跑（覆盖 tfmodisco_export 内相关文件）",
    )
    ap.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="传给 modisco_pattern_rank.py --top-n（默认 10，与单任务脚本一致）",
    )
    ap.add_argument(
        "--timeout-sec",
        type=int,
        default=7200,
        help="单任务超时秒数（默认 7200）",
    )
    ap.add_argument(
        "--meme-pseudocount",
        type=float,
        default=0.001,
        help="传给 modisco_pattern_rank.py：sequence PPM 伪计数（默认 0.001）",
    )
    ap.add_argument("--summary", action="store_true", help="仅打印将处理/跳过的统计后退出")
    ap.add_argument("--task", action="append", dest="tasks", help="仅处理这些任务目录名（可重复传入）")
    args = ap.parse_args()

    out_root = args.outputs_root.resolve()
    if not out_root.is_dir():
        print(f"[error] not a directory: {out_root}", file=sys.stderr)
        raise SystemExit(2)

    rank_script = ROOT / "scripts" / "modisco_pattern_rank.py"
    if not rank_script.is_file():
        print(f"[error] missing {rank_script}", file=sys.stderr)
        raise SystemExit(2)

    skip_existing = not args.no_skip_existing
    all_dirs = _discover_task_dirs(out_root)
    if args.tasks:
        want = set(args.tasks)
        all_dirs = [d for d in all_dirs if d.name in want]
        missing = want - {d.name for d in all_dirs}
        if missing:
            print(f"[warn] 未在 outputs 下找到任务目录: {sorted(missing)}", flush=True)
    queued: list[Path] = []
    n_no_h5 = 0
    n_has_queries = 0
    for td in all_dirs:
        need, reason = _needs_run(td, skip_existing=skip_existing)
        if reason == "no_h5":
            n_no_h5 += 1
        elif reason == "has_queries":
            n_has_queries += 1
        if need:
            queued.append(td)

    n_tasks = len(all_dirs)
    print(
        f"[batch_queries] outputs={out_root} task_dirs={n_tasks} "
        f"queue={len(queued)} no_h5={n_no_h5} skip_has_queries={n_has_queries} jobs={args.jobs}",
        flush=True,
    )

    if args.summary:
        for td in queued[:50]:
            print(f"  queue: {td.name}")
        if len(queued) > 50:
            print(f"  ... +{len(queued) - 50} more")
        return

    if not queued:
        print("[batch_queries] nothing to do.", flush=True)
        return

    if args.top_n <= 0:
        print("[error] --top-n must be > 0 to write queries.meme", file=sys.stderr)
        raise SystemExit(2)
    extra_args = ["--top-n", str(args.top_n), "--meme-pseudocount", str(args.meme_pseudocount)]

    payload_base = {
        "python": sys.executable,
        "rank_script": str(rank_script),
        "extra_args": extra_args,
        "timeout_sec": args.timeout_sec,
    }

    log_dir = ROOT / "outputs" / "batch_motif_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    joblog = log_dir / f"queries_export_joblog_{stamp}.csv"

    rows: list[dict[str, str]] = []
    jobs = min(max(1, args.jobs), len(queued))

    lock_fd = None
    try:
        lock_fd = acquire_lock(ROOT, "stage_b_queries_export", disable=lock_disabled_from_env())
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            futs = {
                ex.submit(_run_one, {**payload_base, "task_dir": str(td)}): td.name
                for td in queued
            }
            for fut in as_completed(futs):
                name = futs[fut]
                try:
                    row = fut.result()
                except Exception as e:
                    row = {"task": name, "exit": "exc", "stdout": "", "stderr": repr(e), "time": ""}
                rows.append(
                    {
                        "task": row["task"],
                        "exit": row["exit"],
                        "time": row.get("time", ""),
                    }
                )
                tag = "OK" if row["exit"] == "0" else f"FAIL({row['exit']})"
                print(f"[batch_queries] {tag} {row['task']}", flush=True)
                if row["exit"] != "0" and (row.get("stderr") or row.get("stdout")):
                    tail = (row.get("stderr") or "") + (row.get("stdout") or "")
                    for ln in tail.splitlines()[-5:]:
                        if ln.strip():
                            print(f"    {ln[:500]}", flush=True)

        rows.sort(key=lambda x: x["task"])
        with open(joblog, "w", encoding="utf-8") as f:
            f.write("task,exit,time\n")
            for r in rows:
                f.write(f"{r['task']},{r['exit']},{r['time']}\n")

        n_ok = sum(1 for r in rows if r["exit"] == "0")
        print(f"[batch_queries] done ok={n_ok}/{len(rows)} joblog={joblog}", flush=True)
        raise SystemExit(0 if n_ok == len(rows) else 1)
    finally:
        release_lock(lock_fd)


if __name__ == "__main__":
    main()
