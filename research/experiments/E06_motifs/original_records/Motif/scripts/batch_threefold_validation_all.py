#!/usr/bin/env python3
"""
对 outputs/<TASK>/（已有 tfmodisco_export/queries.meme）批量运行 run_threefold_validation.py。
可加 --fill-missing-ism：仅补全尚未生成 threefold_validation/ism_profile_compare.png 的任务。

与 batch_figures_and_viz_all 共用阶段 D 单实例锁（outputs/batch_motif_logs/locks/stage_d.lock）；
同机同时只能跑其中一个后台批任务。DISABLE_STAGE_LOCK=1 可跳过（不推荐）。

  cd 可解释性分析/Motif
  mkdir -p outputs/batch_motif_logs   # nohup 重定向前目录须存在，否则 shell 直接报错 Exit 1
  python3 scripts/batch_threefold_validation_all.py --jobs 2 --summary
  nohup python3 scripts/batch_threefold_validation_all.py --jobs 2 >> outputs/batch_motif_logs/threefold_batch.log 2>&1 &

  # 已有三项验证报告、但缺 ISM 曲线 PNG 时，只补跑缺失曲线（推荐）：
  python3 scripts/batch_threefold_validation_all.py --jobs 2 --fill-missing-ism --summary
  nohup python3 scripts/batch_threefold_validation_all.py --jobs 2 --fill-missing-ism >> outputs/batch_motif_logs/threefold_ism_fill.log 2>&1 &
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from motif_analysis.stage_lock import acquire_lock, lock_disabled_from_env, release_lock

MARCP = Path(os.environ.get("MARCP_ROOT", "/data/CARBON_PROJECT/MARCP"))

SKIP_DIR_NAMES = frozenset({"batch_motif_logs", "fixtures", "demo", ".git"})


def _ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _ism_png_ok(task_dir: Path) -> bool:
    p = task_dir / "threefold_validation" / "ism_profile_compare.png"
    return p.is_file() and p.stat().st_size > 0


def _discover(
    outputs_root: Path,
    *,
    skip_existing: bool,
    fill_missing_ism: bool,
) -> list[Path]:
    out: list[Path] = []
    for child in sorted(outputs_root.iterdir()):
        if not child.is_dir() or child.name in SKIP_DIR_NAMES or child.name.startswith("."):
            continue
        qm = child / "tfmodisco_export" / "queries.meme"
        if not qm.is_file() or qm.stat().st_size <= 0:
            continue
        rep = child / "threefold_validation" / "validation_report.json"

        if fill_missing_ism:
            # 仅补全「尚无 ISM 曲线图」的任务（可有 validation_report；以 PNG 为准）
            if _ism_png_ok(child):
                continue
            out.append(child)
            continue

        if skip_existing and rep.is_file():
            continue
        out.append(child)
    return out


def _run_one(task_dir: str) -> dict:
    td = Path(task_dir)
    t0 = datetime.now().isoformat(timespec="seconds")
    env = {
        **os.environ,
        "PYTHONPATH": f"{MARCP}:{os.environ.get('PYTHONPATH', '')}",
        "PYTHONUNBUFFERED": "1",
    }
    cmd = [
        sys.executable,
        "-u",
        str(ROOT / "scripts" / "run_threefold_validation.py"),
        "--task-dir",
        str(td),
    ]
    try:
        r = subprocess.run(cmd, cwd=str(ROOT), env=env, capture_output=True, text=True, timeout=7200)
    except subprocess.TimeoutExpired:
        return {"task": td.name, "exit": "timeout", "time": t0}
    return {
        "task": td.name,
        "exit": str(r.returncode),
        "time": t0,
        "tail": ((r.stderr or "") + (r.stdout or ""))[-2000:],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument("--jobs", type=int, default=2, help="并行数（三项验证用 GPU，宜小，默认 2）")
    ap.add_argument(
        "--no-skip-existing",
        action="store_true",
        help="即使已有 threefold_validation/validation_report.json 也重跑",
    )
    ap.add_argument("--summary", action="store_true", help="只统计将运行的任务数")
    ap.add_argument(
        "--fill-missing-ism",
        action="store_true",
        help="只排队「已有 queries.meme 但缺少 threefold_validation/ism_profile_compare.png」的任务，用于补全 ISM 曲线图；与 --no-skip-existing 互斥时以本选项为准（忽略 validation_report 是否已存在）",
    )
    args = ap.parse_args()

    out_root = args.outputs_root.resolve()
    skip_existing = not args.no_skip_existing
    if args.fill_missing_ism and args.no_skip_existing:
        print(
            f"{_ts()} [threefold_batch] note: --fill-missing-ism ignores --no-skip-existing (only tasks missing ISM PNG are queued)",
            flush=True,
        )
    queued = _discover(
        out_root,
        skip_existing=skip_existing,
        fill_missing_ism=args.fill_missing_ism,
    )
    mode = "fill_missing_ism" if args.fill_missing_ism else "full_threefold"
    print(
        f"{_ts()} [threefold_batch] mode={mode} queue={len(queued)} jobs={args.jobs} skip_existing={skip_existing}",
        flush=True,
    )

    if args.summary:
        for p in queued[:40]:
            print(f"  {p.name}")
        if len(queued) > 40:
            print(f"  ... +{len(queued) - 40}")
        return

    if not queued:
        print(f"{_ts()} [threefold_batch] nothing to do", flush=True)
        return

    log_dir = ROOT / "outputs" / "batch_motif_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    joblog = log_dir / f"threefold_batch_joblog_{stamp}.csv"

    rows: list[dict[str, str]] = []
    j = max(1, min(args.jobs, len(queued)))

    lock_fd = None
    try:
        lock_fd = acquire_lock(ROOT, "stage_d", disable=lock_disabled_from_env())
        with ProcessPoolExecutor(max_workers=j) as ex:
            futs = {ex.submit(_run_one, str(p)): p.name for p in queued}
            for fut in as_completed(futs):
                name = futs[fut]
                try:
                    row = fut.result()
                except Exception as e:
                    row = {"task": name, "exit": "exc", "time": "", "tail": repr(e)}
                rows.append(
                    {
                        "task": row["task"],
                        "exit": row["exit"],
                        "time": row.get("time", ""),
                    }
                )
                tag = "OK" if row["exit"] == "0" else f"FAIL({row['exit']})"
                print(f"{_ts()} [threefold_batch] {tag} {row['task']}", flush=True)
                if row["exit"] not in ("0",) and row.get("tail"):
                    for ln in str(row["tail"]).splitlines()[-4:]:
                        if ln.strip():
                            print(f"    {ln[:400]}", flush=True)

        rows.sort(key=lambda x: x["task"])
        with open(joblog, "w", encoding="utf-8") as f:
            f.write("task,exit,time\n")
            for r in rows:
                f.write(f"{r['task']},{r['exit']},{r['time']}\n")

        n_ok = sum(1 for r in rows if r["exit"] == "0")
        print(f"{_ts()} [threefold_batch] done ok={n_ok}/{len(rows)} joblog={joblog}", flush=True)
        raise SystemExit(0 if n_ok == len(rows) else 1)
    finally:
        release_lock(lock_fd)


if __name__ == "__main__":
    main()
