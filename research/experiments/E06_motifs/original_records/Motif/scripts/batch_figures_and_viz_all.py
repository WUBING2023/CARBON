#!/usr/bin/env python3
"""
对 outputs/<TASK>/ 批量运行 plot_task_summary.py + build_task_visualization.py（与参考任务 figures/、visualization/ 对齐）。

与 batch_threefold_validation_all 共用阶段 D 单实例锁（locks/stage_d.lock）。
DISABLE_STAGE_LOCK=1 可跳过（不推荐）。

  cd 可解释性分析/Motif
  python3 scripts/batch_figures_and_viz_all.py --jobs 8 --summary
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


def _discover(outputs_root: Path) -> list[Path]:
    out: list[Path] = []
    for child in sorted(outputs_root.iterdir()):
        if not child.is_dir() or child.name in SKIP_DIR_NAMES or child.name.startswith("."):
            continue
        if not (child / "run_config.json").is_file():
            continue
        out.append(child)
    return out


def _run_one(task_dir: str) -> dict:
    td = Path(task_dir)
    name = td.name
    env = {
        **os.environ,
        "PYTHONPATH": f"{ROOT}:{MARCP}:{os.environ.get('PYTHONPATH', '')}",
        "PYTHONUNBUFFERED": "1",
    }
    plot_cmd = [sys.executable, "-u", str(ROOT / "scripts" / "plot_task_summary.py"), "--task-dir", str(td)]
    viz_cmd = [sys.executable, "-u", str(ROOT / "scripts" / "build_task_visualization.py"), "--task-dir", str(td)]
    try:
        p1 = subprocess.run(plot_cmd, cwd=str(ROOT), env=env, capture_output=True, text=True, timeout=600)
        if p1.returncode != 0:
            return {
                "task": name,
                "exit": str(p1.returncode),
                "step": "plot",
                "tail": ((p1.stderr or "") + (p1.stdout or ""))[-1500:],
            }
        p2 = subprocess.run(viz_cmd, cwd=str(ROOT), env=env, capture_output=True, text=True, timeout=600)
        return {
            "task": name,
            "exit": str(p2.returncode),
            "step": "viz" if p2.returncode != 0 else "ok",
            "tail": ((p2.stderr or "") + (p2.stdout or ""))[-1500:] if p2.returncode != 0 else "",
        }
    except subprocess.TimeoutExpired:
        return {"task": name, "exit": "timeout", "step": "?", "tail": ""}
    except Exception as e:
        return {"task": name, "exit": "exc", "step": "?", "tail": repr(e)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--summary", action="store_true")
    args = ap.parse_args()

    out_root = args.outputs_root.resolve()
    queued = _discover(out_root)
    print(f"{_ts()} [figviz_batch] queue={len(queued)} jobs={args.jobs}", flush=True)

    if args.summary:
        for p in queued[:30]:
            print(f"  {p.name}")
        if len(queued) > 30:
            print(f"  ... +{len(queued) - 30}")
        return

    log_dir = ROOT / "outputs" / "batch_motif_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    joblog = log_dir / f"figviz_batch_joblog_{stamp}.csv"

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
                    row = {"task": name, "exit": "exc", "step": "?", "tail": repr(e)}
                rows.append(
                    {"task": row["task"], "exit": row["exit"], "step": row.get("step", "")}
                )
                tag = "OK" if row["exit"] == "0" else f"FAIL({row['exit']} {row.get('step','')})"
                print(f"{_ts()} [figviz_batch] {tag} {row['task']}", flush=True)
                if row["exit"] != "0" and row.get("tail"):
                    for ln in str(row["tail"]).splitlines()[-3:]:
                        if ln.strip():
                            print(f"    {ln[:400]}", flush=True)

        rows.sort(key=lambda x: x["task"])
        with open(joblog, "w", encoding="utf-8") as f:
            f.write("task,exit,step\n")
            for r in rows:
                f.write(f"{r['task']},{r['exit']},{r.get('step','')}\n")

        n_ok = sum(1 for r in rows if r["exit"] == "0")
        print(f"{_ts()} [figviz_batch] done ok={n_ok}/{len(rows)} joblog={joblog}", flush=True)
        raise SystemExit(0 if n_ok == len(rows) else 1)
    finally:
        release_lock(lock_fd)


if __name__ == "__main__":
    main()
