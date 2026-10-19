#!/usr/bin/env python3
"""
对所有 outputs/<TASK>/（已有 tfmodisco/modisco_results.h5）批量运行 plot_modisco_ic_logos.py，
生成与 RBFOX2 同结构的 tfmodisco/ic_scaled_logos_mode_<MODE>/（默认 mode=A）。

  cd 可解释性分析/Motif
  python3 scripts/batch_plot_ic_logos.py --mode A --jobs 8
  python3 scripts/batch_plot_ic_logos.py --mode A --jobs 8 --summary
  nohup python3 scripts/batch_plot_ic_logos.py --mode A --jobs 8 >> outputs/batch_motif_logs/ic_logos_batch.log 2>&1 &
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SKIP_DIR_NAMES = frozenset({"batch_motif_logs", "fixtures", "demo", ".git"})


def _discover(outputs_root: Path) -> list[Path]:
    out: list[Path] = []
    for child in sorted(outputs_root.iterdir()):
        if not child.is_dir() or child.name in SKIP_DIR_NAMES or child.name.startswith("."):
            continue
        out.append(child)
    return out


def _needs_run(task_dir: Path, mode: str, skip_existing: bool) -> tuple[bool, str]:
    h5 = task_dir / "tfmodisco" / "modisco_results.h5"
    if not h5.is_file():
        return False, "no_h5"
    out_dir = task_dir / "tfmodisco" / f"ic_scaled_logos_mode_{mode}"
    if skip_existing and out_dir.is_dir() and any(out_dir.glob("*.ic_fwd.png")):
        return False, "has_output"
    return True, "queue"


def _run_one(payload: dict) -> dict:
    task_dir = Path(payload["task_dir"])
    mode = payload["mode"]
    py = payload["python"]
    script = payload["script"]
    cmd = [py, str(script), "--task-dir", str(task_dir), "--mode", mode]
    t0 = datetime.now().isoformat(timespec="seconds")
    try:
        r = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return {"task": task_dir.name, "exit": "timeout", "stderr": "timeout", "time": t0}
    err = (r.stderr or "") + (r.stdout or "")
    return {
        "task": task_dir.name,
        "exit": str(r.returncode),
        "stderr": err[-6000:],
        "time": t0,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument("--mode", type=str.upper, default="A", choices=("A", "B", "C"))
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument(
        "--skip-existing",
        action="store_true",
        default=True,
        help="若已有 *.ic_fwd.png 则跳过（默认开）",
    )
    ap.add_argument(
        "--no-skip-existing",
        action="store_true",
        help="覆盖重跑",
    )
    ap.add_argument("--summary", action="store_true")
    args = ap.parse_args()
    out_root = args.outputs_root.resolve()
    skip_existing = args.skip_existing and not args.no_skip_existing
    mode = args.mode

    script = ROOT / "scripts" / "plot_modisco_ic_logos.py"
    if not script.is_file():
        print(f"[error] missing {script}", file=sys.stderr)
        raise SystemExit(2)

    queued: list[Path] = []
    n_no_h5 = 0
    n_skip = 0
    for td in _discover(out_root):
        need, reason = _needs_run(td, mode, skip_existing)
        if reason == "no_h5":
            n_no_h5 += 1
        elif reason == "has_output":
            n_skip += 1
        if need:
            queued.append(td)

    n_all = len(_discover(out_root))
    print(
        f"[batch_ic_logos] outputs={out_root} mode={mode} task_dirs={n_all} "
        f"queue={len(queued)} no_h5={n_no_h5} skip_has_output={n_skip} jobs={args.jobs}",
        flush=True,
    )

    if args.summary:
        for td in queued[:80]:
            print(f"  queue: {td.name}")
        if len(queued) > 80:
            print(f"  ... +{len(queued) - 80} more")
        return

    if not queued:
        print("[batch_ic_logos] nothing to do.", flush=True)
        return

    payload_base = {"python": sys.executable, "script": str(script), "mode": mode}
    jobs = min(max(1, args.jobs), len(queued))
    log_dir = ROOT / "outputs" / "batch_motif_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    joblog = log_dir / f"ic_logos_mode{mode}_joblog_{stamp}.csv"

    rows: list[dict[str, str]] = []
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(_run_one, {**payload_base, "task_dir": str(td)}): td.name for td in queued}
        for fut in as_completed(futs):
            name = futs[fut]
            try:
                row = fut.result()
            except Exception as e:
                row = {"task": name, "exit": "exc", "stderr": repr(e), "time": ""}
            rows.append({"task": row["task"], "exit": row["exit"], "time": row.get("time", "")})
            tag = "OK" if row["exit"] == "0" else f"FAIL({row['exit']})"
            print(f"[batch_ic_logos] {tag} {row['task']}", flush=True)
            if row["exit"] != "0" and row.get("stderr"):
                for ln in str(row["stderr"]).splitlines()[-4:]:
                    if ln.strip():
                        print(f"    {ln[:400]}", flush=True)

    rows.sort(key=lambda x: x["task"])
    with open(joblog, "w", encoding="utf-8") as f:
        f.write("task,exit,time\n")
        for r in rows:
            f.write(f"{r['task']},{r['exit']},{r['time']}\n")

    n_ok = sum(1 for r in rows if r["exit"] == "0")
    print(f"[batch_ic_logos] done ok={n_ok}/{len(rows)} joblog={joblog}", flush=True)
    raise SystemExit(0 if n_ok == len(rows) else 1)


if __name__ == "__main__":
    main()
