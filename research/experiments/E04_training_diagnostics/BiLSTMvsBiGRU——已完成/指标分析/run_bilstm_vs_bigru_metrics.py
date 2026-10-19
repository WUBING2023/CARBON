#!/usr/bin/env python3
"""
Schedule 20 random RBP tasks × (BiGRU on GPU3 + BiLSTM on GPU2 in parallel queues).
All MARCP outputs go to metrics_root/workspace/<run_id>/ until 40/40 succeed; then
atomic rename to metrics_root/output/<run_id>/ and write summary TSVs.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import shutil
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


TRAIN_CHRS_DEFAULT = "3,4,5,6,7,10,11,12,13,14,17,18,19,20,21,22"
VAL_CHRS_DEFAULT = "2,9,16"
TEST_CHRS_DEFAULT = "1,8,15"

MARCP_TRAIN = "marcp_train.py"


def _script_dir() -> Path:
    return Path(__file__).resolve().parent


def _rbp_root() -> Path:
    return _script_dir().parent.parent


def _marcp_root() -> Path:
    return _rbp_root() / "MARCP"


@dataclass
class TeeStream:
    """Write to multiple text streams (e.g. console + log file)."""

    streams: list

    def write(self, data: str) -> None:
        for s in self.streams:
            s.write(data)
            s.flush()

    def flush(self) -> None:
        for s in self.streams:
            s.flush()


def _file_list(data_root: Path, task: str, chrs_csv: str) -> str:
    parts: list[str] = []
    for c in chrs_csv.split(","):
        c = c.strip()
        if not c:
            continue
        p = data_root / task / "_5ch" / f"{task}_chr{c}_5ch.h5"
        if p.is_file():
            parts.append(str(p))
    return ",".join(parts)


def _task_splits_ok(data_root: Path, task: str, train_c: str, val_c: str, test_c: str) -> bool:
    return all(
        [
            bool(_file_list(data_root, task, train_c)),
            bool(_file_list(data_root, task, val_c)),
            bool(_file_list(data_root, task, test_c)),
        ]
    )


def discover_tasks(data_root: Path) -> list[str]:
    if not data_root.is_dir():
        return []
    out: list[str] = []
    for d in sorted(data_root.iterdir()):
        if d.is_dir() and (d / "_5ch").is_dir():
            out.append(d.name)
    return out


def parse_exclude(s: str | None) -> set[str]:
    if not s:
        return set()
    raw = s.replace(",", " ").split()
    return {x.strip() for x in raw if x.strip()}


def pick_tasks(
    data_root: Path,
    n: int,
    seed: int,
    exclude: set[str],
    train_c: str,
    val_c: str,
    test_c: str,
) -> list[str]:
    candidates: list[str] = []
    for t in discover_tasks(data_root):
        if t in exclude:
            continue
        if not _task_splits_ok(data_root, t, train_c, val_c, test_c):
            continue
        candidates.append(t)
    if len(candidates) < n:
        print(
            f"Need {n} tasks with full chr splits, found {len(candidates)}.",
            file=sys.stderr,
        )
        sys.exit(1)
    rng = random.Random(seed)
    return rng.sample(candidates, n)


def effective_marcp_out(run_base: Path, use_lstm: bool) -> Path:
    """Directory where marcp_train writes (matches ablation subdir rules)."""
    if use_lstm:
        return run_base / "lstm_instead_of_gru"
    return run_base


def build_marcp_cmd(
    *,
    task: str,
    train_csv: str,
    val_csv: str,
    test_csv: str,
    out_dir: Path,
    cuda: str,
    use_lstm: bool,
    epochs: int,
) -> list[str]:
    cmd = [
        sys.executable,
        "-u",
        MARCP_TRAIN,
        "-train",
        train_csv,
        "-validation",
        val_csv,
        "-test",
        test_csv,
        "--output",
        str(out_dir),
        "--task",
        task,
        "--CUDA",
        cuda,
        "--lr",
        "0.000785",
        "--bsz",
        "256",
        "--epochs",
        str(epochs),
        "--input_channels",
        "5",
        "--encode",
        "OneHot",
        "--modelFramework",
        "Treat",
        "--lossM",
        "CombinedProfileLoss",
        "--alpha",
        "1.42",
        "--beta",
        "0.54",
        "--gamma",
        "1.80",
        "--weight_decay",
        "0.000114",
        "--grad_clip",
        "2.0",
        "--max_stopping_step",
        "9",
        "--seqLength",
        "101",
        "--signalMax",
        "0",
    ]
    if use_lstm:
        cmd.append("--ablation_lstm_instead_of_gru")
    return cmd


def run_marcp(
    cmd: list[str],
    *,
    cwd: Path,
    dry_run: bool,
) -> int:
    if dry_run:
        print("DRY_RUN:", " ".join(cmd))
        return 0
    env = os.environ.copy()
    proc = subprocess.run(cmd, cwd=str(cwd), env=env)
    return proc.returncode


def expected_metric_paths(ws_run: Path, tasks: list[str]) -> list[tuple[str, str, Path]]:
    """(task, mode, path_to_test_metric)"""
    out: list[tuple[str, str, Path]] = []
    for t in tasks:
        bigru = effective_marcp_out(ws_run / t / "bigru", False) / "test_metric.txt"
        bilstm = effective_marcp_out(ws_run / t / "bilstm", True) / "test_metric.txt"
        out.append((t, "bigru", bigru))
        out.append((t, "bilstm", bilstm))
    return out


def validate_metrics(paths: list[tuple[str, str, Path]]) -> list[str]:
    missing: list[str] = []
    for task, mode, p in paths:
        if not p.is_file() or p.stat().st_size == 0:
            missing.append(f"{task}/{mode}: missing or empty {p}")
    return missing


def read_test_metric(path: Path) -> dict[str, str] | None:
    if not path.is_file():
        return None
    lines = path.read_text(encoding="utf-8", errors="replace").strip().splitlines()
    if len(lines) < 2:
        return None
    keys = lines[0].strip().split("\t")
    vals = lines[1].strip().split("\t")
    if len(keys) != len(vals):
        return None
    return dict(zip(keys, vals))


def write_summaries(ws_or_out: Path, tasks: list[str]) -> None:
    rows_all: list[dict[str, Any]] = []
    by_task: dict[str, dict[str, dict[str, str]]] = {}
    for t in tasks:
        by_task[t] = {}
    paths = expected_metric_paths(ws_or_out, tasks)
    for task, mode, p in paths:
        row = read_test_metric(p)
        if row is None:
            continue
        flat = {"task": task, "mode": mode, **row}
        rows_all.append(flat)
        by_task[task][mode] = row

    # summary_all_runs.tsv
    if rows_all:
        keys = list(rows_all[0].keys())
        all_path = ws_or_out / "summary_all_runs.tsv"
        with open(all_path, "w", encoding="utf-8") as f:
            f.write("\t".join(keys) + "\n")
            for r in rows_all:
                f.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")

    # summary_by_task_paired.tsv
    pair_path = ws_or_out / "summary_by_task_paired.tsv"
    numeric_keys = [
        "Mean_PCC",
        "Spearman",
        "Signal_AUC",
        "Signal_AP",
        "Peak_AUC",
        "Peak_AP",
        "Training_time_sec",
        "Testing_time_sec",
    ]
    with open(pair_path, "w", encoding="utf-8") as f:
        hdr = [
            "task",
            "bigru_Mean_PCC",
            "bilstm_Mean_PCC",
            "diff_Mean_PCC",
            "bigru_Peak_AP",
            "bilstm_Peak_AP",
            "diff_Peak_AP",
        ]
        f.write("\t".join(hdr) + "\n")
        for t in tasks:
            g = by_task[t].get("bigru", {})
            l = by_task[t].get("bilstm", {})

            def fnum(x: str) -> float | None:
                try:
                    return float(x)
                except (TypeError, ValueError):
                    return None

            g_pcc = fnum(g.get("Mean_PCC", ""))
            l_pcc = fnum(l.get("Mean_PCC", ""))
            g_ap = fnum(g.get("Peak_AP", ""))
            l_ap = fnum(l.get("Peak_AP", ""))
            d_pcc = (l_pcc - g_pcc) if g_pcc is not None and l_pcc is not None else ""
            d_ap = (l_ap - g_ap) if g_ap is not None and l_ap is not None else ""
            f.write(
                "\t".join(
                    [
                        t,
                        g.get("Mean_PCC", ""),
                        l.get("Mean_PCC", ""),
                        str(d_pcc),
                        g.get("Peak_AP", ""),
                        l.get("Peak_AP", ""),
                        str(d_ap),
                    ]
                )
                + "\n"
                )


def maybe_detach() -> None:
    if os.name != "posix":
        raise SystemExit("--detach requires POSIX (fork)")
    pid = os.fork()
    if pid > 0:
        print(f"[detach] child_pid={pid}", flush=True)
        sys.exit(0)


def main() -> None:
    p = argparse.ArgumentParser(description="BiLSTM vs BiGRU MARCP scheduler (workspace -> output on full success)")
    p.add_argument("--metrics_root", type=Path, default=_script_dir(), help="Contains workspace/ and output/")
    p.add_argument("--data_root", type=Path, default=_rbp_root() / "数据" / "h5data")
    p.add_argument("--n_tasks", type=int, default=20)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--exclude_tasks", type=str, default="", help="Space or comma separated task names")
    p.add_argument("--run_id", type=str, default="", help="Default: timestamp_seed")
    p.add_argument("--cuda_lstm", type=str, default="2")
    p.add_argument("--cuda_gru", type=str, default="3")
    p.add_argument("--epochs", type=int, default=100, help="Training epochs (default 100, README default)")
    p.add_argument("--train_chrs", type=str, default=TRAIN_CHRS_DEFAULT)
    p.add_argument("--val_chrs", type=str, default=VAL_CHRS_DEFAULT)
    p.add_argument("--test_chrs", type=str, default=TEST_CHRS_DEFAULT)
    p.add_argument("--dry_run", action="store_true")
    p.add_argument("--detach", action="store_true")
    p.add_argument("--log_file", type=Path, default=None)
    args = p.parse_args()

    metrics_root: Path = args.metrics_root.resolve()
    data_root: Path = args.data_root.resolve()
    marcp_root = _marcp_root()
    if not (marcp_root / MARCP_TRAIN).is_file():
        print(f"Missing {marcp_root / MARCP_TRAIN}", file=sys.stderr)
        sys.exit(1)

    run_id = args.run_id.strip()
    if not run_id:
        ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        run_id = f"{ts}_seed{args.seed}"

    if args.log_file:
        args.log_file.parent.mkdir(parents=True, exist_ok=True)
        logf = open(args.log_file, "a", encoding="utf-8")
        sys.stdout = TeeStream([sys.__stdout__, logf])
        sys.stderr = TeeStream([sys.__stderr__, logf])

    if args.detach:
        maybe_detach()

    exclude = parse_exclude(args.exclude_tasks)
    tasks = pick_tasks(
        data_root,
        args.n_tasks,
        args.seed,
        exclude,
        args.train_chrs,
        args.val_chrs,
        args.test_chrs,
    )

    if args.dry_run:
        train_csv_fn = lambda t: _file_list(data_root, t, args.train_chrs)
        val_csv_fn = lambda t: _file_list(data_root, t, args.val_chrs)
        test_csv_fn = lambda t: _file_list(data_root, t, args.test_chrs)
        for t in tasks:
            b = Path("/tmp/DRYRUN") / t / "bigru"
            run_marcp(
                build_marcp_cmd(
                    task=t,
                    train_csv=train_csv_fn(t),
                    val_csv=val_csv_fn(t),
                    test_csv=test_csv_fn(t),
                    out_dir=b,
                    cuda=args.cuda_gru,
                    use_lstm=False,
                    epochs=args.epochs,
                ),
                cwd=marcp_root,
                dry_run=True,
            )
        for t in tasks:
            b = Path("/tmp/DRYRUN") / t / "bilstm"
            run_marcp(
                build_marcp_cmd(
                    task=t,
                    train_csv=train_csv_fn(t),
                    val_csv=val_csv_fn(t),
                    test_csv=test_csv_fn(t),
                    out_dir=b,
                    cuda=args.cuda_lstm,
                    use_lstm=True,
                    epochs=args.epochs,
                ),
                cwd=marcp_root,
                dry_run=True,
            )
        print("DRY_RUN done (no workspace/output written).", flush=True)
        return

    ws_root = metrics_root / "workspace" / run_id
    out_root = metrics_root / "output" / run_id
    ws_root.parent.mkdir(parents=True, exist_ok=True)
    (metrics_root / "output").mkdir(parents=True, exist_ok=True)

    meta = {
        "run_id": run_id,
        "tasks": tasks,
        "seed": args.seed,
        "cuda_lstm": args.cuda_lstm,
        "cuda_gru": args.cuda_gru,
        "data_root": str(data_root),
        "epochs": args.epochs,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    ws_root.mkdir(parents=True, exist_ok=True)
    (ws_root / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    if args.detach:
        (ws_root / "orchestrator.pid").write_text(str(os.getpid()), encoding="utf-8")

    train_csv_fn = lambda t: _file_list(data_root, t, args.train_chrs)
    val_csv_fn = lambda t: _file_list(data_root, t, args.val_chrs)
    test_csv_fn = lambda t: _file_list(data_root, t, args.test_chrs)

    errors: list[str] = []
    err_lock = threading.Lock()

    def queue_gru() -> None:
        for t in tasks:
            out_base = ws_root / t / "bigru"
            out_base.mkdir(parents=True, exist_ok=True)
            cmd = build_marcp_cmd(
                task=t,
                train_csv=train_csv_fn(t),
                val_csv=val_csv_fn(t),
                test_csv=test_csv_fn(t),
                out_dir=out_base,
                cuda=args.cuda_gru,
                use_lstm=False,
                epochs=args.epochs,
            )
            print(f"\n=== BiGRU GPU{args.cuda_gru} {t} ===", flush=True)
            rc = run_marcp(cmd, cwd=marcp_root, dry_run=args.dry_run)
            if rc != 0:
                with err_lock:
                    errors.append(f"BiGRU {t} exit {rc}")

    def queue_lstm() -> None:
        for t in tasks:
            out_base = ws_root / t / "bilstm"
            out_base.mkdir(parents=True, exist_ok=True)
            cmd = build_marcp_cmd(
                task=t,
                train_csv=train_csv_fn(t),
                val_csv=val_csv_fn(t),
                test_csv=test_csv_fn(t),
                out_dir=out_base,
                cuda=args.cuda_lstm,
                use_lstm=True,
                epochs=args.epochs,
            )
            print(f"\n=== BiLSTM GPU{args.cuda_lstm} {t} ===", flush=True)
            rc = run_marcp(cmd, cwd=marcp_root, dry_run=args.dry_run)
            if rc != 0:
                with err_lock:
                    errors.append(f"BiLSTM {t} exit {rc}")

    with ThreadPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(queue_gru), ex.submit(queue_lstm)]
        for fu in as_completed(futs):
            fu.result()

    if errors:
        print("Some runs failed; leaving results under workspace (no move to output):", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)

    expected = expected_metric_paths(ws_root, tasks)
    miss = validate_metrics(expected)
    if miss:
        print("Metric files missing after successful exit codes; not moving to output:", file=sys.stderr)
        for m in miss:
            print(f"  - {m}", file=sys.stderr)
        sys.exit(1)

    if out_root.exists():
        print(f"Refusing to overwrite existing {out_root}", file=sys.stderr)
        sys.exit(1)

    try:
        os.rename(ws_root, out_root)
    except OSError:
        shutil.move(str(ws_root), str(out_root))

    write_summaries(out_root, tasks)
    (out_root / "RUN_COMPLETE.txt").write_text(
        f"completed_utc={datetime.now(timezone.utc).isoformat()}\n",
        encoding="utf-8",
    )
    print(f"\nDone. Final artifacts: {out_root}", flush=True)


if __name__ == "__main__":
    main()
