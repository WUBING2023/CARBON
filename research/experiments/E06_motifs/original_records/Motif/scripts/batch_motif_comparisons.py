#!/usr/bin/env python3
"""
对多个 outputs/<TASK>/（已有 tfmodisco_export/queries.meme）并行运行 CISBP-RNA Tomtom。

同阶段仅允许一个后台实例：默认占用 outputs/batch_motif_logs/locks/stage_c_cisbp_tomtom.lock
（DISABLE_STAGE_LOCK=1 可跳过）。

  cd 可解释性分析/Motif
  # 并行数 --jobs；每子进程 CPU 线程 --threads-per-job（或环境变量 COMPARISON_THREADS_PER_JOB）
  nohup python3 scripts/batch_motif_comparisons.py --jobs 70 --threads-per-job 2 --tomtom-binary ~/.local/share/mamba/envs/meme-suite/bin/tomtom \\
    >> outputs/batch_motif_logs/comparisons.log 2>&1 &
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from motif_analysis.stage_lock import acquire_lock, lock_disabled_from_env, release_lock

_DEFAULT_HS_PWM = ROOT / "reference_data" / "CISBP-RNA" / "homo_sapiens" / "pwms"
_DEFAULT_HS_RBP = (
    ROOT / "reference_data" / "CISBP-RNA" / "homo_sapiens" / "RBP_Information_all_motifs.txt"
)

SKIP_DIR_NAMES = frozenset(
    {
        "batch_motif_logs",
        "fixtures",
        "demo",
    }
)

DEFAULT_TOMTOM = Path.home() / ".local/share/mamba/envs/meme-suite/bin/tomtom"


def _thread_limited_env(base: dict[str, str], threads: int) -> dict[str, str]:
    t = str(max(1, int(threads)))
    out = {**base, "OMP_NUM_THREADS": t, "MKL_NUM_THREADS": t, "OPENBLAS_NUM_THREADS": t, "NUMEXPR_MAX_THREADS": t}
    return out


def _load_gene_map(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        rd = csv.DictReader(f)
        for row in rd:
            t = (row.get("task") or row.get("TASK") or "").strip()
            g = (row.get("gene") or row.get("Gene_name") or row.get("GENE") or "").strip()
            if t and g:
                out[t] = g
    return out


def _infer_gene(task_name: str) -> str | None:
    if "_" not in task_name:
        return None
    tok = task_name.split("_", 1)[0].strip()
    if len(tok) >= 2 and tok[0].isalpha():
        return tok
    return None


def _skip_task(task_dir: Path, *, skip_existing: bool) -> bool:
    if not skip_existing:
        return False
    cis = task_dir / "tomtom_cisbp" / "tomtom_out" / "tomtom.tsv"
    return cis.is_file()


def _discover_tasks(outputs_root: Path, limit_names: set[str] | None) -> list[Path]:
    out: list[Path] = []
    for child in sorted(outputs_root.iterdir()):
        if not child.is_dir() or child.name in SKIP_DIR_NAMES or child.name.startswith("."):
            continue
        if limit_names is not None and child.name not in limit_names:
            continue
        q = child / "tfmodisco_export" / "queries.meme"
        if q.is_file():
            out.append(child)
    return out


def _run_one(payload: dict) -> dict:
    task_dir = Path(payload["task_dir"])
    name = task_dir.name
    tomtom = payload["tomtom_binary"]
    python = payload["python"]
    gene = payload.get("expected_gene") or None

    r: dict = {"task": name, "cisbp_rc": None, "error": ""}

    threads = max(1, int(payload.get("threads_per_job") or 2))
    env = _thread_limited_env(
        {**os.environ, "PYTHONPATH": f"{ROOT}:{os.environ.get('PYTHONPATH', '')}"},
        threads,
    )
    cmd = [
        python,
        str(ROOT / "scripts" / "tomtom_modisco_vs_cisbp.py"),
        "--task-dir",
        str(task_dir),
        "--tomtom-binary",
        tomtom,
        "--cisbp-pwms",
        str(payload["cisbp_pwms"]),
        "--cisbp-rbp-table",
        str(payload["cisbp_rbp_table"]),
    ]
    if gene:
        cmd += ["--expected-gene", gene]
    p = subprocess.run(cmd, cwd=str(ROOT), env=env)
    r["cisbp_rc"] = p.returncode
    return r


def main() -> None:
    ap = argparse.ArgumentParser(description="Batch CISBP-RNA Tomtom for tasks with queries.meme")
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument(
        "--jobs",
        type=int,
        default=int(os.environ.get("COMPARISON_JOBS", "70")),
        help="并行子进程数（默认 70，或环境变量 COMPARISON_JOBS）",
    )
    ap.add_argument(
        "--threads-per-job",
        type=int,
        default=int(os.environ.get("COMPARISON_THREADS_PER_JOB", "2")),
        help="每个子任务内 tomtom/Python 的 OpenMP/BLAS 线程上限（默认 2）",
    )
    ap.add_argument(
        "--tomtom-binary",
        type=Path,
        default=None,
        help=f"tomtom path (default: {DEFAULT_TOMTOM} if exists, else PATH)",
    )
    ap.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip if tomtom_cisbp/tomtom_out/tomtom.tsv exists",
    )
    ap.add_argument("--gene-map", type=Path, default=None, help="CSV: task,gene")
    ap.add_argument(
        "--infer-gene",
        action="store_true",
        help="First token of task folder name as --expected-gene (e.g. RBFOX2_...)",
    )
    ap.add_argument("--cisbp-pwms", type=Path, default=_DEFAULT_HS_PWM)
    ap.add_argument("--cisbp-rbp-table", type=Path, default=_DEFAULT_HS_RBP)
    ap.add_argument("--task", action="append", dest="tasks", help="Only these task names (repeatable)")

    args, unknown = ap.parse_known_args()
    legacy = {"--attract-only", "--cisbp-only"}
    bad = [x for x in unknown if x in legacy]
    if bad:
        print(
            "[error] Options removed: " + " ".join(bad) + ". ATtRACT 已移除；本脚本仅跑 CISBP Tomtom。",
            file=sys.stderr,
        )
        raise SystemExit(2)
    if unknown:
        print("[error] unknown: " + " ".join(unknown), file=sys.stderr)
        raise SystemExit(2)

    tomtom = args.tomtom_binary
    if tomtom is None:
        tomtom = DEFAULT_TOMTOM if DEFAULT_TOMTOM.is_file() else None
    if tomtom is None:
        import shutil

        w = shutil.which("tomtom")
        if not w:
            print("[error] tomtom not found; pass --tomtom-binary", file=sys.stderr)
            raise SystemExit(2)
        tomtom_str = w
    else:
        tomtom_str = str(tomtom.resolve())
        if not Path(tomtom_str).is_file():
            print(f"[error] not a file: {tomtom_str}", file=sys.stderr)
            raise SystemExit(2)

    outputs_root = args.outputs_root.resolve()
    limit = set(args.tasks) if args.tasks else None
    tasks = _discover_tasks(outputs_root, limit)
    if not tasks:
        print("No tasks with tfmodisco_export/queries.meme found.", flush=True)
        return

    gene_map = _load_gene_map(args.gene_map.resolve()) if args.gene_map else {}

    work: list[dict] = []
    n_preskip = 0
    for td in tasks:
        gene = gene_map.get(td.name)
        if args.infer_gene and not gene:
            gene = _infer_gene(td.name)
        if gene:
            gene = str(gene)

        if _skip_task(td, skip_existing=args.skip_existing):
            n_preskip += 1
            print(f"[skip] {td.name} (already done)", flush=True)
            continue
        work.append(
            {
                "task_dir": str(td),
                "tomtom_binary": tomtom_str,
                "python": sys.executable,
                "expected_gene": gene,
                "cisbp_pwms": str(args.cisbp_pwms.resolve()),
                "cisbp_rbp_table": str(args.cisbp_rbp_table.resolve()),
                "threads_per_job": args.threads_per_job,
            }
        )

    if not work:
        print(f"[batch] nothing to run (pre-skipped={n_preskip})", flush=True)
        return

    log_dir = outputs_root / "batch_motif_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    fail_path = log_dir / f"comparisons_failures_{stamp}.csv"

    print(
        f"[batch] queued={len(work)} pre_skipped={n_preskip} workers={args.jobs} "
        f"threads_per_job={args.threads_per_job} tomtom={tomtom_str}",
        flush=True,
    )

    results: list[dict] = []
    failures: list[dict] = []

    lock_fd = None
    try:
        lock_fd = acquire_lock(ROOT, "stage_c_cisbp_tomtom", disable=lock_disabled_from_env())
        try:
            with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as ex:
                futs = {ex.submit(_run_one, p): p for p in work}
                for fut in as_completed(futs):
                    res = fut.result()
                    results.append(res)
                    name = res["task"]
                    if res.get("error"):
                        print(f"[fail] {name} {res['error']}", flush=True)
                        failures.append({"task": name, "step": "error", "rc": "", "detail": res["error"]})
                        continue
                    cr = res.get("cisbp_rc")
                    if cr is not None and cr != 0:
                        print(f"[fail] {name} cisbp_rc={cr}", flush=True)
                        failures.append({"task": name, "step": "cisbp", "rc": str(cr), "detail": ""})
                    else:
                        print(f"[ok] {name}", flush=True)
        except KeyboardInterrupt:
            print("[batch] interrupted", flush=True)

        if failures:
            with open(fail_path, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=["task", "step", "rc", "detail"])
                w.writeheader()
                w.writerows(failures)
            print(f"[batch] failures: {fail_path}", flush=True)

        n_ok = sum(1 for r in results if (r.get("cisbp_rc") in (None, 0)) and not r.get("error"))
        print(f"[batch] done ok={n_ok} pre_skipped={n_preskip} fail={len(failures)}", flush=True)
        raise SystemExit(0 if not failures else 1)
    finally:
        release_lock(lock_fd)


if __name__ == "__main__":
    main()
