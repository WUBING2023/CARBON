#!/usr/bin/env python3
"""
实验一 Tier A：test 染色体 1/8/15 上，窗内 find_peaks + PhyloP 富集 + 匹配对照。

用法（在目录 保守性分数实验/实验一 下）：
  python run_exp1.py --output-dir ./outputs/run1
  # 默认约 140 个进程池 worker（整段 ex.map 期间常驻），可用环境变量 EXP1_WORKERS 覆盖
  python run_exp1.py --output-dir ./outputs/full
  python run_exp1.py --output-dir ./outputs/full --workers 70 --threads-per-worker 2
  python run_exp1.py --workers 1   # 单进程顺序（与旧 RNG 行为一致）
  python run_exp1.py --max-tasks 5 --no-bed --no-plots
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from exp1.mp_worker import init_thread_env, run_one_task
from exp1.paths import h5data_root, rbp_root
from exp1.pipeline import MetaplotAccumulator, run_task, summarize_task_peaks
from exp1.plots import plot_metaplot, plot_peak_vs_matched, plot_task_delta_heatmap
from exp1.task_ids import list_tasks, parse_rbp_cell_line


def _cross_task_stats(deltas: np.ndarray) -> dict:
    deltas = deltas[np.isfinite(deltas)]
    if deltas.size == 0:
        return {}
    try:
        w = stats.wilcoxon(deltas, zero_method="wilcox", alternative="two-sided", method="auto")
        wp, ws = float(w.pvalue), float(w.statistic)
    except ValueError:
        wp, ws = float("nan"), float("nan")
    return {
        "n_tasks": int(deltas.size),
        "median_delta_across_tasks": float(np.median(deltas)),
        "mean_delta_across_tasks": float(np.mean(deltas)),
        "fraction_delta_gt0": float(np.mean(deltas > 0)),
        "wilcoxon_vs0_statistic": ws,
        "wilcoxon_vs0_pvalue": wp,
    }


def _default_workers() -> int:
    """默认并行进程数；可用环境变量 EXP1_WORKERS 覆盖（默认 140）。"""
    return max(1, int(os.environ.get("EXP1_WORKERS", "140")))


def main() -> None:
    ap = argparse.ArgumentParser(description="实验一 Tier A：峰级保守性富集")
    ap.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "default_run")
    ap.add_argument("--h5data", type=Path, default=None, help="覆盖默认 h5data 根目录")
    ap.add_argument("--max-tasks", type=int, default=0, help="0 表示全部 task")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--no-bed", action="store_true", help="不用 test BED 子集，使用 test 染色体上全部窗")
    ap.add_argument("--sum-threshold", type=float, default=4.0, help="find_peaks signal_sum_threshold（默认 4，敏感性可试 2）")
    ap.add_argument("--n-controls", type=int, default=3, help="每峰匹配对照条数")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--meta-half", type=int, default=25, help="metaplot 半窗宽（bp）")
    ap.add_argument(
        "--workers",
        type=int,
        default=_default_workers(),
        help="进程池最大 worker 数（默认 140，或由 EXP1_WORKERS）；1 表示单进程顺序执行",
    )
    ap.add_argument(
        "--threads-per-worker",
        type=int,
        default=1,
        help="每个子进程内 BLAS/OpenMP 线程上限（默认 1；与 140 worker 搭配约 140 路并行）",
    )
    args = ap.parse_args()

    requested_workers = max(1, int(args.workers))
    if requested_workers > 1:
        init_thread_env(int(args.threads_per_worker))

    h5 = Path(args.h5data) if args.h5data else h5data_root()
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    tasks = list_tasks(h5)
    if args.max_tasks > 0:
        tasks = tasks[: args.max_tasks]

    h5s = str(h5.resolve())
    use_bed = not args.no_bed
    packed_list = [
        (tid, h5s, use_bed, float(args.sum_threshold), int(args.n_controls), int(args.seed), int(args.meta_half))
        for tid in tasks
    ]

    all_peak_rows: list[dict] = []
    task_summaries: list[dict] = []
    qc_tasks: list[dict] = []
    meta_acc = MetaplotAccumulator(half=args.meta_half)

    if requested_workers <= 1:
        rng = np.random.default_rng(args.seed)
        for ti, task_id in enumerate(tasks):
            rows, qc = run_task(
                task_id,
                h5,
                use_bed=use_bed,
                signal_sum_threshold=args.sum_threshold,
                n_controls=args.n_controls,
                rng=rng,
                meta=meta_acc,
            )
            qc["task_id"] = task_id
            qc_tasks.append(qc)
            all_peak_rows.extend(rows)
            pdf = pd.DataFrame(rows)
            summ = summarize_task_peaks(pdf)
            if summ is None:
                summ = {
                    "task_id": task_id,
                    "n_peaks": 0,
                    "delta_mean": float("nan"),
                    "delta_median": float("nan"),
                    "frac_delta_gt0": float("nan"),
                    "wilcoxon_statistic": float("nan"),
                    "wilcoxon_pvalue": float("nan"),
                }
            else:
                summ["task_id"] = task_id
            task_summaries.append(summ)
            if (ti + 1) % 20 == 0 or ti == len(tasks) - 1:
                print(f"[{ti+1}/{len(tasks)}] {task_id} peaks={summ.get('n_peaks', 0)}", flush=True)
    else:
        n = len(packed_list)
        # 不超过 task 数，避免 task 很少时仍 fork 满额空进程
        pool_workers = max(1, min(requested_workers, n))
        print(
            f"进程池: pool_workers={pool_workers}（请求 {requested_workers}，task 数 {n}），"
            f"threads_per_worker={int(args.threads_per_worker)}",
            flush=True,
        )
        # 进程池在 with 块内保持 worker 常驻并复用，直到全部 task 消费完毕
        with ProcessPoolExecutor(
            max_workers=pool_workers,
            initializer=init_thread_env,
            initargs=(int(args.threads_per_worker),),
        ) as ex:
            fut_to_task = {ex.submit(run_one_task, pack): pack[0] for pack in packed_list}
            for k, fut in enumerate(as_completed(fut_to_task), start=1):
                r = fut.result()
                all_peak_rows.extend(r["rows"])
                qc_tasks.append(r["qc"])
                task_summaries.append(r["summary"])
                meta_acc.sum_w += np.asarray(r["meta_sum"], dtype=np.float64)
                meta_acc.cnt_w += np.asarray(r["meta_cnt"], dtype=np.float64)
                step = 1
                if k % step == 0 or k == n:
                    print(
                        f"[{k}/{n}] {r['task_id']} peaks={r['summary'].get('n_peaks', 0)}",
                        flush=True,
                    )

    peak_df = pd.DataFrame(all_peak_rows)
    peak_path = out / "peak_table.csv"
    peak_df.to_csv(peak_path, index=False)

    ts_df = pd.DataFrame(task_summaries)
    ts_df["cell_line"] = ts_df["task_id"].map(lambda t: parse_rbp_cell_line(str(t))[1])
    ts_path = out / "task_summary.csv"
    ts_df.to_csv(ts_path, index=False)

    qc_df = pd.DataFrame(qc_tasks)
    qc_df.to_csv(out / "qc_per_task.csv", index=False)

    deltas = ts_df["delta_median"].to_numpy(dtype=np.float64)
    cross = _cross_task_stats(deltas)
    cross["rbp_root"] = str(rbp_root())
    cross["h5data"] = str(h5)
    cross["use_bed_filter"] = not args.no_bed
    cross["signal_sum_threshold"] = args.sum_threshold
    pool_workers_done = (
        max(1, min(requested_workers, len(packed_list))) if requested_workers > 1 else 1
    )
    cross["workers_requested"] = requested_workers
    cross["pool_workers"] = pool_workers_done
    cross["workers"] = pool_workers_done  # 兼容旧字段：实际进程池大小
    cross["threads_per_worker"] = int(args.threads_per_worker)
    (out / "cross_task_summary.json").write_text(json.dumps(cross, indent=2), encoding="utf-8")

    # 分层 K562 / HepG2
    for cl in ("K562", "HepG2"):
        sub = ts_df[ts_df["cell_line"] == cl]
        if not sub.empty:
            cd = _cross_task_stats(sub["delta_median"].to_numpy(dtype=np.float64))
            (out / f"cross_task_summary_{cl}.json").write_text(json.dumps(cd, indent=2), encoding="utf-8")

    if not args.no_plots and not peak_df.empty:
        plot_peak_vs_matched(peak_df, out / "fig_violin_peak_vs_matched.png")
        mx, my = meta_acc.curve()
        plot_metaplot(mx, my, out / "fig_metaplot_summit.png", title="Summit-centered mean PhyloP (Tier A)")
        plot_task_delta_heatmap(ts_df, out / "fig_heatmap_task_delta.png")

    print("Done.", peak_path, ts_path, flush=True)


if __name__ == "__main__":
    main()
