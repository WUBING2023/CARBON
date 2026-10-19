#!/usr/bin/env python3
"""
实验四：Conservation-guided performance stratification。

用法（在目录 保守性分数实验/实验四 下）：
  python run_exp4.py --output-dir ./outputs/run1 --max-tasks 3
  python run_exp4.py --output-dir ./outputs/full --device cuda:0 --batch-size 4096
  python run_exp4.py --multi-gpu --gpu-ids 0,1,2,3 --workers-gpu0 10 --workers-per-other-gpu 6
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import multiprocessing as mp
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from exp4.conservation import (
    attach_motif_core_means,
    load_motif_core_table,
    peak_core_mean_phylop,
)
from exp4.marcp_infer import load_marcp, predict_profile_batched
from exp4.metrics import (
    assign_quantile_bins,
    flattened_signal_peak_ap,
    row_pearson,
    row_spearman,
)
from exp4.paths import h5data_root, weights_baseline, weights_no_conservative
from exp4.pipeline import collect_task_windows
from exp4.plots import plot_delta_pcc_by_bin, plot_delta_spearman_by_bin
from exp4.stats import cross_task_enrichment_vs_gain, safe_wilcoxon, spearman_trend
from exp4.task_ids import list_tasks, parse_rbp_cell_line


def _set_worker_thread_env(threads_per_worker: int) -> None:
    s = str(max(1, int(threads_per_worker)))
    for k in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[k] = s


def _summarize_bins_one_task(
    task_id: str,
    strat: str,
    bin_id: np.ndarray,
    delta_pcc: np.ndarray,
    delta_spr: np.ndarray,
    Y: np.ndarray,
    pred_m: np.ndarray,
    pred_n: np.ndarray,
    n_bins: int,
) -> tuple[list[dict], list[dict]]:
    """返回 (bin_summary_rows, trend_rows)。"""
    bin_rows: list[dict] = []
    trend_rows: list[dict] = []
    max_b = int(np.max(bin_id)) if np.any(bin_id >= 0) else -1
    med_d_pcc: list[float] = []
    med_d_spr: list[float] = []
    bin_ids_used: list[int] = []

    for b in range(max_b + 1):
        m = bin_id == b
        n_w = int(np.sum(m))
        if n_w == 0:
            continue
        dp = delta_pcc[m]
        ds = delta_spr[m]
        w_stat_p, w_p_p = safe_wilcoxon(dp)
        w_stat_s, w_p_s = safe_wilcoxon(ds)
        yf = Y[m].reshape(-1)
        pm = pred_m[m].reshape(-1)
        pn = pred_n[m].reshape(-1)
        sig_m, pk_m = flattened_signal_peak_ap(yf, pm)
        sig_n, pk_n = flattened_signal_peak_ap(yf, pn)
        md_pcc = float(np.nanmedian(dp))
        md_spr = float(np.nanmedian(ds))
        med_d_pcc.append(md_pcc)
        med_d_spr.append(md_spr)
        bin_ids_used.append(b)
        bin_rows.append(
            {
                "task_id": task_id,
                "stratification": strat,
                "bin_id": int(b),
                "n_bins_requested": int(n_bins),
                "n_windows": n_w,
                "median_delta_pcc": md_pcc,
                "median_delta_spearman": md_spr,
                "signal_ap_marcp": sig_m,
                "signal_ap_noc": sig_n,
                "delta_signal_ap": float(sig_m - sig_n) if np.isfinite(sig_m) and np.isfinite(sig_n) else float("nan"),
                "peak_ap_marcp": pk_m,
                "peak_ap_noc": pk_n,
                "delta_peak_ap": float(pk_m - pk_n) if np.isfinite(pk_m) and np.isfinite(pk_n) else float("nan"),
                "wilcoxon_stat_delta_pcc": w_stat_p,
                "wilcoxon_p_delta_pcc": w_p_p,
                "wilcoxon_stat_delta_spearman": w_stat_s,
                "wilcoxon_p_delta_spearman": w_p_s,
            }
        )

    if bin_ids_used:
        bx = np.asarray(bin_ids_used, dtype=np.float64)
        tr_pcc = spearman_trend(bx, np.asarray(med_d_pcc, dtype=np.float64))
        tr_spr = spearman_trend(bx, np.asarray(med_d_spr, dtype=np.float64))
        trend_rows.append(
            {
                "task_id": task_id,
                "stratification": strat,
                "metric": "median_delta_pcc_vs_bin",
                **tr_pcc,
            }
        )
        trend_rows.append(
            {
                "task_id": task_id,
                "stratification": strat,
                "metric": "median_delta_spearman_vs_bin",
                **tr_spr,
            }
        )
        d_sig = [r["delta_signal_ap"] for r in bin_rows]
        d_pk = [r["delta_peak_ap"] for r in bin_rows]
        if len(d_sig) == len(bin_ids_used):
            tr_sig = spearman_trend(bx, np.asarray(d_sig, dtype=np.float64))
            trend_rows.append({"task_id": task_id, "stratification": strat, "metric": "delta_signal_ap_vs_bin", **tr_sig})
        if len(d_pk) == len(bin_ids_used):
            tr_pk = spearman_trend(bx, np.asarray(d_pk, dtype=np.float64))
            trend_rows.append({"task_id": task_id, "stratification": strat, "metric": "delta_peak_ap_vs_bin", **tr_pk})

    return bin_rows, trend_rows


def _run_single_task(
    task_id: str,
    *,
    h5: Path,
    use_bed: bool,
    n_bins: int,
    sum_threshold: float,
    batch_size: int,
    device: torch.device,
    motif_df: pd.DataFrame,
    windows_tmp_dir: Path,
) -> dict:
    wb = weights_baseline(task_id)
    wn = weights_no_conservative(task_id)
    if not wb.is_file() or not wn.is_file():
        return {"status": "skip", "task_id": task_id, "reason": "missing_weights", "baseline": str(wb), "no_conservative": str(wn)}

    packed = collect_task_windows(task_id, h5, use_bed=use_bed)
    if packed is None:
        return {"status": "skip", "task_id": task_id, "reason": "no_test_windows"}

    meta, X, Y = packed
    n = X.shape[0]
    if n == 0:
        return {"status": "skip", "task_id": task_id, "reason": "empty_windows"}

    m5 = load_marcp(5, wb, device)
    m4 = load_marcp(4, wn, device)
    pred_m = predict_profile_batched(m5, X, device, batch_size=int(batch_size))
    pred_n = predict_profile_batched(m4, X[..., :4], device, batch_size=int(batch_size))
    del m5, m4
    if device.type == "cuda":
        torch.cuda.empty_cache()

    win_phy = np.asarray(X[:, :, 4].mean(axis=1), dtype=np.float64)
    peak_phy = np.array(
        [
            peak_core_mean_phylop(
                Y[i],
                X[i, :, 4],
                signal_max_threshold=2.0,
                signal_sum_threshold=float(sum_threshold),
                signal_length=1,
            )
            for i in range(n)
        ],
        dtype=np.float64,
    )
    if not motif_df.empty:
        motif_phy = attach_motif_core_means(meta, motif_df, X)
    else:
        motif_phy = np.full(n, np.nan, dtype=np.float64)

    pcc_m = row_pearson(pred_m, Y)
    pcc_n = row_pearson(pred_n, Y)
    spr_m = row_spearman(pred_m, Y)
    spr_n = row_spearman(pred_n, Y)
    delta_pcc = pcc_m - pcc_n
    delta_spr = spr_m - spr_n

    rbp, cell_line = parse_rbp_cell_line(task_id)
    wdf = pd.DataFrame(
        {
            "task_id": task_id,
            "rbp": rbp,
            "cell_line": cell_line,
            "window_mean_phyloP": win_phy,
            "peak_core_mean_phyloP": peak_phy,
            "motif_core_mean_phyloP": motif_phy,
            "pcc_marcp": pcc_m,
            "pcc_no_conservation": pcc_n,
            "delta_pcc": delta_pcc,
            "spearman_marcp": spr_m,
            "spearman_no_conservation": spr_n,
            "delta_spearman": delta_spr,
        }
    )
    wdf = pd.concat([meta.reset_index(drop=True), wdf], axis=1)

    all_bin_rows: list[dict] = []
    all_trend_rows: list[dict] = []
    top_mask = np.zeros(n, dtype=bool)
    window_bins_ok = False
    strats: list[tuple[str, np.ndarray]] = [
        ("window_phyloP", win_phy),
        ("peak_core_phyloP", peak_phy),
    ]
    if np.any(np.isfinite(motif_phy)):
        strats.append(("motif_core_phyloP", motif_phy))

    for strat_name, strat_vals in strats:
        bin_id, n_eff = assign_quantile_bins(strat_vals, int(n_bins))
        if n_eff <= 0:
            continue
        if strat_name == "window_phyloP":
            top_mask = top_mask | (bin_id == (n_eff - 1))
            window_bins_ok = True
        br, tr = _summarize_bins_one_task(
            task_id,
            strat_name,
            bin_id,
            delta_pcc,
            delta_spr,
            Y,
            pred_m,
            pred_n,
            int(n_bins),
        )
        all_bin_rows.extend(br)
        all_trend_rows.extend(tr)
        wdf[f"bin_{strat_name}"] = bin_id

    enrich = float(np.mean(top_mask)) if window_bins_ok else float("nan")
    ablation_gain = float(np.nanmean(delta_pcc))
    task_level_row = {
        "task_id": task_id,
        "rbp": rbp,
        "cell_line": cell_line,
        "n_windows": int(n),
        "conservation_enrichment": enrich,
        "ablation_gain_mean_delta_pcc": ablation_gain,
    }

    windows_tmp_dir.mkdir(parents=True, exist_ok=True)
    task_window_path = windows_tmp_dir / f"{task_id}.window_metrics.csv.gz"
    wdf.to_csv(task_window_path, index=False, compression="gzip")

    return {
        "status": "ok",
        "task_id": task_id,
        "n_windows": int(n),
        "window_file": str(task_window_path),
        "bin_rows": all_bin_rows,
        "trend_rows": all_trend_rows,
        "task_level_row": task_level_row,
    }


def _run_task_queue_on_gpu(
    gpu_id: int,
    task_ids: list[str],
    *,
    h5_str: str,
    use_bed: bool,
    n_bins: int,
    sum_threshold: float,
    batch_size: int,
    cpu_threads_per_worker: int,
    motif_core_table: str | None,
    windows_tmp_dir: str,
) -> list[dict]:
    """
    单进程独占一张 GPU，在进程内串行跑完 task_ids。
    避免同卡多进程各自加载模型导致显存爆炸（OOM）。
    """
    _set_worker_thread_env(cpu_threads_per_worker)
    os.environ["CUDA_VISIBLE_DEVICES"] = str(int(gpu_id))
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    motif_df = pd.DataFrame()
    if motif_core_table:
        motif_df = load_motif_core_table(Path(motif_core_table))

    out_list: list[dict] = []
    for task_id in task_ids:
        r = _run_single_task(
            task_id,
            h5=Path(h5_str),
            use_bed=bool(use_bed),
            n_bins=int(n_bins),
            sum_threshold=float(sum_threshold),
            batch_size=int(batch_size),
            device=device,
            motif_df=motif_df,
            windows_tmp_dir=Path(windows_tmp_dir),
        )
        r["gpu_id"] = int(gpu_id)
        out_list.append(r)
    return out_list


def _weighted_gpu_schedule(gpu_ids: list[int], w0: int, w_other: int) -> list[int]:
    """构造循环调度表：含 GPU0 时 GPU0 权重 w0、其余卡 w_other；不含 0 时各卡均 w_other。"""
    w0 = max(1, int(w0))
    w_other = max(1, int(w_other))
    sched: list[int] = []
    if 0 in gpu_ids:
        for g in gpu_ids:
            if int(g) == 0:
                sched.extend([g] * w0)
            else:
                sched.extend([g] * w_other)
    else:
        for g in gpu_ids:
            sched.extend([g] * w_other)
    return sched


def _distribute_tasks_to_gpus(
    tasks: list[str],
    gpu_ids: list[int],
    w0: int,
    w_other: int,
) -> dict[int, list[str]]:
    sched = _weighted_gpu_schedule(gpu_ids, w0, w_other)
    buckets: dict[int, list[str]] = {g: [] for g in gpu_ids}
    for i, tid in enumerate(tasks):
        g = sched[i % len(sched)]
        buckets[g].append(tid)
    return buckets


def _parse_gpu_ids(spec: str) -> list[int]:
    vals: list[int] = []
    for t in str(spec).split(","):
        t = t.strip()
        if not t:
            continue
        vals.append(int(t))
    return vals


def main() -> None:
    ap = argparse.ArgumentParser(description="实验四：conservation 分层与 MARCP 消融差")
    ap.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "default_run")
    ap.add_argument("--h5data", type=Path, default=None)
    ap.add_argument("--device", type=str, default="cuda:0")
    ap.add_argument("--batch-size", type=int, default=4096)
    ap.add_argument(
        "--multi-gpu",
        action="store_true",
        help="多 GPU：每卡一个进程并行，单卡内串行 task（避免同卡多进程 OOM）；跨卡并行",
    )
    ap.add_argument("--gpu-ids", type=str, default="", help="参与并发的 GPU 编号，如 0,1,2,3；为空时自动用全部可见 GPU")
    ap.add_argument(
        "--workers-gpu0",
        type=int,
        default=10,
        help="多 GPU 下任务分配权重：GPU0 相对其他卡的权重（不是并发进程数）",
    )
    ap.add_argument(
        "--workers-per-other-gpu",
        type=int,
        default=6,
        help="多 GPU 下任务分配权重：非 0 号 GPU 的权重（不是并发进程数）",
    )
    ap.add_argument("--cpu-threads-per-worker", type=int, default=4, help="每个进程 CPU 线程数")
    ap.add_argument("--max-tasks", type=int, default=0, help="0 表示全部 task")
    ap.add_argument("--no-bed", action="store_true", help="不用 peak_test BED，使用 chr1/8/15 全部窗")
    ap.add_argument("--n-bins", type=int, default=4)
    ap.add_argument(
        "--sum-threshold",
        type=float,
        default=4.0,
        help="窗内 peak 检测 signal_sum_threshold（与实验一 Tier A 默认一致）",
    )
    ap.add_argument("--motif-core-table", type=Path, default=None, help="可选 motif 窗内坐标表")
    args = ap.parse_args()

    h5 = Path(args.h5data) if args.h5data else h5data_root()
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    fig_dir = out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        print("CUDA 不可用，改用 cpu", flush=True)
        device = torch.device("cpu")

    tasks = list_tasks(h5)
    if int(args.max_tasks) > 0:
        tasks = tasks[: int(args.max_tasks)]

    all_window_frames: list[pd.DataFrame] = []
    all_bin_rows: list[dict] = []
    all_trend_rows: list[dict] = []
    qc_missing: list[dict] = []
    task_level: list[dict] = []
    windows_tmp_dir = out / "tmp_task_windows"
    results: list[dict] = []

    if args.multi_gpu:
        if args.gpu_ids.strip():
            gpu_ids = _parse_gpu_ids(args.gpu_ids)
        else:
            gpu_ids = list(range(torch.cuda.device_count()))
        if not gpu_ids:
            raise RuntimeError("启用 --multi-gpu 但未检测到可用 GPU。请设置 --gpu-ids 或检查 CUDA。")

        print(
            "[multi-gpu] 每块 GPU 仅 1 个 Python 进程，进程内串行处理 task，"
            f"避免同卡多进程重复加载模型导致 OOM；任务按权重轮询分配（GPU0={int(args.workers_gpu0)}, 其他={int(args.workers_per_other_gpu)}）。",
            flush=True,
        )

        task_buckets = _distribute_tasks_to_gpus(
            tasks,
            gpu_ids,
            int(args.workers_gpu0),
            int(args.workers_per_other_gpu),
        )
        n_submit = sum(1 for g in gpu_ids if task_buckets.get(g))
        spawn_ctx = mp.get_context("spawn")
        ex = cf.ProcessPoolExecutor(max_workers=max(1, n_submit), mp_context=spawn_ctx)
        pools = [ex]
        futures: list[cf.Future] = []
        for g in gpu_ids:
            tids = task_buckets.get(g, [])
            if not tids:
                continue
            futures.append(
                ex.submit(
                    _run_task_queue_on_gpu,
                    int(g),
                    tids,
                    h5_str=str(h5.resolve()),
                    use_bed=not args.no_bed,
                    n_bins=int(args.n_bins),
                    sum_threshold=float(args.sum_threshold),
                    batch_size=int(args.batch_size),
                    cpu_threads_per_worker=int(args.cpu_threads_per_worker),
                    motif_core_table=str(args.motif_core_table) if args.motif_core_table else None,
                    windows_tmp_dir=str(windows_tmp_dir),
                )
            )
        done = 0
        total = len(tasks)
        try:
            for fu in cf.as_completed(futures):
                batch = fu.result()
                for r in batch:
                    results.append(r)
                    done += 1
                    if r.get("status") == "ok":
                        print(f"[{done}/{total}] {r['task_id']} gpu={r.get('gpu_id')} n={r['n_windows']}", flush=True)
                    else:
                        print(f"[{done}/{total}] [skip] {r.get('task_id')} reason={r.get('reason')}", flush=True)
        finally:
            for pool in pools:
                pool.shutdown(wait=True)
    else:
        motif_df = pd.DataFrame()
        if args.motif_core_table is not None:
            motif_df = load_motif_core_table(Path(args.motif_core_table))
        _set_worker_thread_env(int(args.cpu_threads_per_worker))
        for ti, task_id in enumerate(tasks, start=1):
            r = _run_single_task(
                task_id,
                h5=h5,
                use_bed=not args.no_bed,
                n_bins=int(args.n_bins),
                sum_threshold=float(args.sum_threshold),
                batch_size=int(args.batch_size),
                device=device,
                motif_df=motif_df,
                windows_tmp_dir=windows_tmp_dir,
            )
            results.append(r)
            if r.get("status") == "ok":
                print(f"[{ti}/{len(tasks)}] {r['task_id']} n={r['n_windows']}", flush=True)
            else:
                print(f"[{ti}/{len(tasks)}] [skip] {r.get('task_id')} reason={r.get('reason')}", flush=True)

    window_files: list[str] = []
    for r in results:
        if r.get("status") != "ok":
            qc_missing.append(
                {
                    "task_id": r.get("task_id"),
                    "reason": r.get("reason", "unknown"),
                    "baseline": r.get("baseline", ""),
                    "no_conservative": r.get("no_conservative", ""),
                }
            )
            continue
        window_files.append(str(r["window_file"]))
        all_bin_rows.extend(r.get("bin_rows", []))
        all_trend_rows.extend(r.get("trend_rows", []))
        task_level.append(r.get("task_level_row", {}))

    if window_files:
        all_window_frames = [pd.read_csv(p, compression="gzip") for p in window_files]

    if all_window_frames:
        window_df = pd.concat(all_window_frames, ignore_index=True)
        window_path = out / "window_metrics.csv.gz"
        window_df.to_csv(window_path, index=False, compression="gzip")
    else:
        window_path = None

    if all_bin_rows:
        bin_df = pd.DataFrame(all_bin_rows)
        bin_df.to_csv(out / "bin_summary.csv", index=False)
    else:
        bin_df = pd.DataFrame()

    trend_path = out / "task_trends.json"
    with open(trend_path, "w", encoding="utf-8") as f:
        json.dump(all_trend_rows, f, indent=2, ensure_ascii=False)

    task_table = pd.DataFrame(task_level)
    if not task_table.empty:
        task_table.to_csv(out / "task_level_scores.csv", index=False)
    if not task_table.empty and "ablation_gain_mean_delta_pcc" in task_table.columns:
        cross = cross_task_enrichment_vs_gain(
            task_table.rename(columns={"ablation_gain_mean_delta_pcc": "ablation_gain"})
        )
    else:
        cross = cross_task_enrichment_vs_gain(
            pd.DataFrame(columns=["conservation_enrichment", "ablation_gain"])
        )

    cross_path = out / "cross_task_correlation.json"
    with open(cross_path, "w", encoding="utf-8") as f:
        json.dump(cross, f, indent=2, ensure_ascii=False)

    if qc_missing:
        pd.DataFrame(qc_missing).to_csv(out / "qc_missing_weights.csv", index=False)

    meta_doc = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "h5data": str(h5.resolve()),
        "device": str(device),
        "multi_gpu": bool(args.multi_gpu),
        "multi_gpu_policy": "每 GPU 单进程串行 task；workers_* 为任务分配权重，不是同卡并发进程数",
        "gpu_ids": args.gpu_ids if args.gpu_ids else "auto",
        "task_weight_gpu0": int(args.workers_gpu0),
        "task_weight_other_gpu": int(args.workers_per_other_gpu),
        "cpu_threads_per_worker": int(args.cpu_threads_per_worker),
        "n_bins": int(args.n_bins),
        "use_bed": not bool(args.no_bed),
        "flattened_ap_note": (
            "bin 内 Signal AP / Peak AP：将该 bin 内所有窗口的 Y 与 pred 展平为一维向量后，"
            "定义与 MARCP/MARCP/operations/train_ops.py calculate_batch_metrics_whole 一致 "
            "(Signal: y>median; Peak: find_peaks on y_flat with thresholds 2/2)。"
        ),
        "peak_core_note": "无峰窗 peak_core_mean_phyloP 为 nan，不参与 peak_core_phyloP 分箱。",
        "conservation_enrichment_note": "window_phyloP 分层下落在最高 conservation bin 的窗口比例。",
        "ablation_gain_note": "全 test 窗上 mean(delta_pcc)。",
    }
    with open(out / "run_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta_doc, f, indent=2, ensure_ascii=False)

    if not bin_df.empty:
        plot_delta_pcc_by_bin(bin_df, fig_dir / "delta_pcc_by_bin_window.png", strat="window_phyloP")
        plot_delta_spearman_by_bin(bin_df, fig_dir / "delta_spearman_by_bin_window.png", strat="window_phyloP")
        if np.any(bin_df["stratification"] == "peak_core_phyloP"):
            plot_delta_pcc_by_bin(bin_df, fig_dir / "delta_pcc_by_bin_peakcore.png", strat="peak_core_phyloP")

    print(f"Done. windows={window_path} bins={len(all_bin_rows)} tasks_scored={len(task_level)}", flush=True)


if __name__ == "__main__":
    main()
