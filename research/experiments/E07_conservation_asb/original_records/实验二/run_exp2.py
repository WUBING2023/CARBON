#!/usr/bin/env python3
"""
实验二：TF-MoDISco seqlet 中心邻域 PhyloP 分析。

用法（在目录 保守性分数实验/实验二 下）：
  python3 run_exp2.py --output-dir ./outputs/run1
  python3 run_exp2.py --output-dir ./outputs/full --workers 70 --threads-per-worker 2
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from exp2.mp_worker import init_thread_env, run_one_task, task_has_modisco
from exp2.paths import h5data_root, motif_outputs_root, motif_task_dir
from exp2.plots import (
    plot_core_vs_flank_violin,
    plot_ic_cwm_vs_phylop_scatter,
    plot_metaprofile,
    plot_motif_position_heatmap,
)
from exp2.stats_agg import (
    cross_task_ci_stats,
    instance_core_vs_flank_wilcoxon,
    motif_ci_summary,
    permutation_spearman_labels,
    position_spearman,
)
from exp2.task_ids import list_tasks, parse_rbp_cell_line


def main() -> None:
    ap = argparse.ArgumentParser(description="实验二：Motif-centered conservation")
    ap.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "default_run")
    ap.add_argument("--h5data", type=Path, default=None)
    ap.add_argument("--motif-root", type=Path, default=None, help="Motif outputs 根目录，默认 RBP/可解释性分析/Motif/outputs")
    ap.add_argument("--max-tasks", type=int, default=0, help="0 表示仅跑含 modisco h5 的全部 task")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--no-bed", action="store_true", help="不按 test BED 过滤窗（仍用 BED 做 overlaps_peak 注释）")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--flank-bp", type=int, default=10)
    ap.add_argument("--meta-half", type=int, default=25)
    ap.add_argument("--min-instances-per-motif", type=int, default=5)
    ap.add_argument("--n-perm", type=int, default=200, help="position 标签 permutation 次数")
    ap.add_argument("--workers", type=int, default=70)
    ap.add_argument("--threads-per-worker", type=int, default=2)
    args = ap.parse_args()

    if int(args.workers) > 1:
        init_thread_env(int(args.threads_per_worker))

    h5 = Path(args.h5data) if args.h5data else h5data_root()
    mroot = Path(args.motif_root) if args.motif_root else motif_outputs_root()
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    tasks_all = list_tasks(h5)
    tasks = [t for t in tasks_all if task_has_modisco(mroot / t)]
    if args.max_tasks > 0:
        tasks = tasks[: args.max_tasks]

    h5s = str(h5.resolve())
    mroots = str(mroot.resolve())
    use_bed = not args.no_bed
    packed_list = [
        (
            tid,
            h5s,
            str((mroot / tid).resolve()),
            use_bed,
            int(args.flank_bp),
            int(args.meta_half),
            int(args.min_instances_per_motif),
            int(args.seed),
        )
        for tid in tasks
    ]

    all_inst: list[pd.DataFrame] = []
    all_pos: list[pd.DataFrame] = []
    qc_rows: list[dict] = []
    meta_sum = None
    meta_cnt = None
    meta_half = int(args.meta_half)

    if args.workers <= 1:
        for ti, pack in enumerate(packed_list):
            r = run_one_task(pack)
            meta_sum, meta_cnt = _merge_meta(meta_sum, meta_cnt, r, meta_half)
            all_inst.append(r["instance_df"])
            all_pos.append(r["position_df"])
            qc_rows.append(r["qc"])
            if (ti + 1) % 10 == 0 or ti == len(packed_list) - 1:
                print(f"[{ti+1}/{len(packed_list)}] {pack[0]}", flush=True)
    else:
        n = len(packed_list)
        with ProcessPoolExecutor(
            max_workers=int(args.workers),
            initializer=init_thread_env,
            initargs=(int(args.threads_per_worker),),
        ) as ex:
            for k, r in enumerate(ex.map(run_one_task, packed_list, chunksize=1), start=1):
                meta_sum, meta_cnt = _merge_meta(meta_sum, meta_cnt, r, meta_half)
                all_inst.append(r["instance_df"])
                all_pos.append(r["position_df"])
                qc_rows.append(r["qc"])
                if k % 10 == 0 or k == n:
                    print(f"[{k}/{n}] {r['qc'].get('task_id', '')}", flush=True)

    inst_df = pd.concat(all_inst, ignore_index=True) if all_inst else pd.DataFrame()
    pos_df = pd.concat(all_pos, ignore_index=True) if all_pos else pd.DataFrame()
    qc_df = pd.DataFrame(qc_rows) if qc_rows else pd.DataFrame()

    inst_path = out / "motif_instance_table.csv"
    inst_df.to_csv(inst_path, index=False)
    pos_df.to_csv(out / "motif_position_table.csv", index=False)
    qc_df.to_csv(out / "qc_per_task.csv", index=False)

    if not inst_df.empty:
        inst_df = inst_df.copy()
        if "cell_line" not in inst_df.columns:
            inst_df["cell_line"] = inst_df["task_id"].map(lambda t: parse_rbp_cell_line(str(t))[1])
        inst_df["CI"] = inst_df["mean_phyloP_core"] - inst_df["mean_phyloP_flank"]
    cross: dict = {
        "h5data": str(h5),
        "motif_root": str(mroot),
        "use_bed_filter_windows": use_bed,
        "workers": int(args.workers),
        "threads_per_worker": int(args.threads_per_worker),
    }
    if not inst_df.empty:
        cross["instance_core_flank_wilcoxon"] = instance_core_vs_flank_wilcoxon(inst_df)
        cross["cross_task_CI"] = cross_task_ci_stats(inst_df["CI"].to_numpy(dtype=np.float64))
        cross["position_spearman_global"] = position_spearman(pos_df)
        cross["position_perm"] = permutation_spearman_labels(pos_df, n_perm=int(args.n_perm), seed=int(args.seed))

    msum = meta_sum if meta_sum is not None else np.zeros(2 * meta_half + 1)
    mcnt = meta_cnt if meta_cnt is not None else np.zeros(2 * meta_half + 1)
    with np.errstate(invalid="ignore"):
        meta_y = np.where(mcnt > 0, msum / np.maximum(mcnt, 1e-9), np.nan)
    meta_x = np.arange(-meta_half, meta_half + 1, dtype=np.int32)
    cross["metaprofile_x"] = meta_x.tolist()
    cross["metaprofile_y"] = [float(v) if v == v else None for v in meta_y]

    for cl in ("K562", "HepG2"):
        sub = inst_df[inst_df["cell_line"] == cl] if "cell_line" in inst_df.columns else pd.DataFrame()
        if not sub.empty:
            cross[f"cross_task_CI_{cl}"] = cross_task_ci_stats(sub["CI"].to_numpy(dtype=np.float64))

    (out / "cross_task_summary.json").write_text(json.dumps(cross, indent=2), encoding="utf-8")

    msum_m = motif_ci_summary(inst_df)
    msum_m.to_csv(out / "motif_level_CI_summary.csv", index=False)

    if not args.no_plots and not inst_df.empty:
        plot_core_vs_flank_violin(inst_df, out / "fig_violin_core_vs_flank.png")
        plot_metaprofile(meta_x, meta_y, out / "fig_metaprofile_motif_centered.png", title="Motif-centered mean PhyloP")
    if not args.no_plots and not pos_df.empty:
        plot_motif_position_heatmap(pos_df, out / "fig_heatmap_motif_position.png")
        plot_ic_cwm_vs_phylop_scatter(pos_df, out / "fig_scatter_cwm_vs_phylop.png")

    print("Done.", inst_path, flush=True)


def _merge_meta(
    meta_sum: np.ndarray | None,
    meta_cnt: np.ndarray | None,
    r: dict,
    meta_half: int,
) -> tuple[np.ndarray, np.ndarray]:
    w = 2 * int(meta_half) + 1
    s = np.asarray(r["meta_sum"], dtype=np.float64).ravel()
    c = np.asarray(r["meta_cnt"], dtype=np.float64).ravel()
    if s.size != w:
        s = np.resize(s, w)
        c = np.resize(c, w)
    if meta_sum is None:
        return s.copy(), c.copy()
    return meta_sum + s, meta_cnt + c


if __name__ == "__main__":
    main()
