#!/usr/bin/env python3
"""
对某任务 attributions 运行 TF-MoDISco（``modisco motifs`` + ``modisco report-simple``）。

优先从 tfmodisco_inputs.h5 读取 hyp_scores；
否则从 ig_hyp_4ch.npz 读 **gxi**（原始 grad×input，与当前主流水线一致，不做每位置去均值）。

用法：
  cd 可解释性分析/Motif
  python3 scripts/run_tfmodisco.py --task-dir outputs/AGGF1_HepG2_ENCSR543TPH

依赖：pip install modisco h5py hdf5plugin
若 import modisco 因 memelite+numba 报错，见 docs/tfmodisco_环境.md

标准 MARCP 输入勿使用 ``--modisco-input-permute-augc-to-acgt``。若此前误开旧版无条件 permute，
请重跑本脚本以刷新 ``modisco_results.h5`` 与 ``report_simple/``。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from motif_analysis.tfmodisco_cli import (
    RELAXED_MOTIF_PRESET,
    choose_modisco_window,
    run_modisco_motifs_and_report,
)


def _load_gxi_ohe(att_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    h5_path = att_dir / "tfmodisco_inputs.h5"
    if h5_path.is_file():
        import h5py

        with h5py.File(h5_path, "r") as f:
            g = f["task0"]
            gxi = g["hyp_scores"][:].astype(np.float32)
            ohe = g["one_hot"][:].astype(np.float32)
        return gxi, ohe

    npz_path = att_dir / "ig_hyp_4ch.npz"
    if not npz_path.is_file():
        raise FileNotFoundError(f"缺少 {h5_path} 与 {npz_path}")
    z = np.load(npz_path)
    gxi = z["gxi"].astype(np.float32)
    ohe = z["ohe"].astype(np.float32)
    if np.all(gxi >= 0):
        print(
            "[warn] hyp 全为非负；MoDISco 通常需要有符号贡献。若结果异常，请检查归因是否为 grad×input。",
            file=sys.stderr,
        )
    return gxi, ohe


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--task-dir",
        type=Path,
        required=True,
        help="任务目录（含 attributions/）",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="默认 task-dir/tfmodisco",
    )
    ap.add_argument(
        "--max-seqlets",
        type=int,
        default=None,
        help="modisco -n；未指定时默认 4000，若与 --relaxed-motifs 联用则为 4500",
    )
    ap.add_argument("--window", type=int, default=None, help="覆盖自动窗口长度")
    ap.add_argument("--n-leiden", type=int, default=None, help="modisco motifs -l")
    ap.add_argument("--sliding-window-size", type=int, default=None, help="modisco motifs -z")
    ap.add_argument("--trim-size", type=int, default=None, help="modisco motifs -t（trim_to_window_size）")
    ap.add_argument("--seqlet-flank-size", type=int, default=None, help="modisco motifs -f")
    ap.add_argument("--initial-flank-to-add", type=int, default=None, help="modisco motifs -g")
    ap.add_argument("--final-flank-to-add", type=int, default=None, help="modisco motifs -j")
    ap.add_argument(
        "--modisco-report-trim-threshold",
        type=float,
        default=0.28,
        help="已保留兼容：report-simple 的 trim 由 MEME 默认决定，本参数不生效",
    )
    ap.add_argument(
        "--relaxed-motifs",
        action="store_true",
        help=(
            "略宽于 MEME 默认的 motifs 参数："
            + ", ".join(f"{k}={v}" for k, v in RELAXED_MOTIF_PRESET.items())
            + "（显式传入的对应参数优先）"
        ),
    )
    ap.add_argument(
        "--strict-motifs",
        action="store_true",
        help="已弃用，等同于 --relaxed-motifs（历史名称）",
    )
    ap.add_argument("--skip-report", action="store_true")
    ap.add_argument(
        "--modisco-input-permute-augc-to-acgt",
        action="store_true",
        help=(
            "仅当 hyp/one-hot 最后一维按 A,U,G,C 排列（非标准 MARCP mapper 顺序）时启用；"
            "标准 MARCP one-hot 为 A,C,G,U，勿开此开关"
        ),
    )
    args = ap.parse_args()

    if args.strict_motifs:
        print("[warn] --strict-motifs 已弃用，请改用 --relaxed-motifs", file=sys.stderr)
        args.relaxed_motifs = True

    if args.relaxed_motifs:
        if args.n_leiden is None:
            args.n_leiden = RELAXED_MOTIF_PRESET["n_leiden"]
        if args.sliding_window_size is None:
            args.sliding_window_size = RELAXED_MOTIF_PRESET["sliding_window_size"]
        if args.trim_size is None:
            args.trim_size = RELAXED_MOTIF_PRESET["trim_size"]
        if args.seqlet_flank_size is None:
            args.seqlet_flank_size = RELAXED_MOTIF_PRESET["seqlet_flank_size"]
        if args.initial_flank_to_add is None:
            args.initial_flank_to_add = RELAXED_MOTIF_PRESET["initial_flank_to_add"]

    if args.max_seqlets is None:
        args.max_seqlets = (
            RELAXED_MOTIF_PRESET["max_seqlets"] if args.relaxed_motifs else 4000
        )

    task_dir = args.task_dir.resolve()
    att_dir = task_dir / "attributions"
    out_dir = (args.out_dir or (task_dir / "tfmodisco")).resolve()

    gxi, ohe = _load_gxi_ohe(att_dir)
    if gxi.shape != ohe.shape or gxi.ndim != 3 or gxi.shape[-1] != 4:
        raise ValueError(f"bad shapes gxi={gxi.shape} ohe={ohe.shape}")

    seq_len = gxi.shape[1]
    window = args.window if args.window is not None else choose_modisco_window(seq_len)

    code = run_modisco_motifs_and_report(
        gxi_nlc4=gxi,
        ohe_nlc4=ohe,
        out_dir=out_dir,
        max_seqlets=args.max_seqlets,
        window=window,
        skip_report=args.skip_report,
        report_trim_threshold=args.modisco_report_trim_threshold,
        n_leiden=args.n_leiden,
        sliding_window_size=args.sliding_window_size,
        trim_size=args.trim_size,
        seqlet_flank_size=args.seqlet_flank_size,
        initial_flank_to_add=args.initial_flank_to_add,
        final_flank_to_add=args.final_flank_to_add,
        modisco_input_permute_augc_to_acgt=args.modisco_input_permute_augc_to_acgt,
    )
    if code != 0:
        sys.exit(code)

    print(f"[done] {out_dir / 'modisco_results.h5'}")
    if not args.skip_report:
        print(f"[done] 打开 {out_dir / 'report_simple'} 下 HTML（视 modisco 版本而定）")


if __name__ == "__main__":
    main()
