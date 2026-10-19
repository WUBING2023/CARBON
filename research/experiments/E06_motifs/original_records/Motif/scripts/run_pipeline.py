#!/usr/bin/env python3
"""
MARCP Motif & interpretability pipeline (all outputs under 可解释性分析/Motif/outputs).

流程：加载数据 → 前向预测与保存 →（可选）Top-K 筛选 → conv → 子集 IG →
TF-MoDISco（hyp 为原始 grad×input，不做每位置去均值）→ 可选 MEME/ISM。

Examples:
  python scripts/run_pipeline.py --demo
  python scripts/run_pipeline.py --checkpoint /path/model_params.pth --h5 /path/test.h5 --max_samples 256
  python scripts/run_pipeline.py --task MYTASK --h5_mode all_chr --checkpoint ... --select_top 10000
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import torch
import torch.nn.functional as F

from motif_analysis.paths import MARCP_ROOT, OUTPUTS, ensure_marcp_on_path
from motif_analysis.tfmodisco_cli import RELAXED_MOTIF_PRESET
from motif_analysis.wrapper import MARCPBindingIntensityWrapper, MARCPProfileWrapper


def _make_attr_wrapper(model, args):
    """归因标量：默认与训练解码一致（binding 强度），避免仅对 logits 归因。"""
    t = getattr(args, "attr_target", "binding_topk")
    if t == "binding_topk":
        return MARCPBindingIntensityWrapper(model, agg="topk", topk=args.topk)
    if t == "binding_sum":
        return MARCPBindingIntensityWrapper(model, agg="sum", topk=args.topk)
    if t == "binding_max":
        return MARCPBindingIntensityWrapper(model, agg="max", topk=args.topk)
    if t == "logits_topk":
        return MARCPProfileWrapper(model, agg="sum_topk", topk=args.topk)
    if t == "logits_sum":
        return MARCPProfileWrapper(model, agg="sum", topk=args.topk)
    if t == "logits_max":
        return MARCPProfileWrapper(model, agg="max", topk=args.topk)
    raise ValueError(f"Unknown attr_target={t}")


def _build_model(args, device: torch.device):
    ensure_marcp_on_path()
    from MARCP.models.marcp_model import MARCP

    eff_ch = 4 if args.ablation_no_conservative else args.input_channels
    ab_kw = {
        "ablation_no_se": args.ablation_no_se,
        "ablation_no_multiscale": args.ablation_no_multiscale,
        "ablation_no_gru": args.ablation_no_gru,
        "ablation_lstm_instead_of_gru": args.ablation_lstm_instead_of_gru,
    }
    m = MARCP(
        dropOut=0.0,
        control=False,
        in_channels=eff_ch,
        **ab_kw,
    ).to(device)
    if args.checkpoint:
        try:
            state = torch.load(args.checkpoint, map_location=device, weights_only=True)
        except TypeError:
            state = torch.load(args.checkpoint, map_location=device)
        try:
            m.load_state_dict(state, strict=True)
        except RuntimeError:
            m.load_state_dict(state, strict=False)
    else:
        print("[warn] No --checkpoint: using random init (demo only).")
    m.eval()
    return m


def cmd_conv(args, device, model, out_root: Path):
    from motif_analysis.conv_motifs import plot_conv_logos, save_conv_motif_report

    d = out_root / "conv_motifs"
    save_conv_motif_report(model, d)
    plot_conv_logos(model, d, max_filters_per_branch=args.max_filters_logo)
    print(f"[conv] wrote {d}")


def run_forward_profiles_batch(
    model,
    X: np.ndarray,
    device: torch.device,
    batch: int = 64,
) -> tuple[np.ndarray, np.ndarray]:
    """
    返回 softmax_profile (N,L), binding_intensity (N,L)，与 wrapper / 训练解码一致。
    """
    softmax_list = []
    bind_list = []
    for start in range(0, X.shape[0], batch):
        sl = slice(start, min(start + batch, X.shape[0]))
        xb = torch.from_numpy(X[sl]).float().to(device)
        with torch.no_grad():
            o = model(xb)
            lp = o["logits_profile"]
            ltc = o["log_total_count"]
            prob = F.softmax(lp, dim=-1)
            intensity = prob * torch.exp(ltc)
        softmax_list.append(prob.detach().cpu().numpy())
        bind_list.append(intensity.detach().cpu().numpy())
    return np.concatenate(softmax_list, axis=0), np.concatenate(bind_list, axis=0)


def cmd_predictions_and_rank(
    args,
    model,
    X: np.ndarray,
    Y: np.ndarray,
    chr_id: np.ndarray | None,
    source_idx: np.ndarray | None,
    out_root: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """写出 predictions/，返回用于下游的 X_sel, Y_sel, indices。"""
    from motif_analysis.sample_ranking import rank_samples, save_ranking

    device = torch.device(args.device)
    pred_dir = out_root / "predictions"
    pred_dir.mkdir(parents=True, exist_ok=True)

    softmax_profile, binding_intensity = run_forward_profiles_batch(
        model, X, device, batch=max(32, args.batch_size)
    )
    np.savez_compressed(
        pred_dir / "pred_profile.npz",
        softmax_profile=softmax_profile.astype(np.float32),
        binding_intensity=binding_intensity.astype(np.float32),
    )

    n = X.shape[0]
    select_top = args.select_top
    eff_k = n if (select_top <= 0 or select_top >= n) else select_top
    indices, scores, meta = rank_samples(
        softmax_profile=softmax_profile,
        binding_intensity=binding_intensity,
        Y=Y,
        rank_by=args.rank_by,
        select_top=eff_k,
    )
    if select_top <= 0 or select_top >= n:
        meta = {
            **meta,
            "requested_select_top": select_top,
            "note": "select_top<=0 或 >=N：使用全样本",
        }
    save_ranking(pred_dir, indices, meta, scores=scores, chr_id=chr_id, source_idx=source_idx)

    X_sel = X[indices]
    Y_sel = Y[indices]
    return X_sel, Y_sel, indices


def cmd_predictions_and_rank_chunked(
    args,
    model,
    paths: list[str],
    out_root: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """全染色体多文件：分块前向 + 流式 Top-K，再按全局下标加载子集并写 predictions/。"""
    from motif_analysis.chunked_ranking import stream_forward_topk
    from motif_analysis.data_io import load_h5_rows_by_global_indices
    from motif_analysis.sample_ranking import save_ranking

    device = torch.device(args.device)
    pred_dir = out_root / "predictions"
    pred_dir.mkdir(parents=True, exist_ok=True)

    max_total = args.max_samples
    indices, n_scanned, entries = stream_forward_topk(
        paths=paths,
        model=model,
        device=device,
        rank_by=args.rank_by,
        select_top=args.select_top,
        forward_batch=max(32, args.batch_size),
        chunk_rows=args.forward_chunk_rows,
        max_total_samples=max_total,
        ablation_no_conservative=args.ablation_no_conservative,
        run_forward_profiles_batch=run_forward_profiles_batch,
    )

    from motif_analysis.data_io import h5_concat_cumsum

    cum, _ = h5_concat_cumsum(paths)
    cum_n = int(cum[-1])
    if max_total is not None:
        cum_n = min(cum_n, max_total)

    if indices.size == 0:
        raise SystemExit("[chunked_rank] Top-K 为空；检查数据与 select_top")

    indices = indices[indices < cum_n]
    scores_topk = np.array([e[0] for e in entries], dtype=np.float64)[: indices.size]

    X_sel, Y_sel, _, chr_id, source_idx = load_h5_rows_by_global_indices(
        paths,
        indices,
        max_global=cum_n,
        ablation_no_conservative=args.ablation_no_conservative,
    )

    softmax_profile, binding_intensity = run_forward_profiles_batch(
        model, X_sel, device, batch=max(32, args.batch_size)
    )
    np.savez_compressed(
        pred_dir / "pred_profile.npz",
        softmax_profile=softmax_profile.astype(np.float32),
        binding_intensity=binding_intensity.astype(np.float32),
    )

    meta = {
        "rank_by": args.rank_by,
        "select_top": args.select_top,
        "n_total": cum_n,
        "n_scanned": n_scanned,
        "n_selected": int(indices.size),
        "chunked_rank": True,
        "forward_chunk_rows": args.forward_chunk_rows,
        "score_mean_selected": float(np.mean(scores_topk)) if indices.size else None,
    }
    save_ranking(pred_dir, indices, meta, scores=scores_topk, chr_id=chr_id, source_idx=source_idx)
    return X_sel, Y_sel, indices


def cmd_attrib_modisco(
    args,
    device: torch.device,
    model,
    X: np.ndarray,
    out_root: Path,
):
    from motif_analysis.attribution import grad_x_input_4ch, integrated_gradients_batch
    from motif_analysis.export_tfmodisco_h5 import save_tfmodisco_inputs
    from motif_analysis.seqlet_motifs import cluster_motifs, extract_seqlets, plot_cluster_logos, save_seqlet_report
    from motif_analysis.tfmodisco_cli import run_modisco_motifs_and_report

    wrap = _make_attr_wrapper(model, args)
    batch = min(args.batch_size, X.shape[0])
    hyps = []
    gxis = []
    ohes = []
    for start in range(0, X.shape[0], batch):
        sl = slice(start, min(start + batch, X.shape[0]))
        xb = torch.from_numpy(X[sl]).float().to(device)
        attr = integrated_gradients_batch(wrap, xb, n_steps=args.ig_steps)
        gxi = grad_x_input_4ch(attr, xb)
        hyps.append(attr[..., :4].detach().cpu().numpy())
        gxis.append(gxi.detach().cpu().numpy())
        ohes.append(xb[..., :4].detach().cpu().numpy())
    hyp = np.concatenate(hyps, axis=0)
    gxi_all = np.concatenate(gxis, axis=0)
    ohe = np.concatenate(ohes, axis=0)

    d_att = out_root / "attributions"
    d_att.mkdir(parents=True, exist_ok=True)
    gxi_abs = np.abs(gxi_all)
    hyp_for_h5 = gxi_all

    save_kw = dict(hyp=hyp, gxi=gxi_all, gxi_abs=gxi_abs, ohe=ohe)
    np.savez_compressed(d_att / "ig_hyp_4ch.npz", **save_kw)
    save_tfmodisco_inputs(d_att / "tfmodisco_inputs.h5", hyp_for_h5, ohe)

    backend = args.motif_backend
    if backend == "modisco" and not args.skip_tfmodisco:
        d_tf = out_root / "tfmodisco"
        run_modisco_motifs_and_report(
            gxi_nlc4=hyp_for_h5,
            ohe_nlc4=ohe,
            out_dir=d_tf,
            max_seqlets=args.max_seqlets,
            window=args.modisco_window,
            skip_report=args.skip_modisco_report,
            report_trim_threshold=args.modisco_report_trim_threshold,
            n_leiden=args.modisco_n_leiden,
            sliding_window_size=args.modisco_sliding_window_size,
            trim_size=args.modisco_trim_size,
            seqlet_flank_size=args.modisco_seqlet_flank_size,
            initial_flank_to_add=args.modisco_initial_flank_to_add,
            final_flank_to_add=args.modisco_final_flank_to_add,
            modisco_input_permute_augc_to_acgt=args.modisco_input_permute_augc_to_acgt,
        )
        print(f"[tfmodisco] {d_tf}")

    if backend == "kmeans":
        wo, _, _ = extract_seqlets(gxi_abs, ohe, window=args.seqlet_window)
        lab, ppms = cluster_motifs(wo, n_clusters=args.n_motif_clusters)
        d_seq = out_root / "seqlet_motifs"
        save_seqlet_report(
            d_seq,
            lab,
            ppms,
            meta={
                "window": args.seqlet_window,
                "n_clusters": args.n_motif_clusters,
                "note": "KMeans on one-hot seqlets from max |grad*input| window (legacy)",
            },
        )
        plot_cluster_logos(ppms, d_seq, prefix="motif")
        print(f"[seqlet kmeans] {d_seq}")

    if backend == "modisco" and getattr(args, "legacy_kmeans", False):
        wo, _, _ = extract_seqlets(gxi_abs, ohe, window=args.seqlet_window)
        lab, ppms = cluster_motifs(wo, n_clusters=args.n_motif_clusters)
        d_seq = out_root / "seqlet_motifs_legacy"
        save_seqlet_report(
            d_seq,
            lab,
            ppms,
            meta={"window": args.seqlet_window, "n_clusters": args.n_motif_clusters, "note": "optional legacy KMeans"},
        )
        plot_cluster_logos(ppms, d_seq, prefix="motif")
        print(f"[seqlet legacy] {d_seq}")

    print(f"[attrib] {d_att}")


def cmd_meme(args, X: np.ndarray, out_root: Path, pred_profile: np.ndarray | None):
    from motif_analysis.meme_export import peaks_to_fasta, write_streme_shell_script

    d = out_root / "meme_streme"
    d.mkdir(parents=True, exist_ok=True)
    ohe = X[..., :4]
    if pred_profile is None:
        pred_profile = X[..., 0].copy()
        print("[meme] no preds passed; using dummy peak proxy")
    n = peaks_to_fasta(ohe, pred_profile, d / "peaks_for_streme.fa", window=args.meme_window, max_seqs=args.max_streme_seqs)
    write_streme_shell_script(d / "run_streme.sh", d / "peaks_for_streme.fa", d / "streme_output")
    print(f"[meme] wrote {n} sequences to {d / 'peaks_for_streme.fa'}; stub {d / 'run_streme.sh'}")


def cmd_ism(args, device, model, X: np.ndarray, out_root: Path):
    from motif_analysis.attribution import grad_x_input_4ch, integrated_gradients_batch
    from motif_analysis.ism import compare_ig_ism, ism_delta_at_positions, save_ism_report

    wrap = _make_attr_wrapper(model, args)
    x = torch.from_numpy(X[:1]).float()
    attr = integrated_gradients_batch(wrap, x.to(device), n_steps=args.ig_steps)
    gxi = np.abs(grad_x_input_4ch(attr, x.to(device)).detach().cpu().numpy())
    L = x.size(1)
    rng = np.random.default_rng(0)
    positions = sorted(rng.choice(np.arange(10, L - 10), size=min(args.ism_positions, L - 20), replace=False).tolist())
    rows, _ = ism_delta_at_positions(model, x, positions, device)
    cmp = compare_ig_ism(gxi, positions, rows)
    d = out_root / "ism_check"
    save_ism_report(d, rows, cmp)
    print(f"[ism] {cmp} -> {d}")


def run_predict_profile_batch(model, X: np.ndarray, device, batch: int = 64) -> np.ndarray:
    preds = []
    for start in range(0, X.shape[0], batch):
        sl = slice(start, min(start + batch, X.shape[0]))
        xb = torch.from_numpy(X[sl]).float().to(device)
        with torch.no_grad():
            lp = model(xb)["logits_profile"]
            preds.append(torch.softmax(lp, dim=-1).detach().cpu().numpy())
    return np.concatenate(preds, axis=0)


def resolve_h5_paths(args) -> tuple[list[str], bool]:
    """返回 (paths, used_concat_meta)。"""
    if getattr(args, "demo", False):
        return [], False
    if args.h5_mode == "all_chr":
        if not args.task:
            raise SystemExit("--h5_mode all_chr 需要 --task")
        from motif_analysis.user_paths import list_chr5ch

        paths = [str(p) for p in list_chr5ch(args.task)]
        if not paths:
            raise FileNotFoundError(f"未找到 {args.task} 的 chr*_5ch.h5")
        return paths, True
    if not args.h5:
        raise SystemExit("--h5 required unless --demo 或 --h5_mode all_chr")
    parts = [p.strip() for p in args.h5.split(",") if p.strip()]
    return parts, len(parts) > 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--demo", action="store_true", help="Random model + synthetic data; no files required.")
    p.add_argument("--task", type=str, default=None, help="任务名；用于 h5_mode=all_chr 或默认输出目录名")
    p.add_argument("--checkpoint", type=str, default=None)
    p.add_argument("--h5", type=str, default=None, help="单个 5ch h5，或多个逗号分隔路径")
    p.add_argument(
        "--h5_mode",
        type=str,
        default="single",
        choices=["single", "all_chr"],
        help="single：用 --h5；all_chr：用 user_paths.list_chr5ch(--task) 合并全部染色体",
    )
    p.add_argument("--out", type=str, default=None, help="默认 outputs/<task> 或 outputs/run_default")
    p.add_argument("--max_samples", type=int, default=None, help="合并后最多使用样本数；默认不截断")
    p.add_argument(
        "--forward_chunk_rows",
        type=int,
        default=512,
        help="h5_mode=all_chr 且分块 Top-K 时，每次从 h5 读入的行数",
    )
    p.add_argument(
        "--no_chunked_rank",
        action="store_true",
        help="all_chr 时仍整表 concat 再排序（占内存；调试用）",
    )
    p.add_argument("--select_top", type=int, default=10000, help="归因前保留的 Top-K；<=0 表示全样本")
    p.add_argument(
        "--rank_by",
        type=str,
        default="pred_mean",
        choices=["pred_mean", "softmax_mean", "spearman_y"],
        help="筛选指标：binding 强度均值 / softmax 均值 / pred 与 Y 的 Spearman",
    )
    p.add_argument(
        "--motif_backend",
        type=str,
        default="modisco",
        choices=["modisco", "kmeans", "none"],
        help="motif 主路径：TF-MoDISco（默认）或旧版 KMeans seqlet，或仅写归因不写 motif",
    )
    p.add_argument("--skip_tfmodisco", action="store_true", help="仍写 npz/h5，但不调用 modisco CLI")
    p.add_argument("--skip_modisco_report", action="store_true", help="modisco motifs 后跳过 report-simple")
    p.add_argument(
        "--modisco-report-trim-threshold",
        type=float,
        default=0.28,
        help="已保留兼容：当前由 modisco report-simple 决定 trim，本参数不生效",
    )
    p.add_argument(
        "--modisco-relaxed-motifs",
        action="store_true",
        help=(
            "略宽于 MEME 默认的 motifs 参数（-z/-l/-t/-f/-g 与 max_seqlets），"
            "见 RELAXED_MOTIF_PRESET；显式传入的对应项优先"
        ),
    )
    p.add_argument(
        "--modisco-strict-motifs",
        action="store_true",
        help="已弃用，等同于 --modisco-relaxed-motifs",
    )
    p.add_argument(
        "--max_seqlets",
        type=int,
        default=None,
        help="modisco -n；未指定时默认 4000，若启用 --modisco-relaxed-motifs 则为 4500",
    )
    p.add_argument("--modisco_window", type=int, default=None, help="覆盖 modisco 窗口长度，默认自动")
    p.add_argument("--modisco-n-leiden", type=int, default=None, help="modisco motifs -l")
    p.add_argument("--modisco-sliding-window-size", type=int, default=None, help="modisco motifs -z（seqlet core）")
    p.add_argument("--modisco-trim-size", type=int, default=None, help="modisco motifs -t（trim_to_window_size）")
    p.add_argument("--modisco-seqlet-flank-size", type=int, default=None, help="modisco motifs -f")
    p.add_argument("--modisco-initial-flank-to-add", type=int, default=None, help="modisco motifs -g")
    p.add_argument("--modisco-final-flank-to-add", type=int, default=None, help="modisco motifs -j")
    p.add_argument(
        "--modisco-input-permute-augc-to-acgt",
        action="store_true",
        help=(
            "仅当输入最后一维按 A,U,G,C 排列时启用；标准 MARCP one-hot 为 A,C,G,U，默认勿开"
        ),
    )
    p.add_argument(
        "--ideepb-motif-settings",
        action="store_true",
        help=(
            "对齐 iDeepB 论文中可复现部分：select_top=1000、rank_by=spearman_y；"
            "MoDISco -t 10 -g 3 -j 3（trim_to_window_size / flank）。"
            "不含 FDR/null 打乱（MEME modisco CLI 无同名项）。"
        ),
    )
    p.add_argument("--legacy_kmeans", action="store_true", help="在 modisco 之外额外跑 KMeans seqlet 到 seqlet_motifs_legacy/")
    p.add_argument("--input_channels", type=int, default=5)
    p.add_argument("--ablation_no_conservative", action="store_true")
    p.add_argument("--ablation_no_se", action="store_true")
    p.add_argument("--ablation_no_multiscale", action="store_true")
    p.add_argument("--ablation_no_gru", action="store_true")
    p.add_argument("--ablation_lstm_instead_of_gru", action="store_true")
    p.add_argument("--batch_size", type=int, default=16)
    p.add_argument("--ig_steps", type=int, default=32)
    p.add_argument(
        "--attr_target",
        type=str,
        default="binding_topk",
        choices=[
            "binding_topk",
            "binding_sum",
            "binding_max",
            "logits_topk",
            "logits_sum",
            "logits_max",
        ],
        help="IG 归因标量：binding_* 与训练解码 prob*exp(log_count) 一致；logits_* 为旧版（易与训练不一致）",
    )
    p.add_argument(
        "--agg",
        type=str,
        default="sum_topk",
        choices=["sum", "sum_topk", "max"],
        help="已弃用：请用 --attr_target。保留仅为兼容旧脚本。",
    )
    p.add_argument("--topk", type=int, default=3)
    p.add_argument("--seqlet_window", type=int, default=11)
    p.add_argument("--n_motif_clusters", type=int, default=5)
    p.add_argument("--max_filters_logo", type=int, default=8)
    p.add_argument("--meme_window", type=int, default=15)
    p.add_argument("--max_streme_seqs", type=int, default=300)
    p.add_argument("--ism_positions", type=int, default=8)
    p.add_argument("--skip_meme", action="store_true")
    p.add_argument("--skip_ism", action="store_true")
    p.add_argument("--skip_conv", action="store_true", help="跳过第一层卷积核 Motif（conv_motifs/）")
    args = p.parse_args()

    if getattr(args, "ideepb_motif_settings", False):
        args.select_top = 1000
        args.rank_by = "spearman_y"
        if args.modisco_trim_size is None:
            args.modisco_trim_size = 10
        if args.modisco_initial_flank_to_add is None:
            args.modisco_initial_flank_to_add = 3
        if args.modisco_final_flank_to_add is None:
            args.modisco_final_flank_to_add = 3

    if getattr(args, "modisco_strict_motifs", False):
        print(
            "[warn] --modisco-strict-motifs 已弃用，请改用 --modisco-relaxed-motifs",
            file=sys.stderr,
        )
        args.modisco_relaxed_motifs = True

    if getattr(args, "modisco_relaxed_motifs", False):
        if args.modisco_n_leiden is None:
            args.modisco_n_leiden = RELAXED_MOTIF_PRESET["n_leiden"]
        if args.modisco_sliding_window_size is None:
            args.modisco_sliding_window_size = RELAXED_MOTIF_PRESET["sliding_window_size"]
        if args.modisco_trim_size is None:
            args.modisco_trim_size = RELAXED_MOTIF_PRESET["trim_size"]
        if args.modisco_seqlet_flank_size is None:
            args.modisco_seqlet_flank_size = RELAXED_MOTIF_PRESET["seqlet_flank_size"]
        if args.modisco_initial_flank_to_add is None:
            args.modisco_initial_flank_to_add = RELAXED_MOTIF_PRESET["initial_flank_to_add"]

    if args.max_seqlets is None:
        args.max_seqlets = (
            RELAXED_MOTIF_PRESET["max_seqlets"]
            if getattr(args, "modisco_relaxed_motifs", False)
            else 4000
        )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    args.device = str(device)

    if args.out:
        out_root = Path(args.out)
    elif args.task:
        out_root = OUTPUTS / args.task
    else:
        out_root = OUTPUTS / "run_default"
    out_root.mkdir(parents=True, exist_ok=True)

    run_cfg = {
        **{k: v for k, v in vars(args).items() if k != "device"},
        "marcp_root": str(MARCP_ROOT),
        "device": str(device),
        "motif_backend": args.motif_backend,
    }
    with open(out_root / "run_config.json", "w", encoding="utf-8") as f:
        json.dump(run_cfg, f, indent=2, ensure_ascii=False)

    paths: list[str] = []
    multi = False
    use_chunked = False

    if args.demo:
        from motif_analysis.data_io import synthetic_batch

        X = synthetic_batch(batch=min(256, 48), length=101, channels=5 if not args.ablation_no_conservative else 4)
        Y = np.random.rand(X.shape[0], X.shape[1]).astype(np.float32)
        args.checkpoint = None
    else:
        paths, multi = resolve_h5_paths(args)
        use_chunked = (
            not args.no_chunked_rank
            and args.h5_mode == "all_chr"
            and args.select_top > 0
            and len(paths) >= 1
        )
        if args.h5_mode == "all_chr" and args.select_top <= 0:
            print(
                "[warn] all_chr 且 select_top<=0 将尝试整表加载，可能 OOM；建议设 select_top>0 以启用分块 Top-K。",
                flush=True,
            )
        if not use_chunked:
            from motif_analysis.data_io import load_h5_5ch, load_h5_5ch_concat

            max_s = args.max_samples
            if multi or len(paths) > 1:
                X, Y, Y_ctl, chr_id, source_idx = load_h5_5ch_concat(paths, max_samples=max_s)
            else:
                X, Y, Y_ctl = load_h5_5ch(paths[0], max_samples=max_s)
                chr_id = np.full(X.shape[0], -1, dtype=np.int64)
                source_idx = np.zeros(X.shape[0], dtype=np.int64)
            if args.ablation_no_conservative:
                X = X[..., :4]

    model = _build_model(args, device)

    if args.demo:
        X_sel, Y_sel, _indices = cmd_predictions_and_rank(
            args, model, X, Y, np.full(X.shape[0], -1, dtype=np.int64), np.zeros(X.shape[0], dtype=np.int64), out_root
        )
    elif use_chunked:
        X_sel, Y_sel, _indices = cmd_predictions_and_rank_chunked(args, model, paths, out_root)
    else:
        X_sel, Y_sel, _indices = cmd_predictions_and_rank(args, model, X, Y, chr_id, source_idx, out_root)

    if not args.skip_conv:
        cmd_conv(args, device, model, out_root)

    cmd_attrib_modisco(args, device, model, X_sel, out_root)

    pred = run_predict_profile_batch(model, X_sel, device)
    if not args.skip_meme:
        cmd_meme(args, X_sel, out_root, pred)
    if not args.skip_ism:
        cmd_ism(args, device, model, X_sel, out_root)

    print(f"Done. Outputs under {out_root}")


if __name__ == "__main__":
    main()
