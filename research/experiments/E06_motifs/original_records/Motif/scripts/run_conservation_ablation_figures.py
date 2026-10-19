#!/usr/bin/env python3
"""
Motif 区域保守性消融图谱对照（非碱基 ISM）。

核心思路：
- 同一条序列 one-hot 不变
- 仅修改第 5 维保守性通道（zero / shuffle）
- 使用同一 5ch 模型前向，比较 binding intensity 曲线
- 可选叠加 4ch no_conservative 模型曲线
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from motif_analysis.motif_scan import best_hit_window
from motif_analysis.query_motifs import load_all_query_motifs
from motif_analysis.user_paths import MARCP_WEIGHTS_ROOT, list_chr5ch
from motif_analysis.validation_ism import ISM_COLOR_MOTIF_REGION, ISM_COLOR_MUTATED, ISM_COLOR_REFERENCE


def _load_run_config(task_dir: Path) -> dict:
    p = task_dir / "run_config.json"
    if not p.is_file():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def _task_name_for_weights(task_dir: Path, cfg: dict) -> str:
    if cfg.get("task"):
        return str(cfg["task"])
    return task_dir.name


def _build_marcp(ckpt: Path, *, in_channels: int, cfg_template: dict, device: torch.device) -> torch.nn.Module:
    from motif_analysis.paths import ensure_marcp_on_path

    ensure_marcp_on_path()
    from MARCP.models.marcp_model import MARCP

    m = MARCP(
        dropOut=0.0,
        control=False,
        in_channels=in_channels,
        ablation_no_se=bool(cfg_template.get("ablation_no_se", False)),
        ablation_no_multiscale=bool(cfg_template.get("ablation_no_multiscale", False)),
        ablation_no_gru=bool(cfg_template.get("ablation_no_gru", False)),
        ablation_lstm_instead_of_gru=bool(cfg_template.get("ablation_lstm_instead_of_gru", False)),
    ).to(device)
    try:
        state = torch.load(str(ckpt), map_location=device, weights_only=True)
    except TypeError:
        state = torch.load(str(ckpt), map_location=device)
    m.load_state_dict(state, strict=False)
    m.eval()
    return m


def _binding_profile(model: torch.nn.Module, x_nlc: np.ndarray, device: torch.device) -> np.ndarray:
    with torch.no_grad():
        xb = torch.from_numpy(x_nlc).float().to(device)
        out = model(xb)
        lp = out["logits_profile"]
        ltc = out["log_total_count"]
        prob = F.softmax(lp, dim=-1)
        intensity = (prob * torch.exp(ltc)).squeeze(0).detach().cpu().numpy()
    return np.asarray(intensity, dtype=np.float64)


def _seq_code_from_ohe(ohe_l4: np.ndarray) -> bytes:
    idx = np.argmax(ohe_l4, axis=-1).astype(np.uint8)
    return idx.tobytes()


def _resolve_x5_candidates(
    task_name: str,
    ohe_sel_nl4: np.ndarray,
    global_indices: np.ndarray | None,
) -> tuple[list[np.ndarray], list[dict]]:
    """返回与 ohe_sel 对齐的 x5 列表；无法匹配的位置返回 None。"""
    h5_paths = list_chr5ch(task_name)
    n_sel = ohe_sel_nl4.shape[0]
    xs: list[np.ndarray | None] = [None] * n_sel
    meta: list[dict] = [{"method": "unresolved"} for _ in range(n_sel)]

    # 方案 A：有全局索引时直接按顺序取
    if global_indices is not None and global_indices.size == n_sel:
        offsets = []
        total = 0
        for p in h5_paths:
            with h5py.File(p, "r") as hf:
                n = int(hf["X"].shape[0])
            offsets.append((p, total, total + n))
            total += n
        for i, gidx in enumerate(global_indices.tolist()):
            gi = int(gidx)
            for p, lo, hi in offsets:
                if lo <= gi < hi:
                    loc = gi - lo
                    with h5py.File(p, "r") as hf:
                        row = np.asarray(hf["X"][loc], dtype=np.float32)
                    xs[i] = row[None, ...]
                    meta[i] = {"method": "global_index", "global_index": gi, "h5": str(p), "local_index": int(loc)}
                    break

    # 方案 B：按 one-hot 码匹配（用于缺少 sample_indices.npy 的历史任务）
    unresolved = [i for i, x in enumerate(xs) if x is None]
    if unresolved:
        wanted = {_seq_code_from_ohe(ohe_sel_nl4[i]): i for i in unresolved}
        for p in h5_paths:
            if not wanted:
                break
            with h5py.File(p, "r") as hf:
                X = np.asarray(hf["X"], dtype=np.float32)
            idx4 = np.argmax(X[..., :4], axis=-1).astype(np.uint8)
            for j in range(idx4.shape[0]):
                code = idx4[j].tobytes()
                i = wanted.get(code)
                if i is None:
                    continue
                xs[i] = X[j : j + 1]
                meta[i] = {"method": "ohe_match", "h5": str(p), "local_index": int(j)}
                del wanted[code]
                if not wanted:
                    break

    x_ok: list[np.ndarray] = []
    m_ok: list[dict] = []
    for i, x in enumerate(xs):
        if x is None:
            continue
        x_ok.append(x)
        m = dict(meta[i])
        m["selected_order"] = i
        m_ok.append(m)
    return x_ok, m_ok


def _calc_metrics(ref: np.ndarray, abl: np.ndarray, st: int, en: int) -> dict:
    sw = slice(st, en)
    ref_w = ref[sw]
    abl_w = abl[sw]
    return {
        "motif_start": int(st),
        "motif_end": int(en),
        "sum_ref_full": float(np.sum(ref)),
        "sum_abl_full": float(np.sum(abl)),
        "delta_sum_full": float(np.sum(ref) - np.sum(abl)),
        "peak_ref_full": float(np.max(ref)),
        "peak_abl_full": float(np.max(abl)),
        "delta_peak_full": float(np.max(ref) - np.max(abl)),
        "l2_full": float(np.linalg.norm(ref - abl)),
        "sum_ref_motif": float(np.sum(ref_w)),
        "sum_abl_motif": float(np.sum(abl_w)),
        "delta_sum_motif": float(np.sum(ref_w) - np.sum(abl_w)),
        "peak_ref_motif": float(np.max(ref_w)),
        "peak_abl_motif": float(np.max(abl_w)),
        "delta_peak_motif": float(np.max(ref_w) - np.max(abl_w)),
        "l2_motif": float(np.linalg.norm(ref_w - abl_w)),
    }


def _plot_main(
    out_png: Path,
    ref: np.ndarray,
    abl: np.ndarray,
    phy: np.ndarray,
    st: int,
    en: int,
    label_abl: str,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    L = ref.shape[0]
    fs = float(plt.rcParams.get("font.size", 10.0)) * 1.5
    fig, ax = plt.subplots(figsize=(12.0, 3.8))
    x = np.arange(L)
    ax.plot(x, ref, label="with conservation", color=ISM_COLOR_REFERENCE, linewidth=1.5)
    ax.plot(x, abl, label=label_abl, color=ISM_COLOR_MUTATED, linewidth=1.4)
    ax.axvspan(st, en - 1, facecolor=ISM_COLOR_MOTIF_REGION, edgecolor="none", alpha=0.2, label="motif region")
    ax2 = ax.twinx()
    ax2.plot(x, phy, color="#888888", linewidth=1.0, alpha=0.5, label="PhyloP")
    ax.set_xlabel("position", fontsize=fs)
    ax.set_ylabel("binding intensity", fontsize=fs)
    ax2.set_ylabel("conservation score", fontsize=fs)
    ax.tick_params(axis="both", labelsize=fs * 0.85)
    ax2.tick_params(axis="y", labelsize=fs * 0.75)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=fs * 0.8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_png.with_suffix(".pdf"), format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _plot_diff(out_png: Path, ref: np.ndarray, abl: np.ndarray, st: int, en: int, title: str) -> None:
    L = ref.shape[0]
    d = ref - abl
    fs = float(plt.rcParams.get("font.size", 10.0)) * 1.4
    fig, ax = plt.subplots(figsize=(12.0, 3.2))
    ax.plot(np.arange(L), d, color="#3C5488", linewidth=1.4, label="ref - ablation")
    ax.axhline(0.0, color="#999999", linewidth=0.8, linestyle="--")
    ax.axvspan(st, en - 1, facecolor=ISM_COLOR_MOTIF_REGION, edgecolor="none", alpha=0.2, label="motif region")
    ax.set_xlabel("position", fontsize=fs)
    ax.set_ylabel("delta intensity", fontsize=fs)
    ax.set_title(title, fontsize=fs)
    ax.tick_params(axis="both", labelsize=fs * 0.85)
    ax.legend(fontsize=fs * 0.75)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_png.with_suffix(".pdf"), format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _plot_multiline(
    out_png: Path,
    curves: list[tuple[str, np.ndarray, str]],
    st: int,
    en: int,
    title: str,
) -> None:
    L = curves[0][1].shape[0]
    fs = float(plt.rcParams.get("font.size", 10.0)) * 1.35
    fig, ax = plt.subplots(figsize=(12.0, 3.5))
    x = np.arange(L)
    for label, y, color in curves:
        ax.plot(x, y, label=label, linewidth=1.3, color=color)
    ax.axvspan(st, en - 1, facecolor=ISM_COLOR_MOTIF_REGION, edgecolor="none", alpha=0.2, label="motif region")
    ax.set_xlabel("position", fontsize=fs)
    ax.set_ylabel("binding intensity", fontsize=fs)
    ax.set_title(title, fontsize=fs)
    ax.tick_params(axis="both", labelsize=fs * 0.82)
    ax.legend(fontsize=fs * 0.72, ncol=2)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_png.with_suffix(".pdf"), format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _plot_summary(out_png: Path, metric_rows: list[dict], key: str, title: str) -> None:
    labels = []
    vals = []
    for r in metric_rows:
        labels.append(r["sample_id"])
        vals.append(float(r[key]))
    fs = float(plt.rcParams.get("font.size", 10.0)) * 1.25
    fig, ax = plt.subplots(figsize=(max(8.0, 0.42 * len(vals) + 4), 3.8))
    ax.bar(np.arange(len(vals)), vals, color="#4DBBD5", edgecolor="#3C5488", linewidth=0.5)
    ax.axhline(0.0, color="#999999", linewidth=0.8, linestyle="--")
    ax.set_xticks(np.arange(len(vals)))
    ax.set_xticklabels(labels, rotation=65, ha="right", fontsize=fs * 0.6)
    ax.set_ylabel(key, fontsize=fs)
    ax.set_title(title, fontsize=fs)
    ax.tick_params(axis="y", labelsize=fs * 0.8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_png.with_suffix(".pdf"), format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _resolve_global_indices(pred_dir: Path, pick_idx: np.ndarray) -> np.ndarray | None:
    p = pred_dir / "sample_indices.npy"
    if not p.is_file():
        return None
    all_idx = np.load(p)
    if all_idx.ndim != 1:
        return None
    if int(np.max(pick_idx)) >= all_idx.shape[0]:
        return None
    return np.asarray(all_idx[pick_idx], dtype=np.int64)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-dir", type=Path, required=True)
    ap.add_argument("--ckpt-5ch", type=Path, default=None)
    ap.add_argument("--ckpt-4ch", type=Path, default=None)
    ap.add_argument("--n-samples", type=int, default=12)
    ap.add_argument("--ablation", type=str, default="both", choices=["zero", "shuffle", "both"])
    ap.add_argument("--query-index", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--skip-4ch", action="store_true")
    args = ap.parse_args()

    task_dir = args.task_dir.resolve()
    cfg = _load_run_config(task_dir)
    task_name = _task_name_for_weights(task_dir, cfg)
    out_dir = task_dir / "conservation_ablation_viz"
    out_dir.mkdir(parents=True, exist_ok=True)

    att_path = task_dir / "attributions" / "ig_hyp_4ch.npz"
    if not att_path.is_file():
        raise FileNotFoundError(f"缺少 attributions 文件: {att_path}")
    pred_path = task_dir / "predictions" / "pred_profile.npz"
    if not pred_path.is_file():
        raise FileNotFoundError(f"缺少 predictions 文件: {pred_path}")

    z = np.load(att_path)
    ohe = np.asarray(z["ohe"], dtype=np.float32)
    pz = np.load(pred_path)
    bind = np.asarray(
        pz["binding_intensity"] if "binding_intensity" in pz.files else pz["softmax_profile"],
        dtype=np.float32,
    )

    scores = np.max(bind, axis=1)
    order = np.argsort(-scores)
    pick = order[: min(args.n_samples, order.size)]
    ohe_sel = ohe[pick]
    bind_sel = bind[pick]

    pred_dir = task_dir / "predictions"
    global_indices = _resolve_global_indices(pred_dir, pick)
    x5_list, x5_meta = _resolve_x5_candidates(task_name, ohe_sel, global_indices)
    if not x5_list:
        raise RuntimeError("未能从 5ch h5 匹配到样本，请确认 h5data 与 task 一致。")

    # 与 x5 对齐的 ohe / pred 子集
    kept_orders = np.array([m["selected_order"] for m in x5_meta], dtype=np.int64)
    ohe_keep = ohe_sel[kept_orders]
    bind_keep = bind_sel[kept_orders]

    queries = load_all_query_motifs(task_dir)
    if not queries:
        raise RuntimeError(f"未找到 query motifs: {task_dir}")
    qid, ppm = queries[min(args.query_index, len(queries) - 1)]

    ckpt5 = args.ckpt_5ch or Path(cfg.get("checkpoint", ""))
    if not ckpt5 or not Path(ckpt5).is_file():
        ckpt5 = MARCP_WEIGHTS_ROOT / task_name / "baseline" / "model_params.pth"
    ckpt5 = Path(ckpt5)
    if not ckpt5.is_file():
        raise FileNotFoundError(f"5ch checkpoint 不存在: {ckpt5}")

    ckpt4 = args.ckpt_4ch or (MARCP_WEIGHTS_ROOT / task_name / "no_conservative" / "model_params.pth")
    ckpt4 = Path(ckpt4)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model5 = _build_marcp(ckpt5, in_channels=5, cfg_template=cfg, device=device)
    model4 = None
    if (not args.skip_4ch) and ckpt4.is_file():
        model4 = _build_marcp(ckpt4, in_channels=4, cfg_template=cfg, device=device)

    rng = np.random.default_rng(args.seed)
    rows = []
    for k, (x5, meta_row) in enumerate(zip(x5_list, x5_meta)):
        ohe_one = ohe_keep[k]
        _, st = best_hit_window(ohe_one, ppm)
        en = st + ppm.shape[0]
        sample_id = f"s{k:03d}"
        fig_dir = out_dir / "figures" / sample_id
        fig_dir.mkdir(parents=True, exist_ok=True)

        ref = _binding_profile(model5, x5, device)
        x5_zero = x5.copy()
        x5_zero[..., 4] = 0.0
        zero = _binding_profile(model5, x5_zero, device)

        x5_shuffle = x5.copy()
        perm = rng.permutation(x5.shape[1])
        x5_shuffle[0, :, 4] = x5_shuffle[0, perm, 4]
        shuf = _binding_profile(model5, x5_shuffle, device)

        prof4 = None
        if model4 is not None:
            prof4 = _binding_profile(model4, x5[..., :4], device)

        if args.ablation == "zero":
            abl = zero
            abl_name = "without conservation (zero)"
        elif args.ablation == "shuffle":
            abl = shuf
            abl_name = "without conservation (shuffle)"
        else:
            abl = zero
            abl_name = "without conservation (zero)"

        m_main = _calc_metrics(ref, abl, st, en)
        m_zero = _calc_metrics(ref, zero, st, en)
        m_shuf = _calc_metrics(ref, shuf, st, en)
        m4 = _calc_metrics(ref, prof4, st, en) if prof4 is not None else None

        _plot_main(fig_dir / "figA_profile_compare.png", ref, abl, x5[0, :, 4], st, en, abl_name)
        _plot_diff(
            fig_dir / "figB_diff_curve.png",
            ref,
            abl,
            st,
            en,
            title="reference - ablation",
        )
        curves = [
            ("with conservation", ref, ISM_COLOR_REFERENCE),
            ("zero-PhyloP", zero, "#4DBBD5"),
            ("shuffle-PhyloP", shuf, "#00A087"),
        ]
        if prof4 is not None:
            curves.append(("4ch no_conservative", prof4, "#F39B7F"))
        _plot_multiline(
            fig_dir / "figC_multi_curve.png",
            curves,
            st,
            en,
            title="Conservation ablation profile curves",
        )

        row = {
            "sample_id": sample_id,
            "query_motif_id": qid,
            "picked_rank": int(k),
            "source_meta": meta_row,
            "motif_start": int(st),
            "motif_end": int(en),
            "main_ablation": args.ablation,
            "metrics_main": m_main,
            "metrics_zero": m_zero,
            "metrics_shuffle": m_shuf,
            "metrics_4ch": m4,
            "figures": {
                "figA": str(fig_dir / "figA_profile_compare.png"),
                "figB": str(fig_dir / "figB_diff_curve.png"),
                "figC": str(fig_dir / "figC_multi_curve.png"),
            },
            "pred_peak_selected_subset": float(np.max(bind_keep[k])),
        }
        rows.append(row)

    # 跨样本汇总图
    rows_zero = [{"sample_id": r["sample_id"], "delta_peak_motif": r["metrics_zero"]["delta_peak_motif"]} for r in rows]
    rows_sum = [{"sample_id": r["sample_id"], "delta_sum_motif": r["metrics_zero"]["delta_sum_motif"]} for r in rows]
    _plot_summary(
        out_dir / "figD_delta_peak_motif_zero.png",
        rows_zero,
        key="delta_peak_motif",
        title="Motif-window delta peak (with conservation - zero-PhyloP)",
    )
    _plot_summary(
        out_dir / "figE_delta_sum_motif_zero.png",
        rows_sum,
        key="delta_sum_motif",
        title="Motif-window delta sum (with conservation - zero-PhyloP)",
    )

    summary = {
        "task_dir": str(task_dir),
        "task_name": task_name,
        "device": str(device),
        "query_motif_id": qid,
        "n_requested": int(args.n_samples),
        "n_matched_x5": int(len(rows)),
        "ckpt_5ch": str(ckpt5),
        "ckpt_4ch": str(ckpt4) if model4 is not None else None,
        "ablation_mode": args.ablation,
        "samples": rows,
        "aggregate": {
            "delta_peak_motif_zero_mean": float(np.mean([r["metrics_zero"]["delta_peak_motif"] for r in rows])),
            "delta_sum_motif_zero_mean": float(np.mean([r["metrics_zero"]["delta_sum_motif"] for r in rows])),
            "delta_peak_motif_shuffle_mean": float(np.mean([r["metrics_shuffle"]["delta_peak_motif"] for r in rows])),
            "delta_sum_motif_shuffle_mean": float(np.mean([r["metrics_shuffle"]["delta_sum_motif"] for r in rows])),
        },
    }
    out_json = out_dir / "conservation_ablation_report.json"
    out_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[ok] {out_json}")
    print(f"[ok] figures in {out_dir / 'figures'}")


if __name__ == "__main__":
    main()
