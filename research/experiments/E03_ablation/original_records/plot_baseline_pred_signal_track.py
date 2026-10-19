#!/usr/bin/env python3
"""
同一测试窗口上绘制多轨「强度」图谱（默认四轨）：eCLIP signal、MARCP、iDeepB、RBPNet；
`--figure2-export` 时为六轨并在 eCLIP 下增加两条倒置展示轨；
`--figure1b` 时仅两轨（MARCP + eCLIP signal），画布更扁（默认约 3.3:1 宽:高），配色对齐 Figure2，且不标注 Pearson r。

- 样本与 Signal：与全量消融一致（5ch hdf + peak_test ∩ chr1/8/15）；默认按 **max(Y)** 降序取若干窗口，各输出一张同风格四轨图。
- iDeepB / RBPNet 权重：`比较实验/output/模型参数/`（RBPNet 优先 `model.h5`，与正式对比流程一致）。
- MARCP：`比较实验/output/模型参数/MARCP` 若存在则用，否则默认 `消融实验/全量消融运行/模型权重/.../baseline/model_params.pth`
  （当前对比实验 output 下通常未单独存 MARCP 权重）。

默认画布宽高比 2:1（宽:高）；`--figure1b` 为约 3.3:1。六轨时 eCLIP/rep 填充为 `图/参考配色` NPG 7/6/4，模型轨为 NPG 1–3；填充 alpha=1、savefig 白底便于取色器对齐。
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt

plt.rcParams.setdefault("pdf.fonttype", 42)  # PDF 内文字为 TrueType，可选中
import numpy as np
import pandas as pd
import torch

BASE = Path(__file__).resolve().parents[1]
MARCP_ROOT = BASE / "MARCP"
IDEEP_ROOT = BASE / "比较实验/iDeepB"
RBP_ROOT = BASE / "比较实验/RBPNet"
COMPARE_OUTPUT = BASE / "比较实验/output"
WEIGHTS_OUT = COMPARE_OUTPUT / "模型参数"

# 图/参考配色 NPG 1–3：三模型主色（与 图/figure2、compare_three_models、circos 统一）
# 轨标签展示名为 CARBON（与论文模型名一致；权重路径仍见 resolve_marcp_weights）
FIG01_VIOLIN_COLOR = {
    "CARBON": "#E64B35",
    "iDeepB": "#4DBBD5",
    "RBPNet": "#00A087",
}
FIGURE2_MODEL_COLOR = dict(FIG01_VIOLIN_COLOR)
# 图/参考配色：Signal 与 rep 轨仅使用 NPG 4/6/7（与模型 1–3、散点细胞系可区分）
_ECLIP_MAIN_FILL = "#91D1C2"
_ENCODE_REP_FILL = "#8491B4"
_ENCODE_REP2_FILL = "#3C5488"
# figure2 模型轨右侧：名称与 "r = …" 之间严格为两个半角空格（U+0020）
FIGURE2_TRACK_R_GAP = "  "
# Figure 1b：扁长画布（figsize 为 matplotlib 的 width×height 英寸）；约 3.3:1 宽:高
FIGURE1B_FIGSIZE_INCH = (10.0, 3.0)

H5DATA = BASE / "数据/h5data"
BED_ROOT = BASE / "数据/new_beddata"
TEST_CHRS = ("1", "8", "15")

DEFAULT_TASK = "FMR1_K562_ENCSR331VNX"


def _read_test_keys(bed_path: Path) -> set[tuple[str, int, int, str]]:
    keys: set[tuple[str, int, int, str]] = set()
    with open(bed_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom = str(parts[0]).replace("chr", "")
            start, end = int(parts[1]), int(parts[2])
            strand = parts[5]
            keys.add((chrom, start, end, strand))
    return keys


def _norm_chrom(c: str) -> str:
    return str(c).replace("chr", "")


def _figure_display_task_name(task: str) -> str:
    """子图标题用名：去掉 `_ENC…`（ENCODE accession），与 Figure2 散点注释一致。"""
    s = str(task)
    if "_ENC" in s:
        return s.split("_ENC", 1)[0]
    return s


def collect_test_rows_sorted_by_peak(
    task: str,
    test_keys: set[tuple[str, int, int, str]],
    chrs: tuple[str, ...],
) -> list[tuple[float, pd.Series, np.ndarray, np.ndarray, np.ndarray]]:
    """peak_test ∩ 各 chr 的窗口；按 max(Y) 降序。每项 (peak_max, meta_row, X_5ch, Y, Y_ctl)。"""
    rows: list[tuple[float, pd.Series, np.ndarray, np.ndarray, np.ndarray]] = []

    for chr_n in chrs:
        c5 = H5DATA / task / "_5ch" / f"{task}_chr{chr_n}_5ch.h5"
        meta_path = H5DATA / task / f"{task}_chr{chr_n}.h5"
        if not c5.is_file() or not meta_path.is_file():
            continue
        df = pd.read_hdf(meta_path, key="df")
        df["_ck"] = list(
            zip(
                df["chr"].map(_norm_chrom),
                df["chromStart"].astype(int),
                df["chromEnd"].astype(int),
                df["strand"].astype(str),
            )
        )
        sub_mask = df["_ck"].apply(lambda t: t in test_keys)
        idxs = np.flatnonzero(sub_mask.to_numpy())
        if idxs.size == 0:
            continue
        with h5py.File(c5, "r") as hf:
            Yall = np.asarray(hf["Y"], dtype=np.float32)
            Yctl_all = np.asarray(hf["Y_ctl"], dtype=np.float32)
            Xall = np.asarray(hf["X"], dtype=np.float32)
        for i in idxs:
            ii = int(i)
            y = Yall[ii]
            yctl = Yctl_all[ii]
            peak = float(np.max(y))
            rows.append((peak, df.iloc[ii], Xall[ii].copy(), y.copy(), yctl.copy()))

    rows.sort(key=lambda t: -t[0])
    if not rows:
        raise RuntimeError(
            "未找到测试集交集窗口：请检查 test BED、5ch 与 h5data；TEST_CHRS 需与训练一致 (1,8,15)。"
        )
    return rows


def pick_spread_by_rank(
    rows: list[tuple[float, pd.Series, np.ndarray, np.ndarray, np.ndarray]],
    k: int,
) -> list[tuple[float, pd.Series, np.ndarray, np.ndarray, np.ndarray]]:
    """在已按 peak 降序的列表中，按排名近似均匀抽样（强峰、中等、较弱…）。"""
    n = len(rows)
    if k >= n:
        return list(rows)
    if k <= 1:
        return [rows[0]]
    pos = np.linspace(0, n - 1, k)
    idxs: list[int] = sorted({int(min(n - 1, round(float(p)))) for p in pos})
    t = 0
    while len(idxs) < k and t < n:
        if t not in idxs:
            idxs.append(t)
        t += 1
    idxs = sorted(idxs)[:k]
    return [rows[i] for i in idxs]


def resolve_marcp_weights(task: str, override: Path | None) -> Path:
    if override is not None:
        return override
    p = WEIGHTS_OUT / "MARCP" / task / "model_params.pth"
    if p.is_file():
        return p
    p = (
        BASE
        / "消融实验/全量消融运行/模型权重/MARCP"
        / task
        / "baseline/model_params.pth"
    )
    if p.is_file():
        return p
    raise FileNotFoundError(
        "未找到 MARCP 权重：可在 比较实验/output/模型参数/MARCP 放置 model_params.pth，"
        f"或确认消融 baseline 存在：{p}"
    )


def resolve_ideep_weights(task: str, override: Path | None) -> Path:
    if override is not None:
        return override
    p = WEIGHTS_OUT / "iDeepB" / task / "model_params.pth"
    if not p.is_file():
        raise FileNotFoundError(p)
    return p


def resolve_rbpnet_weights(task: str, override: Path | None) -> Path:
    if override is not None:
        return override
    d = WEIGHTS_OUT / "RBPNet" / task
    h5 = d / "model.h5"
    pth = d / "model_params.pth"
    if h5.is_file():
        return h5
    if pth.is_file():
        return pth
    raise FileNotFoundError(f"未找到 RBPNet 模型：{h5} 或 {pth}")


def load_marcp_model(weights: Path, device: torch.device):
    sys.path.insert(0, str(MARCP_ROOT))
    from MARCP.models.marcp_model import MARCP  # noqa: E402

    model = MARCP(dropOut=0.2, control=False, in_channels=5).to(device)
    try:
        state = torch.load(weights, map_location=device, weights_only=True)
    except TypeError:
        state = torch.load(weights, map_location=device)
    model.load_state_dict(state)
    model.eval()
    return model


def forward_marcp(model, X_5ch: np.ndarray, device: torch.device) -> np.ndarray:
    t = torch.from_numpy(X_5ch[np.newaxis, ...]).float().to(device)
    with torch.no_grad():
        out = model(t)
        logits = out["logits_profile"]
        log_tc = out["log_total_count"]
        prob = torch.softmax(logits, dim=-1)
        pred = prob * torch.exp(log_tc)
    return pred.squeeze(0).detach().cpu().numpy()


def load_ideepb_model(weights: Path, device: torch.device):
    sys.path.insert(0, str(IDEEP_ROOT))
    from iDeepB.models.iDeepBModel import iDeepB  # noqa: E402

    model = iDeepB(dropOut=0.2, control=False).to(device)
    try:
        state = torch.load(weights, map_location=device, weights_only=True)
    except TypeError:
        state = torch.load(weights, map_location=device)
    model.load_state_dict(state)
    model.eval()
    return model


def forward_ideepb(model, onehot_L4: np.ndarray, device: torch.device) -> np.ndarray:
    t = torch.from_numpy(onehot_L4[np.newaxis, ...]).float().to(device)
    with torch.no_grad():
        pred = model(t)
    return pred.squeeze(0).detach().cpu().numpy()


def seq_to_ideep_onehot(seq: str) -> np.ndarray:
    sys.path.insert(0, str(IDEEP_ROOT))
    from iDeepB.utils.functions import onehot_encode  # noqa: E402

    seq_u = str(seq).upper().replace("T", "U")
    arr = onehot_encode([seq_u], list("AUGC"), 4)
    return np.asarray(arr[0], dtype=np.float32)


def predict_rbpnet_profile(seq_dna: str, weights: Path, rbp_python: Path) -> np.ndarray:
    """RBPNet 依赖 TensorFlow，使用 比较实验/RBPNet/venv/bin/python3 子进程推理。"""
    helper = Path(__file__).resolve().parent / "rbpnet_predict_profile_once.py"
    if not helper.is_file():
        raise FileNotFoundError(helper)
    if not rbp_python.is_file():
        raise FileNotFoundError(
            f"未找到 RBPNet 解释器: {rbp_python}（请使用 比较实验/RBPNet/venv 或 --rbp-python）"
        )
    env = os.environ.copy()
    env.setdefault("CUDA_VISIBLE_DEVICES", "")
    proc = subprocess.run(
        [str(rbp_python), str(helper), seq_dna, str(weights)],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(f"RBPNet 子进程失败 (code {proc.returncode}): {err[:2000]}")
    line = proc.stdout.strip().splitlines()[-1] if proc.stdout else ""
    return np.asarray(json.loads(line), dtype=np.float32)


def _pearson_r_vs_signal(pred: np.ndarray, y_signal: np.ndarray) -> float:
    """与 eCLIP 真值剖面 Y 的 Pearson r，作预测轨「准确率」展示。"""
    a = np.asarray(pred, dtype=np.float64).ravel()
    b = np.asarray(y_signal, dtype=np.float64).ravel()
    if a.size != b.size or a.size < 2:
        return float("nan")
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    r = float(np.corrcoef(a, b)[0, 1])
    return r if np.isfinite(r) else float("nan")


def _track_name_for_color(label: str) -> str:
    """配色键：去掉换行；若在 FIGURE2_TRACK_R_GAP 后接 r = …，则只保留模型名。"""
    base = label.split("\n", 1)[0].strip()
    if FIGURE2_TRACK_R_GAP in base:
        return base.split(FIGURE2_TRACK_R_GAP, 1)[0].strip()
    return base


def plot_profile_tracks(
    tracks: list[tuple[str, np.ndarray, bool]],
    title: str,
    out_png: Path,
    *,
    model_colors: dict[str, str],
    dpi: int = 200,
    figsize_inch: tuple[float, float] | None = None,
    vector_pdf: Path | None = None,
    font_scale: float = 1.0,
    figsize_scale: float | None = None,
    ytick_rotation: float = 0.0,
) -> None:
    """tracks: (label, y_array, invert)。invert=True 时绘制向下填充（y 应 ≤0，常用 −|Y|）。"""
    assert len(tracks) >= 1
    L = int(tracks[0][1].shape[0])
    x = np.arange(L)

    n = len(tracks)
    # 默认宽:高 = 2:1；可由 figsize_inch 覆盖。字号用 font_scale；画布勿同比放大否则整图缩放预览时字相对图区几乎不变
    _bw, _bh = figsize_inch if figsize_inch is not None else (10.0, 5.0)
    if figsize_scale is None:
        figsize_scale = min(2.2, 1.0 + 0.16 * max(0.0, float(font_scale) - 1.0))
    fig_w, fig_h = _bw * figsize_scale, _bh * figsize_scale
    # 扁长画布时略收紧轨间距；Y 轴刻度斜排时略加大，避免相邻轨左侧数字叠影
    _ar = fig_w / fig_h
    _hspace = 0.06 if _ar >= 3.0 else (0.08 if _ar >= 2.4 else 0.12)
    if ytick_rotation and n > 1:
        _hspace += 0.04
    fig, axes = plt.subplots(
        n,
        1,
        figsize=(fig_w, fig_h),
        sharex=True,
        gridspec_kw={"hspace": _hspace, "height_ratios": [1] * n},
    )
    if n == 1:
        axes = [axes]

    def _fill_color(lab: str) -> str:
        key = _track_name_for_color(lab)
        if key in ("Signal", "eCLIP signal"):  # Signal 兼容旧调用
            return _ECLIP_MAIN_FILL
        if key.startswith("Encode replicate"):
            return _ENCODE_REP2_FILL if key.endswith(" 2") else _ENCODE_REP_FILL
        return model_colors.get(key, FIG01_VIOLIN_COLOR.get(key, "#cccccc"))

    fs = 9 * font_scale
    for ax, (label, y, invert) in zip(axes, tracks):
        c = _fill_color(label)
        y = np.asarray(y, dtype=np.float64).ravel()
        if invert:
            y = np.minimum(y, 0.0)
            ax.fill_between(x, 0.0, y, where=(y <= 0.0), color=c, alpha=1.0, linewidth=0)
            ymin = float(np.min(y)) if y.size else -1.0
            pad = max((0.0 - ymin) * 0.06, 1e-6)
            ax.set_ylim(ymin - pad, pad * 0.35)
        else:
            y = np.maximum(y, 0.0)
            ax.fill_between(x, 0.0, y, where=(y >= 0.0), color=c, alpha=1.0, linewidth=0)
            ymax = float(np.max(y)) if y.size else 1.0
            pad = max(ymax * 0.06, 1e-6)
            ax.set_ylim(-pad * 0.02, ymax + pad)
        ax.tick_params(axis="x", labelsize=fs)
        ax.tick_params(axis="y", labelsize=fs)
        if ytick_rotation:
            for _yl in ax.get_yticklabels():
                _yl.set_rotation(ytick_rotation)
                _yl.set_horizontalalignment("right")
                _yl.set_rotation_mode("anchor")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.margins(x=0)
        ax.axhline(0.0, color="#888888", linewidth=0.6, zorder=1)
        # 轨名在子图右侧；模型轨为「名称 + 两个空格 + r = …」同一字符串
        ax.text(
            1.02,
            0.5,
            label,
            transform=ax.transAxes,
            fontsize=fs,
            ha="left",
            va="center",
            clip_on=False,
        )

    axes[-1].set_xlabel("Position", fontsize=fs)
    step = max(L // 5, 1)
    axes[-1].set_xticks(np.arange(0, L, step))
    axes[0].set_title(title, fontsize=fs)
    _pad_left = 0.10
    if ytick_rotation:
        _pad_left = min(0.26, 0.10 + abs(ytick_rotation) / 45.0 * 0.14)
    _h_adj = 0.22
    if ytick_rotation and n > 1:
        _h_adj = 0.30
    fig.subplots_adjust(
        left=_pad_left, right=0.64, top=0.90, bottom=0.10, hspace=_h_adj
    )
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out_png,
        dpi=dpi,
        bbox_inches=None,
        pad_inches=0,
        facecolor="white",
        edgecolor="none",
    )
    if vector_pdf is not None:
        vector_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            vector_pdf,
            format="pdf",
            bbox_inches=None,
            pad_inches=0,
            facecolor="white",
            edgecolor="none",
        )
    plt.close(fig)


def plot_four_tracks(
    tracks: list[tuple[str, np.ndarray]],
    title: str,
    out_png: Path,
    dpi: int = 200,
    vector_pdf: Path | None = None,
    font_scale: float = 1.0,
    figsize_scale: float | None = None,
    ytick_rotation: float = 0.0,
) -> None:
    """兼容旧接口：四轨、无倒置；Signal 轨仍用原浅色名。"""
    model_colors = dict(FIG01_VIOLIN_COLOR)
    spec = [(lab, y, False) for lab, y in tracks]
    # 将旧标签 Signal 映射为 eCLIP signal 的可选统一命名在调用处处理
    plot_profile_tracks(
        spec,
        title,
        out_png,
        model_colors=model_colors,
        dpi=dpi,
        vector_pdf=vector_pdf,
        font_scale=font_scale,
        figsize_scale=figsize_scale,
        ytick_rotation=ytick_rotation,
    )


def _resolve_out_path(base_out: Path, batch_index: int, n_outputs: int) -> Path:
    """多图时在文件名中插入序号：foo.png → foo_01.png；仅输出 1 张时仍为 foo.png。"""
    if n_outputs <= 1:
        return base_out
    return base_out.parent / f"{base_out.stem}_{batch_index:02d}{base_out.suffix}"


def main() -> None:
    parser = argparse.ArgumentParser(description="四轨 Signal + MARCP + iDeepB + RBPNet 强度图。")
    parser.add_argument("--task", default=DEFAULT_TASK)
    parser.add_argument("--marcp-model", type=Path, default=None)
    parser.add_argument("--ideep-model", type=Path, default=None)
    parser.add_argument("--rbpnet-model", type=Path, default=None)
    parser.add_argument("--test-bed", type=Path, default=None)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "figures/baseline_pred_signal_track.png",
    )
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument(
        "--rbp-python",
        type=Path,
        default=RBP_ROOT / "venv" / "bin" / "python3",
        help="含 TensorFlow/gin 的 Python（默认同比较实验 RBPNet venv）",
    )
    parser.add_argument(
        "--count",
        type=int,
        default=1,
        metavar="N",
        help="生成 N 张图（不同测试窗口）；默认 1。",
    )
    parser.add_argument(
        "--selection",
        choices=("top", "spread"),
        default="top",
        help="top=按 max(signal) 最强的 N 个；spread=在排名上均匀取 N 个，兼顾强/中/弱峰。",
    )
    parser.add_argument(
        "--min-peak",
        type=float,
        default=0.0,
        help="仅保留 max(Y) ≥ 此值的窗口（默认 0，不过滤）。",
    )
    parser.add_argument(
        "--figure2-export",
        action="store_true",
        help=(
            "输出论文 Figure2 样式：六轨、模型色为 图/参考配色 NPG(1–3)；"
            "eCLIP 下方两条倒置轨为 −|Y| 与 −|Y_ctl|（Y_ctl 为预处理中的 control 剖面，"
            "非 ENCODE 生物学重复，仅版式对齐示意图）。整图仍保持宽:高=2:1。"
        ),
    )
    parser.add_argument(
        "--figure1b",
        action="store_true",
        help=(
            "论文 Figure 1b：仅两轨 MARCP 预测 + eCLIP signal（与 --figure2-export 同色与版式比例），"
            "不加载 iDeepB/RBPNet，轨标签不含 Pearson r。"
        ),
    )
    parser.add_argument(
        "--vector-pdf-dir",
        type=Path,
        default=None,
        help="若指定，则每张 PNG 同步输出同名 .pdf 矢量到此目录（如 图/figure2/结果）。",
    )
    parser.add_argument(
        "--font-scale",
        type=float,
        default=1.0,
        help="轴标签/轨名/刻度字号与画布（英寸）相对默认的缩放；Figure2 当前约 1.625。",
    )
    parser.add_argument(
        "--ytick-rotation",
        type=float,
        default=0.0,
        help="各子图左侧 Y 轴刻度数值旋转角度（度）；Figure2 常用 45。",
    )
    args = parser.parse_args()
    if args.figure1b and args.figure2_export:
        parser.error("--figure1b 与 --figure2-export 互斥，请只选其一。")

    task = args.task
    n_batch = max(1, int(args.count))
    bed_path = args.test_bed or (BED_ROOT / task / f"{task}.101.peak_test.seqBw.treat_ctl.bed")
    if not bed_path.is_file():
        raise FileNotFoundError(bed_path)

    test_keys = _read_test_keys(bed_path)
    all_rows = collect_test_rows_sorted_by_peak(task, test_keys, TEST_CHRS)
    if args.min_peak > 0:
        all_rows = [r for r in all_rows if r[0] >= args.min_peak]
        if not all_rows:
            raise RuntimeError(f"--min-peak {args.min_peak} 过滤后无可用测试窗口。")

    if args.selection == "top":
        chosen = all_rows[: min(n_batch, len(all_rows))]
    else:
        chosen = pick_spread_by_rank(all_rows, n_batch)

    if len(chosen) < n_batch:
        print(
            f"[warn] 测试交集仅 {len(all_rows)} 条，少于 --count {n_batch}，实际输出 {len(chosen)} 张。",
            file=sys.stderr,
        )

    marcp_p = resolve_marcp_weights(task, args.marcp_model)
    if args.figure1b:
        ideep_p = None
        rbpn_p = None
    else:
        ideep_p = resolve_ideep_weights(task, args.ideep_model)
        rbpn_p = resolve_rbpnet_weights(task, args.rbpnet_model)

    dev = torch.device(args.device)
    marcp_model = load_marcp_model(marcp_p, dev)
    ideep_model = None if args.figure1b else load_ideepb_model(ideep_p, dev)

    model_colors = (
        dict(FIGURE2_MODEL_COLOR)
        if (args.figure2_export or args.figure1b)
        else dict(FIG01_VIOLIN_COLOR)
    )

    written: list[Path] = []
    for bi, (_peak, meta, X_5ch, Y, Y_ctl) in enumerate(chosen, start=1):
        chrom = _norm_chrom(meta["chr"])
        start = int(meta["chromStart"])
        end = int(meta["chromEnd"])
        strand = str(meta["strand"])
        seq_table = str(meta["seq"])

        pred_m = forward_marcp(marcp_model, X_5ch, dev)
        if not args.figure1b:
            assert ideep_model is not None and ideep_p is not None and rbpn_p is not None
            oh = seq_to_ideep_onehot(seq_table)
            pred_i = forward_ideepb(ideep_model, oh, dev)
            seq_rbp = seq_table.upper().replace("U", "T")
            pred_r = predict_rbpnet_profile(seq_rbp, rbpn_p, args.rbp_python)
            if pred_m.shape[0] != Y.shape[0] or pred_i.shape[0] != Y.shape[0]:
                raise ValueError(
                    f"长度不一致: MARCP {pred_m.shape}, iDeepB {pred_i.shape}, signal {Y.shape}"
                )
            if pred_r.shape[0] != Y.shape[0]:
                raise ValueError(f"RBPNet 输出长度 {pred_r.shape} 与窗口 {Y.shape} 不一致")
        elif pred_m.shape[0] != Y.shape[0]:
            raise ValueError(f"长度不一致: MARCP {pred_m.shape}, signal {Y.shape}")

        title = f"{_figure_display_task_name(task)} chr{chrom}:{start}-{end}:{strand}"
        out_png = _resolve_out_path(args.out, bi, len(chosen))
        vpdf = None
        if args.vector_pdf_dir is not None:
            args.vector_pdf_dir.mkdir(parents=True, exist_ok=True)
            vpdf = args.vector_pdf_dir / f"{out_png.stem}.pdf"
        if args.figure1b:
            plot_profile_tracks(
                [
                    ("CARBON", pred_m, False),
                    ("eCLIP signal", Y, False),
                ],
                title,
                out_png,
                model_colors=model_colors,
                dpi=args.dpi,
                figsize_inch=FIGURE1B_FIGSIZE_INCH,
                vector_pdf=vpdf,
                font_scale=args.font_scale,
                ytick_rotation=args.ytick_rotation,
            )
        elif args.figure2_export:
            y_rep1 = -np.abs(np.asarray(Y, dtype=np.float64))
            y_rep2 = -np.abs(np.asarray(Y_ctl, dtype=np.float64))
            rm = _pearson_r_vs_signal(pred_m, Y)
            ri = _pearson_r_vs_signal(pred_i, Y)
            rr = _pearson_r_vs_signal(pred_r, Y)
            def _fmt_r(name: str, r: float) -> str:
                rhs = f"r = {r:.2f}" if np.isfinite(r) else "r = NA"
                return FIGURE2_TRACK_R_GAP.join((name.strip(), rhs))
            # 自上向下：CARBON → eCLIP → rep → iDeepB（倒数第二）→ RBPNet（最下）
            plot_profile_tracks(
                [
                    (_fmt_r("CARBON", rm), pred_m, False),
                    ("eCLIP signal", Y, False),
                    ("Encode replicate 1", y_rep1, True),
                    ("Encode replicate 2", y_rep2, True),
                    (_fmt_r("iDeepB", ri), pred_i, False),
                    (_fmt_r("RBPNet", rr), pred_r, False),
                ],
                title,
                out_png,
                model_colors=model_colors,
                dpi=args.dpi,
                vector_pdf=vpdf,
                font_scale=args.font_scale,
                ytick_rotation=args.ytick_rotation,
            )
        else:
            plot_four_tracks(
                [
                    ("eCLIP signal", Y),
                    ("CARBON", pred_m),
                    ("iDeepB", pred_i),
                    ("RBPNet", pred_r),
                ],
                title,
                out_png,
                dpi=args.dpi,
                vector_pdf=vpdf,
                font_scale=args.font_scale,
                ytick_rotation=args.ytick_rotation,
            )
        written.append(out_png)
        print(out_png)
        print(f"  {title}  max_signal={_peak:.3g}")

    print(f"Wrote {len(written)} file(s).")
    print(f"  MARCP 权重: {marcp_p}")
    if not args.figure1b:
        print(f"  iDeepB 权重: {ideep_p}")
        print(f"  RBPNet 权重: {rbpn_p}")


if __name__ == "__main__":
    main()
