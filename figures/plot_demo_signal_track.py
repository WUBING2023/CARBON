#!/usr/bin/env python3
"""
Demo: load a window from this repo's 5ch h5 and plot the same layout as the first four
tracks in paper Figure 2:
  1) CARBON prediction  2) eCLIP signal  3) Encode replicate 1 (−|Y|, inverted)
  4) Encode replicate 2 (−|Y_ctl|, inverted)

No iDeepB / RBPNet. Run scripts/run_demo_train.sh first to produce model_params.pth,
or pass --weights to a checkpoint.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from CARBON.models.carbon_model import CARBON  # noqa: E402

# NPG-style colors (aligned with paper Figure 2)
MODEL_COLOR = {"CARBON": "#E64B35"}
_ECLIP_MAIN_FILL = "#91D1C2"
_ENCODE_REP_FILL = "#8491B4"
_ENCODE_REP2_FILL = "#3C5488"
FIGURE2_TRACK_R_GAP = "  "


def _pearson_r(pred: np.ndarray, y_signal: np.ndarray) -> float:
    a = np.asarray(pred, dtype=np.float64).ravel()
    b = np.asarray(y_signal, dtype=np.float64).ravel()
    if a.size != b.size or a.size < 2:
        return float("nan")
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    r = float(np.corrcoef(a, b)[0, 1])
    return r if np.isfinite(r) else float("nan")


def _track_name_for_color(label: str) -> str:
    base = label.split("\n", 1)[0].strip()
    if FIGURE2_TRACK_R_GAP in base:
        return base.split(FIGURE2_TRACK_R_GAP, 1)[0].strip()
    return base


def plot_profile_tracks(
    tracks: list[tuple[str, np.ndarray, bool]],
    title: str,
    out_png: Path,
    *,
    dpi: int = 300,
    font_scale: float = 1.625,
    ytick_rotation: float = 45.0,
) -> None:
    assert len(tracks) >= 1
    L = int(tracks[0][1].shape[0])
    x = np.arange(L)
    n = len(tracks)
    _bw, _bh = 10.0, 5.0
    figsize_scale = min(2.2, 1.0 + 0.16 * max(0.0, float(font_scale) - 1.0))
    fig_w, fig_h = _bw * figsize_scale, _bh * figsize_scale
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
        if key in ("Signal", "eCLIP signal"):
            return _ECLIP_MAIN_FILL
        if key.startswith("Encode replicate"):
            return _ENCODE_REP2_FILL if key.endswith(" 2") else _ENCODE_REP_FILL
        return MODEL_COLOR.get(key, "#cccccc")

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
    fig.subplots_adjust(left=_pad_left, right=0.64, top=0.90, bottom=0.10, hspace=_h_adj)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi, bbox_inches=None, pad_inches=0, facecolor="white", edgecolor="none")
    plt.close(fig)


def load_h5_window(h5_path: Path, index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(h5_path, "r") as hf:
        X = np.asarray(hf["X"][index], dtype=np.float32)
        Y = np.asarray(hf["Y"][index], dtype=np.float32)
        Y_ctl = np.asarray(hf["Y_ctl"][index], dtype=np.float32)
    return X, Y, Y_ctl


def pick_best_index(h5_path: Path) -> int:
    """Index of the window with largest max(Y) (strongest peak)."""
    with h5py.File(h5_path, "r") as hf:
        Yall = np.asarray(hf["Y"], dtype=np.float32)
    peaks = np.max(Yall, axis=1)
    return int(np.argmax(peaks))


def forward_carbon(model: torch.nn.Module, X_5ch: np.ndarray, device: torch.device) -> np.ndarray:
    t = torch.from_numpy(X_5ch[np.newaxis, ...]).float().to(device)
    with torch.no_grad():
        out = model(t)
        logits = out["logits_profile"]
        log_tc = out["log_total_count"]
        prob = torch.softmax(logits, dim=-1)
        pred = prob * torch.exp(log_tc)
    return pred.squeeze(0).detach().cpu().numpy()


def main() -> None:
    p = argparse.ArgumentParser(
        description="Demo: four-track prediction intensity figure (matches Figure 2 first four tracks)"
    )
    p.add_argument(
        "--weights",
        type=Path,
        default=REPO_ROOT / "carbon_trained_model/AARS_K562_ENCSR825SVO_demo/model_params.pth",
        help="CARBON model_params.pth",
    )
    p.add_argument(
        "--h5",
        type=Path,
        default=REPO_ROOT / "data/h5data/AARS_K562_ENCSR825SVO/_5ch/AARS_K562_ENCSR825SVO_chr1_5ch.h5",
        help="5ch h5 to sample a window from (default chr1)",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "figures/outputs/demo_signal_track_4rows.png",
    )
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    p.add_argument(
        "--index",
        type=int,
        default=-1,
        help="Row index in h5; default -1 picks strongest peak by max(Y)",
    )
    args = p.parse_args()

    if not args.weights.is_file():
        raise SystemExit(
            f"Weights not found: {args.weights}\nRun first: bash scripts/run_demo_train.sh"
        )
    if not args.h5.is_file():
        raise SystemExit(f"h5 not found: {args.h5}")

    idx = pick_best_index(args.h5) if args.index < 0 else args.index
    X, Y, Y_ctl = load_h5_window(args.h5, idx)
    dev = torch.device(args.device)
    model = CARBON(dropOut=0.2, control=False, in_channels=5).to(dev)
    try:
        state = torch.load(args.weights, map_location=dev, weights_only=True)
    except TypeError:
        state = torch.load(args.weights, map_location=dev)
    model.load_state_dict(state)
    model.eval()
    pred_m = forward_carbon(model, X, dev)

    rm = _pearson_r(pred_m, Y)
    rhs = f"r = {rm:.2f}" if np.isfinite(rm) else "r = NA"
    label_carbon = FIGURE2_TRACK_R_GAP.join(("CARBON", rhs))

    y_rep1 = -np.abs(np.asarray(Y, dtype=np.float64))
    y_rep2 = -np.abs(np.asarray(Y_ctl, dtype=np.float64))

    title = f"AARS_K562 chr1 (demo window idx={idx})"
    plot_profile_tracks(
        [
            (label_carbon, pred_m, False),
            ("eCLIP signal", Y, False),
            ("Encode replicate 1", y_rep1, True),
            ("Encode replicate 2", y_rep2, True),
        ],
        title,
        args.out,
        dpi=300,
        font_scale=1.625,
        ytick_rotation=45.0,
    )
    print(args.out)


if __name__ == "__main__":
    main()
