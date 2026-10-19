#!/usr/bin/env python3
"""
从 modisco_results.h5 按 --mode A|B|C 读取矩阵，经 MARCP ic_scale 得 IC 加权矩阵，再用 logomaker 绘制序列 logo
（非 report_simple/trimmed_logos 的 CWM；碱基字形由 logomaker 统一渲染，含 U）。

在 Motif 根目录执行示例：
  python3 scripts/plot_modisco_ic_logos.py --task-dir outputs/RBFOX2_K562_ENCSR756CKJ --mode C --pattern pattern_0

依赖：pip install logomaker（见 requirements-motif.txt）
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import logomaker

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from motif_analysis.modisco_h5_patterns import (
    ICLogoMode,
    collect_pattern_records,
    load_ppm_for_ic_logo,
    reverse_complement_ppm_acgu,
)
from motif_analysis.paths import ensure_marcp_on_path

RNA_BG = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)

# 与原先 MARCP/tfmodisco 风格一致：A 绿、C 蓝、G 橙、U 红
LOGO_COLOR_SCHEME = {"A": "green", "C": "blue", "G": "orange", "U": "red"}


def group_path_to_filename_stem(group_path: str) -> str:
    """pos_patterns/pattern_0 -> pos_patterns.pattern_0"""
    parts = group_path.strip("/").split("/")
    if not parts:
        return "pattern"
    return parts[0] + "." + ".".join(parts[1:]) if len(parts) > 1 else parts[0]


def pattern_filter_matches(record_pattern_id: str, h5_group_path: str, needle: str) -> bool:
    """needle 可为 pattern_0、pos_pattern_0、或子串。"""
    n = needle.strip()
    if not n:
        return True
    if record_pattern_id == n:
        return True
    if h5_group_path == n or h5_group_path.endswith("/" + n):
        return True
    if n in record_pattern_id or n in h5_group_path:
        return True
    return False


def plot_ic_logo_png(
    out_path: Path,
    ppm_pre_ic: np.ndarray,
    *,
    title: str | None = None,
    figsize: tuple[float, float] | None = None,
    dpi: int = 300,
    ylabel: str | None = "IC bits",
    show_axes: bool = True,
) -> None:
    ensure_marcp_on_path()
    from MARCP.visualization.viz_sequence import ic_scale

    w = ic_scale(np.asarray(ppm_pre_ic, dtype=np.float64), RNA_BG)
    w = np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)
    w = np.maximum(w, 0.0)

    npos = int(w.shape[0])
    df = pd.DataFrame(w, columns=list("ACGU"))

    if figsize is None:
        fig_w = max(8.0, npos * 0.35)
        fig, ax = plt.subplots(figsize=(fig_w, 2.8))
    else:
        fig, ax = plt.subplots(figsize=figsize)
    logomaker.Logo(
        df,
        ax=ax,
        color_scheme=LOGO_COLOR_SCHEME,
        font_name="sans",
        stack_order="big_on_top",
        vsep=0.0,
    )
    if show_axes:
        if ylabel:
            ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title, fontsize=11)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="both", which="both", labelbottom=False, labelleft=False, length=0)
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=int(dpi), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="IC-scaled PPM sequence logos (MARCP ic_scale + logomaker) from modisco_results.h5; modes A/B/C."
    )
    ap.add_argument("--mode", type=str.upper, choices=("A", "B", "C"), required=True)
    ap.add_argument("--task-dir", type=Path, default=None, help="任务目录（将使用 tfmodisco/modisco_results.h5）")
    ap.add_argument("--h5", type=Path, default=None, help="直接指定 modisco_results.h5")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="输出目录，默认 <task-dir>/tfmodisco/ic_scaled_logos_mode_<MODE>/",
    )
    ap.add_argument(
        "--pattern",
        type=str,
        default="",
        help="只处理匹配的 pattern（如 pattern_0 或 pos_pattern_0）",
    )
    ap.add_argument("--include-neg", action="store_true")
    ap.add_argument("--include-subpatterns", action="store_true")
    ap.add_argument("--include-rev", action="store_true", help="额外输出反向互补链 logo")
    args = ap.parse_args()
    mode: ICLogoMode = args.mode  # type: ignore[assignment]

    if args.h5 is not None:
        h5 = args.h5.resolve()
    elif args.task_dir is not None:
        h5 = args.task_dir.resolve() / "tfmodisco" / "modisco_results.h5"
    else:
        print("需要 --task-dir 或 --h5", file=sys.stderr)
        raise SystemExit(2)

    if not h5.is_file():
        print(f"[error] missing {h5}", file=sys.stderr)
        raise SystemExit(1)

    task_root = h5.parent.parent if h5.parent.name == "tfmodisco" else h5.parent
    out_dir = args.out_dir
    if out_dir is None:
        if args.task_dir is None:
            out_dir = h5.parent / f"ic_scaled_logos_mode_{mode}"
        else:
            out_dir = args.task_dir.resolve() / "tfmodisco" / f"ic_scaled_logos_mode_{mode}"
    out_dir = out_dir.resolve()

    records = collect_pattern_records(
        h5,
        include_neg=args.include_neg,
        include_subpatterns=args.include_subpatterns,
    )
    needle = args.pattern.strip()
    if needle:
        records = [r for r in records if pattern_filter_matches(r.pattern_id, r.h5_group_path, needle)]
    if not records:
        print("[error] no patterns after filter", file=sys.stderr)
        raise SystemExit(1)

    safe_mode = re.sub(r"[^\w]+", "", mode)
    for rec in records:
        stem = group_path_to_filename_stem(rec.h5_group_path)
        try:
            ppm_fwd = load_ppm_for_ic_logo(h5, rec.h5_group_path, mode)
        except (KeyError, ValueError, TypeError) as e:
            print(f"[skip] {rec.h5_group_path}: {e}", file=sys.stderr)
            continue

        title = f"{stem} mode={safe_mode} n_seqlets={rec.n_seqlets}"
        plot_ic_logo_png(out_dir / f"{stem}.ic_fwd.png", ppm_fwd, title=title)
        print(f"[ok] {out_dir / f'{stem}.ic_fwd.png'}")

        if args.include_rev:
            ppm_rev = reverse_complement_ppm_acgu(ppm_fwd)
            plot_ic_logo_png(out_dir / f"{stem}.ic_rev.png", ppm_rev, title=title + " rev")
            print(f"[ok] {out_dir / f'{stem}.ic_rev.png'}")


if __name__ == "__main__":
    main()
