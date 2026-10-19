#!/usr/bin/env python3
"""
对 match_motif_p_sorted.csv 中 p < 0.05 的每一行生成四列对照图：
蛋白质名 | TF-Modisco pattern CWM logo (trimmed) | CISBP 标准 motif（PWM→logomaker，**纵向**高度为第 3 列单元格的 40%，水平宽度与第 2 列一致）| p / E / q 三行文本。

  python3 scripts/plot_match_motif_panels_lt005.py
  python3 scripts/plot_match_motif_panels_lt005.py --csv reference_data/motif_p_value/match_motif_p_sorted.csv \\
      --out-dir reference_data/motif_p_value/match_p_lt005_panels
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import logomaker

from motif_analysis.cisbp_compare import load_cisbp_pwm_txt

LOGO_COLOR_SCHEME = {"A": "green", "C": "blue", "G": "orange", "U": "red"}
_DEFAULT_PWM_DIR = ROOT / "reference_data" / "CISBP-RNA" / "homo_sapiens" / "pwms"

# 第 3 列内 CISBP logo 子轴高度占该列单元格高度的比例（缩短上下方向）
_CISBP_LOGO_HEIGHT_FRAC = 0.4

# 单张面板画布：宽:高 = 2.5 : 1（英寸）
_PANEL_FIG_W_INCH = 12.5
_PANEL_FIG_H_INCH = _PANEL_FIG_W_INCH / 2.5


def _fmt_stat(x: str) -> str:
    v = float(x)
    if v == 0.0:
        return "0"
    av = abs(v)
    if av >= 0.0001 and av < 100000:
        return f"{v:.6g}"
    return f"{v:.3e}"


def query_id_to_trimmed_fwd_png(task: str, query_id: str, outputs_root: Path) -> Path:
    """pos_pattern_N -> pos_patterns.pattern_N.cwm.fwd.png"""
    num = query_id.replace("pos_pattern_", "pattern_")
    stem = f"pos_patterns.{num}.cwm.fwd.png"
    return outputs_root / task / "tfmodisco" / "report_simple" / "trimmed_logos" / stem


def safe_stem(s: str, max_len: int = 120) -> str:
    t = re.sub(r"[^\w.\-]+", "_", s.strip())
    return t[:max_len] if len(t) > max_len else t


def plot_cisbp_ppm_logo(ax, ppm: np.ndarray) -> None:
    df = pd.DataFrame(np.asarray(ppm, dtype=np.float64), columns=list("ACGU"))
    df.index.name = "pos"
    row_sums = df.sum(axis=1)
    logo = logomaker.Logo(
        df,
        ax=ax,
        color_scheme=LOGO_COLOR_SCHEME,
        font_name="sans",
        stack_order="big_on_top",
        vsep=0.0,
    )
    logo.style_spines(visible=False)
    ax.set_ylim(min(row_sums.min(), 0), row_sums.max())
    ax.set_xticks([])
    ax.set_yticks([])


def plot_one_panel(
    *,
    protein: str,
    pattern_png: Path,
    pwm_path: Path,
    p_val: str,
    e_val: str,
    q_val: str,
    out_png: Path,
) -> None:
    img = plt.imread(str(pattern_png))
    ppm = load_cisbp_pwm_txt(pwm_path)
    if ppm is None or ppm.shape[0] == 0:
        raise ValueError(f"empty PWM: {pwm_path}")

    fig = plt.figure(figsize=(_PANEL_FIG_W_INCH, _PANEL_FIG_H_INCH))
    gs = fig.add_gridspec(1, 4, width_ratios=[0.72, 1.15, 1.15, 0.82], wspace=0.22)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.axis("off")
    ax0.text(0.5, 0.5, protein, ha="center", va="center", fontsize=12, fontweight="600")

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.imshow(img)
    ax1.axis("off")

    ax2_cell = fig.add_subplot(gs[0, 2])
    ax2_cell.axis("off")
    ax2 = inset_axes(
        ax2_cell,
        width="100%",
        height=f"{100 * _CISBP_LOGO_HEIGHT_FRAC:.0f}%",
        loc="center",
        borderpad=0,
    )
    plot_cisbp_ppm_logo(ax2, ppm)

    ax3 = fig.add_subplot(gs[0, 3])
    ax3.axis("off")
    lines = (
        f"p = {_fmt_stat(p_val)}",
        f"E = {_fmt_stat(e_val)}",
        f"q = {_fmt_stat(q_val)}",
    )
    for i, line in enumerate(lines):
        ax3.text(0.02, 0.78 - i * 0.28, line, transform=ax3.transAxes, va="top", fontsize=10.5, family="monospace")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    # 不用 bbox_inches="tight"，否则裁边后像素宽高比会偏离 figsize 的 2.5:1
    fig.savefig(out_png, dpi=160, facecolor="white")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=ROOT / "reference_data/motif_p_value/match_motif_p_sorted.csv")
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "reference_data/motif_p_value/match_p_lt005_panels",
    )
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument("--pwm-dir", type=Path, default=_DEFAULT_PWM_DIR)
    ap.add_argument("--threshold", type=float, default=0.05)
    ap.add_argument("--max-rows", type=int, default=0, help="0=全部")
    args = ap.parse_args()

    rows: list[dict[str, str]] = []
    with args.csv.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if float(row["p_value"]) < args.threshold:
                rows.append(row)

    if args.max_rows > 0:
        rows = rows[: args.max_rows]

    manifest = args.out_dir / "manifest.csv"
    args.out_dir.mkdir(parents=True, exist_ok=True)

    ok = 0
    with manifest.open("w", newline="", encoding="utf-8") as mf:
        w = csv.writer(mf)
        w.writerow(["index", "task", "Query_ID", "Target_ID", "p_value", "filename", "status"])
        for i, row in enumerate(rows, start=1):
            task = row["task"]
            qid = row["Query_ID"]
            tid = row["Target_ID"]
            prot = row["cisbp_RBP_Name"] or task.split("_")[0]
            pat_png = query_id_to_trimmed_fwd_png(task, qid, args.outputs_root)
            pwm_path = args.pwm_dir / f"{tid}.txt"
            fname = f"{i:04d}_{safe_stem(prot)}_{safe_stem(qid)}_{safe_stem(tid)}.png"
            out_png = args.out_dir / fname
            try:
                plot_one_panel(
                    protein=prot,
                    pattern_png=pat_png,
                    pwm_path=pwm_path,
                    p_val=row["p_value"],
                    e_val=row["E-value"],
                    q_val=row["q-value"],
                    out_png=out_png,
                )
                status = "ok"
                ok += 1
            except Exception as e:
                status = f"error: {e}"
                print(f"[{i}/{len(rows)}] {status} ({task} {qid} {tid})", file=sys.stderr)
            w.writerow([i, task, qid, tid, row["p_value"], fname, status])

    print(f"Wrote {ok}/{len(rows)} panels under {args.out_dir}")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
