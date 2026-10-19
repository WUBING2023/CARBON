#!/usr/bin/env python3
"""
将某任务 outputs 下 conv_motifs + seqlet_motifs 的 PPM 与 CISBP-RNA pwms/ 全库扫描比对，
生成 cisbp_compare/comparison_results.json 与 comparison_report.html，并为每个 query 画「模型 vs 最佳 CISBP」并排图。

用法（在 Motif 根目录）：
  python3 scripts/compare_outputs_to_cisbp.py \\
    --task-dir outputs/AGGF1_HepG2_ENCSR543TPH
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

# Motif 包根目录
_MOTIF_ROOT = Path(__file__).resolve().parents[1]
if str(_MOTIF_ROOT) not in sys.path:
    sys.path.insert(0, str(_MOTIF_ROOT))

from motif_analysis.cisbp_compare import (  # noqa: E402
    extract_aligned_segments,
    load_all_cisbp_pwms,
    rank_queries_against_cisbp,
    save_compare_report_json,
)


def _load_queries_from_task(task_dir: Path, only_existing_png: bool = True) -> list[tuple[str, np.ndarray]]:
    """
    默认仅包含 conv_motifs 中已导出 logo_*.png 的卷积核（与「目录里有哪些图」一致），
    另加全部 seqlet motif。
    """
    out: list[tuple[str, np.ndarray]] = []
    conv_json = task_dir / "conv_motifs" / "conv_kernel_ppms.json"
    logo_keys: set[tuple[int, int, int]] | None = None
    if only_existing_png:
        pat = re.compile(r"logo_branch(\d+)_k(\d+)_f(\d+)\.png$")
        logo_keys = set()
        for png in (task_dir / "conv_motifs").glob("logo_*.png"):
            m = pat.search(png.name)
            if m:
                logo_keys.add((int(m.group(1)), int(m.group(2)), int(m.group(3))))
        if not logo_keys:
            logo_keys = None
    if conv_json.is_file():
        with open(conv_json, encoding="utf-8") as f:
            branches = json.load(f)
        for b in branches:
            bi = b["branch_index"]
            ks = b["kernel_size"]
            for filt in b.get("filters", []):
                fi = filt["filter_index"]
                if logo_keys is not None and (bi, ks, fi) not in logo_keys:
                    continue
                ppm = np.asarray(filt["ppm"], dtype=np.float64)
                out.append((f"conv_branch{bi}_k{ks}_f{fi}", ppm))
    seq_npz = task_dir / "seqlet_motifs" / "seqlet_clusters.npz"
    if seq_npz.is_file():
        z = np.load(seq_npz)
        ppms = z["ppms"]
        for i in range(ppms.shape[0]):
            out.append((f"seqlet_motif_{i}", ppms[i].astype(np.float64)))
    return out


def _plot_side_by_side(
    left: np.ndarray,
    right: np.ndarray,
    left_title: str,
    right_title: str,
    out_png: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from motif_analysis.paths import ensure_marcp_on_path

    ensure_marcp_on_path()
    from MARCP.visualization.viz_sequence import plot_weights_given_ax

    fig, axes = plt.subplots(1, 2, figsize=(max(10, left.shape[0] + right.shape[0]), 2.8))
    plot_weights_given_ax(axes[0], left.astype(np.float32), ylabel=left_title)
    plot_weights_given_ax(axes[1], right.astype(np.float32), ylabel=right_title)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def _build_html(task_name: str, rows: list[dict]) -> str:
    lines = [
        "<!DOCTYPE html><html><head><meta charset='utf-8'/>",
        "<title>CISBP-RNA 扫描比对</title>",
        "<style>body{font-family:sans-serif;margin:24px;} table{border-collapse:collapse;} th,td{border:1px solid #ccc;padding:6px;} img{max-width:100%;height:auto;}</style>",
        "</head><body>",
        f"<h1>CISBP-RNA 全库扫描比对：{task_name}</h1>",
        "<p>说明：对每条模型 PPM 在 CISBP Build 2.00 全 PWM 上做<strong>滑动窗口最小对称 KL</strong>（越小越相似）。",
        "CISBP 中未必有同名 RBP；本表为<strong>全库最佳匹配</strong>，供形态对照，非官方 TOMTOM p 值。</p>",
        "<table><tr><th>Query</th><th>Top-1 CISBP ID</th><th>mean sym-KL</th><th>并排图</th></tr>",
    ]
    for r in rows:
        qid = r["query_id"]
        top = r["top"][0] if r["top"] else None
        if not top:
            continue
        img_rel = f"figures/{qid.replace('/', '_')}_vs_{top['cisbp_id']}.png"
        lines.append(
            "<tr>"
            f"<td>{qid}</td>"
            f"<td>{top['cisbp_id']}</td>"
            f"<td>{top['mean_sym_kl']:.4f}</td>"
            f"<td><img src='{img_rel}' alt=''/></td>"
            "</tr>"
        )
    lines.append("</table></body></html>")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--task-dir",
        type=Path,
        default=_MOTIF_ROOT / "outputs/AGGF1_HepG2_ENCSR543TPH",
        help="含 conv_motifs / seqlet_motifs 的任务目录",
    )
    ap.add_argument(
        "--cisbp-pwms",
        type=Path,
        default=_MOTIF_ROOT / "reference_data/CISBP-RNA/pwms",
        help="CISBP 解压后的 pwms 目录",
    )
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument(
        "--all-conv-filters",
        action="store_true",
        help="包含 json 中全部卷积核（不仅限于已导出的 logo 图）",
    )
    args = ap.parse_args()

    task_dir = args.task_dir.resolve()
    cisbp_dir = args.cisbp_pwms.resolve()
    if not cisbp_dir.is_dir():
        print(f"[error] CISBP pwms 目录不存在: {cisbp_dir}", file=sys.stderr)
        sys.exit(1)

    queries = _load_queries_from_task(task_dir, only_existing_png=not args.all_conv_filters)
    if not queries:
        print(f"[error] 未找到 conv_kernel_ppms.json 或 seqlet_clusters.npz: {task_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"[info] 加载 CISBP PWM: {cisbp_dir} ...")
    cisbp = load_all_cisbp_pwms(cisbp_dir)
    print(f"[info] CISBP 条数: {len(cisbp)}；模型 query 数: {len(queries)}")

    print("[info] 扫描比对（可能需数十秒）...")
    ranked = rank_queries_against_cisbp(queries, cisbp, top_k=args.top_k)

    out_dir = task_dir / "cisbp_compare"
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    cisbp_map = dict(cisbp)
    qmap = dict(queries)
    for block in ranked:
        qid = block["query_id"]
        top = block["top"][0]
        cid = top["cisbp_id"]
        off = top["offset_on_longer"]
        q_ppm = qmap[qid]
        mat = cisbp_map[cid]
        a, b = extract_aligned_segments(q_ppm, mat, off)
        safe = f"{qid}_vs_{cid}".replace("/", "_")
        png = fig_dir / f"{safe}.png"
        _plot_side_by_side(
            a,
            b,
            f"model {qid}",
            f"CISBP {cid}",
            png,
        )

    payload = {
        "task_dir": str(task_dir),
        "cisbp_pwms_dir": str(cisbp_dir),
        "n_cisbp_pwms": len(cisbp),
        "metric": "min sliding-window mean symmetric KL (lower=better)",
        "results": ranked,
    }
    save_compare_report_json(out_dir / "comparison_results.json", payload)

    html = _build_html(task_dir.name, ranked)
    (out_dir / "comparison_report.html").write_text(html, encoding="utf-8")
    print(f"[done] 写入 {out_dir / 'comparison_report.html'}")


if __name__ == "__main__":
    main()
