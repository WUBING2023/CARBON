#!/usr/bin/env python3
"""
为单个任务目录生成综合可视化 HTML（排名 / TF-MoDISco / 流水线 ISM / 三项验证）。

在 Motif 根目录执行：
  python3 scripts/build_task_visualization.py --task-dir outputs/AGGF1_HepG2_ENCSR543TPH

输出：<task-dir>/visualization/summary.html（用浏览器本地打开）。
"""
from __future__ import annotations

import argparse
import html as html_mod
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _count_meme_motifs(meme_path: Path) -> int:
    if not meme_path.is_file():
        return 0
    return len(re.findall(r"^MOTIF\s+", meme_path.read_text(encoding="utf-8"), re.MULTILINE))


def _count_cisbp_txt(dir_path: Path) -> int:
    if not dir_path.is_dir():
        return 0
    return len(list(dir_path.glob("*.txt")))


def _count_pattern_ranking_entries(task_dir: Path) -> int:
    p = task_dir / "tfmodisco_export" / "pattern_ranking.json"
    if not p.is_file():
        return 0
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return 0
    arr = data.get("patterns")
    return len(arr) if isinstance(arr, list) else 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-dir", type=Path, required=True)
    args = ap.parse_args()
    task_dir = args.task_dir.resolve()
    if not task_dir.is_dir():
        print(f"[error] not a directory: {task_dir}", file=sys.stderr)
        raise SystemExit(1)

    out_html = task_dir / "visualization" / "summary.html"
    out_html.parent.mkdir(parents=True, exist_ok=True)

    ranking = {}
    p_rank = task_dir / "predictions" / "ranking_meta.json"
    if p_rank.is_file():
        ranking = json.loads(p_rank.read_text(encoding="utf-8"))

    ism_pipe = {}
    p_ism = task_dir / "ism_check" / "ism_report.json"
    if p_ism.is_file():
        ism_pipe = json.loads(p_ism.read_text(encoding="utf-8"))

    val_path = task_dir / "threefold_validation" / "validation_report.json"
    val = {}
    if val_path.is_file():
        val = json.loads(val_path.read_text(encoding="utf-8"))
    sec = val.get("sections", {})
    cisbp_sec = sec.get("cisbp_alignment", {})
    ranked = cisbp_sec.get("ranked", [])
    tt = sec.get("tomtom", {})
    ism3 = sec.get("ism", {})

    n_meme = _count_meme_motifs(task_dir / "tfmodisco_export" / "queries.meme")
    n_txt = _count_cisbp_txt(task_dir / "tfmodisco_export" / "cisbp_style_txt")
    n_rank = _count_pattern_ranking_entries(task_dir)
    n_cisbp_tgt = _count_meme_motifs(task_dir / "tomtom_cisbp" / "cisbp_targets.meme")

    task_name = task_dir.name
    summary_ism = ism_pipe.get("summary", {})

    # CISBP 表：每条 query 的 top-1
    rows_cisbp = []
    for block in ranked[:12]:
        qid = block.get("query_id", "")
        top = (block.get("top") or [{}])[0]
        rows_cisbp.append(
            (
                qid[:56] + ("…" if len(str(qid)) > 56 else ""),
                top.get("cisbp_id", "—"),
                f"{float(top.get('mean_sym_kl', 0)):.4f}" if top else "—",
            )
        )

    # 流水线 ISM：每位点 |Δlogit| 最大
    pos_rows = []
    bases = ["A", "C", "G", "U"]
    for item in ism_pipe.get("per_position", []):
        pos = item.get("pos")
        d = item.get("ism_delta_sumlogits") or []
        if not d:
            continue
        arr = [float(x) for x in d[:4]]
        j = max(range(4), key=lambda k: abs(arr[k]))
        pos_rows.append((pos, arr, j, bases[j], abs(arr[j])))

    max_eff = max((p[4] for p in pos_rows), default=1.0)

    def esc(s: object) -> str:
        return html_mod.escape(str(s))

    # --- HTML ---
    parts = [
        "<!DOCTYPE html>",
        '<html lang="zh-CN">',
        "<head>",
        '<meta charset="utf-8"/>',
        '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
        f"<title>{esc(task_name)} · Motif 可视化总览</title>",
        '<link rel="preconnect" href="https://fonts.googleapis.com"/>',
        '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin/>',
        '<link href="https://fonts.googleapis.com/css2?family=DM+Sans:ital,opsz,wght@0,9..40,400;0,9..40,600;1,9..40,400&family=Fraunces:opsz,wght@9..144,600&display=swap" rel="stylesheet"/>',
        "<style>",
        ":root{--bg:#0f1419;--card:#1a222d;--text:#e8eaed;--muted:#8b949e;--accent:#3ecf8e;--warn:#f0b429;--line:#2d3848;}",
        "body{margin:0;font-family:'DM Sans',system-ui,sans-serif;background:var(--bg);color:var(--text);min-height:100vh;}",
        "body::before{content:'';position:fixed;inset:0;background:radial-gradient(ellipse 80% 50% at 20% -10%,rgba(62,207,142,.12),transparent),radial-gradient(ellipse 60% 40% at 100% 0%,rgba(240,180,41,.06),transparent);pointer-events:none;z-index:0;}",
        ".wrap{max-width:1100px;margin:0 auto;padding:2.5rem 1.5rem 4rem;position:relative;z-index:1;}",
        "h1{font-family:'Fraunces',Georgia,serif;font-weight:600;font-size:1.85rem;margin:0 0 .35rem;letter-spacing:-.02em;}",
        ".sub{color:var(--muted);font-size:.95rem;margin-bottom:2rem;}",
        ".grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(200px,1fr));gap:1rem;margin-bottom:2rem;}",
        ".card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:1rem 1.1rem;}",
        ".card label{display:block;font-size:.72rem;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin-bottom:.35rem;}",
        ".card .v{font-size:1.35rem;font-weight:600;font-variant-numeric:tabular-nums;}",
        "h2{font-family:'Fraunces',Georgia,serif;font-size:1.15rem;margin:2rem 0 1rem;padding-bottom:.5rem;border-bottom:1px solid var(--line);}",
        "table{width:100%;border-collapse:collapse;font-size:.85rem;}",
        "th,td{padding:.55rem .65rem;text-align:left;border-bottom:1px solid var(--line);}",
        "th{color:var(--muted);font-weight:600;font-size:.72rem;text-transform:uppercase;letter-spacing:.04em;}",
        ".bar-wrap{display:flex;align-items:flex-end;gap:.4rem;flex-wrap:wrap;margin-top:.75rem;}",
        ".bar{flex:1;min-width:36px;max-width:64px;background:var(--line);border-radius:4px 4px 0 0;position:relative;}",
        ".bar i{display:block;background:linear-gradient(180deg,var(--accent),#2a8f63);border-radius:4px 4px 0 0;min-height:4px;}",
        ".bar span{display:block;text-align:center;font-size:.65rem;color:var(--muted);margin-top:.35rem;}",
        ".imgbox{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:1rem;text-align:center;}",
        ".imgbox img{max-width:100%;height:auto;border-radius:8px;}",
        "iframe{width:100%;min-height:420px;border:1px solid var(--line);border-radius:12px;background:#fff;}",
        ".links a{color:var(--accent);text-decoration:none;}",
        ".links a:hover{text-decoration:underline;}",
        ".note{font-size:.8rem;color:var(--muted);margin-top:.5rem;}",
        ".pill{display:inline-block;padding:.2rem .55rem;border-radius:999px;background:rgba(62,207,142,.15);color:var(--accent);font-size:.75rem;margin-left:.5rem;}",
        "</style>",
        "</head>",
        "<body>",
        '<div class="wrap">',
        f"<h1>{esc(task_name)}</h1>",
        '<p class="sub">全染色体 Top-K · TF-MoDISco · 流水线 ISM · 三项验证（对称 KL / 剖面 ISM）<span class="pill">本地汇总页</span></p>',
        '<div class="grid">',
    ]

    def card(label: str, value: str) -> str:
        return f'<div class="card"><label>{esc(label)}</label><div class="v">{esc(value)}</div></div>'

    parts.append(card("全库扫描样本数", f"{ranking.get('n_total', '—'):,}" if isinstance(ranking.get("n_total"), int) else str(ranking.get("n_total", "—"))))
    parts.append(card("Top-K", str(ranking.get("select_top", "—"))))
    parts.append(card("已选入榜", str(ranking.get("n_selected", "—"))))
    parts.append(card("排序指标", str(ranking.get("rank_by", "—"))))
    parts.append(card("选中集均分", f"{float(ranking.get('score_mean_selected', 0)):.4f}" if ranking.get("score_mean_selected") is not None else "—"))
    parts.append(card("导出 MEME motif 数", str(n_meme or "—")))
    if n_rank > 0:
        parts.append(card("pattern 排名条目（seqlet）", str(n_rank)))
    else:
        parts.append(card("CISBP 风格 txt 数", str(n_txt or "—")))
    parts.append(card("Tomtom 靶库 MEME 条数", str(n_cisbp_tgt or "—")))

    parts.append("</div>")

    # TF-MoDISco report iframe
    motifs_rel = "../tfmodisco/report_simple/motifs.html"
    p_motifs = task_dir / "tfmodisco" / "report_simple" / "motifs.html"
    parts.append("<h2>TF-MoDISco（report_simple）</h2>")
    if p_motifs.is_file():
        parts.append(f'<iframe src="{esc(motifs_rel)}" title="modisco motifs"></iframe>')
        parts.append('<p class="note">若 iframe 空白，请确认用「文件」方式打开时允许加载子页面；或直接用浏览器打开 tfmodisco/report_simple/motifs.html。</p>')
    else:
        parts.append("<p class=\"note\">未找到 tfmodisco/report_simple/motifs.html。</p>")

    # Pipeline ISM
    parts.append("<h2>流水线内 ISM（抽样位点 · 与 IG argmax 一致率）</h2>")
    parts.append('<div class="grid">')
    parts.append(card("抽样位点数", str(summary_ism.get("n", "—"))))
    parts.append(card("argmax 一致数", str(summary_ism.get("argmax_match", "—"))))
    parts.append(card("一致率", f"{float(summary_ism.get('rate', 0))*100:.1f}%" if summary_ism.get("rate") is not None else "—"))
    parts.append("</div>")
    if pos_rows:
        parts.append('<p class="note">每位点单碱基替换 |Δ sum logits| 最大者（柱高 ∝ 效应幅度）</p>')
        parts.append('<div class="bar-wrap">')
        for pos, arr, j, b, eff in sorted(pos_rows, key=lambda x: x[0]):
            h_pct = max(8, int(100 * eff / max_eff))
            parts.append(f'<div class="bar"><i style="height:{h_pct}px"></i><span>pos {esc(pos)}<br/>→{esc(b)}</span></div>')
        parts.append("</div>")

    # Threefold ISM image
    parts.append("<h2>三项验证 · 剖面 ISM（motif 区间 shuffle）</h2>")
    png_rel = "../threefold_validation/ism_profile_compare.png"
    p_png = task_dir / "threefold_validation" / "ism_profile_compare.png"
    if p_png.is_file():
        parts.append('<div class="imgbox">')
        parts.append(f'<img src="{esc(png_rel)}" alt="ISM profile compare"/>')
        parts.append("</div>")
    if ism3 and not ism3.get("skipped"):
        parts.append('<div class="grid" style="margin-top:1rem;">')
        parts.append(card("query motif", str(ism3.get("query_motif_id", "—"))[:40] + "…"))
        parts.append(card("样本 index", str(ism3.get("sample_index", "—"))))
        parts.append(card("peak 下降比例", f"{float(ism3.get('peak_drop_ratio', 0))*100:.1f}%"))
        parts.append(card("Δ sum binding", f"{float(ism3.get('delta_sum', 0)):.3f}"))
        parts.append("</div>")

    # CISBP table
    parts.append("<h2>CISBP-RNA 对称 KL（Top-1，粗筛）</h2>")
    if rows_cisbp:
        parts.append("<table><thead><tr><th>Query（modisco 解析名）</th><th>Top-1 CISBP ID</th><th>mean sym-KL</th></tr></thead><tbody>")
        for q, cid, kl in rows_cisbp:
            parts.append(f"<tr><td>{esc(q)}</td><td>{esc(cid)}</td><td>{esc(kl)}</td></tr>")
        parts.append("</tbody></table>")
    else:
        parts.append('<p class="note">无 threefold_validation/validation_report.json 或 ranked 为空。</p>')

    # Tomtom status
    parts.append("<h2>Tomtom（MEME Suite）</h2>")
    if isinstance(tt, dict) and tt.get("skipped"):
        parts.append(f'<p class="note">当前状态：{esc(tt.get("skipped", tt))}。安装 tomtom 并运行 scripts/tomtom_modisco_vs_cisbp.py 后可在此目录生成 tomtom_out/。</p>')
    elif isinstance(tt, dict) and tt.get("tomtom_dir"):
        parts.append(f'<p class="links">结果目录：<a href="../threefold_validation/tomtom/tomtom_out/tomtom.tsv">tomtom.tsv</a>（若存在）</p>')

    # Links
    parts.append("<h2>原始文件与报告</h2>")
    parts.append('<ul class="links">')
    parts.append(f'<li><a href="../threefold_validation/validation_report.html">三项验证完整 HTML</a></li>')
    parts.append(f'<li><a href="../tfmodisco_export/pattern_ranking.json">tfmodisco_export/pattern_ranking.json</a>（按 seqlet 排序）</li>')
    parts.append(f'<li><a href="../tfmodisco_export/top_patterns_for_benchmark.meme">tfmodisco_export/top_patterns_for_benchmark.meme</a></li>')
    parts.append(f'<li><a href="../tfmodisco_export/queries.meme">tfmodisco_export/queries.meme</a>（IC Span 截断版，供 Tomtom；全长见 top_patterns_for_benchmark.meme）</li>')
    parts.append(f'<li><a href="../tomtom_cisbp/cisbp_targets.meme">tomtom_cisbp/cisbp_targets.meme</a>（全库靶）</li>')
    parts.append(f'<li><a href="../predictions/ranking_meta.json">predictions/ranking_meta.json</a></li>')
    parts.append("</ul>")

    parts.append('<p class="note">生成：scripts/build_task_visualization.py · 打开本文件请使用「相对于本 HTML 的路径」，建议将整个任务文件夹置于本地后双击 summary.html。</p>')
    parts.append("</div></body></html>")

    out_html.write_text("\n".join(parts), encoding="utf-8")
    print(f"[ok] {out_html}")


if __name__ == "__main__":
    main()
