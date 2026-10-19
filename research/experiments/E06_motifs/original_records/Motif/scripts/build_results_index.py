#!/usr/bin/env python3
"""
生成 outputs/VISUALIZATION_INDEX.html：任务列表 + conv/seqlet 缩略图 + 搜索框。
在 Motif 根目录执行: python scripts/build_results_index.py
"""
from __future__ import annotations

import html
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _first_png(d: Path) -> str | None:
    if not d.is_dir():
        return None
    for p in sorted(d.glob("*.png")):
        return p.name
    return None


def _discover_run_dirs(out_dir: Path) -> list[str]:
    """存在 run_config.json 的实验目录（outputs/<TASK_NAME>/）。"""
    skip = {"batch_motif_logs", "fixtures", "e2e_full_run", "demo_run", "smoke_fix_v2"}
    names: list[str] = []
    for p in out_dir.iterdir():
        if not p.is_dir() or p.name.startswith(".") or p.name in skip:
            continue
        if (p / "run_config.json").is_file():
            names.append(p.name)
    return sorted(names)


def main():
    out_dir = ROOT / "outputs"
    tasks = _discover_run_dirs(out_dir)
    rows = []
    for name in tasks:
        rel = name  # folder name
        conv_dir = out_dir / rel / "conv_motifs"
        seq_dir = out_dir / rel / "seqlet_motifs"
        c1 = _first_png(conv_dir)
        s0 = "motif_0.png" if (seq_dir / "motif_0.png").is_file() else _first_png(seq_dir)
        conv_src = f"{rel}/conv_motifs/{c1}" if c1 else ""
        seq_src = f"{rel}/seqlet_motifs/{s0}" if s0 else ""
        display_name = rel.replace("task_", "", 1).replace("_chr22", "") if rel.startswith("task_") else rel
        rows.append(
            {
                "name": display_name,
                "folder": rel,
                "conv_img": conv_src,
                "seq_img": seq_src,
            }
        )

    # optional: merge batch summary
    summary_note = ""
    for p in sorted((out_dir / "batch_motif_logs").glob("batch_summary_*.json"), reverse=True):
        try:
            j = json.loads(p.read_text(encoding="utf-8"))
            summary_note = (
                f"汇总：n_ok={j.get('n_ok')}, n_fail={j.get('n_fail')}, "
                f"skipped_no_chr_h5={j.get('skipped_no_chr_h5', '—')} "
                f"（文件 {p.name}）"
            )
            break
        except Exception:
            pass

    parts = [
        "<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'/>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'/>",
        "<title>MARCP Motif 结果可视化索引</title>",
        "<style>",
        ":root { --bg:#0f1419; --card:#1a2332; --text:#e6edf3; --muted:#8b9cb3; --accent:#58a6ff; }",
        "body { font-family: 'Segoe UI', 'PingFang SC', sans-serif; background: var(--bg); color: var(--text); margin:0; padding:24px; line-height:1.5; }",
        "h1 { font-weight:600; font-size:1.35rem; margin-bottom:8px; }",
        ".sub { color: var(--muted); font-size:0.9rem; margin-bottom:20px; }",
        "input#q { width:min(480px,100%); padding:10px 14px; border-radius:8px; border:1px solid #30363d; background:#161b22; color:var(--text); font-size:15px; }",
        "table { border-collapse:collapse; width:100%; margin-top:16px; font-size:13px; }",
        "th, td { border-bottom:1px solid #30363d; padding:10px 8px; text-align:left; vertical-align:middle; }",
        "th { color: var(--muted); font-weight:600; }",
        "tr:hover td { background:#161b22; }",
        ".thumb { max-height:72px; max-width:220px; object-fit:contain; background:#0d1117; border-radius:4px; }",
        "a.path { color: var(--accent); text-decoration:none; font-size:12px; word-break:break-all; }",
        "a.path:hover { text-decoration:underline; }",
        ".hidden { display:none !important; }",
        "</style></head><body>",
        "<h1>MARCP Motif 批量结果可视化索引</h1>",
        f"<p class='sub'>共 <strong>{len(rows)}</strong> 个任务目录（含 run_config.json 的流水线输出）。{html.escape(summary_note)}</p>",
        "<p class='sub'>说明：图片为相对路径；请用浏览器「打开本地文件」打开本 HTML（位于 outputs 目录下）。详细解读见 docs/Motif结果解读完全指南.md</p>",
        "<label for='q'>筛选任务名：</label><br/><input id='q' type='search' placeholder='例如 AGGF1 或 HepG2' autocomplete='off'/>",
        "<table id='tbl'><thead><tr>",
        "<th>#</th><th>任务（RBP_细胞系_ENCODE）</th>",
        "<th>卷积 Motif（首张 logo）</th><th>Seqlet 聚类（motif_0 等）</th><th>目录</th>",
        "</tr></thead><tbody>",
    ]

    for i, r in enumerate(rows, 1):
        fn = r["folder"]
        task_label = html.escape(r["name"])
        conv_img = r["conv_img"]
        seq_img = r["seq_img"]
        ccell = f"<img class='thumb' src='{html.escape(conv_img)}' alt='conv' loading='lazy'/>" if conv_img else "—"
        scell = f"<img class='thumb' src='{html.escape(seq_img)}' alt='seqlet' loading='lazy'/>" if seq_img else "—"
        link = f"{fn}/"
        parts.append(
            f"<tr data-name='{html.escape(r['name'].lower())}'><td>{i}</td>"
            f"<td><strong>{task_label}</strong></td><td>{ccell}</td><td>{scell}</td>"
            f"<td><a class='path' href='{html.escape(link)}'>{html.escape(fn)}/</a></td></tr>"
        )

    parts.extend(
        [
            "</tbody></table>",
            "<script>",
            "const q=document.getElementById('q'); const tr=[...document.querySelectorAll('#tbl tbody tr')];",
            "function filt(){ const v=q.value.trim().toLowerCase(); tr.forEach(t=>{",
            "const n=t.getAttribute('data-name')||''; t.classList.toggle('hidden', v && !n.includes(v)); }); }",
            "q.addEventListener('input', filt);",
            "</script>",
            "</body></html>",
        ]
    )

    html_path = out_dir / "VISUALIZATION_INDEX.html"
    html_path.write_text("".join(parts), encoding="utf-8")
    print(f"Wrote {html_path} ({len(rows)} tasks)")


if __name__ == "__main__":
    main()
