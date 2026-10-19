#!/usr/bin/env python3
"""
多源参考 Motif 与模型输出对比：CISBP-RNA、用户 MEME、DEWSeq/catRAPID 导出、可选 STREME（峰区序列）。

生成 task_dir/motif_benchmark/benchmark_full.json 与 benchmark_report.html。

参考库目录（可空）：
  reference_data/motif_libraries/meme_rna/*.meme
  reference_data/motif_libraries/dewseq/*.meme 或 *.txt（CISBP 制表符格式）

DEWSeq 本身不发布独立 PWM 包；论文使用 catRAPID 等外部注释，请从 R/文献导出后放入 dewseq/。

用法：
  cd 可解释性分析/Motif
  python3 scripts/build_motif_benchmark.py --task-dir outputs/AGGF1_HepG2_ENCSR543TPH
"""
from __future__ import annotations

import argparse
import glob
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

_MOTIF_ROOT = Path(__file__).resolve().parents[1]
if str(_MOTIF_ROOT) not in sys.path:
    sys.path.insert(0, str(_MOTIF_ROOT))

from motif_analysis.cisbp_compare import (  # noqa: E402
    extract_aligned_segments,
    load_all_cisbp_pwms,
    rank_queries_against_cisbp,
)
from motif_analysis.meme_pwm_io import load_cisbp_style_txt_dir, load_meme_dir, parse_meme_file  # noqa: E402


def _load_queries(task_dir: Path, only_png_conv: bool = True):
    import json as json_mod
    import re

    out: list[tuple[str, np.ndarray]] = []
    conv_json = task_dir / "conv_motifs" / "conv_kernel_ppms.json"
    logo_keys = None
    if only_png_conv:
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
            branches = json_mod.load(f)
        for b in branches:
            bi, ks = b["branch_index"], b["kernel_size"]
            for filt in b.get("filters", []):
                fi = filt["filter_index"]
                if logo_keys is not None and (bi, ks, fi) not in logo_keys:
                    continue
                out.append((f"conv_branch{bi}_k{ks}_f{fi}", np.asarray(filt["ppm"], dtype=np.float64)))
    for sub, name_fmt in (
        ("seqlet_motifs", "seqlet_motif_{}"),
        ("seqlet_motifs_legacy", "seqlet_legacy_{}"),
    ):
        seq_npz = task_dir / sub / "seqlet_clusters.npz"
        if seq_npz.is_file():
            z = np.load(seq_npz)
            ppms = z["ppms"]
            for i in range(ppms.shape[0]):
                out.append((name_fmt.format(i), ppms[i].astype(np.float64)))
    return out


def _plot_side_by_side(left, right, left_title, right_title, out_png: Path) -> None:
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


def _try_run_streme(task_dir: Path) -> Path | None:
    sh = task_dir / "meme_streme" / "run_streme.sh"
    if not sh.is_file():
        return None
    if not shutil.which("streme"):
        return None
    try:
        subprocess.run(["bash", str(sh)], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError:
        return None
    out_base = task_dir / "meme_streme" / "streme_output" / "streme_out"
    for pat in ("meme.txt", "streme.txt", "*.txt"):
        for p in glob.glob(str(out_base / pat)):
            if Path(p).stat().st_size > 50:
                return Path(p)
    return None


def _load_streme_motifs(streme_txt: Path) -> list[tuple[str, np.ndarray]]:
    return parse_meme_file(streme_txt)


def _dedupe_lib(entries: list[tuple[str, np.ndarray]]) -> list[tuple[str, np.ndarray]]:
    seen: dict[str, int] = {}
    out: list[tuple[str, np.ndarray]] = []
    for k, v in entries:
        if k not in seen:
            seen[k] = 0
            out.append((k, v))
        else:
            seen[k] += 1
            out.append((f"{k}__dup{seen[k]}", v))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-dir", type=Path, required=True)
    ap.add_argument("--cisbp-pwms", type=Path, default=_MOTIF_ROOT / "reference_data/CISBP-RNA/pwms")
    ap.add_argument("--meme-lib", type=Path, default=_MOTIF_ROOT / "reference_data/motif_libraries/meme_rna")
    ap.add_argument("--dewseq-lib", type=Path, default=_MOTIF_ROOT / "reference_data/motif_libraries/dewseq")
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument("--skip-streme", action="store_true")
    ap.add_argument(
        "--html-thumbnail-queries",
        type=str,
        default="seqlet_motif_",
        help="仅对此前缀的 query 在 HTML 中嵌入缩略图（JSON 仍含全部）；默认仅 seqlet",
    )
    args = ap.parse_args()

    task_dir = args.task_dir.resolve()
    queries = _load_queries(task_dir)
    if not queries:
        print("[error] 无 conv/seqlet query", file=sys.stderr)
        sys.exit(1)

    libraries: dict[str, dict] = {}
    cisbp_dir = args.cisbp_pwms.resolve()
    cisbp_lib_full: list[tuple[str, np.ndarray]] = []

    # 1) CISBP
    if cisbp_dir.is_dir():
        cisbp_lib_full = load_all_cisbp_pwms(cisbp_dir)
        libraries["CISBP-RNA"] = {
            "n_motifs": len(cisbp_lib_full),
            "results": rank_queries_against_cisbp(queries, cisbp_lib_full, top_k=args.top_k),
        }
    else:
        libraries["CISBP-RNA"] = {"error": f"missing {cisbp_dir}"}

    # 2) MEME 用户库
    meme_path = args.meme_lib.resolve()
    ml = _dedupe_lib(load_meme_dir(meme_path) + load_cisbp_style_txt_dir(meme_path))
    if ml:
        libraries["MEME_or_literature"] = {"n_motifs": len(ml), "results": rank_queries_against_cisbp(queries, ml, top_k=args.top_k)}
    else:
        libraries["MEME_or_literature"] = {"n_motifs": 0, "note": f"未在 {meme_path} 找到 .meme 或 .txt", "results": []}

    # 3) DEWSeq / catRAPID 导出目录
    dew_path = args.dewseq_lib.resolve()
    dl = _dedupe_lib(load_meme_dir(dew_path) + load_cisbp_style_txt_dir(dew_path))
    if dl:
        libraries["DEWSeq_catRAPID_export"] = {"n_motifs": len(dl), "results": rank_queries_against_cisbp(queries, dl, top_k=args.top_k)}
    else:
        libraries["DEWSeq_catRAPID_export"] = {
            "n_motifs": 0,
            "note": "DEWSeq 无独立全库 PWM；请将 Bioconductor/文献导出的 .meme 或 CISBP 风格 .txt 放入 dewseq/。",
            "results": [],
        }

    # 4) STREME（数据驱动 MEME）
    streme_motifs: list[tuple[str, np.ndarray]] = []
    streme_note = ""
    if not args.skip_streme:
        st_txt = _try_run_streme(task_dir)
        if st_txt:
            streme_motifs = _load_streme_motifs(st_txt)
            streme_note = str(st_txt)
        elif not shutil.which("streme"):
            streme_note = "streme 不在 PATH，已跳过（conda install -c bioconda meme）"
        else:
            streme_note = "STREME 运行失败或无输出文件"
    else:
        streme_note = "用户指定 --skip-streme"

    if streme_motifs:
        streme_motifs = _dedupe_lib(streme_motifs)
        libraries["STREME_on_peaks"] = {
            "n_motifs": len(streme_motifs),
            "source_file": streme_note,
            "results": rank_queries_against_cisbp(queries, streme_motifs, top_k=args.top_k),
        }
    else:
        libraries["STREME_on_peaks"] = {"n_motifs": 0, "note": streme_note, "results": []}

    out_dir = task_dir / "motif_benchmark"
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 为每个 query × 每个有结果的库画 top-1 并排图
    qmap = dict(queries)
    thumb_pref: str = args.html_thumbnail_queries or ""
    lib_maps: dict[str, dict] = {}
    for section, payload in libraries.items():
        res = payload.get("results") or []
        if not res:
            continue
        # 重建该库的 id->matrix（仅 top 命中需要）
        if section == "CISBP-RNA" and cisbp_lib_full:
            lib_maps[section] = dict(cisbp_lib_full)
        elif section == "MEME_or_literature":
            lib_maps[section] = dict(ml)
        elif section == "DEWSeq_catRAPID_export":
            lib_maps[section] = dict(dl)
        elif section == "STREME_on_peaks":
            lib_maps[section] = dict(streme_motifs)

    for section, res in libraries.items():
        if not res.get("results"):
            continue
        lmap = lib_maps.get(section, {})
        for block in res["results"]:
            qid = block["query_id"]
            if thumb_pref and not qid.startswith(thumb_pref):
                continue
            top = block["top"][0] if block.get("top") else None
            if not top:
                continue
            cid = top["cisbp_id"]
            off = top["offset_on_longer"]
            q_ppm = qmap[qid]
            mat = lmap.get(cid)
            if mat is None:
                continue
            a, b = extract_aligned_segments(q_ppm, mat, off)
            safe = f"{qid}__{section[:24]}__{cid}".replace("/", "_").replace(":", "_")
            _plot_side_by_side(a, b, f"model {qid}", f"{section} {cid}", fig_dir / f"{safe}.png")

    payload = {
        "task_dir": str(task_dir),
        "n_queries": len(queries),
        "libraries": libraries,
    }
    (out_dir / "benchmark_full.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    # HTML
    rows = []
    for qid, _ in queries:
        row = [f"<td>{qid}</td>"]
        for sec in ("CISBP-RNA", "MEME_or_literature", "DEWSeq_catRAPID_export", "STREME_on_peaks"):
            lib = libraries.get(sec, {})
            br = next((x for x in (lib.get("results") or []) if x["query_id"] == qid), None)
            if not br or not br.get("top"):
                row.append("<td>—</td>")
                continue
            t = br["top"][0]
            kl = t["mean_sym_kl"]
            cid = t["cisbp_id"]
            safe = f"{qid}__{sec[:24]}__{cid}".replace("/", "_").replace(":", "_")
            img = f"figures/{safe}.png"
            if thumb_pref and qid.startswith(thumb_pref) and (fig_dir / f"{safe}.png").is_file():
                cell = f"{cid}<br/>KL={kl:.4f}<br/><img src='{img}' width='380'/>"
            else:
                cell = f"{cid}<br/>KL={kl:.4f}"
            row.append(f"<td>{cell}</td>")
        rows.append("<tr>" + "".join(row) + "</tr>")

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"/><title>Motif 多源对比 {task_dir.name}</title>
<style>
body {{ font-family: sans-serif; margin: 16px; }}
table {{ border-collapse: collapse; }}
th, td {{ border: 1px solid #ccc; padding: 6px; vertical-align: top; }}
h2 {{ margin-top: 1.2em; }}
</style></head><body>
<h1>Motif 多源对比（CISBP / MEME库 / DEWSeq导出 / STREME）</h1>
<p>指标：滑动窗口最小对称 KL（越小越相似）。STREME 需本机安装 MEME Suite 并成功运行 meme_streme/run_streme.sh。</p>
<p>DEWSeq 列为用户放入 <code>reference_data/motif_libraries/dewseq/</code> 的 motif；若为空则仅显示说明。</p>
<table>
<tr><th>Query</th><th>CISBP-RNA</th><th>MEME/文献库</th><th>DEWSeq/catRAPID 导出</th><th>STREME（峰区）</th></tr>
{"".join(rows)}
</table>
<h2>库状态</h2>
<pre>{json.dumps({k: {x: v for x, v in val.items() if x != "results"} for k, val in libraries.items()}, indent=2, ensure_ascii=False)}</pre>
</body></html>"""
    (out_dir / "benchmark_report.html").write_text(html, encoding="utf-8")
    print(f"[done] {out_dir / 'benchmark_report.html'}")


if __name__ == "__main__":
    main()
