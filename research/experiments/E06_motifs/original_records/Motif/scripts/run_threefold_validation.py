#!/usr/bin/env python3
"""
三项验证：① Tomtom（或 CISBP 全库对齐作备选）② ISM 结合图谱 ③ Motif 富集（Mann–Whitney）

在 Motif 根目录执行：
  export PYTHONPATH="/data/CARBON_PROJECT/MARCP:${PYTHONPATH}"
  python3 scripts/run_threefold_validation.py --task-dir outputs/AGGF1_HepG2_ENCSR543TPH

依赖：numpy scipy torch；可选 MEME suite 的 `tomtom`（否则跳过 Tomtom 并提示）。
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from motif_analysis.cisbp_compare import load_all_cisbp_pwms, rank_queries_against_cisbp
from motif_analysis.meme_write import write_meme_file
from motif_analysis.motif_scan import best_hit_window, max_loglik_scores
from motif_analysis.query_motifs import load_all_query_motifs
from motif_analysis.validation_ism import run_ism_profile_compare


def _load_run_config(task_dir: Path) -> dict:
    p = task_dir / "run_config.json"
    if not p.is_file():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def _task_name_for_weights(task_dir: Path, cfg: dict) -> str:
    if cfg.get("task"):
        return str(cfg["task"])
    name = task_dir.name
    if name.startswith("task_"):
        name = name[5:]
    for marker in ("_chr22", "_chrX"):
        if marker in name:
            return name.split(marker)[0]
    return name


def _build_marcp(cfg: dict, ckpt: str, device):
    import torch
    from motif_analysis.paths import ensure_marcp_on_path

    ensure_marcp_on_path()
    from MARCP.models.marcp_model import MARCP

    eff_ch = 4 if cfg.get("ablation_no_conservative") else int(cfg.get("input_channels", 5))
    ab_kw = {
        "ablation_no_se": cfg.get("ablation_no_se", False),
        "ablation_no_multiscale": cfg.get("ablation_no_multiscale", False),
        "ablation_no_gru": cfg.get("ablation_no_gru", False),
        "ablation_lstm_instead_of_gru": cfg.get("ablation_lstm_instead_of_gru", False),
    }
    m = MARCP(dropOut=0.0, control=False, in_channels=eff_ch, **ab_kw).to(device)
    try:
        state = torch.load(ckpt, map_location=device, weights_only=True)
    except TypeError:
        state = torch.load(ckpt, map_location=device)
    m.load_state_dict(state, strict=False)
    m.eval()
    return m


def _build_cisbp_queries_ranking(
    task_dir: Path,
    cisbp_dir: Path,
    top_k: int,
    cisbp_max: int,
) -> tuple[list[dict], list[tuple[str, np.ndarray]]]:
    queries = load_all_query_motifs(task_dir)
    if not queries:
        return [], []
    cisbp = load_all_cisbp_pwms(cisbp_dir)
    if cisbp_max > 0 and len(cisbp) > cisbp_max:
        cisbp = cisbp[:cisbp_max]
    ranked = rank_queries_against_cisbp(queries, cisbp, top_k=top_k)
    return ranked, queries


def _run_tomtom(
    queries: list[tuple[str, np.ndarray]],
    cisbp_dir: Path,
    out_dir: Path,
    max_targets: int,
) -> dict | None:
    exe = shutil.which("tomtom")
    if not exe or not queries:
        return None
    out_dir.mkdir(parents=True, exist_ok=True)
    q_meme = out_dir / "queries.meme"
    write_meme_file(q_meme, queries)
    cisbp = load_all_cisbp_pwms(cisbp_dir)[:max_targets]
    if not cisbp:
        return {"error": "no CISBP pwms"}
    write_meme_file(out_dir / "cisbp_subset.meme", cisbp)
    t_out = out_dir / "tomtom_out"
    t_out.mkdir(exist_ok=True)
    cmd = [
        exe,
        "-no-ssc",
        "-dist",
        "pearson",
        "-min-overlap",
        "4",
        str(q_meme),
        str(out_dir / "cisbp_subset.meme"),
        "-o",
        str(t_out),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return {"error": r.stderr or r.stdout or "tomtom failed", "returncode": r.returncode}
    tsv = t_out / "tomtom.tsv"
    rows = []
    if tsv.is_file():
        with open(tsv, encoding="utf-8", errors="replace") as f:
            rd = csv.DictReader(f, delimiter="\t")
            for row in rd:
                rows.append(row)
    return {"tomtom_dir": str(t_out), "n_rows": len(rows), "rows": rows[:200]}


def _enrichment_one_motif(
    ohe_fg: np.ndarray,
    ohe_bg: np.ndarray,
    ppm: np.ndarray,
) -> dict:
    from scipy.stats import mannwhitneyu

    s_fg = max_loglik_scores(ohe_fg, ppm)
    s_bg = max_loglik_scores(ohe_bg, ppm)
    try:
        stat, pval = mannwhitneyu(s_fg, s_bg, alternative="greater")
    except Exception as e:
        return {"error": str(e)}
    return {
        "mean_fg": float(np.mean(s_fg)),
        "mean_bg": float(np.mean(s_bg)),
        "mannwhitney_pvalue": float(pval),
        "n_fg": int(s_fg.size),
        "n_bg": int(s_bg.size),
    }


def _shuffle_bg(ohe: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """每条序列独立打乱碱基顺序（保长度与组成近似）。 """
    n, L, c = ohe.shape
    out = ohe.copy()
    for i in range(n):
        perm = rng.permutation(L)
        out[i] = out[i, perm]
    return out


def _pick_sequence_for_ism(
    ohe: np.ndarray,
    pred_intensity: np.ndarray | None,
    gxi_abs: np.ndarray | None,
    ppm: np.ndarray,
) -> tuple[int, int, int]:
    """返回 (sample_idx, motif_start, motif_end)。"""
    n = ohe.shape[0]
    if pred_intensity is not None:
        best_i = int(np.argmax(np.max(pred_intensity, axis=1)))
    elif gxi_abs is not None:
        best_i = int(np.argmax(np.max(np.mean(gxi_abs, axis=2), axis=1)))
    else:
        best_i = 0
    sc, st = best_hit_window(ohe[best_i], ppm)
    w = ppm.shape[0]
    return best_i, st, st + w


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-dir", type=Path, required=True)
    ap.add_argument(
        "--cisbp-pwms",
        type=Path,
        default=ROOT / "reference_data" / "CISBP-RNA" / "pwms",
        help="CISBP-RNA pwm .txt 目录",
    )
    ap.add_argument("--cisbp-top-k", type=int, default=5)
    ap.add_argument(
        "--cisbp-max-align",
        type=int,
        default=1200,
        help="与 CISBP 全库对齐时最多使用条目前 N 条（过大会很慢）",
    )
    ap.add_argument("--tomtom-max-targets", type=int, default=1500)
    ap.add_argument("--enrichment-motifs", type=int, default=10, help="对前 N 个 query 做富集")
    ap.add_argument("--skip-tomtom", action="store_true")
    ap.add_argument("--skip-ism", action="store_true")
    ap.add_argument("--skip-enrichment", action="store_true")
    args = ap.parse_args()

    task_dir = args.task_dir.resolve()
    out_dir = task_dir / "threefold_validation"
    out_dir.mkdir(parents=True, exist_ok=True)

    report: dict = {"task_dir": str(task_dir), "sections": {}}

    # --- 1) 数据库比对 ---
    ranked, queries = _build_cisbp_queries_ranking(
        task_dir, args.cisbp_pwms, args.cisbp_top_k, args.cisbp_max_align
    )
    report["sections"]["cisbp_alignment"] = {
        "n_queries": len(queries),
        "ranked": ranked,
        "note": "对称 KL 滑动对齐（与 Tomtom 的 E-value 不同，可作粗筛）。",
    }

    if not args.skip_tomtom and queries:
        tt = _run_tomtom(
            queries[:50],
            args.cisbp_pwms,
            out_dir / "tomtom",
            args.tomtom_max_targets,
        )
        report["sections"]["tomtom"] = tt or {"skipped": "tomtom not in PATH or failed"}
    else:
        report["sections"]["tomtom"] = {"skipped": True}

    # --- 2) ISM ---
    cfg = _load_run_config(task_dir)
    att = task_dir / "attributions" / "ig_hyp_4ch.npz"
    pred_npz = task_dir / "predictions" / "pred_profile.npz"
    if not args.skip_ism and att.is_file() and queries:
        import torch

        z = np.load(att)
        ohe = z["ohe"].astype(np.float32)
        bind = None
        gxa = None
        if pred_npz.is_file():
            pr = np.load(pred_npz)
            bind = pr["binding_intensity"] if "binding_intensity" in pr.files else pr["softmax_profile"]
        elif "gxi_abs" in z.files:
            gxa = z["gxi_abs"].astype(np.float32)
        qid, ppm = queries[0]
        si, st, en = _pick_sequence_for_ism(ohe, bind, gxa, ppm)
        x_one = ohe[si : si + 1]
        if x_one.shape[-1] == 4 and not cfg.get("ablation_no_conservative", False):
            xc = np.concatenate([x_one, np.zeros((1, x_one.shape[1], 1), dtype=np.float32)], axis=-1)
        else:
            xc = x_one
        ckpt = cfg.get("checkpoint")
        if not ckpt:
            from motif_analysis.user_paths import baseline_checkpoint

            ckpt = str(baseline_checkpoint(_task_name_for_weights(task_dir, cfg)))
        if not Path(ckpt).is_file():
            report["sections"]["ism"] = {"error": f"checkpoint not found: {ckpt}"}
        else:
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            m = _build_marcp(cfg, ckpt, device)
            ism_json = run_ism_profile_compare(
                m,
                xc,
                st,
                en,
                device,
                out_dir / "ism_profile_compare.png",
                mode="single_point",
                ppm=ppm,
            )
            ism_json["query_motif_id"] = qid
            ism_json["sample_index"] = si
            report["sections"]["ism"] = ism_json
    else:
        report["sections"]["ism"] = {
            "skipped": True,
            "reason": "missing attributions/predictions/queries or --skip-ism",
        }

    # --- 3) 富集 ---
    if not args.skip_enrichment and att.is_file() and queries:
        z = np.load(att)
        ohe = z["ohe"].astype(np.float32)
        rng = np.random.default_rng(0)
        ohe_bg = _shuffle_bg(ohe, rng)
        enr = []
        for qid, ppm in queries[: args.enrichment_motifs]:
            r = _enrichment_one_motif(ohe, ohe_bg, ppm)
            r["query_id"] = qid
            enr.append(r)
        report["sections"]["enrichment"] = {
            "foreground": "attributions ohe (全部归因样本)",
            "background": "每条序列位置随机打乱",
            "per_motif": enr,
            "criterion": "Mann–Whitney U，alternative=greater；p 小且 mean_fg>mean_bg 表示富集。",
        }
    else:
        report["sections"]["enrichment"] = {"skipped": True}

    json_path = out_dir / "validation_report.json"
    json_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    # HTML
    html = _render_html(report)
    html_path = out_dir / "validation_report.html"
    html_path.write_text(html, encoding="utf-8")
    print(f"[ok] {json_path}")
    print(f"[ok] {html_path}")


def _render_html(report: dict) -> str:
    import html as html_mod

    parts = [
        "<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'/>",
        "<title>三项验证报告</title>",
        "<style>body{font-family:system-ui,sans-serif;max-width:960px;margin:24px;line-height:1.5;}",
        "table{border-collapse:collapse;width:100%;} td,th{border:1px solid #ccc;padding:6px;font-size:13px;}",
        "h2{border-bottom:1px solid #eee;} .note{color:#555;font-size:14px;}</style></head><body>",
        "<h1>Motif 三项验证报告</h1>",
        f"<p class='note'>任务目录：{html_mod.escape(report.get('task_dir',''))}</p>",
    ]
    sec = report.get("sections", {})
    parts.append("<h2>1. 外部数据库比对（CISBP + 可选 Tomtom）</h2>")
    parts.append("<p class='note'>Tomtom 需安装 MEME Suite；若跳过则仅展示 CISBP 滑动对称 KL 对齐。</p>")
    ca = sec.get("cisbp_alignment", {})
    parts.append(f"<p>Query 数量：{ca.get('n_queries', 0)}</p>")
    parts.append("<pre style='white-space:pre-wrap;background:#f6f8fa;padding:12px;'>")
    parts.append(html_mod.escape(json.dumps(ca.get("ranked", [])[:3], indent=2, ensure_ascii=False)[:8000]))
    parts.append("</pre>")
    tt = sec.get("tomtom")
    if tt and not tt.get("skipped"):
        parts.append("<h3>Tomtom</h3><pre>")
        parts.append(html_mod.escape(json.dumps(tt, indent=2, ensure_ascii=False)[:12000]))
        parts.append("</pre>")
    parts.append("<h2>2. ISM（突变前后 binding profile）</h2>")
    parts.append("<pre>")
    parts.append(html_mod.escape(json.dumps(sec.get("ism", {}), indent=2, ensure_ascii=False)))
    parts.append("</pre>")
    if sec.get("ism", {}).get("png"):
        parts.append(f"<p><img src='ism_profile_compare.png' style='max-width:100%;' alt='ism'/></p>")
    parts.append("<h2>3. Motif 富集（Mann–Whitney）</h2>")
    parts.append("<pre>")
    parts.append(html_mod.escape(json.dumps(sec.get("enrichment", {}), indent=2, ensure_ascii=False)))
    parts.append("</pre>")
    parts.append("</body></html>")
    return "".join(parts)


if __name__ == "__main__":
    main()
