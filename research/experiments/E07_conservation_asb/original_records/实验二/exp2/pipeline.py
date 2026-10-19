from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from exp2.data_bridge import (
    MetaTableCache,
    global_index_to_file_local,
    h5_concat_cumsum,
    load_x_phylop_row,
    parse_chr_from_5ch_path,
    read_test_bed_keys,
    resolve_5ch_paths,
    row_key,
)
from exp2.features import (
    MetaprofileAccumulator,
    flank_mean_outside_core,
    mean_abs_cwm_on_core,
    mean_ic_on_core,
    mean_phylop_slice,
    random_non_core_segment_mean_phylop,
)
from exp2.genome_coords import norm_chrom, seqlet_genomic_span, seqlet_overlaps_half_open_interval
from exp2.modisco_seqlets import load_all_seqlets, load_pattern_contrib_and_sequence, load_pattern_ranking
from exp2.paths import TEST_CHROMOSOMES, default_bed_path
from exp2.task_ids import parse_rbp_cell_line
from exp2.tomtom_status import best_tomtom_q_by_query, motif_match_status


def _ensure_motif_syspath() -> Path:
    from exp2.paths import motif_marcp_root

    root = motif_marcp_root()
    s = str(root)
    if s not in sys.path:
        sys.path.insert(0, s)
    return root


def _canonical_cwm(ppm_raw: np.ndarray, modisco_h5: Path) -> np.ndarray | None:
    _ensure_motif_syspath()
    try:
        from motif_analysis.channel_order import contrib_scores_to_canonical_acgt, read_contrib_channel_order

        order = read_contrib_channel_order(modisco_h5)
        return contrib_scores_to_canonical_acgt(np.asarray(ppm_raw, dtype=np.float64), order)
    except Exception:
        return None


def _canonical_ppm(seq_raw: np.ndarray, modisco_h5: Path) -> np.ndarray | None:
    _ensure_motif_syspath()
    try:
        from motif_analysis.channel_order import contrib_scores_to_canonical_acgt, read_contrib_channel_order

        order = read_contrib_channel_order(modisco_h5)
        return contrib_scores_to_canonical_acgt(np.asarray(seq_raw, dtype=np.float64), order)
    except Exception:
        return None


def run_task_exp2(
    task_id: str,
    motif_task_dir: Path,
    h5data: Path,
    *,
    use_bed: bool,
    flank_bp: int,
    meta_half: int,
    min_instances_per_motif: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    motif_task_dir = Path(motif_task_dir)
    h5data = Path(h5data)
    qc: dict[str, Any] = {
        "task_id": task_id,
        "status": "ok",
        "n_seqlets_raw": 0,
        "n_instances_out": 0,
        "n_skipped_example": 0,
        "n_skipped_chrom": 0,
        "n_skipped_bed": 0,
        "n_skipped_motif_count": 0,
        "message": "",
    }

    modisco_h5 = motif_task_dir / "tfmodisco" / "modisco_results.h5"
    pred_idx = motif_task_dir / "predictions" / "sample_indices.npy"
    att_npz = motif_task_dir / "attributions" / "ig_hyp_4ch.npz"
    rc_path = motif_task_dir / "run_config.json"

    if not modisco_h5.is_file():
        qc["status"] = "skip_no_modisco_h5"
        qc["message"] = str(modisco_h5)
        return _empty_result(qc, meta_half=int(meta_half))
    if not pred_idx.is_file() or not att_npz.is_file():
        qc["status"] = "skip_no_predictions_or_attr"
        return _empty_result(qc, meta_half=int(meta_half))

    sample_indices = np.load(pred_idx).astype(np.int64)
    z = np.load(att_npz)
    ohe_all = np.asarray(z["ohe"], dtype=np.float32)
    n_attr = int(ohe_all.shape[0])
    if sample_indices.shape[0] != n_attr:
        qc["status"] = "skip_indices_mismatch"
        qc["message"] = f"sample_indices {sample_indices.shape} vs ohe {ohe_all.shape}"
        return _empty_result(qc, meta_half=int(meta_half))

    h5_mode = "all_chr"
    single_h5: str | None = None
    if rc_path.is_file():
        try:
            rc = json.loads(rc_path.read_text(encoding="utf-8"))
            h5_mode = str(rc.get("h5_mode", "all_chr"))
            h = rc.get("h5")
            single_h5 = str(h) if h else None
        except (json.JSONDecodeError, OSError):
            pass

    paths = resolve_5ch_paths(h5data, task_id, h5_mode, single_h5)
    if not paths:
        qc["status"] = "skip_no_5ch_paths"
        return _empty_result(qc, meta_half=int(meta_half))

    cum, path_list = h5_concat_cumsum(paths)
    meta_cache = MetaTableCache(h5data, task_id)
    bed_path = default_bed_path(task_id)
    bed_keys_overlap = read_test_bed_keys(bed_path)

    patterns = load_pattern_ranking(motif_task_dir)
    if not patterns:
        _ensure_motif_syspath()
        from motif_analysis.modisco_h5_patterns import collect_pattern_records

        recs = collect_pattern_records(modisco_h5, include_neg=False, include_subpatterns=False)
        patterns = [
            {"pattern_id": r.pattern_id, "h5_group_path": r.h5_group_path, "track": r.track}
            for r in recs
        ]
    if not patterns:
        qc["status"] = "skip_no_patterns"
        return _empty_result(qc, meta_half=int(meta_half))

    qmap = best_tomtom_q_by_query(motif_task_dir)
    seqlets = load_all_seqlets(modisco_h5, patterns)
    qc["n_seqlets_raw"] = int(len(seqlets))

    rbp, cell_line = parse_rbp_cell_line(task_id)
    instance_rows: list[dict[str, Any]] = []

    motif_cwm: dict[str, np.ndarray] = {}
    motif_ppm: dict[str, np.ndarray] = {}
    motif_len: dict[str, int] = {}
    for row in patterns:
        pid = str(row.get("pattern_id", ""))
        gp = str(row.get("h5_group_path", ""))
        c_raw, s_raw = load_pattern_contrib_and_sequence(modisco_h5, gp)
        if c_raw is not None:
            cw = _canonical_cwm(c_raw, modisco_h5)
            if cw is not None:
                motif_cwm[pid] = cw
                motif_len[pid] = int(cw.shape[0])
        if s_raw is not None:
            pm = _canonical_ppm(s_raw, modisco_h5)
            if pm is not None:
                motif_ppm[pid] = pm
                if pid not in motif_len:
                    motif_len[pid] = int(pm.shape[0])

    for sl in seqlets:
        ex = int(sl.example_idx)
        if ex < 0 or ex >= n_attr:
            qc["n_skipped_example"] += 1
            continue
        gidx = int(sample_indices[ex])
        fi, loc = global_index_to_file_local(cum, gidx)
        chr5_path = path_list[fi]
        ch = parse_chr_from_5ch_path(chr5_path)
        if ch not in TEST_CHROMOSOMES:
            qc["n_skipped_chrom"] += 1
            continue

        srow = meta_cache.row_for_5ch_row(chr5_path, loc)
        if srow is None:
            qc["n_skipped_example"] += 1
            continue
        chrom = str(srow["chr"])
        cs = int(srow["chromStart"])
        ce = int(srow["chromEnd"])
        strand = str(srow["strand"])
        rk = row_key(chrom, cs, ce, strand)
        if use_bed and rk not in bed_keys_overlap:
            qc["n_skipped_bed"] += 1
            continue

        phylop = load_x_phylop_row(chr5_path, loc)
        if phylop is None:
            qc["n_skipped_example"] += 1
            continue

        L = int(phylop.shape[0])
        st, en = int(sl.start), int(sl.end)
        if en <= st or st < 0 or en > L:
            qc["n_skipped_example"] += 1
            continue

        core_mean = mean_phylop_slice(phylop, st, en)
        flank_mean = flank_mean_outside_core(phylop, st, en, flank_bp)
        null_mean = random_non_core_segment_mean_phylop(phylop, st, en, rng)

        g_lo, g_hi = seqlet_genomic_span(cs, ce, strand, st, en)
        overlaps = False
        for bk in bed_keys_overlap:
            bch, bs, be, bst = bk
            if norm_chrom(chrom) != bch or strand != bst:
                continue
            if seqlet_overlaps_half_open_interval(g_lo, g_hi, bs, be):
                overlaps = True
                break

        pid = sl.pattern_id
        Lp = motif_len.get(pid, en - st)
        cwm_c = mean_abs_cwm_on_core(motif_cwm.get(pid), st, en, Lp) if pid in motif_cwm else float("nan")
        ic_c = mean_ic_on_core(motif_ppm.get(pid), st, en, Lp) if pid in motif_ppm else float("nan")

        qv = qmap.get(pid)
        inst_id = f"{task_id}_{pid}_ex{ex}_i{sl.instance_index}"
        instance_rows.append(
            {
                "task_id": task_id,
                "rbp": rbp,
                "cell_line": cell_line,
                "motif_id": pid,
                "motif_match_status": motif_match_status(qv),
                "tomtom_min_q": qv,
                "instance_id": inst_id,
                "chr": norm_chrom(chrom),
                "chromStart": cs,
                "chromEnd": ce,
                "strand": strand,
                "overlaps_peak": bool(overlaps),
                "mean_phyloP_core": core_mean,
                "mean_phyloP_flank": flank_mean,
                "mean_phyloP_null_seg": null_mean,
                "cwm_core_score": cwm_c,
                "ppm_ic_core": ic_c,
                "seqlet_score": sl.seqlet_score,
                "seqlet_start": st,
                "seqlet_end": en,
                "is_revcomp": bool(sl.is_revcomp),
                "example_idx": ex,
                "h5_group_path": sl.h5_group_path,
            }
        )

    pdf = pd.DataFrame(instance_rows)
    n_before_motif = int(len(pdf))
    if not pdf.empty:
        cnt = pdf.groupby("motif_id").size()
        keep = cnt[cnt >= int(min_instances_per_motif)].index
        pdf = pdf[pdf["motif_id"].isin(keep)].reset_index(drop=True)
        qc["n_skipped_motif_count"] = int(n_before_motif - len(pdf))

    meta_acc = MetaprofileAccumulator(half=int(meta_half))
    position_acc: dict[str, dict[str, list[Any]]] = {}
    for _, row in pdf.iterrows():
        ex = int(row["example_idx"])
        gidx = int(sample_indices[ex])
        fi, loc = global_index_to_file_local(cum, gidx)
        chr5_path = path_list[fi]
        phylop = load_x_phylop_row(chr5_path, loc)
        if phylop is None:
            continue
        st = int(row["seqlet_start"])
        en = int(row["seqlet_end"])
        pid = str(row["motif_id"])
        meta_acc.add_instance(phylop, st, en)
        center = 0.5 * (st + en - 1)
        cwm_mat = motif_cwm.get(pid)
        if cwm_mat is not None and cwm_mat.shape[0] == (en - st):
            for j in range(st, en):
                pj = j - st
                rel = int(round(j - center))
                abs_cwm = float(np.sum(np.abs(cwm_mat[pj])))
                ppm_row = motif_ppm.get(pid)
                icp = float("nan")
                if ppm_row is not None and pj < ppm_row.shape[0]:
                    from exp2.features import ppm_information_content_row

                    pr = ppm_row[pj]
                    s = float(pr.sum())
                    prn = pr / s if s > 1e-12 else pr
                    icp = ppm_information_content_row(prn)
                phy_j = float(phylop[j])
                if pid not in position_acc:
                    position_acc[pid] = {"rels": [], "phy": [], "cwm": [], "ic": []}
                position_acc[pid]["rels"].append(rel)
                position_acc[pid]["phy"].append(phy_j)
                position_acc[pid]["cwm"].append(abs_cwm)
                position_acc[pid]["ic"].append(icp)

    position_rows = _aggregate_position_table(task_id, position_acc, pdf)

    qc["n_instances_out"] = int(len(pdf))
    return {
        "qc": qc,
        "instance_df": pdf,
        "position_df": position_rows,
        "meta_sum": meta_acc.sum_w.copy(),
        "meta_cnt": meta_acc.cnt_w.copy(),
        "meta_half": int(meta_half),
    }


def _aggregate_position_table(
    task_id: str,
    position_acc: dict[str, dict[str, list[Any]]],
    instance_df: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if instance_df.empty:
        return pd.DataFrame(columns=["task_id", "motif_id", "rel_pos", "mean_phyloP", "abs_cwm", "ppm_ic"])
    for pid, acc in position_acc.items():
        if not instance_df[instance_df["motif_id"] == pid].shape[0]:
            continue
        rels = np.asarray(acc["rels"], dtype=np.int64)
        phy = np.asarray(acc["phy"], dtype=np.float64)
        cwm = np.asarray(acc["cwm"], dtype=np.float64)
        icv = np.asarray(acc["ic"], dtype=np.float64)
        for r in sorted(set(rels.tolist())):
            m = rels == r
            rows.append(
                {
                    "task_id": task_id,
                    "motif_id": pid,
                    "rel_pos": int(r),
                    "mean_phyloP": float(np.nanmean(phy[m])),
                    "abs_cwm": float(np.nanmean(cwm[m])),
                    "ppm_ic": float(np.nanmean(icv[m])),
                }
            )
    return pd.DataFrame(rows)


def _empty_result(qc: dict[str, Any], *, meta_half: int = 25) -> dict[str, Any]:
    w = 2 * int(meta_half) + 1
    return {
        "qc": qc,
        "instance_df": pd.DataFrame(),
        "position_df": pd.DataFrame(),
        "meta_sum": np.zeros(w, dtype=np.float64),
        "meta_cnt": np.zeros(w, dtype=np.float64),
        "meta_half": int(meta_half),
    }
