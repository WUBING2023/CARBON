#!/usr/bin/env python3
"""1) 每个蛋白质一张对比图；2) 全局热图/条形汇总（matplotlib，无 seaborn 依赖）。

ordered_heatmap_variants/（默认生成，可用 --no-heatmap-variants 关闭）：
variant_A_*（cividis+网格）、variant_B_*（环状热图+环内聚类树+叶端蛋白名）、merge_M1/M2/M3（合并版式）；
行序由 Mean_PCC 与 Spearman 矩阵列拼接后统一聚类得到。
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import leaves_list, linkage, to_tree


def _midpoint_shortest_arc(th1: float, th2: float) -> float:
    """两方向角在圆上的最短弧中点（弧度）。"""
    d = (th2 - th1 + np.pi) % (2 * np.pi) - np.pi
    return th1 + d / 2.0


def _radial_hierarchy_tree_polar(
    ax,
    Z,
    *,
    r_inner: float,
    r_outer: float,
    color: str = "#2a2a2a",
    lw: float = 0.55,
) -> None:
    """
    使用 scipy to_tree 在极坐标下画径向聚类树（非直角 dendrogram 映射），减少枝干交叉重叠。
    叶节点等角分布（与 leaves_list 顺序一致）；合并点角度为两子结点在圆上最短弧中点；
    半径由 linkage 的 merge 高度归一化到 [r_inner, r_outer]。
    """
    n = Z.shape[0] + 1
    if n < 2:
        return
    root = to_tree(Z, rd=False)
    leaves_order = leaves_list(Z)
    theta_by_obs = {int(leaves_order[i]): 2.0 * np.pi * (i + 0.5) / n for i in range(n)}
    hmax = float(np.max(Z[:, 2]))
    if hmax <= 1e-15:
        hmax = 1.0

    def visit(node):
        if node.is_leaf():
            oid = int(node.get_id())
            th = theta_by_obs[oid]
            return float(th), float(r_outer)
        left = node.get_left()
        right = node.get_right()
        tl, rl = visit(left)
        tr, rr = visit(right)
        h = float(node.dist) if getattr(node, "dist", None) is not None else 0.0
        h = min(max(h, 0.0), hmax)
        rm = r_inner + (1.0 - h / hmax) * (r_outer - r_inner)
        thm = _midpoint_shortest_arc(tl, tr)
        ax.plot([thm, tl], [rm, rl], color=color, lw=lw, zorder=3, clip_on=False, alpha=0.92)
        ax.plot([thm, tr], [rm, rr], color=color, lw=lw, zorder=3, clip_on=False, alpha=0.92)
        return float(thm), float(rm)

    th_root, r_root = visit(root)
    r_stub = max(r_inner * 0.45, 0.02)
    ax.plot(
        [th_root, th_root],
        [r_stub, r_root],
        color=color,
        lw=lw,
        zorder=3,
        clip_on=False,
        alpha=0.85,
    )


def _polar_ring_heatmap_with_tree(
    mat: pd.DataFrame,
    Z,
    title: str,
    vmax: float,
    *,
    fname: str,
    out_dir: Path,
    dpi: int,
) -> None:
    """环状热图（外圈）+ 环内层次聚类树；蛋白名标在叶端（树与热图之间）。"""
    n_r, n_c = mat.shape
    if n_r < 1:
        return

    r_tree_inner = 0.06
    r_tree_outer = 0.38
    r_heat_inner = 0.41
    r_heat_outer = 1.0

    fig = plt.figure(figsize=(11, 11))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    theta_edges = np.linspace(0.0, 2 * np.pi, n_r + 1)
    r_edges = np.linspace(r_heat_inner, r_heat_outer, n_c + 1)
    TH, RR = np.meshgrid(theta_edges, r_edges, indexing="ij")
    C = np.ma.masked_invalid(mat.to_numpy(dtype=float))
    cmap = plt.get_cmap("cividis").copy()
    cmap.set_bad("#e8e8e8")
    im = ax.pcolormesh(
        TH,
        RR,
        C,
        cmap=cmap,
        vmin=0.0,
        vmax=vmax,
        shading="flat",
        zorder=2,
    )

    _radial_hierarchy_tree_polar(ax, Z, r_inner=r_tree_inner, r_outer=r_tree_outer)

    if n_r >= 2:
        fontsize = max(3.5, min(7.0, 220.0 / n_r))
        # 标签放在树外缘与热图之间；文字沿叶端径向（由圆心指向外），与枝末端方向一致
        gap_th = r_heat_inner - r_tree_outer
        r_label = r_tree_outer + max(0.01, min(0.06, 0.35 * gap_th))
        for i, name in enumerate(mat.index):
            th = 2 * np.pi * (i + 0.5) / n_r
            short = str(name) if len(str(name)) <= 14 else str(name)[:12] + "…"
            rot_deg = np.degrees(th)
            ax.text(
                th,
                r_label,
                short,
                rotation=rot_deg,
                ha="center",
                va="center",
                fontsize=fontsize,
                color="#222222",
                zorder=4,
                clip_on=False,
            )

    ax.set_ylim(0.0, r_heat_outer * 1.02)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title(
        f"{title}\n(ring heatmap + radial hierarchy tree via to_tree; unified row order)",
        pad=16,
    )
    ax.grid(False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.12, shrink=0.55)
    cbar.ax.tick_params(labelsize=8)
    try:
        cbar.set_label(f"max = {vmax:.4g}", loc="top", fontsize=9)
    except TypeError:
        cbar.ax.set_title(f"max = {vmax:.4g}", fontsize=9, pad=6)

    fig.subplots_adjust(left=0.06, right=0.9, top=0.9, bottom=0.06)
    fig.savefig(out_dir / fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

# 默认不纳入任何图的蛋白质（可 --exclude-symbols 覆盖）
DEFAULT_EXCLUDE_SYMBOLS = frozenset({"HNRNPU"})


def parse_exclude_symbols(arg: str) -> frozenset[str]:
    if arg is None or not str(arg).strip():
        return frozenset()
    return frozenset(x.strip() for x in str(arg).split(",") if x.strip())


def parse_test_metrics_summary(path: Path) -> dict[str, float]:
    d: dict[str, float] = {}
    if not path.is_file():
        return d
    for line in path.read_text(encoding="utf-8").splitlines():
        if "\t" not in line or line.startswith("Metric"):
            continue
        k, v = line.split("\t", 1)
        key = k.strip()
        try:
            d[key] = float(v.strip())
        except ValueError:
            pass
    return d


def collect_symbol_conditions(
    bundle_root: Path,
    excluded: frozenset[str],
) -> dict[str, list[tuple[str, dict[str, float]]]]:
    out: dict[str, list[tuple[str, dict[str, float]]]] = {}
    for sym_dir in sorted(bundle_root.iterdir()):
        if not sym_dir.is_dir():
            continue
        symbol = sym_dir.name
        if symbol in excluded:
            continue
        rows: list[tuple[str, dict[str, float]]] = []
        for sub in sorted(sym_dir.iterdir()):
            if not sub.is_dir():
                continue
            m = parse_test_metrics_summary(sub / "test_metrics_summary.txt")
            if not m:
                continue
            label = sub.name
            rows.append((label, m))
        if rows:
            out[symbol] = rows
    return out


# 全局 ordered 热图横轴固定顺序（不含 gland / adrenalgland self 列）
ORDERED_COL_KEYS = [
    "HepG2",
    "K562_to_HepG2",
    "HepG2_to_K562",
    "K562",
]

# per_symbol 柱状图：先 within-cell（self），细胞顺序与跨细胞热图一致；再 cross，按 train→test 排序
PER_SYMBOL_CELL_ORDER = ["HepG2", "K562", "adrenalgland"]


def _per_symbol_cell_rank(cell_slug: str) -> int:
    try:
        return PER_SYMBOL_CELL_ORDER.index(cell_slug)
    except ValueError:
        return len(PER_SYMBOL_CELL_ORDER) + (sum(ord(c) for c in cell_slug) % 997)


def per_symbol_condition_sort_key(label: str) -> tuple:
    """
    1) 所有 *self 在前，按 HepG2 → K562 → adrenalgland。
    2) cross_train_to_test 在后，先按 train 细胞顺序，再按 test 细胞顺序。
    3) 其它目录名排在最后，按字符串排序。
    """
    if label.endswith("self"):
        slug = label[: -len("self")]
        return (0, _per_symbol_cell_rank(slug), label)
    if label.startswith("cross_"):
        rest = label[len("cross_") :]
        if "_to_" in rest:
            train, test = rest.split("_to_", 1)
            return (
                1,
                _per_symbol_cell_rank(train),
                _per_symbol_cell_rank(test),
                label,
            )
    return (2, 0, 0, label)


def sort_per_symbol_rows(
    rows: list[tuple[str, dict[str, float]]],
) -> list[tuple[str, dict[str, float]]]:
    return sorted(rows, key=lambda r: per_symbol_condition_sort_key(r[0]))


def _self_slug_for_ordered_col(col: str) -> str | None:
    if col == "HepG2":
        return "HepG2"
    if col == "K562":
        return "K562"
    if col == "gland":
        return "adrenalgland"
    return None


def build_ordered_five_column_matrix(
    ok: pd.DataFrame,
    bundle_root: Path,
    *,
    mode: str,
) -> pd.DataFrame:
    """
    mode: mean_pcc | spearman | auc_ratio | pcc_ratio | spearman_ratio
    self 列从 bundle/<symbol>/<cell>self/ 读取；交叉两列仅 K562↔HepG2（来自 ok）。
    *_ratio：self 列 = 1.0（有该指标时）；交叉列 = cross / 目标细胞 self（与 AUC ratio 对齐）。
    """
    symbols = sorted(ok["symbol"].astype(str).unique())
    rows: list[dict[str, float | str]] = []
    ratio_mode = mode == "auc_ratio"
    pcc_ratio_mode = mode == "pcc_ratio"
    sp_ratio_mode = mode == "spearman_ratio"
    cross_key = (
        "cross_Mean_PCC"
        if mode in ("mean_pcc", "pcc_ratio")
        else "cross_Spearman"
        if mode in ("spearman", "spearman_ratio")
        else None
    )
    self_key = (
        "Mean_PCC"
        if mode in ("mean_pcc", "pcc_ratio")
        else "Spearman"
        if mode in ("spearman", "spearman_ratio")
        else None
    )

    def self_metric(symbol: str, slug: str, key: str) -> float:
        d = parse_test_metrics_summary(
            bundle_root
            / symbol
            / f"{slug}self"
            / "test_metrics_summary.txt"
        )
        v = d.get(key, np.nan)
        return float(v) if not pd.isna(v) else float("nan")

    for sym in symbols:
        sub = ok[ok["symbol"].astype(str) == sym]
        row: dict[str, float | str] = {"symbol": sym}

        def cross_val(train_cell: str, test_cell: str) -> float:
            tr = sub["train_term_name"].astype(str).str.strip()
            te = sub["test_term_name"].astype(str).str.strip()
            tr_m = tr == train_cell
            te_m = te == test_cell
            hit = sub[tr_m & te_m]
            if hit.empty:
                return float("nan")
            r0 = hit.iloc[0]
            if ratio_mode:
                c = r0.get("cross_Signal_AUC", np.nan)
                b = r0.get("baseline_test_Signal_AUC", np.nan)
                if pd.isna(c) or pd.isna(b) or b == 0:
                    return float("nan")
                return float(c) / float(b)
            if pcc_ratio_mode or sp_ratio_mode:
                c = r0.get(cross_key, np.nan)
                if test_cell == "HepG2":
                    denom_slug = "HepG2"
                elif test_cell == "K562":
                    denom_slug = "K562"
                else:
                    return float("nan")
                den = self_metric(sym, denom_slug, self_key)
                if pd.isna(c) or pd.isna(den) or den == 0:
                    return float("nan")
                return float(c) / float(den)
            v = r0.get(cross_key, np.nan)
            return float(v) if not pd.isna(v) else float("nan")

        for col in ORDERED_COL_KEYS:
            slug = _self_slug_for_ordered_col(col)
            if slug is not None:
                d = parse_test_metrics_summary(
                    bundle_root
                    / sym
                    / f"{slug}self"
                    / "test_metrics_summary.txt"
                )
                if ratio_mode:
                    row[col] = (
                        1.0 if d.get("Signal_AUC") is not None else float("nan")
                    )
                elif pcc_ratio_mode or sp_ratio_mode:
                    row[col] = (
                        1.0 if d.get(self_key) is not None else float("nan")
                    )
                else:
                    row[col] = float(d.get(self_key, np.nan))
            elif col == "K562_to_HepG2":
                row[col] = cross_val("K562", "HepG2")
            elif col == "HepG2_to_K562":
                row[col] = cross_val("HepG2", "K562")
            else:
                row[col] = float("nan")

        rows.append(row)

    out = pd.DataFrame(rows).set_index("symbol")
    return out.reindex(columns=ORDERED_COL_KEYS)


def order_rows_by_clustering(mat: pd.DataFrame) -> pd.DataFrame:
    """按行（蛋白质）对数值矩阵做层次聚类排序；聚类前 NaN 用列/行中位数再 0 填充。"""
    if mat.shape[0] < 2:
        return mat
    filled = (
        mat.fillna(mat.median(axis=0))
        .fillna(mat.median(axis=1))
        .fillna(0.0)
    )
    try:
        Z = linkage(filled.to_numpy(dtype=float), method="average", metric="euclidean")
        order = leaves_list(Z)
    except Exception:
        order = np.arange(mat.shape[0], dtype=int)
    return mat.iloc[order]


def _fill_matrix_for_cluster(mat: pd.DataFrame) -> np.ndarray:
    m = mat.fillna(mat.median(axis=0)).fillna(mat.median(axis=1)).fillna(0.0)
    return m.to_numpy(dtype=float)


def unified_row_order(mat_a: pd.DataFrame, mat_b: pd.DataFrame) -> np.ndarray:
    """对 [Mean_PCC 块 | Spearman 块] 行向拼接后做层次聚类，返回行索引顺序（与两矩阵 index 对齐）。"""
    if not mat_a.index.equals(mat_b.index):
        raise ValueError("unified_row_order: index mismatch")
    n = mat_a.shape[0]
    if n < 2:
        return np.arange(n, dtype=int)
    a = mat_a.copy()
    b = mat_b.copy()
    a.columns = [f"{c}_pcc" for c in a.columns]
    b.columns = [f"{c}_sp" for c in b.columns]
    comb = pd.concat([a, b], axis=1)
    filled = _fill_matrix_for_cluster(comb)
    try:
        Z = linkage(filled, method="average", metric="euclidean")
        return leaves_list(Z)
    except Exception:
        return np.arange(n, dtype=int)


def linkage_concatenated(mat_a: pd.DataFrame, mat_b: pd.DataFrame):
    """与 unified_row_order 相同的填充与拼接，返回 linkage 矩阵供 dendrogram。"""
    a = mat_a.copy()
    b = mat_b.copy()
    a.columns = [f"{c}_pcc" for c in a.columns]
    b.columns = [f"{c}_sp" for c in b.columns]
    comb = pd.concat([a, b], axis=1)
    filled = _fill_matrix_for_cluster(comb)
    return linkage(filled, method="average", metric="euclidean")


def _apply_heatmap_grid(ax, nrows: int, ncols: int) -> None:
    ax.set_xticks(np.arange(ncols) - 0.5, minor=True)
    ax.set_yticks(np.arange(nrows) - 0.5, minor=True)
    ax.grid(which="minor", color="0.85", linewidth=0.35)
    ax.tick_params(which="minor", bottom=False, left=False)


def _colorbar_with_max_label(
    fig,
    im,
    ax,
    vmin: float,
    vmax: float,
    *,
    fraction: float = 0.035,
    pad: float = 0.02,
) -> None:
    cbar = fig.colorbar(im, ax=ax, fraction=fraction, pad=pad)
    cbar.ax.tick_params(labelsize=8)
    try:
        cbar.set_label(f"color scale max = {vmax:.4g}", loc="top", fontsize=9)
    except TypeError:
        cbar.ax.set_title(f"color scale max = {vmax:.4g}", fontsize=9, pad=6)
    ticks = [vmin, vmax]
    if vmin < vmax:
        mid = 0.5 * (vmin + vmax)
        if abs(mid - vmin) > 1e-6 and abs(vmax - mid) > 1e-6:
            ticks = [vmin, mid, vmax]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:.4g}" for t in ticks])


def plot_ordered_heatmap_variants(
    mat_pcc: pd.DataFrame,
    mat_sp: pd.DataFrame,
    vmax_pcc: float,
    vmax_sp: float,
    out_dir: Path,
    dpi: int,
) -> None:
    """
    审稿改版：统一行序（[PCC|Spearman] 拼接后聚类）、cividis、网格、树状图、合并版式。
    输出至 out_dir（通常为 global/ordered_heatmap_variants/）。
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    order = unified_row_order(mat_pcc, mat_sp)
    pcc_o = mat_pcc.iloc[order]
    sp_o = mat_sp.iloc[order]
    Z = linkage_concatenated(mat_pcc, mat_sp)
    n_r, n_c = pcc_o.shape
    h = max(6.0, n_r * 0.12)
    w_base = max(7.0, n_c * 1.05)

    def draw_acd_single(
        mat: pd.DataFrame, fname: str, title: str, vmax: float
    ) -> None:
        fig, ax = plt.subplots(figsize=(w_base, h))
        arr = np.ma.masked_invalid(mat.to_numpy(dtype=float))
        im = ax.imshow(arr, aspect="auto", cmap="cividis", vmin=0, vmax=vmax)
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_xticklabels(mat.columns, rotation=25, ha="right", fontsize=9)
        ax.set_yticklabels(mat.index, fontsize=6)
        _apply_heatmap_grid(ax, mat.shape[0], mat.shape[1])
        ax.set_facecolor("white")
        ax.set_title(
            f"{title}\n(unified row order: [Mean_PCC | Spearman] concat clustering; cividis)"
        )
        _colorbar_with_max_label(fig, im, ax, 0.0, vmax)
        fig.subplots_adjust(left=0.1, right=0.88, top=0.9, bottom=0.12)
        fig.savefig(out_dir / fname, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    draw_acd_single(
        pcc_o,
        "variant_A_Mean_PCC_cividis_grid.png",
        "Mean_PCC",
        vmax_pcc,
    )
    draw_acd_single(
        sp_o,
        "variant_A_Spearman_cividis_grid.png",
        "Spearman",
        vmax_sp,
    )

    _polar_ring_heatmap_with_tree(
        pcc_o,
        Z,
        "Mean_PCC",
        vmax_pcc,
        fname="variant_B_Mean_PCC_with_dendrogram.png",
        out_dir=out_dir,
        dpi=dpi,
    )
    _polar_ring_heatmap_with_tree(
        sp_o,
        Z,
        "Spearman",
        vmax_sp,
        fname="variant_B_Spearman_with_dendrogram.png",
        out_dir=out_dir,
        dpi=dpi,
    )

    # M1: 1×2 共享行序与 y
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(w_base * 2.0, h),
        sharey=True,
        gridspec_kw={"wspace": 0.14},
    )
    for ax, mat, vmax, ttl in (
        (axes[0], pcc_o, vmax_pcc, "Mean_PCC"),
        (axes[1], sp_o, vmax_sp, "Spearman"),
    ):
        arr = np.ma.masked_invalid(mat.to_numpy(dtype=float))
        im = ax.imshow(arr, aspect="auto", cmap="cividis", vmin=0, vmax=vmax)
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_xticklabels(mat.columns, rotation=25, ha="right", fontsize=9)
        ax.set_yticks(np.arange(mat.shape[0]))
        if ax is axes[1]:
            ax.set_yticklabels([])
        else:
            ax.set_yticklabels(mat.index, fontsize=6)
        _apply_heatmap_grid(ax, mat.shape[0], mat.shape[1])
        ax.set_facecolor("white")
        ax.set_title(ttl)
        _colorbar_with_max_label(fig, im, ax, 0.0, vmax, fraction=0.046, pad=0.04)
    fig.suptitle(
        "merge_M1: Mean_PCC | Spearman (unified row order, separate scales)",
        fontsize=10,
        y=1.01,
    )
    fig.subplots_adjust(left=0.06, right=0.97, top=0.88, bottom=0.12)
    fig.savefig(out_dir / "merge_M1_side_by_side.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    # M2: 8 列（PCC×4 | Spearman×4），左右块各自除以 vmax 缩放到 [0,1] 共用色条
    pcc_n = (pcc_o.to_numpy(dtype=float) / vmax_pcc).clip(0.0, 1.0)
    sp_n = (sp_o.to_numpy(dtype=float) / vmax_sp).clip(0.0, 1.0)
    wide = np.hstack([pcc_n, sp_n])
    fig, ax = plt.subplots(figsize=(w_base * 1.85, h))
    im = ax.imshow(wide, aspect="auto", cmap="cividis", vmin=0.0, vmax=1.0)
    ax.set_xticks(np.arange(wide.shape[1]))
    xlabs = [f"{c}\n(PCC)" for c in ORDERED_COL_KEYS] + [
        f"{c}\n(Sp)" for c in ORDERED_COL_KEYS
    ]
    ax.set_xticklabels(xlabs, rotation=30, ha="right", fontsize=7)
    ax.set_yticks(np.arange(n_r))
    ax.set_yticklabels(pcc_o.index, fontsize=6)
    _apply_heatmap_grid(ax, wide.shape[0], wide.shape[1])
    ax.axvline(x=n_c - 0.5, color="white", lw=2.0)
    ax.set_facecolor("white")
    ax.set_title(
        "merge_M2: PCC (÷vmax_pcc) | Spearman (÷vmax_sp), 4+4 cols, shared [0,1] colormap"
    )
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="scaled value [0,1]")
    fig.subplots_adjust(left=0.12, right=0.92, top=0.9, bottom=0.2)
    fig.savefig(out_dir / "merge_M2_wide_8col_scaled.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    # M3: 上下，共享 x
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(w_base, h * 1.65),
        sharex=True,
        gridspec_kw={"hspace": 0.18},
    )
    for ax, mat, vmax, ttl in (
        (axes[0], pcc_o, vmax_pcc, "Mean_PCC"),
        (axes[1], sp_o, vmax_sp, "Spearman"),
    ):
        arr = np.ma.masked_invalid(mat.to_numpy(dtype=float))
        im = ax.imshow(arr, aspect="auto", cmap="cividis", vmin=0, vmax=vmax)
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_yticklabels(mat.index, fontsize=6)
        _apply_heatmap_grid(ax, mat.shape[0], mat.shape[1])
        ax.set_facecolor("white")
        ax.set_title(ttl)
        _colorbar_with_max_label(fig, im, ax, 0.0, vmax, fraction=0.03, pad=0.02)
    axes[0].tick_params(axis="x", labelbottom=False)
    axes[1].set_xticklabels(pcc_o.columns, rotation=25, ha="right", fontsize=9)
    fig.suptitle("merge_M3: vertical stack (unified row order)", fontsize=10, y=0.995)
    fig.subplots_adjust(left=0.1, right=0.93, top=0.93, bottom=0.08)
    fig.savefig(out_dir / "merge_M3_vertical_stack.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def plot_ordered_heatmaps(
    ok: pd.DataFrame,
    bundle_root: Path,
    plot_dir: Path,
    dpi: int = 150,
    *,
    heatmap_variants: bool = True,
) -> None:
    """固定横轴 4 列：HepG2, K562_to_HepG2, HepG2_to_K562, K562（不含 gland）。"""

    def draw(mat: pd.DataFrame, title: str, fname: str, vmin, vmax, cmap: str):
        mat_o = order_rows_by_clustering(mat)
        fig, ax = plt.subplots(
            figsize=(max(7, mat_o.shape[1] * 1.1), max(6, mat_o.shape[0] * 0.14))
        )
        arr = np.ma.masked_invalid(mat_o.to_numpy(dtype=float))
        im = ax.imshow(arr, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks(np.arange(mat_o.shape[1]))
        ax.set_yticks(np.arange(mat_o.shape[0]))
        ax.set_xticklabels(mat_o.columns, rotation=25, ha="right", fontsize=9)
        ax.set_yticklabels(mat_o.index, fontsize=6)
        ax.set_title(f"{title}\n(rows: hierarchical clustering)")
        _colorbar_with_max_label(fig, im, ax, vmin, vmax)
        fig.tight_layout()
        fig.savefig(plot_dir / fname, dpi=dpi)
        plt.close(fig)

    mat_pcc = build_ordered_five_column_matrix(ok, bundle_root, mode="mean_pcc")
    v = np.nanpercentile(mat_pcc.to_numpy(), 98)
    vmax_pcc = float(min(0.85, max(0.35, v * 1.05))) if np.isfinite(v) else 0.6
    draw(
        mat_pcc,
        "Ordered panel: Mean_PCC (within-cell *self + K562↔HepG2 cross)",
        "global_ordered_Mean_PCC_heatmap.png",
        0,
        vmax_pcc,
        "viridis",
    )

    mat_sp = build_ordered_five_column_matrix(ok, bundle_root, mode="spearman")
    v2 = np.nanpercentile(mat_sp.to_numpy(), 98)
    vmax_sp = float(min(0.9, max(0.4, v2 * 1.05))) if np.isfinite(v2) else 0.65
    draw(
        mat_sp,
        "Ordered panel: Spearman (within-cell *self + K562↔HepG2 cross)",
        "global_ordered_Spearman_heatmap.png",
        0,
        vmax_sp,
        "plasma",
    )

    mat_ratio = build_ordered_five_column_matrix(ok, bundle_root, mode="auc_ratio")
    draw(
        mat_ratio,
        "Ordered panel: cross/baseline Signal_AUC (self columns = 1.0)",
        "global_cross_over_baseline_Signal_AUC_ratio.png",
        0,
        1.2,
        "RdYlGn",
    )

    mat_pcc_r = build_ordered_five_column_matrix(ok, bundle_root, mode="pcc_ratio")
    draw(
        mat_pcc_r,
        "Ordered panel: cross/target-self Mean_PCC ratio (self columns = 1.0)",
        "global_ordered_Mean_PCC_ratio_heatmap.png",
        0,
        1.2,
        "RdYlGn",
    )

    mat_sp_r = build_ordered_five_column_matrix(ok, bundle_root, mode="spearman_ratio")
    draw(
        mat_sp_r,
        "Ordered panel: cross/target-self Spearman ratio (self columns = 1.0)",
        "global_ordered_Spearman_ratio_heatmap.png",
        0,
        1.2,
        "RdYlGn",
    )

    if heatmap_variants:
        variant_dir = plot_dir / "ordered_heatmap_variants"
        plot_ordered_heatmap_variants(
            mat_pcc, mat_sp, vmax_pcc, vmax_sp, variant_dir, dpi
        )
        print(f"改版热图（统一行序）-> {variant_dir}")


def _symbol_from_per_symbol_png(path: Path) -> str:
    stem = path.stem
    suf = "_metrics"
    return stem[: -len(suf)] if stem.endswith(suf) else stem


def plot_per_symbol(
    bundle_root: Path,
    plot_dir: Path,
    dpi: int = 120,
    *,
    excluded: frozenset[str] = frozenset(),
) -> int:
    """仅输出至少 2 个条件（可比较）的蛋白图，并删除不再满足条件的旧 PNG。"""
    plot_dir.mkdir(parents=True, exist_ok=True)
    data = collect_symbol_conditions(bundle_root, excluded)
    compare_data = {sym: rows for sym, rows in data.items() if len(rows) >= 2}
    allowed = set(compare_data.keys())

    for png in plot_dir.glob("*_metrics.png"):
        if _symbol_from_per_symbol_png(png) not in allowed:
            try:
                png.unlink()
            except OSError:
                pass

    n = 0
    for symbol, rows in sorted(compare_data.items()):
        rows = sort_per_symbol_rows(rows)
        labels = [r[0] for r in rows]
        aucs = [r[1].get("Signal_AUC", np.nan) for r in rows]
        aps = [r[1].get("Signal_AP", np.nan) for r in rows]
        pccs = [r[1].get("Mean_PCC", np.nan) for r in rows]
        spears = [r[1].get("Spearman", np.nan) for r in rows]

        x = np.arange(len(labels))
        bw = 0.18
        offs = np.array([-1.5, -0.5, 0.5, 1.5]) * bw
        fig, ax = plt.subplots(figsize=(max(7, len(labels) * 0.55), 5))
        ax.bar(x + offs[0], aucs, width=bw, label="Signal_AUC", color="#4C72B0")
        ax.bar(x + offs[1], aps, width=bw, label="Signal_AP", color="#DD8452")
        ax.bar(x + offs[2], pccs, width=bw, label="Mean_PCC", color="#8172B3")
        ax.bar(x + offs[3], spears, width=bw, label="Spearman", color="#55A868")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Score")
        ax.set_title(f"{symbol}: within-cell (×self) vs cross-cell")
        ax.legend(loc="upper right", fontsize=7, ncol=2)
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(plot_dir / f"{symbol}_metrics.png", dpi=dpi)
        plt.close(fig)
        n += 1
    if len(data) > len(compare_data):
        print(
            f"per_symbol: 仅保留 ≥2 个条件的图，跳过 {len(data) - len(compare_data)} 个 symbol；"
            f"已清理不在此集合内的旧 PNG"
        )
    return n


def plot_global_from_summary(
    summary_tsv: Path,
    plot_dir: Path,
    bundle_root: Path,
    dpi: int = 150,
    *,
    excluded: frozenset[str] = frozenset(),
    heatmap_variants: bool = True,
) -> None:
    plot_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(summary_tsv, sep="\t")
    ok = df[df["status"] == "ok"].copy()
    if excluded:
        ok = ok[~ok["symbol"].astype(str).isin(excluded)]
    if ok.empty:
        print("cross_cell_summary 无 status==ok 行，跳过全局交叉热图")
        return

    ok["col"] = (
        ok["train_term_name"].astype(str) + "→" + ok["test_term_name"].astype(str)
    )
    pivot_auc = ok.pivot_table(
        index="symbol", columns="col", values="cross_Signal_AUC", aggfunc="first"
    )
    pivot_ap = ok.pivot_table(
        index="symbol", columns="col", values="cross_Signal_AP", aggfunc="first"
    )

    ratio = ok.copy()
    ratio["r_auc"] = ratio["cross_Signal_AUC"] / ratio["baseline_test_Signal_AUC"]
    ratio["r_ap"] = ratio["cross_Signal_AP"] / ratio["baseline_test_Signal_AP"]
    pivot_r = ratio.pivot_table(
        index="symbol", columns="col", values="r_auc", aggfunc="first"
    )

    def heatmap(mat: pd.DataFrame, title: str, fname: str, vmin=None, vmax=None):
        mat_o = order_rows_by_clustering(mat)
        fig, ax = plt.subplots(figsize=(max(8, mat_o.shape[1] * 0.35), max(6, mat_o.shape[0] * 0.18)))
        arr = np.ma.masked_invalid(mat_o.to_numpy(dtype=float))
        im = ax.imshow(arr, aspect="auto", cmap="viridis", vmin=vmin, vmax=vmax)
        ax.set_xticks(np.arange(mat_o.shape[1]))
        ax.set_yticks(np.arange(mat_o.shape[0]))
        ax.set_xticklabels(mat_o.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(mat_o.index, fontsize=7)
        ax.set_title(f"{title}\n(rows: hierarchical clustering)")
        if vmin is not None and vmax is not None:
            _colorbar_with_max_label(fig, im, ax, float(vmin), float(vmax), fraction=0.02, pad=0.02)
        else:
            fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
        fig.tight_layout()
        fig.savefig(plot_dir / fname, dpi=dpi)
        plt.close(fig)

    heatmap(pivot_auc, "Cross-cell Signal_AUC (chr3 test)", "global_cross_Signal_AUC_heatmap.png", vmin=0, vmax=1)
    heatmap(pivot_ap, "Cross-cell Signal_AP (AUPRC)", "global_cross_Signal_AP_heatmap.png", vmin=0, vmax=1)
    heatmap(
        pivot_r,
        "Cross / target within-cell baseline (Signal_AUC ratio, all directions)",
        "global_cross_over_baseline_Signal_AUC_ratio_all_directions.png",
        vmin=0,
        vmax=1.2,
    )

    if bundle_root.is_dir():
        plot_ordered_heatmaps(
            ok, bundle_root, plot_dir, dpi=dpi, heatmap_variants=heatmap_variants
        )

    # 全蛋白均值条形：各 cross 方向
    mean_auc = pivot_auc.mean(axis=0, skipna=True).sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(max(6, len(mean_auc) * 0.4), 4))
    mean_auc.plot(kind="bar", ax=ax, color="#4C72B0")
    ax.set_ylabel("Mean Signal_AUC")
    ax.set_title("Mean cross-cell AUC by train→test direction (all symbols with ok)")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(plot_dir / "global_mean_AUC_by_direction.png", dpi=dpi)
    plt.close(fig)

    # 散点：目标细胞 within-cell baseline AUC vs 交叉预测 AUC（全点汇总）
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    x = ok["baseline_test_Signal_AUC"].to_numpy()
    y = ok["cross_Signal_AUC"].to_numpy()
    ax.scatter(x, y, alpha=0.35, s=18, c="#4C72B0", edgecolors="none")
    lim = (0.45, 1.0)
    ax.plot(lim, lim, "k--", lw=1, alpha=0.7, label="y=x")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Target within-cell baseline Signal_AUC")
    ax.set_ylabel("Cross-cell Signal_AUC")
    ax.set_title("All cross pairs: cross vs target self baseline")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(plot_dir / "global_scatter_cross_vs_target_baseline_AUC.png", dpi=dpi)
    plt.close(fig)

    # 组合大图：热图缩略 + 条形 + 散点 + ratio 分布
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    # (0,0) cross AUC heatmap 缩小字体（行聚类）
    ax = axes[0, 0]
    p_auc_o = order_rows_by_clustering(pivot_auc)
    arr = np.ma.masked_invalid(p_auc_o.to_numpy(dtype=float))
    im = ax.imshow(arr, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_title("Cross Signal_AUC (rows clustered)")
    ax.set_xticks(np.arange(p_auc_o.shape[1]))
    ax.set_yticks(np.arange(p_auc_o.shape[0]))
    ax.set_xticklabels(p_auc_o.columns, rotation=45, ha="right", fontsize=6)
    ax.set_yticklabels(p_auc_o.index, fontsize=5)
    _colorbar_with_max_label(fig, im, ax, 0.0, 1.0, fraction=0.046, pad=0.04)

    ax = axes[0, 1]
    mean_auc2 = pivot_auc.mean(axis=0, skipna=True).sort_values(ascending=False)
    mean_auc2.plot(kind="bar", ax=ax, color="#55A868")
    ax.set_title("Mean AUC by direction")
    ax.tick_params(axis="x", rotation=45, labelsize=7)

    ax = axes[1, 0]
    ax.scatter(x, y, alpha=0.3, s=14, c="#4C72B0")
    ax.plot(lim, lim, "k--", lw=0.8, alpha=0.7)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Target baseline AUC")
    ax.set_ylabel("Cross AUC")
    ax.set_title("Cross vs baseline")
    ax.grid(alpha=0.25)

    ax = axes[1, 1]
    r = ratio["r_auc"].dropna()
    ax.hist(r, bins=30, color="#C44E52", alpha=0.85, edgecolor="white")
    ax.axvline(1.0, color="k", linestyle="--", lw=1)
    ax.set_title("Distribution of cross/baseline AUC ratio")
    ax.set_xlabel("Ratio")
    ax.set_ylabel("Count")

    fig.suptitle("MARCP cross-cell overview (all ok pairs)", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(plot_dir / "global_four_panel_overview.png", dpi=dpi)
    plt.close(fig)


def plot_global_self_heatmap(
    bundle_root: Path,
    plot_dir: Path,
    dpi: int = 150,
    *,
    excluded: frozenset[str] = frozenset(),
) -> None:
    """由 metrics_viz_bundle 中 *self 文件夹汇总 within-cell Signal_AUC 热图。"""
    rows = []
    for sym_dir in sorted(bundle_root.iterdir()):
        if not sym_dir.is_dir():
            continue
        sym = sym_dir.name
        if sym in excluded:
            continue
        for sub in sym_dir.iterdir():
            if not sub.is_dir() or not sub.name.endswith("self"):
                continue
            cell = sub.name[: -len("self")]
            m = parse_test_metrics_summary(sub / "test_metrics_summary.txt")
            if not m:
                continue
            rows.append(
                {
                    "symbol": sym,
                    "cell": cell,
                    "Signal_AUC": m.get("Signal_AUC", np.nan),
                    "Signal_AP": m.get("Signal_AP", np.nan),
                }
            )
    if not rows:
        print("bundle 中无 *self 指标，跳过 within-cell 热图")
        return
    df = pd.DataFrame(rows)
    pivot = df.pivot_table(
        index="symbol", columns="cell", values="Signal_AUC", aggfunc="first"
    )
    pivot_o = order_rows_by_clustering(pivot)

    fig, ax = plt.subplots(
        figsize=(max(6, pivot_o.shape[1] * 0.5), max(6, pivot_o.shape[0] * 0.15))
    )
    arr = np.ma.masked_invalid(pivot_o.to_numpy(dtype=float))
    im = ax.imshow(arr, aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(np.arange(pivot_o.shape[1]))
    ax.set_yticks(np.arange(pivot_o.shape[0]))
    ax.set_xticklabels(pivot_o.columns, rotation=45, ha="right")
    ax.set_yticklabels(pivot_o.index, fontsize=6)
    ax.set_title(
        "Within-cell baseline Signal_AUC (*self, from MARCP metrics)\n"
        "(rows: hierarchical clustering)"
    )
    _colorbar_with_max_label(fig, im, ax, 0.0, 1.0, fraction=0.02, pad=0.02)
    fig.tight_layout()
    plot_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_dir / "global_within_cell_Signal_AUC_heatmap.png", dpi=dpi)
    plt.close(fig)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--bundle",
        type=Path,
        default=root / "metrics_viz_bundle",
    )
    ap.add_argument(
        "--per-symbol-dir",
        type=Path,
        default=root / "metrics_viz_plots" / "per_symbol",
    )
    ap.add_argument(
        "--global-dir",
        type=Path,
        default=root / "metrics_viz_plots" / "global",
    )
    ap.add_argument(
        "--summary",
        type=Path,
        default=root / "cross_cell_summary.tsv",
    )
    ap.add_argument("--dpi", type=int, default=120)
    ap.add_argument(
        "--exclude-symbols",
        type=str,
        default=None,
        metavar="LIST",
        help=(
            "逗号分隔、不绘制的 symbol。省略本参数时默认排除 HNRNPU；"
            "传空字符串 '' 表示不排除任何 symbol；否则仅排除列表中的名称。"
        ),
    )
    ap.add_argument(
        "--no-heatmap-variants",
        action="store_true",
        help=(
            "不生成 global/ordered_heatmap_variants/ 下的统一行序改版热图 "
            "(variant_A/B、merge_M1/M2/M3)。"
        ),
    )
    args = ap.parse_args()

    if args.exclude_symbols is None:
        excluded = DEFAULT_EXCLUDE_SYMBOLS
    elif str(args.exclude_symbols).strip() == "":
        excluded = frozenset()
    else:
        excluded = parse_exclude_symbols(args.exclude_symbols)

    if excluded:
        print(f"排除不绘制的 symbol: {', '.join(sorted(excluded))}")

    if not args.bundle.is_dir():
        raise SystemExit(f"请先运行 prepare_symbol_metric_folders.py；缺少: {args.bundle}")

    n = plot_per_symbol(
        args.bundle, args.per_symbol_dir, dpi=args.dpi, excluded=excluded
    )
    print(f"单蛋白图: {n} 张 -> {args.per_symbol_dir}")

    plot_global_self_heatmap(
        args.bundle, args.global_dir, dpi=max(args.dpi, 150), excluded=excluded
    )
    if args.summary.is_file():
        plot_global_from_summary(
            args.summary,
            args.global_dir,
            args.bundle,
            dpi=max(args.dpi, 150),
            excluded=excluded,
            heatmap_variants=not args.no_heatmap_variants,
        )
    else:
        print(f"未找到 {args.summary}，跳过 cross_cell 汇总热图")
    print(f"全局图 -> {args.global_dir}")


if __name__ == "__main__":
    main()
