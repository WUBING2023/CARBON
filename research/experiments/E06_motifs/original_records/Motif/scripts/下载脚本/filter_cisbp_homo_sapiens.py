#!/usr/bin/env python3
"""
从已下载的全库 CISBP-RNA（RBP_Information*.txt + pwms/）生成 Homo_sapiens 子集。

- 元数据：筛 RBP_Species == Homo_sapiens。
- PWM：对 all_motifs 中 Motif_ID 非「.」的唯一 Motif_ID，从全库 pwms/ 复制或硬链接到 homo_sapiens/pwms/。

全库表中人类唯一 Motif_ID 数量通常多于官网「物种 zip」内 pwms_all_motifs 文件数（服务器侧可能再筛过），
本脚本与全库 pwms 命名 100% 一致。

用法（在 Motif 目录）：
  python3 scripts/filter_cisbp_homo_sapiens.py
  python3 scripts/filter_cisbp_homo_sapiens.py --out-dir reference_data/CISBP-RNA/homo_sapiens --copy
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

_MOTIF_ROOT = Path(__file__).resolve().parents[1]
_DEFAULT_CISBP = _MOTIF_ROOT / "reference_data/CISBP-RNA"
_DEFAULT_OUT = _DEFAULT_CISBP / "homo_sapiens"


def main() -> int:
    ap = argparse.ArgumentParser(description="从全库 CISBP 生成 Homo_sapiens 子集")
    ap.add_argument("--cisbp-dir", type=Path, default=_DEFAULT_CISBP, help="含 RBP_Information*.txt 与 pwms/ 的目录")
    ap.add_argument("--out-dir", type=Path, default=_DEFAULT_OUT, help="输出目录（将创建 pwms/ 与两个 txt）")
    ap.add_argument(
        "--copy",
        action="store_true",
        help="用复制而非硬链接写入 pwms（默认硬链接以省空间）",
    )
    args = ap.parse_args()

    try:
        import pandas as pd
    except ImportError:
        print("[error] 需要 pandas", file=sys.stderr)
        return 1

    cis = args.cisbp_dir.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    inf = cis / "RBP_Information.txt"
    inf_all = cis / "RBP_Information_all_motifs.txt"
    pwms_src = cis / "pwms"

    for p in (inf, inf_all, pwms_src):
        if not p.exists():
            print(f"[error] 缺少: {p}", file=sys.stderr)
            return 1

    print(f"[info] 读取 {inf_all} ...")
    df_all = pd.read_csv(inf_all, sep="\t", dtype=str, low_memory=False)
    if "RBP_Species" not in df_all.columns:
        print("[error] 无 RBP_Species 列", file=sys.stderr)
        return 1

    hom_all = df_all[df_all["RBP_Species"] == "Homo_sapiens"].copy()
    hom_all.to_csv(out / "RBP_Information_all_motifs.txt", sep="\t", index=False)
    print(f"[ok] Homo_sapiens 行（all_motifs）: {len(hom_all)} -> {out / 'RBP_Information_all_motifs.txt'}")

    print(f"[info] 读取 {inf} ...")
    df_inf = pd.read_csv(inf, sep="\t", dtype=str, low_memory=False)
    hom_inf = df_inf[df_inf["RBP_Species"] == "Homo_sapiens"].copy()
    hom_inf.to_csv(out / "RBP_Information.txt", sep="\t", index=False)
    print(f"[ok] Homo_sapiens 行（RBP_Information）: {len(hom_inf)} -> {out / 'RBP_Information.txt'}")

    sub_m = hom_all[hom_all["Motif_ID"].notna() & (hom_all["Motif_ID"] != ".")]
    motif_ids = sorted(sub_m["Motif_ID"].unique())
    print(f"[info] 人类且 Motif_ID 非「.」的唯一 Motif_ID: {len(motif_ids)}")

    pwms_dst = out / "pwms"
    pwms_dst.mkdir(parents=True, exist_ok=True)

    n_ok = 0
    missing: list[str] = []
    for mid in motif_ids:
        name = f"{mid}.txt"
        src = pwms_src / name
        dst = pwms_dst / name
        if not src.is_file():
            missing.append(mid)
            continue
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        if args.copy:
            shutil.copy2(src, dst)
        else:
            try:
                if hasattr(dst, "hardlink_to"):
                    dst.hardlink_to(src)
                else:
                    os.link(src, dst)
            except OSError:
                shutil.copy2(src, dst)
        n_ok += 1

    print(f"[ok] pwms/ 写入 {n_ok} 个文件 -> {pwms_dst}")
    if missing:
        print(f"[warn] 全库 pwms 中缺失 {len(missing)} 个 Motif_ID 文件（前 10 个）: {missing[:10]}", file=sys.stderr)

    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
