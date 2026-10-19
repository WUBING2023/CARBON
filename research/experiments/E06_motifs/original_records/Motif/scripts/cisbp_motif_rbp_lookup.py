#!/usr/bin/env python3
"""
CISBP-RNA：Motif_ID（M*，与 pwms/M*_2.00.txt 对应）与 RBP_ID（T*）对照与任务目录解析。

为何「对不上」：
  - PWM 文件名 = Motif_ID（例如 M00634_2.00.txt），不是 RBP_ID。
  - 表中用列 Motif_ID 与 pwms 对应；RBP_ID（T 开头）标识一条 RBP 记录。
  - RBP_Information.txt 按行排序，前几行常为 Motif_ID 为「.」的无 motif 占位，看起来像「只有 T」。
  - RBP_Information_all_motifs.txt 展开 motif–RBP 关系，查 M 应用此文件（或 RBP_Information 中 Motif_ID 非空的行）。

为何 Top hit 的 RBP_Name 常不是你的实验蛋白：
  - Tomtom / 全库 KL 是在**全库**里找**最相似 PWM**；该 Motif_ID 可能来自其它物种或「共享」文献 motif（同一 M 多行、多 RBP）。
  - 要与「文件夹里的基因名」对照，应限制 Homo_sapiens 且 RBP_Name = 任务基因符号。

用法示例：
  python3 scripts/cisbp_motif_rbp_lookup.py motif --id M00634_2.00
  python3 scripts/cisbp_motif_rbp_lookup.py human-gene --symbol RBFOX2
  python3 scripts/cisbp_motif_rbp_lookup.py scan-outputs --outputs-dir outputs
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

_MOTIF_ROOT = Path(__file__).resolve().parents[1]
_DEFAULT_ALL = _MOTIF_ROOT / "reference_data/CISBP-RNA/RBP_Information_all_motifs.txt"

# 任务目录：<Gene>_<HepG2|K562>_ENCSR...
_TASK_RE = re.compile(
    r"^(?P<gene>[A-Za-z0-9]+)_(?P<cell>HepG2|K562)_(?P<encsr>ENCSR[A-Za-z0-9]+)$"
)


def load_all_motifs(path: Path):
    import pandas as pd

    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def cmd_motif(args: argparse.Namespace) -> int:
    df = load_all_motifs(Path(args.table))
    mid = args.id.strip()
    if not mid.endswith(".txt"):
        pass
    else:
        mid = mid.replace(".txt", "")
    sub = df[df["Motif_ID"] == mid]
    if sub.empty:
        print(f"[无匹配] Motif_ID={mid!r} 在 {args.table} 中不存在。", file=sys.stderr)
        return 1
    cols = [
        "RBP_ID",
        "Motif_ID",
        "RBP_Name",
        "RBP_Species",
        "RBP_Status",
        "Motif_Type",
        "MSource_Type",
    ]
    cols = [c for c in cols if c in sub.columns]
    print(sub[cols].to_string(index=False))
    n = len(sub)
    n_hsap = (sub["RBP_Species"] == "Homo_sapiens").sum() if "RBP_Species" in sub.columns else 0
    print(f"\n[统计] 该 Motif_ID 共 {n} 行；其中 Homo_sapiens: {n_hsap} 行。")
    if n > 1 and n_hsap == 0:
        print(
            "[提示] 全库命中该 motif 时，若未过滤物种，RBP_Name 可能是非人物种或共享文献矩阵。",
            file=sys.stderr,
        )
    return 0


def cmd_human_gene(args: argparse.Namespace) -> int:
    df = load_all_motifs(Path(args.table))
    sym = args.symbol.strip()
    sub = df[(df["RBP_Species"] == "Homo_sapiens") & (df["RBP_Name"].str.upper() == sym.upper())]
    if sub.empty:
        print(f"[无匹配] Homo_sapiens 且 RBP_Name={sym!r} 无记录。", file=sys.stderr)
        return 1
    cols = ["RBP_ID", "Motif_ID", "RBP_Name", "RBP_Status", "Motif_Type", "MSource_Type"]
    cols = [c for c in cols if c in sub.columns]
    print(sub[cols].to_string(index=False))
    return 0


def parse_task_name(name: str) -> dict | None:
    m = _TASK_RE.match(name)
    if not m:
        return None
    return {"gene": m.group("gene"), "cell": m.group("cell"), "encsr": m.group("encsr")}


def read_tomtom_top_target(task_dir: Path) -> str | None:
    tsv = task_dir / "tomtom_cisbp" / "tomtom_out" / "tomtom.tsv"
    if not tsv.is_file() or tsv.stat().st_size == 0:
        return None
    with open(tsv, newline="", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            tid = row.get("Target_ID", "").strip()
            if tid:
                return tid
    return None


def cmd_scan_outputs(args: argparse.Namespace) -> int:
    df = load_all_motifs(Path(args.table))
    out = Path(args.outputs_dir).resolve()
    if not out.is_dir():
        print(f"[error] 不是目录: {out}", file=sys.stderr)
        return 1

    rows_out: list[tuple[str, str, str, str, str]] = []
    nlim = int(args.max_tasks)

    # 允许直接传入单个任务目录（例如 .../outputs/RBFOX2_K562_ENCSR756CKJ）
    if parse_task_name(out.name):
        children = [out]
    else:
        children = sorted(p for p in out.iterdir() if p.is_dir())

    for child in children:
        if nlim > 0 and len(rows_out) >= nlim:
            break
        info = parse_task_name(child.name)
        if not info:
            continue
        gene = info["gene"]
        tt = read_tomtom_top_target(child)
        match_note = ""
        if tt:
            sub = df[(df["Motif_ID"] == tt) & (df["RBP_Species"] == "Homo_sapiens")]
            if sub.empty:
                match_note = "top_motif无人类行"
            else:
                names = sub["RBP_Name"].dropna().unique().tolist()
                if gene.upper() in {n.upper() for n in names}:
                    match_note = "人类行含任务基因"
                else:
                    match_note = f"人类: {','.join(names[:5])}{'…' if len(names) > 5 else ''}"
        else:
            match_note = "无tomtom.tsv或为空"

        rows_out.append((child.name, gene, info["cell"], tt or "—", match_note))

    if not rows_out:
        print("[无输出] 未解析到 <Gene>_<Cell>_ENCSR* 任务目录。", file=sys.stderr)
        return 1

    w = max(len(r[0]) for r in rows_out)
    print(f"{'task':<{w}}  gene        cell   top_Target_ID    cisbp人类备注")
    for name, gene, cell, tt, note in rows_out:
        print(f"{name:<{w}}  {gene:<10}  {cell:<5}  {tt:<15}  {note}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="CISBP Motif_ID / RBP 与 outputs 任务对照")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("motif", help="按 Motif_ID（如 M00634_2.00）查所有物种行")
    p1.add_argument("--id", required=True, help="Motif_ID，与 pwms 中文件名一致（可带或不带 .txt）")
    p1.add_argument("--table", type=Path, default=_DEFAULT_ALL)
    p1.set_defaults(func=cmd_motif)

    p2 = sub.add_parser("human-gene", help="人类某基因符号在 CISBP 中的 motif 行")
    p2.add_argument("--symbol", required=True, help="如 RBFOX2、PRPF8")
    p2.add_argument("--table", type=Path, default=_DEFAULT_ALL)
    p2.set_defaults(func=cmd_human_gene)

    p3 = sub.add_parser("scan-outputs", help="扫描 outputs 下任务目录，解析基因并核对 tomtop-1 与人类注释")
    p3.add_argument("--outputs-dir", type=Path, default=_MOTIF_ROOT / "outputs")
    p3.add_argument("--table", type=Path, default=_DEFAULT_ALL)
    p3.add_argument(
        "--max-tasks",
        type=int,
        default=0,
        help="最多处理多少个任务目录（0=不限制；目录很多时建议设 50 试跑）",
    )
    p3.set_defaults(func=cmd_scan_outputs)

    args = ap.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
