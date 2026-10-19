#!/usr/bin/env python3
"""
从已下载的 CISBP pwms/ 抽取若干条写入 reference_data/motif_libraries/meme_rna/aggf1_cisbp_refs.meme，
作为「MEME 格式文献子集」示例（与全库 CISBP 列互补）。
"""
from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from motif_analysis.cisbp_compare import load_cisbp_pwm_txt  # noqa: E402

# 与 AGGF1 任务 cisbp_compare 中常见 top hit 相关的 CISBP ID（可改）
DEFAULT_IDS = ("M00590_2.00", "M00624_2.00", "M00787_2.00", "M00581_2.00", "M07652_2.00")


def _block(name: str, m) -> str:
    w = m.shape[0]
    lines = [
        f"MOTIF {name}",
        f"letter-probability matrix: alength= 4 w= {w} nsites= 1",
    ]
    for row in m:
        lines.append(" ".join(f"{float(x):.8f}" for x in row))
    return "\n".join(lines)


def main() -> None:
    pwms_dir = _ROOT / "reference_data/CISBP-RNA/pwms"
    out = _ROOT / "reference_data/motif_libraries/meme_rna/aggf1_cisbp_refs.meme"
    chunks = ["MEME version 4", "", "ALPHABET= ACGU", ""]
    for stem in DEFAULT_IDS:
        p = pwms_dir / f"{stem}.txt"
        if not p.is_file():
            print("[skip] missing", p, file=sys.stderr)
            continue
        m = load_cisbp_pwm_txt(p)
        if m is None:
            continue
        chunks.append(_block(stem.replace(".", "_"), m))
        chunks.append("")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(chunks).strip() + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
