#!/usr/bin/env python3
"""将单个 CISBP 风格 PWM 文本转为最小 MEME 4 文件（内部 ACGU，写出标准 ALPHABET= ACGU）。"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from motif_analysis.cisbp_compare import load_cisbp_pwm_txt  # noqa: E402
from motif_analysis.meme_write import write_meme_file  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cisbp_txt", type=Path)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("-n", "--motif-name", type=str, default="motif1")
    args = ap.parse_args()
    m = load_cisbp_pwm_txt(args.cisbp_txt)
    if m is None:
        raise SystemExit("failed to parse")
    write_meme_file(args.out, [(args.motif_name, m)], alphabet_rna=True, nsites=1)
    print(args.out)


if __name__ == "__main__":
    main()
