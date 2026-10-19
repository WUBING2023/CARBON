#!/usr/bin/env python3
"""
生成最小可复现的 5ch h5 + MARCP state_dict，用于在无外部数据时跑通全流程。

  python scripts/prepare_fixtures.py --out outputs/fixtures
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import h5py
import numpy as np
import torch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=str, default=str(ROOT / "outputs" / "fixtures"))
    ap.add_argument("--n", type=int, default=48, help="number of sequences")
    ap.add_argument("--length", type=int, default=101)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    n, L = args.n, args.length
    X = np.zeros((n, L, 5), dtype=np.float32)
    for i in range(n):
        for j in range(L):
            b = int(rng.integers(0, 4))
            X[i, j, b] = 1.0
            X[i, j, 4] = float(rng.normal(0.0, 0.2))
    # pseudo eCLIP-like profiles (non-negative)
    Y = rng.random((n, L)).astype(np.float32) * 5.0
    Y_ctl = rng.random((n, L)).astype(np.float32) * 3.0

    h5_path = out / "minimal_5ch.h5"
    with h5py.File(h5_path, "w") as hf:
        hf.create_dataset("X", data=X, compression="gzip")
        hf.create_dataset("Y", data=Y, compression="gzip")
        hf.create_dataset("Y_ctl", data=Y_ctl, compression="gzip")

    sys.path.insert(0, str(Path("/data/CARBON_PROJECT/MARCP")))
    from MARCP.models.marcp_model import MARCP

    torch.manual_seed(args.seed)
    model = MARCP(dropOut=0.0, control=False, in_channels=5)
    ckpt_path = out / "tiny_marcp.pth"
    torch.save(model.state_dict(), ckpt_path)

    meta = {
        "h5": str(h5_path.resolve()),
        "checkpoint": str(ckpt_path.resolve()),
        "shapes": {"X": list(X.shape), "Y": list(Y.shape)},
    }
    import json

    (out / "fixture_manifest.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("Wrote:", h5_path)
    print("Wrote:", ckpt_path)
    print("Manifest:", out / "fixture_manifest.json")


if __name__ == "__main__":
    main()
