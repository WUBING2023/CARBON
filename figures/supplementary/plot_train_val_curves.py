#!/usr/bin/env python3
"""
Supplementary-style figure: plot loss and Pearson vs epoch from train_metric.txt /
val_metric.txt (tab-separated, header loss\\tpearson) in a training output directory.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--train-dir", type=Path, required=True,
                   help="Directory containing train_metric.txt and val_metric.txt")
    p.add_argument("--out", type=Path, default=None,
                   help="Output PNG (default: <train-dir>/train_val_curves.png)")
    args = p.parse_args()
    d = args.train_dir
    tr = d / "train_metric.txt"
    va = d / "val_metric.txt"
    if not tr.is_file() or not va.is_file():
        raise SystemExit(f"Required: {tr} and {va}")
    tdf = pd.read_csv(tr, sep="\t")
    vdf = pd.read_csv(va, sep="\t")
    ep = range(1, len(tdf) + 1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(ep, tdf["loss"], label="train")
    axes[0].plot(ep, vdf["loss"], label="val")
    axes[0].set_xlabel("Epoch")
    axes[0].set_ylabel("Loss")
    axes[0].legend()
    axes[0].set_title("Loss")
    axes[1].plot(ep, tdf["pearson"], label="train")
    axes[1].plot(ep, vdf["pearson"], label="val")
    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Pearson")
    axes[1].legend()
    axes[1].set_title("Pearson")
    fig.tight_layout()
    out = args.out or (d / "train_val_curves.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
