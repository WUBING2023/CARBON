#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert a 17-column eCLIP training BED file to pandas HDF5 for per-chr splitting or preprocess_5ch.

Output columns:
    chr, chromStart, chromEnd, peak, score, strand, signal, signal_ctl, seq

Default HDF5 key is **df** (matches carbon_train.py / MARCP-style pipelines).
"""
import argparse
import os

import pandas as pd


def bed_to_h5(bed_file: str, h5_file: str, key: str = "df"):
    df = pd.read_table(bed_file, sep="\t", header=None)
    if df.shape[1] != 17:
        raise ValueError(f"Expected 17 columns in BED, got {df.shape[1]}.")
    df.columns = [
        "chrom_peak", "start_peak", "stop_peak",
        "name_peak", "score_peak", "strand_peak",
        "max_peak", "sum_peak",
        "signal_CL", "leftPadding", "rightPadding",
        "start_CL", "stop_CL",
        "signal", "signal_ctl", "signal_ctl_CL",
        "seq",
    ]
    bed_fa = pd.DataFrame({
        "chr": df["chrom_peak"],
        "chromStart": df["start_peak"],
        "chromEnd": df["stop_peak"],
        "peak": df["name_peak"],
        "score": df["score_peak"],
        "strand": df["strand_peak"],
        "signal": df["signal"],
        "signal_ctl": df["signal_ctl"],
        "seq": df["seq"],
    })
    bed_fa.to_hdf(
        h5_file,
        key=key,
        mode="w",
        format="table",
        complevel=5,
        complib="zlib",
    )
    print(f"Saved: {h5_file}  (rows={len(bed_fa)}, cols={bed_fa.shape[1]})")


def main():
    parser = argparse.ArgumentParser(
        description="Convert 17-column BED to h5 with 9 columns (key default: df)."
    )
    parser.add_argument("-i", "--input", required=True, help="Input BED path")
    parser.add_argument("-o", "--output", required=True, help="Output HDF5 path")
    parser.add_argument("-k", "--key", default="df", help="HDF5 key (default df, matches training script)")
    args = parser.parse_args()
    bed_to_h5(os.path.abspath(args.input), os.path.abspath(args.output), key=args.key)


if __name__ == "__main__":
    main()
