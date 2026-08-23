#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert h5data/<TASK>/*_chr*.h5 (pandas key=df) to 5-channel (OneHot + PhyloP),
writing h5data/<TASK>/_5ch/<TASK>_chr*_5ch.h5 (h5py: X, Y, Y_ctl).
Only chr-split h5 files are processed; not .101.peak_* style files.

Requires: pyBigWig, numpy, pandas, h5py; supply hg38.phyloP100way.bw yourself.
"""
import os
import sys
import re
import glob
import argparse
import numpy as np
import h5py
import pandas as pd

try:
    import pyBigWig
except ImportError:
    pyBigWig = None


def extract_phylop_for_sample(bw, chrom: str, start: int, end: int,
                              strand: str, seq_len: int) -> np.ndarray:
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    try:
        vals = bw.values(chrom, int(start), int(end))
    except (RuntimeError, Exception):
        vals = [np.nan] * seq_len
    phylop = np.array(vals, dtype=np.float32)
    phylop = np.nan_to_num(phylop, nan=0.0)
    if len(phylop) != seq_len:
        padded = np.zeros(seq_len, dtype=np.float32)
        n = min(len(phylop), seq_len)
        padded[:n] = phylop[:n]
        phylop = padded
    if strand == '-':
        phylop = phylop[::-1].copy()
    return phylop


def onehot_encode_single(seq: str) -> np.ndarray:
    mapper = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}
    seq = seq.upper().replace('T', 'U')
    one_hot = np.zeros((len(seq), 4), dtype=np.float32)
    for i, base in enumerate(seq):
        if base in mapper:
            one_hot[i, mapper[base]] = 1.0
    return one_hot


def _to_signal_array(val, seq_len: int) -> np.ndarray:
    if isinstance(val, (list, np.ndarray)):
        arr = np.array(val, dtype=np.float32).ravel()
    elif isinstance(val, str):
        nums = re.findall(r"\(([^)]+)\)", val)
        out = []
        for x in nums:
            try:
                out.append(float(x.strip()))
            except ValueError:
                continue
        if not out:
            for x in re.findall(r"-?\d+\.?\d*|-?\d*\.\d+", val):
                try:
                    out.append(float(x))
                except ValueError:
                    continue
        arr = np.array(out, dtype=np.float32) if out else np.zeros(seq_len, dtype=np.float32)
    else:
        arr = np.array(val, dtype=np.float32).ravel()
    if len(arr) >= seq_len:
        return arr[:seq_len]
    return np.pad(arr, (0, seq_len - len(arr)), constant_values=0.0)


def process_single_h5(input_path: str, output_path: str, bw_path: str):
    if pyBigWig is None:
        raise RuntimeError("pyBigWig is required: pip install pyBigWig")
    df = pd.read_hdf(input_path, key='df')
    df.columns = ["chr", "chromStart", "chromEnd", "peak", "score",
                  "strand", "signal", "signal_ctl", "seq"]
    df = df[~df['seq'].str.contains('N', case=False)].dropna().reset_index(drop=True)
    n_samples = len(df)
    seq_len = len(df.iloc[0]['seq'])
    n_channels = 5
    bw = pyBigWig.open(bw_path)
    with h5py.File(output_path, 'w') as hf:
        ds_x = hf.create_dataset(
            'X', shape=(n_samples, seq_len, n_channels),
            dtype='float32', chunks=(min(256, n_samples), seq_len, n_channels),
            compression='gzip', compression_opts=4
        )
        ds_y = hf.create_dataset(
            'Y', shape=(n_samples, seq_len),
            dtype='float32', chunks=(min(256, n_samples), seq_len),
            compression='gzip', compression_opts=4
        )
        ds_y_ctl = hf.create_dataset(
            'Y_ctl', shape=(n_samples, seq_len),
            dtype='float32', chunks=(min(256, n_samples), seq_len),
            compression='gzip', compression_opts=4
        )
        for i in range(n_samples):
            row = df.iloc[i]
            onehot = onehot_encode_single(row['seq'])
            phylop = extract_phylop_for_sample(
                bw, row['chr'], row['chromStart'], row['chromEnd'],
                row['strand'], seq_len
            )
            feature = np.concatenate([onehot, phylop[:, np.newaxis]], axis=1)
            ds_x[i] = feature
            ds_y[i] = _to_signal_array(row['signal'], seq_len)
            ds_y_ctl[i] = _to_signal_array(row['signal_ctl'], seq_len)
    bw.close()


def main():
    parser = argparse.ArgumentParser(
        description="One task directory under h5data: write _5ch files to <input_dir>/_5ch/"
    )
    parser.add_argument("--input_dir", type=str, required=True,
                        help="Task directory under h5data, e.g. .../AGGF1_HepG2_ENCSR543TPH")
    parser.add_argument("--bw_path", type=str, required=True,
                        help="Path to hg38.phyloP100way.bw")
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Default: <input_dir>/_5ch")
    parser.add_argument("--skip_existing", action="store_true", default=True,
                        help="Skip if output _5ch.h5 already exists")
    args = parser.parse_args()
    out_dir = args.output_dir or os.path.join(args.input_dir, "_5ch")
    os.makedirs(out_dir, exist_ok=True)
    all_h5 = sorted(glob.glob(os.path.join(args.input_dir, "*.h5")))
    h5_files = [f for f in all_h5 if "_chr" in os.path.basename(f) and "peak" not in os.path.basename(f)]
    for h5_path in h5_files:
        base = os.path.basename(h5_path).replace(".h5", "_5ch.h5")
        output_path = os.path.join(out_dir, base)
        if args.skip_existing and os.path.exists(output_path):
            continue
        process_single_h5(h5_path, output_path, args.bw_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
