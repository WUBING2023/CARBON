#!/usr/bin/env python3
"""单次序列 RBPNet profile 预测，向 stdout 打印 JSON 数组（供主绘图脚本子进程调用）。"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

RBP_ROOT = Path(__file__).resolve().parents[1] / "比较实验/RBPNet"


def _prepend_tf_cuda_libs() -> None:
    pyver = f"python{sys.version_info.major}.{sys.version_info.minor}"
    site_pkg = Path(sys.prefix) / "lib" / pyver / "site-packages"
    candidates = [
        site_pkg / "nvidia" / "cudnn" / "lib",
        site_pkg / "nvidia" / "cublas" / "lib",
    ]
    existing = [str(p) for p in candidates if p.exists()]
    if not existing:
        return
    ld_old = os.environ.get("LD_LIBRARY_PATH", "")
    parts = [p for p in ld_old.split(":") if p]
    for p in reversed(existing):
        if p not in parts:
            parts.insert(0, p)
    os.environ["LD_LIBRARY_PATH"] = ":".join(parts)


def main() -> None:
    if len(sys.argv) != 3:
        print("usage: rbpnet_predict_profile_once.py <DNA_SEQ> <model.h5>", file=sys.stderr)
        sys.exit(2)
    seq = sys.argv[1].strip().upper().replace("U", "T")
    wpath = Path(sys.argv[2]).resolve()
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
    _prepend_tf_cuda_libs()
    sys.path.insert(0, str(RBP_ROOT))
    from rbpnet.load import load_model  # noqa: E402

    model = load_model(str(wpath))
    out = model.predict_from_sequence(seq)
    task = model.tasks[0]
    key = f"{task}_profile"
    if key not in out:
        raise KeyError(f"missing {key}, keys={list(out.keys())}")
    prof = out[key]
    if hasattr(prof, "numpy"):
        prof = prof.numpy()
    arr = np.asarray(prof, dtype=np.float32).reshape(-1)
    print(json.dumps(arr.tolist()))


if __name__ == "__main__":
    main()
