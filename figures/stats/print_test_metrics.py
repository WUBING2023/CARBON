#!/usr/bin/env python3
"""Print carbon_train's test_metrics_summary.txt (or same format) to stdout."""
from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("summary", type=Path, nargs="?", default=Path("test_metrics_summary.txt"))
    args = p.parse_args()
    path = args.summary
    if not path.is_file():
        raise SystemExit(f"Not found: {path}")
    print(path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
