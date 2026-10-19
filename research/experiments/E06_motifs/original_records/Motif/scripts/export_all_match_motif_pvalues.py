#!/usr/bin/env python3
"""
扫描 outputs/<TASK>/tomtom_cisbp/tomtom_out/tomtom_with_rbp.tsv，
收集 cisbp_priority == match_task_gene 的所有行，按 Tomtom p-value 升序写出 CSV。

  python3 scripts/export_all_match_motif_pvalues.py
  python3 scripts/export_all_match_motif_pvalues.py --out reference_data/motif_p_value/match_motif_p_sorted.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SKIP_DIRS = frozenset({"batch_motif_logs", "fixtures", "demo", "workspace", ".git"})
MATCH_PRIORITY = "match_task_gene"


def _parse_float(s: str) -> float:
    s = (s or "").strip()
    if not s or s.lower() in ("nan", "na", ""):
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def _p_from_row(row: dict[str, str]) -> float:
    for k in ("_p_num", "p-value", "P-value", "p_value", "_p"):
        if k in row and row[k] is not None:
            v = _parse_float(str(row[k]))
            if v == v:  # not nan
                return v
    return float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outputs-root", type=Path, default=ROOT / "outputs")
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "reference_data" / "motif_p_value" / "match_motif_p_sorted.csv",
    )
    args = ap.parse_args()

    out_root = args.outputs_root.resolve()
    out_path = args.out.resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows_out: list[dict[str, str | float]] = []
    n_files = 0

    for child in sorted(out_root.iterdir()):
        if not child.is_dir() or child.name in SKIP_DIRS or child.name.startswith("."):
            continue
        tsv = child / "tomtom_cisbp" / "tomtom_out" / "tomtom_with_rbp.tsv"
        if not tsv.is_file():
            continue
        n_files += 1
        with open(tsv, encoding="utf-8", errors="replace", newline="") as f:
            rd = csv.DictReader(f, delimiter="\t")
            fieldnames = rd.fieldnames or []
            pri_col = next((h for h in fieldnames if h and str(h).strip() == "cisbp_priority"), None)
            if pri_col is None:
                print(f"[warn] skip {tsv}: no cisbp_priority column", file=sys.stderr)
                continue

            qcol = "Query_ID" if "Query_ID" in fieldnames else None
            tcol = "Target_ID" if "Target_ID" in fieldnames else None

            for row in rd:
                if (row.get(pri_col) or "").strip() != MATCH_PRIORITY:
                    continue
                p = _p_from_row(row)
                rows_out.append(
                    {
                        "task": child.name,
                        "Query_ID": (row.get(qcol) or "") if qcol else "",
                        "Target_ID": (row.get(tcol) or "") if tcol else "",
                        "p_value": p,
                        "E-value": row.get("E-value", row.get("E-value ", "")) or "",
                        "q-value": row.get("q-value", row.get("Q-value", "")) or "",
                        "cisbp_RBP_Name": row.get("cisbp_RBP_Name", "") or "",
                        "cisbp_Motif_ID": row.get("cisbp_Motif_ID", "") or "",
                    }
                )

    def _sort_key(r: dict) -> tuple[float, str, str, str]:
        p = r["p_value"]
        pv = float(p) if isinstance(p, float) and p == p else 1e300
        return (pv, str(r.get("task", "")), str(r.get("Query_ID", "")), str(r.get("Target_ID", "")))

    rows_out.sort(key=_sort_key)

    fieldnames = [
        "task",
        "Query_ID",
        "Target_ID",
        "p_value",
        "E-value",
        "q-value",
        "cisbp_RBP_Name",
        "cisbp_Motif_ID",
    ]
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows_out:
            line = {k: r[k] for k in fieldnames}
            if isinstance(line["p_value"], float):
                line["p_value"] = f"{line['p_value']:.17g}"
            w.writerow(line)

    print(
        f"[export_match_p] tasks_with_tsv={n_files} match_rows={len(rows_out)} -> {out_path}",
        flush=True,
    )


if __name__ == "__main__":
    main()
