"""Recalculate CARBON Figure 2c from the centered sums in Data S40."""

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import fmean


MODELS = ("CARBON", "iDeepB", "RBPNet")
DIRECTIONS = {("K562", "HepG2"), ("HepG2", "K562")}


def calculate(path):
    tasks = defaultdict(dict)
    with path.open(encoding="utf-8-sig", newline="") as table:
        for row in csv.DictReader(table):
            key = (row["rbp"], row["source_cell"], row["target_cell"])
            model = row["model"]
            if model not in MODELS or model in tasks[key]:
                raise ValueError(f"Unexpected or duplicate model {model} in {key}")
            tasks[key][model] = row

    scores = defaultdict(lambda: defaultdict(list))
    directions = defaultdict(set)
    windows = 0
    for (protein, source, target), rows in tasks.items():
        if rows.keys() != set(MODELS):
            raise ValueError(f"Incomplete comparison for {protein}, {source} to {target}")
        identifiers = {
            (r["source_task"], r["target_task"], r["cohort_sha256"], int(r["n_windows"]))
            for r in rows.values()
        }
        if len(identifiers) != 1:
            raise ValueError(f"Model inputs disagree for {protein}, {source} to {target}")
        window_count = int(rows["CARBON"]["n_windows"])
        if window_count <= 0:
            raise ValueError(f"Empty task for {protein}, {source} to {target}")
        windows += window_count
        directions[protein].add((source, target))
        for model, row in rows.items():
            denominator = math.sqrt(float(row["sum_pp"]) * float(row["sum_qq"]))
            if not denominator > 0:
                raise ValueError(f"Undefined PCC for {protein}, {source} to {target}, {model}")
            pcc = float(row["sum_pq"]) / denominator
            reported = float(row["pcc_unit_mass_task_pooled"])
            if not math.isclose(pcc, reported, rel_tol=0, abs_tol=1e-12):
                raise ValueError(f"Stored PCC differs for {protein}, {source} to {target}, {model}")
            scores[protein][model].append(pcc)
    if not tasks or any(pair != DIRECTIONS for pair in directions.values()):
        raise ValueError("Each protein must have both transfer directions")

    by_protein = {
        protein: {model: fmean(values) for model, values in models.items()}
        for protein, models in scores.items()
    }
    means = {model: fmean(row[model] for row in by_protein.values()) for model in MODELS}
    contrasts = {}
    for comparator in ("iDeepB", "RBPNet"):
        differences = [row["CARBON"] - row[comparator] for row in by_protein.values()]
        contrasts[comparator] = {
            "mean_paired_difference": fmean(differences),
            "proteins_with_positive_difference": sum(value > 0 for value in differences),
        }
    return {
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "proteins": len(by_protein),
        "directed_tasks": len(tasks),
        "target_windows": windows,
        "protein_equal_mean_pcc": means,
        "carbon_minus_comparator": contrasts,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", type=Path, help="Path to cross_model_task_metrics.csv")
    print(json.dumps(calculate(parser.parse_args().table), indent=2))
