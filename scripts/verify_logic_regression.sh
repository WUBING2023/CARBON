#!/usr/bin/env bash
# Deterministic regression: fixed seed + same data/hyperparameters, short run, comparable signature file.
# Usage: from CARBON repo root:  bash scripts/verify_logic_regression.sh
# Then run the baseline (e.g. ablation/MARCP) with the same args and diff the two regression_signature.txt files.

set -e
CARBON_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$CARBON_ROOT"

OUT_DIR="${CARBON_ROOT}/exampledata/out_regression_refactor"
TRAIN="${CARBON_ROOT}/exampledata/AGGF1_HepG2_ENCSR543TPH_chr1_5ch.h5"
VAL="${CARBON_ROOT}/exampledata/AGGF1_HepG2_ENCSR543TPH_chr2_5ch.h5"
TEST="${CARBON_ROOT}/exampledata/AGGF1_HepG2_ENCSR543TPH_chr3_5ch.h5"

echo "CARBON regression test (refactored version)"
echo "  output dir: $OUT_DIR"
echo "  epochs: 2 (fixed seed in carbon_train.py)"

python3 carbon_train.py \
  -train "$TRAIN" \
  -validation "$VAL" \
  -test "$TEST" \
  --output "$OUT_DIR" \
  --modelFramework Treat \
  --lossM CombinedProfileLoss \
  --input_channels 5 \
  --epochs 2 \
  --bsz 256

SIG="${OUT_DIR}/regression_signature.txt"
{
  echo "=== train_metric.txt ==="
  cat "${OUT_DIR}/train_metric.txt" 2>/dev/null || true
  echo ""
  echo "=== test_metric.txt ==="
  cat "${OUT_DIR}/test_metric.txt" 2>/dev/null || true
  echo ""
  echo "=== test_metrics_summary.txt ==="
  cat "${OUT_DIR}/test_metrics_summary.txt" 2>/dev/null || true
} > "$SIG"

echo "Regression signature written to: $SIG"
echo "Compare this file with the baseline MARCP run (same args) to verify logic equivalence."
