#!/usr/bin/env bash
# 确定性回归测试：固定种子 + 相同数据/超参，跑短训练并生成可对比的签名文件。
# 用法：在 MARCP 项目根目录执行  bash scripts/verify_logic_regression.sh
# 然后用基准版本（如 消融实验/MARCP）用相同参数跑一遍，对比两份 regression_signature.txt。

set -e
MARCP_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$MARCP_ROOT"

OUT_DIR="${MARCP_ROOT}/exampledata/out_regression_refactor"
TRAIN="${MARCP_ROOT}/exampledata/AGGF1_HepG2_ENCSR543TPH_chr1_5ch.h5"
VAL="${MARCP_ROOT}/exampledata/AGGF1_HepG2_ENCSR543TPH_chr2_5ch.h5"
TEST="${MARCP_ROOT}/exampledata/AGGF1_HepG2_ENCSR543TPH_chr3_5ch.h5"

echo "MARCP regression test (refactored version)"
echo "  output dir: $OUT_DIR"
echo "  epochs: 2 (fixed seed in marcp_train.py)"

python marcp_train.py \
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
echo "Compare this file with the same from 消融实验/MARCP (same args) to verify logic equivalence."
