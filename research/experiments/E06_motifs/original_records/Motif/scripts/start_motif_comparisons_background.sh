#!/usr/bin/env bash
# 后台启动 batch_motif_comparisons.py（CPU 并行 Tomtom），日志带时间戳。
#
# 用法：
#   ./scripts/start_motif_comparisons_background.sh
#   NUM_JOBS=70 THREADS_PER_JOB=2 TOMTOM="$HOME/.../tomtom" ./scripts/start_motif_comparisons_background.sh --infer-gene
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MOTIF_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$MOTIF_ROOT"

export PYTHONPATH="${MARCP_ROOT:-/data/CARBON_PROJECT/MARCP}:${PYTHONPATH:-}"

LOGDIR="$MOTIF_ROOT/outputs/batch_motif_logs"
mkdir -p "$LOGDIR"
STAMP=$(date +%Y%m%d_%H%M%S)
OUT="$LOGDIR/comparisons_nohup_${STAMP}.log"

NUM_JOBS="${NUM_JOBS:-70}"
THREADS_PER_JOB="${THREADS_PER_JOB:-2}"
TOMTOM="${TOMTOM:-$HOME/.local/share/mamba/envs/meme-suite/bin/tomtom}"
EXTRA=("$@")

if [[ ! -f "$TOMTOM" ]]; then
  echo "WARN: $TOMTOM not found; batch script will try PATH. Set TOMTOM=..." >&2
  TOMTOM_ARGS=()
else
  TOMTOM_ARGS=(--tomtom-binary "$TOMTOM")
fi

echo "nohup log: $OUT"
echo "jobs=$NUM_JOBS threads_per_job=$THREADS_PER_JOB tomtom=${TOMTOM:-PATH}"

nohup python3 "$MOTIF_ROOT/scripts/batch_motif_comparisons.py" \
  --jobs "$NUM_JOBS" \
  --threads-per-job "$THREADS_PER_JOB" \
  "${TOMTOM_ARGS[@]}" \
  "${EXTRA[@]}" \
  >>"$OUT" 2>&1 &

echo "PID $!"
