#!/usr/bin/env bash
# 串联：① sequence PPM 重导出 queries.meme（--no-skip-existing 覆盖 tfmodisco_export）
#       ② 并行 CISBP-RNA Tomtom（batch_motif_comparisons → tomtom_modisco_vs_cisbp.py，含 -no-ssc 与 CISBP 联表）
#
# 用法（在 可解释性分析/Motif 下）：
#   ./scripts/run_full_motif_recompare_pipeline.sh
#   ./scripts/run_full_motif_recompare_pipeline.sh RBFOX2_K562_ENCSR756CKJ
#   EXPORT_JOBS=8 COMPARE_JOBS=2 ./scripts/run_full_motif_recompare_pipeline.sh TASK_A TASK_B
#
# 环境变量：
#   EXPORT_JOBS   默认 40
#   COMPARE_JOBS  默认 4
#   TOMTOM        默认 ~/.local/share/mamba/envs/meme-suite/bin/tomtom
#   TOP_N         传给 modisco_pattern_rank --top-n，默认 10
#   MEME_PSEUDOCOUNT  默认 0.001
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MOTIF_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$MOTIF_ROOT"

export PYTHONPATH="${MOTIF_ROOT}:${PYTHONPATH:-}"

EXPORT_JOBS="${EXPORT_JOBS:-40}"
COMPARE_JOBS="${COMPARE_JOBS:-4}"
TOMTOM="${TOMTOM:-$HOME/.local/share/mamba/envs/meme-suite/bin/tomtom}"
TOP_N="${TOP_N:-10}"
MEME_PSEUDOCOUNT="${MEME_PSEUDOCOUNT:-0.001}"

TASK_FILTER=()
if [[ $# -gt 0 ]]; then
  for t in "$@"; do
    TASK_FILTER+=(--task "$t")
  done
fi

TOMTOM_ARGS=()
if [[ -f "$TOMTOM" ]]; then
  TOMTOM_ARGS=(--tomtom-binary "$TOMTOM")
else
  echo "WARN: $TOMTOM 不存在，batch_motif_comparisons 将尝试 PATH 中的 tomtom" >&2
fi

echo "[pipeline] step 1: batch_modisco_queries_export.py --no-skip-existing jobs=$EXPORT_JOBS top_n=$TOP_N"
python3 "$MOTIF_ROOT/scripts/batch_modisco_queries_export.py" \
  --no-skip-existing \
  --jobs "$EXPORT_JOBS" \
  --top-n "$TOP_N" \
  --meme-pseudocount "$MEME_PSEUDOCOUNT" \
  "${TASK_FILTER[@]}"

echo "[pipeline] step 2: batch_motif_comparisons.py --infer-gene jobs=$COMPARE_JOBS"
python3 "$MOTIF_ROOT/scripts/batch_motif_comparisons.py" \
  --jobs "$COMPARE_JOBS" \
  --infer-gene \
  "${TOMTOM_ARGS[@]}" \
  "${TASK_FILTER[@]}"

echo "[pipeline] done."
