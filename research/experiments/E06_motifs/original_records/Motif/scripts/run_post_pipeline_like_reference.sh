#!/usr/bin/env bash
# 将各任务 outputs/<TASK>/ 补全到与 RBFOX2_K562_ENCSR756CKJ 同类的产物：
#   1) tfmodisco_export（queries.meme 等）— batch_modisco_queries_export.py
#   2) tomtom_cisbp — batch_motif_comparisons.py（子进程调用 tomtom_modisco_vs_cisbp.py，RNA 用 -no-ssc）
#   3) threefold_validation/ — batch_threefold_validation_all.py（用 GPU，勿与大规模 run_pipeline 同开）
#   4) figures/ + visualization/ — batch_figures_and_viz_all.py
#
# 与 run_full_motif_recompare_pipeline.sh 的分工：
#   - 默认（不设 FORCE_*）：跳过已有 queries.meme、已有 tomtom.tsv（增量补全）。
#   - 仅重导出 + 重 CISBP Tomtom、不跑 threefold/fig：用 run_full_motif_recompare_pipeline.sh。
#   - 在本脚本内全量重刷 queries/Tomtom：FORCE_EXPORT=1 FORCE_TOMTOM=1（见下）。
#
# 环境变量（可选全量重刷，与 PIPELINE.md 阶段 B / CISBP 一致）：
#   FORCE_EXPORT=1   → batch_modisco_queries_export 加 --no-skip-existing（覆盖已有 MEME）
#   FORCE_TOMTOM=1   → batch_motif_comparisons 不加 --skip-existing（重跑已有 tomtom）
#
# 用法（Motif 根目录）：
#   export MARCP_ROOT=/data/CARBON_PROJECT/MARCP
#   nohup ./scripts/run_post_pipeline_like_reference.sh >> outputs/batch_motif_logs/post_pipeline_nohup_$(date +%Y%m%d_%H%M%S).log 2>&1 &
#
set -uo pipefail
# 各步可能部分失败，不因 set -e 中断后续步骤

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MOTIF_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$MOTIF_ROOT"

export PYTHONPATH="${MARCP_ROOT:-/data/CARBON_PROJECT/MARCP}:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

TOMTOM="${TOMTOM:-$HOME/.local/share/mamba/envs/meme-suite/bin/tomtom}"
JOBS_QUERIES="${JOBS_QUERIES:-40}"
JOBS_COMP="${JOBS_COMP:-6}"
JOBS_THREE="${JOBS_THREE:-2}"
JOBS_FIG="${JOBS_FIG:-8}"
DISABLE_LOCK="${DISABLE_LOCK:-0}"
FORCE_EXPORT="${FORCE_EXPORT:-0}"
FORCE_TOMTOM="${FORCE_TOMTOM:-0}"

ts() { date '+%Y-%m-%d %H:%M:%S'; }

echo "$(ts) [post] MOTIF_ROOT=$MOTIF_ROOT"
echo "$(ts) [post] TOMTOM=$TOMTOM JOBS queries=$JOBS_QUERIES comp=$JOBS_COMP threefold=$JOBS_THREE fig=$JOBS_FIG"
echo "$(ts) [post] FORCE_EXPORT=$FORCE_EXPORT FORCE_TOMTOM=$FORCE_TOMTOM"

# 防止误启动多份编排同时写 outputs（可通过 DISABLE_LOCK=1 关闭）
if [[ "$DISABLE_LOCK" != "1" ]]; then
  mkdir -p "$MOTIF_ROOT/outputs/batch_motif_logs"
  LOCKFILE="$MOTIF_ROOT/outputs/batch_motif_logs/post_pipeline.lock"
  exec 200>"$LOCKFILE"
  if ! flock -n 200; then
    echo "$(ts) [post] Another run_post_pipeline_like_reference.sh is already running (lock: $LOCKFILE)." >&2
    echo "$(ts) [post] If you really want parallel runs, export DISABLE_LOCK=1 (NOT recommended)." >&2
    exit 2
  fi
fi

if [[ ! -f "$TOMTOM" ]]; then
  echo "$(ts) [post] WARN: tomtom not at $TOMTOM — batch_motif_comparisons will try PATH" >&2
  TOMTOM_ARGS=()
else
  TOMTOM_ARGS=(--tomtom-binary "$TOMTOM")
fi

ERR=0
QUERIES_EXTRA=()
if [[ "$FORCE_EXPORT" == "1" ]]; then
  QUERIES_EXTRA+=(--no-skip-existing)
fi

echo "$(ts) [post] === 1/4 queries export (modisco_pattern_rank) ==="
python3 scripts/batch_modisco_queries_export.py --jobs "$JOBS_QUERIES" "${QUERIES_EXTRA[@]}" || ERR=1

echo "$(ts) [post] === 2/4 CISBP tomtom ==="
COMP_ARGS=(
  python3 scripts/batch_motif_comparisons.py
  --jobs "$JOBS_COMP"
  "${TOMTOM_ARGS[@]}"
  --infer-gene
)
if [[ "$FORCE_TOMTOM" != "1" ]]; then
  COMP_ARGS+=(--skip-existing)
fi
"${COMP_ARGS[@]}" || ERR=1

echo "$(ts) [post] === 3/4 threefold validation (GPU) ==="
python3 scripts/batch_threefold_validation_all.py --jobs "$JOBS_THREE" || ERR=1

echo "$(ts) [post] === 4/4 figures + visualization HTML ==="
python3 scripts/batch_figures_and_viz_all.py --jobs "$JOBS_FIG" || ERR=1

echo "$(ts) [post] ALL STEPS FINISHED (ERR=$ERR, 0=all subprocesses reported success)"
exit "$ERR"
