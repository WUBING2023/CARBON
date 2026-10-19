#!/usr/bin/env bash
# 多 GPU 并行跑 Motif 全量：与 run_all_tasks.py 默认参数对齐。
#
# 默认走 Python 线程池（scripts/run_all_tasks_pool.py）：常驻队列按槽位填充。
# 未设置 POOL_SLOTS 时：第 0 号 GPU 20 槽，其余每卡 9 槽（与 GPU 列表同序），例如 4 卡 → 20,9,9,9。
# 每槽 CPU 线程数：POOL_THREADS_PER_SLOT（默认 3）→ --threads-per-slot。
# 若 POOL_NO_RESUME=1，传给 pool --no-resume（全量重跑，不跳过已有 modisco_results.h5）。
# 可选：run-parallel 子命令在已安装 parallel 时使用 GNU parallel。
#
# 用法：
#   ./scripts/run_all_tasks_parallel.sh count         # 待跑任务数（--resume）
#   ./scripts/run_all_tasks_parallel.sh list          # 列出待跑任务名
#   ./scripts/run_all_tasks_parallel.sh run           # 前台 Python 多 GPU 池
#   ./scripts/run_all_tasks_parallel.sh start         # nohup 后台同上
#   ./scripts/run_all_tasks_parallel.sh run-parallel  # 需 GNU parallel
#
# 环境变量：
#   NUM_GPU   并发工作线程数 = GPU 数（默认 nvidia-smi -L 行数，失败则 4）；用于生成 --gpus 0,1,...
#   MARCP_ROOT 默认 /data/CARBON_PROJECT/MARCP
#   POOL_SLOTS            可选；未设时默认首卡 20 槽、其余各 9 槽（与 NUM_GPU 一致）
#   POOL_THREADS_PER_SLOT 可选，默认 3；传给 pool --threads-per-slot
#   POOL_NO_RESUME        设为 1 时加 --no-resume（忽略 modisco h5，强制全量重跑）

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MOTIF_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MARCP_ROOT="${MARCP_ROOT:-/data/CARBON_PROJECT/MARCP}"
cd "$MOTIF_ROOT"

export PYTHONPATH="${MARCP_ROOT}:${PYTHONPATH:-}"

if [[ -z "${NUM_GPU:-}" ]]; then
  if command -v nvidia-smi &>/dev/null; then
    NUM_GPU=$(nvidia-smi -L 2>/dev/null | wc -l)
  else
    NUM_GPU=4
  fi
fi
# trim
NUM_GPU="${NUM_GPU// /}"
[[ "$NUM_GPU" =~ ^[0-9]+$ ]] || NUM_GPU=4
if (( NUM_GPU < 1 )); then NUM_GPU=1; fi

LOGDIR="$MOTIF_ROOT/outputs/batch_motif_logs"
mkdir -p "$LOGDIR"

usage() {
  echo "Usage: $0 {count|list|run|start|run-parallel}" >&2
  echo "  count         tasks remaining (run_all_tasks.py --dry-run --resume)" >&2
  echo "  list          print those task names" >&2
  echo "  run           Python pool: one thread per GPU (default)" >&2
  echo "  start         nohup '$0 run' -> $LOGDIR/pool_nohup_<stamp>.log" >&2
  echo "  run-parallel  GNU parallel (must: apt install parallel)" >&2
  echo "Env: NUM_GPU (default $NUM_GPU), GPU_LIST=0,2, POOL_SLOTS (default 20,9,9,...), POOL_THREADS_PER_SLOT (default 3), POOL_NO_RESUME=1, MARCP_ROOT" >&2
}

cmd_count() {
  python3 scripts/run_all_tasks.py --dry-run --resume 2>/dev/null | wc -l
}

cmd_list() {
  python3 scripts/run_all_tasks.py --dry-run --resume 2>/dev/null
}

ensure_parallel() {
  if ! command -v parallel &>/dev/null; then
    echo "ERROR: GNU parallel not found. Install: sudo apt-get install -y parallel" >&2
    exit 1
  fi
}

gpus_csv_from_num() {
  local n="$1"
  local i s=""
  (( n < 1 )) && n=1
  for ((i=0; i<n; i++)); do
    [[ -n "$s" ]] && s+=","
    s+="$i"
  done
  echo "$s"
}

cmd_run() {
  if [[ -n "${GPU_LIST:-}" ]]; then
    GPU_CSV="$GPU_LIST"
  else
    GPU_CSV="$(gpus_csv_from_num "$NUM_GPU")"
  fi
  export MARCP_ROOT

  # 默认槽位：GPU0=20，其余每卡=9（与 GPU 数量一致）
  if [[ -z "${POOL_SLOTS:-}" ]]; then
    POOL_SLOTS="20"
    local _i
    for ((_i = 1; _i < NUM_GPU; _i++)); do
      POOL_SLOTS+=",9"
    done
  fi
  export POOL_THREADS_PER_SLOT="${POOL_THREADS_PER_SLOT:-3}"

  POOL_ARGS=(--slots "$POOL_SLOTS" --threads-per-slot "$POOL_THREADS_PER_SLOT")
  if [[ "${POOL_NO_RESUME:-}" == "1" ]]; then
    POOL_ARGS+=(--no-resume)
  fi
  exec python3 scripts/run_all_tasks_pool.py --gpus "$GPU_CSV" --marcp-root "$MARCP_ROOT" "${POOL_ARGS[@]}"
}

cmd_run_parallel() {
  ensure_parallel
  mapfile -t TASKS < <(python3 scripts/run_all_tasks.py --dry-run --resume 2>/dev/null)
  if [[ ${#TASKS[@]} -eq 0 ]]; then
    echo "No tasks to run (all have outputs/<task>/run_config.json or none eligible)."
    exit 0
  fi
  echo "Tasks: ${#TASKS[@]}, parallel jobs: $NUM_GPU, joblog under $LOGDIR"
  STAMP=$(date +%Y%m%d_%H%M%S)
  JOBLOG="$LOGDIR/parallel_joblog_${STAMP}.txt"
  parallel -j "$NUM_GPU" --joblog "$JOBLOG" --linebuffer \
    'cd '"$MOTIF_ROOT"' && export CUDA_VISIBLE_DEVICES=$(( ({%}-1) % '"$NUM_GPU"' )) && export PYTHONPATH="'"$MARCP_ROOT"':${PYTHONPATH}" && python3 scripts/run_from_ablation_task.py --task {} --all-chromosomes --select_top 10000 --rank_by pred_mean --forward-chunk-rows 512 --batch_size 16 --ig_steps 24 --out outputs/{}' \
    ::: "${TASKS[@]}"
  echo "Joblog: $JOBLOG"
}

cmd_start() {
  STAMP=$(date +%Y%m%d_%H%M%S)
  OUT="$LOGDIR/pool_nohup_${STAMP}.log"
  echo "$(date '+%Y-%m-%d %H:%M:%S') Starting nohup (Python GPU pool); log: $OUT"
  nohup bash -c "
    echo \"\$(date '+%Y-%m-%d %H:%M:%S') [launcher] pool_nohup -> $OUT\"
    cd \"$MOTIF_ROOT\" && exec \"$SCRIPT_DIR/run_all_tasks_parallel.sh\" run
  " >"$OUT" 2>&1 &
  echo "$(date '+%Y-%m-%d %H:%M:%S') PID $!"
}

case "${1:-}" in
  count) cmd_count ;;
  list) cmd_list ;;
  run) cmd_run ;;
  start) cmd_start ;;
  run-parallel) cmd_run_parallel ;;
  *)
    usage
    exit 1
    ;;
esac
