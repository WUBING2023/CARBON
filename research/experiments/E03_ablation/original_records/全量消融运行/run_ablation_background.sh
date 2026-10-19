#!/bin/bash
# =============================================================================
# MARCP 全量消融：后台启动训练（调用 run_ablation_full_scale.sh，使用 消融实验/MARCP/marcp_train.py）
# 用法：bash run_ablation_background.sh
# 查询进度：bash query_ablation_progress.sh  停止：bash stop_ablation_training.sh
# =============================================================================
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
OUT_ROOT="${SCRIPT_DIR}"
PID_FILE="${OUT_ROOT}/.run_ablation.pid"
LOG_FILE="${OUT_ROOT}/run_ablation.log"
# 固定本批 11 条（batch11_tasks.txt，含 SERBP1 六消融 + 另 5 任务）：1=启用；0=恢复扫描全部数据集×消融
USE_FIXED_BATCH11="${USE_FIXED_BATCH11:-1}"
if [ "$USE_FIXED_BATCH11" = "1" ] && [ -f "${SCRIPT_DIR}/batch11_tasks.txt" ]; then
    export ABLATION_FIXED_LIST_FILE="${SCRIPT_DIR}/batch11_tasks.txt"
elif [ "$USE_FIXED_BATCH11" = "1" ] && [ ! -f "${SCRIPT_DIR}/batch11_tasks.txt" ]; then
    echo "注意: USE_FIXED_BATCH11=1 但未找到 ${SCRIPT_DIR}/batch11_tasks.txt，将使用全量数据扫描。" >&2
fi

if [ -z "${N_GPU}" ]; then
    N_GPU=$(nvidia-smi -L 2>/dev/null | wc -l)
    [ "$N_GPU" -lt 1 ] && N_GPU=1
fi
# 默认 GPU 槽位：GPU0=15，其它 GPU=20（即 15 20 20 20，共 75 槽）
if [ -z "${SLOTS_GPU0}" ]; then
    SLOTS_GPU0=15
fi
if [ -z "${SLOTS_OTHER}" ]; then
    SLOTS_OTHER=20
fi
export N_GPU
export SLOTS_GPU0
export SLOTS_OTHER
[ -n "${BATCH_SIZE}" ] && export BATCH_SIZE
[ -n "${TASK_EXCLUDE}" ] && export TASK_EXCLUDE
[ -n "${FORCE_RERUN}" ] && export FORCE_RERUN

# 限制每进程的数值库/BLAS/OpenMP 线程数，避免多进程并发导致线程过度并发
: "${OMP_NUM_THREADS:=1}"
: "${MKL_NUM_THREADS:=1}"
: "${OPENBLAS_NUM_THREADS:=1}"
: "${NUMEXPR_NUM_THREADS:=1}"
export OMP_NUM_THREADS MKL_NUM_THREADS OPENBLAS_NUM_THREADS NUMEXPR_NUM_THREADS

if [ -f "$PID_FILE" ]; then
    pid=$(cat "$PID_FILE")
    if kill -0 "$pid" 2>/dev/null; then
        echo "全量消融训练已在运行 (PID=$pid)。查看进度: bash query_ablation_progress.sh  退出: bash stop_ablation_training.sh"
        exit 1
    fi
    rm -f "$PID_FILE"
fi

echo "在后台启动全量消融训练 (N_GPU=$N_GPU)，日志: $LOG_FILE"
nohup env N_GPU="$N_GPU" SLOTS_GPU0="${SLOTS_GPU0:-}" SLOTS_OTHER="${SLOTS_OTHER:-}" BATCH_SIZE="${BATCH_SIZE:-}" TASK_EXCLUDE="${TASK_EXCLUDE:-}" FORCE_RERUN="${FORCE_RERUN:-}" ABLATION_FIXED_LIST_FILE="${ABLATION_FIXED_LIST_FILE:-}" OMP_NUM_THREADS="$OMP_NUM_THREADS" MKL_NUM_THREADS="$MKL_NUM_THREADS" OPENBLAS_NUM_THREADS="$OPENBLAS_NUM_THREADS" NUMEXPR_NUM_THREADS="$NUMEXPR_NUM_THREADS" bash ./run_ablation_full_scale.sh >> "$LOG_FILE" 2>&1 &
echo $! > "$PID_FILE"
echo "已启动 PID=$(cat "$PID_FILE")。查看进度: bash query_ablation_progress.sh  退出: bash stop_ablation_training.sh"
