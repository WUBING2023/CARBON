#!/bin/bash
# =============================================================================
# 本批 STATE 任务：全部并行（每个任务一个后台进程，最后对每个 $! 显式 wait，不用 wait -n）
#
# 用法:
#   bash stop_ablation_training.sh
#   bash run_ablation_batch_one_by_one.sh
#
# 环境变量:
#   STATE_FILE, FORCE_RERUN, GPU_ROTATE, BATCH_SIZE（同前）
#
# 任务行来源（按优先级）:
#   1) 第一个参数：任意任务列表文件（每行格式同 STATE 的 jobs 后一行）
#   2) $STATE_FILE（若存在）
#   3) ${OUT_ROOT}/.run_ablation.state
#   4) ${OUT_ROOT}/batch11_tasks.txt
#   5) 若皆无：自动运行 generate_batch11_task_lines.sh 生成 batch11_tasks.txt（不依赖 STATE）
# =============================================================================
set +e
BASE="/data/CARBON_PROJECT"
OUT_ROOT="${BASE}/消融实验/全量消融运行"
MODEL_12="${BASE}/MARCP"
STATE_FILE="${STATE_FILE:-${OUT_ROOT}/.run_ablation.state}"
EPOCHS=100
BATCH_SIZE="${BATCH_SIZE:-256}"
DATA_WORKERS_PER_JOB="${DATA_WORKERS_PER_JOB:-0}"
RESOURCE_FRAC_PER_JOB="${RESOURCE_FRAC_PER_JOB:-0.9}"
FORCE_RERUN="${FORCE_RERUN:-0}"
GPU_ROTATE="${GPU_ROTATE:-1}"

VENV_PYTHON=""
[ -x "${OUT_ROOT}/venv/bin/python3" ] && VENV_PYTHON="${OUT_ROOT}/venv/bin/python3"

COMPLETED_LOG="${OUT_ROOT}/.run_ablation.completed.standalone"
SUMMARY_LOCK="${OUT_ROOT}/.summary.standalone.lock"
mkdir -p "${OUT_ROOT}/logs" "${OUT_ROOT}/模型权重/MARCP" "${OUT_ROOT}/各个具体实验的指标/MARCP"

TASK_SOURCE=""
if [ -n "${1:-}" ] && [ -f "$1" ]; then
  TASK_SOURCE="$1"
elif [ -f "${STATE_FILE}" ]; then
  TASK_SOURCE="${STATE_FILE}"
elif [ -f "${OUT_ROOT}/.run_ablation.state" ]; then
  TASK_SOURCE="${OUT_ROOT}/.run_ablation.state"
elif [ -f "${OUT_ROOT}/batch11_tasks.txt" ]; then
  TASK_SOURCE="${OUT_ROOT}/batch11_tasks.txt"
else
  echo "未找到任务列表，正在根据本批默认 11 项生成 ${OUT_ROOT}/batch11_tasks.txt …"
  if ! bash "${OUT_ROOT}/generate_batch11_task_lines.sh" > "${OUT_ROOT}/batch11_tasks.txt"; then
    echo "生成失败。请手动指定: $0 /path/to/task_lines.txt"
    exit 1
  fi
  TASK_SOURCE="${OUT_ROOT}/batch11_tasks.txt"
fi

TASK_LINES=()
if grep -q '^jobs$' "$TASK_SOURCE" 2>/dev/null; then
  in_jobs=0
  while IFS= read -r line || [ -n "$line" ]; do
    case "$line" in
      jobs) in_jobs=1; continue ;;
      start*|total*|n_gpu*|total_slots*|gpu_slots_desc*) continue ;;
    esac
    [ "$in_jobs" = 1 ] && [ -n "$line" ] && TASK_LINES+=("$line")
  done < "$TASK_SOURCE"
else
  while IFS= read -r line || [ -n "$line" ]; do
    case "$line" in
      \#*|'') continue ;;
    esac
    [ -n "$line" ] && TASK_LINES+=("$line")
  done < "$TASK_SOURCE"
fi

n=${#TASK_LINES[@]}
if [ "$n" -lt 1 ]; then
  echo "任务列表为空: $TASK_SOURCE"
  exit 1
fi

echo "========== 本批全并行补跑（${n} 路同时，末尾 wait 各 PID）=========="
echo "来源: $TASK_SOURCE  FORCE_RERUN=$FORCE_RERUN  GPU_ROTATE=$GPU_ROTATE  BATCH_SIZE=$BATCH_SIZE"
echo "================================================================"

run_single_task() {
  local line="$1" run_idx="$2" total_n="$3" gpu_id="$4"
  local task train_5ch val_5ch test_5ch ablation_short ablation_flag idx
  IFS='|' read -r task train_5ch val_5ch test_5ch ablation_short ablation_flag idx <<< "$line"
  if [ -z "$task" ] || [ -z "$train_5ch" ]; then
    echo "[SKIP] 第 $run_idx 行解析失败"
    return 0
  fi

  local weight_dir="${OUT_ROOT}/模型权重/MARCP/${task}/${ablation_short}"
  local metric_dir="${OUT_ROOT}/各个具体实验的指标/MARCP/${task}/${ablation_short}"
  if [ "$FORCE_RERUN" != "1" ] && [ -f "${weight_dir}/model_params.pth" ] && [ -f "${metric_dir}/test_metrics_summary.txt" ]; then
    echo "[$run_idx/$total_n] [SKIP] $task $ablation_short 已完成"
    return 0
  fi

  local base_out="${OUT_ROOT}/workspace/MARCP/${task}_${ablation_short}_standalone"
  mkdir -p "$base_out" "$weight_dir" "$metric_dir"

  echo "========== [$run_idx/$total_n] 启动 $task | $ablation_short | GPU $gpu_id =========="

  (
    cd "${MODEL_12}" || exit 1
    export CUDA_VISIBLE_DEVICES="$gpu_id"
    ${VENV_PYTHON:-python3} -u marcp_train.py \
      --lr 0.000785 \
      --CUDA "$gpu_id" \
      --signalMax 0 \
      --seqLength 101 \
      --bsz "${BATCH_SIZE}" \
      --encode OneHot \
      --modelFramework Treat \
      --epochs ${EPOCHS} \
      --max_stopping_step 9 \
      --input_channels 5 \
      --num_workers "${DATA_WORKERS_PER_JOB}" \
      --prefetch_factor "${PREFETCH_FACTOR:-4}" \
      --resource_fraction "${RESOURCE_FRAC_PER_JOB}" \
      -train "$train_5ch" \
      -validation "$val_5ch" \
      -test "$test_5ch" \
      --task "$task" \
      --lossM CombinedProfileLoss \
      --alpha 1.42 --beta 0.54 --gamma 1.80 \
      --weight_decay 0.000114 \
      --grad_clip 2.0 \
      --output "$base_out" \
      $ablation_flag \
      2>&1 | tee "${OUT_ROOT}/logs/MARCP_${task}_${ablation_short}_standalone.log"
  )
  local train_ret=$?

  local actual_out="$base_out"
  [ "$ablation_short" != "baseline" ] && actual_out="$base_out/$ablation_short"
  if [ -f "${actual_out}/model_params.pth" ] && [ -f "${actual_out}/test_metrics_summary.txt" ]; then
    cp "${actual_out}/model_params.pth" "${weight_dir}/"
    local f
    for f in test_metrics_summary.txt test_metric.txt run_info.txt val_metric.txt train_metric.txt ablation_config.txt; do
      [ -f "${actual_out}/$f" ] && cp "${actual_out}/$f" "${metric_dir}/"
    done
    rm -rf "${base_out}"
    echo "[OK] $task $ablation_short"
    ( flock 9; echo "$task $ablation_short $(date +%s)" >> "$COMPLETED_LOG"; ) 9>"$SUMMARY_LOCK" 2>/dev/null || echo "$task $ablation_short $(date +%s)" >> "$COMPLETED_LOG"
  else
    rm -rf "${base_out}" 2>/dev/null || true
    echo "[FAIL] $task $ablation_short exit=$train_ret 详见 logs/MARCP_${task}_${ablation_short}_standalone.log"
  fi
  return "$train_ret"
}

pids=()
run_idx=0
for line in "${TASK_LINES[@]}"; do
  run_idx=$((run_idx + 1))
  if [ "$GPU_ROTATE" = "1" ]; then
    gpu_id=$(( (run_idx - 1) % 4 ))
  else
    gpu_id=0
  fi
  run_single_task "$line" "$run_idx" "$n" "$gpu_id" &
  pids+=($!)
done

echo "已全部启动 ${#pids[@]} 个后台任务，等待全部结束…"
for pid in "${pids[@]}"; do
  wait "$pid" 2>/dev/null
done

echo ""
echo "========== 本批并行补跑结束 =========="
echo "独立完成记录: $COMPLETED_LOG"
exit 0
