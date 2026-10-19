#!/bin/bash
# =============================================================================
# MARCP 自动全量消融守护脚本
# 与 run_ablation_background.sh / query_ablation_progress.sh / stop_ablation_training.sh 配合，
# 使用 消融实验/MARCP/marcp_train.py 进行全量消融（baseline、no_se、no_multiscale、no_gru、no_combined_loss、no_conservative）。
#
# 每 5 分钟查询一次进度:
# - 若「状态: 未运行」且本批未完成>0，则自动以
#   SLOTS_GPU0=15 SLOTS_OTHER=20 BATCH_SIZE=256（GPU 槽位 15 20 20 20）启动一轮
# - 若「状态: 未运行」且本批未完成=0，则认为全部完成并退出
# - 若正在运行，则只打印状态并等待下一轮
# - 若正在运行但「其中 N 个占用 GPU」的 N < 7，且本批仍有大量任务待跑（未完成>16），则先 stop 再重启，
#   用于 OOM 后只剩少数进程仍占着调度的情况；小批量（如仅补跑十几个任务）启动初期进程数常<7，不得误杀。
#
# 用法:
#   前台运行: cd 全量消融运行 && bash auto_ablation_runner.sh
#   后台持续运行: bash run_auto_ablation_runner_background.sh（关终端后仍运行，日志见 auto_ablation_runner.log）
# =============================================================================

BASE="/data/CARBON_PROJECT"
OUT_ROOT="${BASE}/消融实验/全量消融运行"

cd "${OUT_ROOT}" || {
  echo "切换到 ${OUT_ROOT} 失败，退出。"
  exit 1
}

INTERVAL_SEC=300  # 每 5 分钟检查一次

echo "========== 自动全量消融守护启动 =========="
echo "工作目录: ${OUT_ROOT}"
echo "检查间隔: ${INTERVAL_SEC} 秒"
echo "启动参数: SLOTS_GPU0=15 SLOTS_OTHER=20 BATCH_SIZE=256（GPU 槽位 15 20 20 20）"

while true; do
  echo ""
  echo "---- $(date '+%F %T') 查询进度 ----"

  # 运行查询脚本并缓存输出，后续多次解析
  PROG_OUT=$(bash query_ablation_progress.sh 2>&1)
  echo "${PROG_OUT}"

  # 提取关键信息
  status_line=$(printf '%s\n' "${PROG_OUT}" | grep -E '^状态:' || true)
  batch_total=$(printf '%s\n' "${PROG_OUT}" | grep -E '^本批消融任务数:' | awk -F':' '{gsub(/ /,"",$2); print $2+0}' || echo 0)
  this_run_completed=$(printf '%s\n' "${PROG_OUT}" | grep -E '^本轮已完成:' | awk -F':' '{gsub(/ /,"",$2); print $2+0}' || echo 0)
  pending_this_run=$(printf '%s\n' "${PROG_OUT}" | grep -E '^本批未完成（尚在跑/队列）:' | awk -F':' '{gsub(/ /,"",$2); print $2+0}' || echo 0)
  # 解析「其中 N 个占用 GPU」的 N
  gpu_occupied=$(printf '%s\n' "${PROG_OUT}" | grep -oE '其中 [0-9]+ 个占用 GPU' | grep -oE '[0-9]+' | head -1)
  gpu_occupied=$((gpu_occupied + 0))

  # 判断是否在运行中
  if printf '%s\n' "${status_line}" | grep -q '运行中'; then
    # 仅当「待完成任务仍很多」时才认为可能 OOM 后僵死；小批量启动阶段 gpu 进程常少于 7，不得 stop
    if [ "${gpu_occupied}" -lt 7 ] 2>/dev/null && [ "${pending_this_run}" -gt 16 ] 2>/dev/null; then
      echo "状态: 主调度运行中，但占用 GPU 的进程数=${gpu_occupied} < 7 且本批未完成=${pending_this_run} > 16，将先退出再重新加载（OOM 恢复）。"
      echo "执行: bash stop_ablation_training.sh"
      bash stop_ablation_training.sh
      sleep 10
      echo "执行: SLOTS_GPU0=15 SLOTS_OTHER=20 BATCH_SIZE=256 bash run_ablation_background.sh"
      SLOTS_GPU0=15 SLOTS_OTHER=20 BATCH_SIZE=256 bash run_ablation_background.sh
      echo "已重新启动，等待 ${INTERVAL_SEC} 秒后再次检查。"
      sleep "${INTERVAL_SEC}"
      continue
    fi
    echo "状态: 主调度正在运行，本轮已完成=${this_run_completed}，未完成=${pending_this_run}，占用 GPU 进程数=${gpu_occupied}。等待下一次检查。"
    sleep "${INTERVAL_SEC}"
    continue
  fi

  # 未运行：根据是否还有未完成任务决定是补跑还是结束
  if [ "${pending_this_run}" -gt 0 ] && [ "${batch_total}" -gt 0 ]; then
    echo "状态: 未运行，且本批共有 ${batch_total} 个任务、其中 ${pending_this_run} 个未完成 → 自动补跑一轮。"
    echo "执行: SLOTS_GPU0=15 SLOTS_OTHER=20 BATCH_SIZE=256 bash run_ablation_background.sh"
    SLOTS_GPU0=15 SLOTS_OTHER=20 BATCH_SIZE=256 bash run_ablation_background.sh
    echo "已触发补跑，等待 ${INTERVAL_SEC} 秒后再次检查。"
    sleep "${INTERVAL_SEC}"
    continue
  fi

  # 未运行且本批未完成=0：认为全部完成
  if [ "${batch_total}" -eq 0 ] || [ "${pending_this_run}" -eq 0 ]; then
    echo "状态: 未运行，且本批任务已全部完成 (本批=${batch_total}，本轮已完成=${this_run_completed}，未完成=${pending_this_run})。"
    echo "自动守护脚本结束。"
    exit 0
  fi

  # 兜底：其他情况（极少见），再次等待
  echo "状态: 未运行，但无法明确是否全部完成 (本批=${batch_total}，未完成=${pending_this_run})，保守起见暂不自动启动。"
  echo "请手动检查日志后再运行本脚本。"
  sleep "${INTERVAL_SEC}"
done

