#!/bin/bash
# =============================================================================
# SERBP1_K562_ENCSR121GQH 的 6 种消融（与 batch11_tasks.txt 第 6～11 行一致）删除本地权重与指标，
# 以便下一次 run_ablation_background.sh 固定 11 条时真正重训（固定列表模式不跳过「已完成」，
# 但若仍残留旧权重可能影响 marcp 续训行为，建议先清再跑）。
# 用法: bash prepare_serbp1_batch11_rerun.sh
# 然后: bash run_ablation_background.sh
# =============================================================================
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_ROOT="${SCRIPT_DIR}"
TASK="SERBP1_K562_ENCSR121GQH"
for ab in no_se no_multiscale no_gru no_combined_loss no_conservative baseline; do
    rm -rf "${OUT_ROOT}/模型权重/MARCP/${TASK}/${ab}"
    rm -rf "${OUT_ROOT}/各个具体实验的指标/MARCP/${TASK}/${ab}"
done
echo "已清除 ${TASK} 的 6 种消融目录（模型权重/指标，若存在）。"
echo "下一步: bash run_ablation_background.sh"
echo "  （默认 USE_FIXED_BATCH11=1，使用 batch11_tasks.txt 共 11 条同批调度；query 请用进程数判断并发）"
