#!/bin/bash
# =============================================================================
# MARCP 全量消融：在后台持续运行 auto_ablation_runner.sh（监督/守护），关闭终端后仍会运行
# 守护脚本会定期查询进度，未完成时自动启动 run_ablation_background.sh，异常时先 stop 再重起
# 用法：bash run_auto_ablation_runner_background.sh
# 查看日志：tail -f auto_ablation_runner.log
# 结束守护：kill $(cat .run_auto_ablation_runner.pid)
# =============================================================================
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
PID_FILE="${SCRIPT_DIR}/.run_auto_ablation_runner.pid"
LOG_FILE="${SCRIPT_DIR}/auto_ablation_runner.log"

if [ -f "$PID_FILE" ]; then
    pid=$(cat "$PID_FILE")
    if kill -0 "$pid" 2>/dev/null; then
        echo "MARCP 全量消融守护已在后台运行 (PID=$pid)。日志: $LOG_FILE  结束: kill $pid"
        exit 0
    fi
    rm -f "$PID_FILE"
fi

echo "在后台启动 MARCP 全量消融守护，日志: $LOG_FILE"
nohup bash ./auto_ablation_runner.sh >> "$LOG_FILE" 2>&1 &
echo $! > "$PID_FILE"
echo "已启动 PID=$(cat "$PID_FILE")。查看日志: tail -f auto_ablation_runner.log  结束守护: kill \$(cat .run_auto_ablation_runner.pid)"
