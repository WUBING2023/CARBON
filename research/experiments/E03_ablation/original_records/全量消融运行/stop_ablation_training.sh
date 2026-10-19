#!/bin/bash
# =============================================================================
# MARCP 全量消融：退出训练（结束主调度及 marcp_train 子进程，清理 workspace/MARCP）
# 用法：bash stop_ablation_training.sh
# 仅结束 cmdline 含「全量消融运行」的进程，不杀 论文data/全量运行 等其他训练
# =============================================================================
BASE="/data/CARBON_PROJECT"
OUT_ROOT="${BASE}/消融实验/全量消融运行"
PID_FILE="${OUT_ROOT}/.run_ablation.pid"
STATE_FILE="${OUT_ROOT}/.run_ablation.state"
WORKSPACE="${OUT_ROOT}/workspace/MARCP"

OUR_RUN_MARKER="全量消融运行"
kill_marcp_residual() {
    local count=0
    for c in $(pgrep -f "marcp_train" 2>/dev/null); do
        [ "$c" = "$$" ] && continue
        [ -r "/proc/$c/cmdline" ] || continue
        cmd=$(tr '\0' ' ' < "/proc/$c/cmdline" 2>/dev/null)
        echo "$cmd" | grep -q "$OUR_RUN_MARKER" || continue
        kill -TERM "$c" 2>/dev/null || true
        count=$((count+1))
    done
    sleep 2
    for c in $(pgrep -f "marcp_train" 2>/dev/null); do
        [ "$c" = "$$" ] && continue
        [ -r "/proc/$c/cmdline" ] || continue
        cmd=$(tr '\0' ' ' < "/proc/$c/cmdline" 2>/dev/null)
        echo "$cmd" | grep -q "$OUR_RUN_MARKER" || continue
        kill -KILL "$c" 2>/dev/null || true
    done
    if command -v nvidia-smi >/dev/null 2>&1; then
        for c in $(nvidia-smi --query-compute-apps=pid --format=csv,noheader,nounits 2>/dev/null); do
            c=$(echo "$c" | tr -d ' ')
            [ -z "$c" ] || [ "$c" = "$$" ] && continue
            [ ! -r "/proc/$c/cmdline" ] && continue
            cmd=$(tr '\0' ' ' < "/proc/$c/cmdline" 2>/dev/null)
            if echo "$cmd" | grep -q "marcp_train" && echo "$cmd" | grep -q "$OUR_RUN_MARKER"; then
                kill -KILL "$c" 2>/dev/null || true
                count=$((count+1))
            fi
        done
    fi
    [ "$count" -gt 0 ] && echo "已结束 $count 个本运行残留训练进程。"
}

echo "========== 退出全量消融训练 =========="
if [ ! -f "$PID_FILE" ]; then
    echo "未发现运行中的消融训练 (无 PID 文件)。"
    kill_marcp_residual
    if [ -d "$WORKSPACE" ]; then
        n=$(find "$WORKSPACE" -maxdepth 1 -type d ! -path "$WORKSPACE" 2>/dev/null | wc -l)
        if [ "$n" -gt 0 ]; then
            echo "清理 workspace 中未完整训练数据 ($n 个)..."
            rm -rf "${WORKSPACE:?}"/*
            echo "已清理。"
        fi
    fi
    rm -f "$PID_FILE" "$STATE_FILE" "${OUT_ROOT}/.gpu_slot_counts" "${OUT_ROOT}/.gpu_done" "${OUT_ROOT}/.gpu_slot_counts.lock" "${OUT_ROOT}/.gpu_done.lock"
    echo "===================================="
    exit 0
fi

pid=$(cat "$PID_FILE")
kill_tree() {
    local p=$1
    local sig=${2:-TERM}
    local c
    for c in $(pgrep -P "$p" 2>/dev/null); do
        kill_tree "$c" "$sig"
        kill -"$sig" "$c" 2>/dev/null || true
    done
}
if kill -0 "$pid" 2>/dev/null; then
    echo "正在结束消融训练进程及其子进程 (主 PID=$pid)..."
    kill_tree "$pid" "TERM"
    kill -TERM "$pid" 2>/dev/null || true
    for _ in 1 2 3 4 5 6 7 8 9 10; do
        kill -0 "$pid" 2>/dev/null || break
        sleep 1
    done
    if kill -0 "$pid" 2>/dev/null; then
        echo "主进程未退出，发送 SIGKILL..."
        kill_tree "$pid" "KILL"
        kill -KILL "$pid" 2>/dev/null || true
    fi
    sleep 2
    kill_marcp_residual
    echo "进程已结束。"
else
    echo "主进程 (PID=$pid) 已不存在；正在检查并清理可能残留的训练子进程..."
    kill_marcp_residual
fi

echo "再次检查并清理占用 GPU 的 marcp_train 进程..."
kill_marcp_residual
echo "清除 workspace 中未完整训练的数据..."
[ -d "$WORKSPACE" ] && rm -rf "${WORKSPACE:?}"/*
echo "清除运行状态文件..."
rm -f "$PID_FILE" "$STATE_FILE" "${OUT_ROOT}/.gpu_slot_counts" "${OUT_ROOT}/.gpu_done" "${OUT_ROOT}/.gpu_slot_counts.lock" "${OUT_ROOT}/.gpu_done.lock"
echo "退出完成。已保存的权重与指标（正式目录）未改动。"
echo "===================================="
