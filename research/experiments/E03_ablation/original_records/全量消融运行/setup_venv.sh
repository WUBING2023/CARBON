#!/bin/bash
# =============================================================================
# 在全量消融运行目录下创建 venv 并安装依赖，供 run_ablation_*.sh 使用
# 用法：bash setup_venv.sh
# =============================================================================
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
VENV_DIR="${SCRIPT_DIR}/venv"
REQ_FILE="${SCRIPT_DIR}/requirements_venv.txt"

if [ ! -f "$REQ_FILE" ]; then
    echo "错误: 未找到 $REQ_FILE"
    exit 1
fi

echo "========== 创建并配置 venv =========="
echo "目录: $VENV_DIR"

if [ -d "$VENV_DIR" ]; then
    echo "venv 已存在，将仅执行 pip install -r（升级/补装依赖）。若要重建请先删除 venv 目录。"
else
    echo "创建 venv ..."
    if ! python3 -m venv "$VENV_DIR" 2>/dev/null; then
        echo "python3 -m venv 不可用，尝试使用 virtualenv ..."
        if command -v virtualenv >/dev/null 2>&1; then
            virtualenv -p python3 "$VENV_DIR" || { echo "virtualenv 失败"; exit 1; }
        else
            echo "请先安装 python3-venv 或 virtualenv："
            echo "  Debian/Ubuntu: sudo apt install python3-venv"
            echo "  或: pip install --user virtualenv 后重新执行本脚本"
            exit 1
        fi
    fi
fi

echo "激活 venv 并安装依赖 ..."
# shellcheck source=/dev/null
source "${VENV_DIR}/bin/activate"
pip install -q --upgrade pip
pip install -r "$REQ_FILE"

echo ""
echo "venv 已就绪。运行全量消融时将自动使用此环境。"
echo "手动激活: source ${VENV_DIR}/bin/activate"
