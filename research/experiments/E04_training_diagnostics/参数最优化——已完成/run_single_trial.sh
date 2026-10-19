#!/bin/bash
#
# MARCP单次超参数试验训练脚本
# 用于手动测试特定超参数组合
#
# 用法:
#   ./run_single_trial.sh <task_name> [options]
#
# 示例:
#   ./run_single_trial.sh AGGF1_HepG2_ENCSR543TPH --lr 0.0005 --bsz 256
#   ./run_single_trial.sh HNRNPC_K562_ENCSR249ROI --lr 0.0002 --weight_decay 0.01
#

set -e

# ============== 默认配置 ==============
MARCP_DIR="/data/CARBON_PROJECT/MARCP"
DATA_ROOT="/data/CARBON_PROJECT/数据/h5data"
OUTPUT_ROOT="/data/CARBON_PROJECT/参数最优化/results/manual_trials"

# 默认超参数
DEFAULT_LR="0.0002"
DEFAULT_BSZ="256"
DEFAULT_EPOCHS="100"
DEFAULT_WEIGHT_DECAY="0.01"
DEFAULT_GRAD_CLIP="1.0"
DEFAULT_MAX_STOPPING_STEP="8"
DEFAULT_ALPHA="1.0"
DEFAULT_BETA="0.5"
DEFAULT_GAMMA="1.0"
DEFAULT_CUDA="0"

# ============== 解析参数 ==============
if [ $# -lt 1 ]; then
    echo "用法: $0 <task_name> [options]"
    echo ""
    echo "必需参数:"
    echo "  task_name    数据集名称，如 AGGF1_HepG2_ENCSR543TPH"
    echo ""
    echo "可选参数:"
    echo "  --lr <float>              学习率 (默认: $DEFAULT_LR)"
    echo "  --bsz <int>               Batch size (默认: $DEFAULT_BSZ)"
    echo "  --epochs <int>            训练轮数 (默认: $DEFAULT_EPOCHS)"
    echo "  --weight_decay <float>    权重衰减 (默认: $DEFAULT_WEIGHT_DECAY)"
    echo "  --grad_clip <float>       梯度裁剪 (默认: $DEFAULT_GRAD_CLIP)"
    echo "  --max_stopping_step <int> 早停耐心 (默认: $DEFAULT_MAX_STOPPING_STEP)"
    echo "  --alpha <float>           Profile Shape Loss权重 (默认: $DEFAULT_ALPHA)"
    echo "  --beta <float>            Jitter Loss权重 (默认: $DEFAULT_BETA)"
    echo "  --gamma <float>           Count Loss权重 (默认: $DEFAULT_GAMMA)"
    echo "  --cuda <str>              CUDA设备ID (默认: $DEFAULT_CUDA)"
    echo "  --output <dir>            输出目录 (默认: 自动生成)"
    echo ""
    echo "示例:"
    echo "  $0 AGGF1_HepG2_ENCSR543TPH --lr 0.0005 --bsz 256"
    exit 1
fi

TASK_NAME=$1
shift

# 初始化参数为默认值
LR=$DEFAULT_LR
BSZ=$DEFAULT_BSZ
EPOCHS=$DEFAULT_EPOCHS
WEIGHT_DECAY=$DEFAULT_WEIGHT_DECAY
GRAD_CLIP=$DEFAULT_GRAD_CLIP
MAX_STOPPING_STEP=$DEFAULT_MAX_STOPPING_STEP
ALPHA=$DEFAULT_ALPHA
BETA=$DEFAULT_BETA
GAMMA=$DEFAULT_GAMMA
CUDA=$DEFAULT_CUDA
OUTPUT_DIR=""

# 解析可选参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --lr)
            LR="$2"
            shift 2
            ;;
        --bsz)
            BSZ="$2"
            shift 2
            ;;
        --epochs)
            EPOCHS="$2"
            shift 2
            ;;
        --weight_decay)
            WEIGHT_DECAY="$2"
            shift 2
            ;;
        --grad_clip)
            GRAD_CLIP="$2"
            shift 2
            ;;
        --max_stopping_step)
            MAX_STOPPING_STEP="$2"
            shift 2
            ;;
        --alpha)
            ALPHA="$2"
            shift 2
            ;;
        --beta)
            BETA="$2"
            shift 2
            ;;
        --gamma)
            GAMMA="$2"
            shift 2
            ;;
        --cuda)
            CUDA="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "未知参数: $1"
            exit 1
            ;;
    esac
done

# ============== 设置数据路径 ==============
TASK_DIR="${DATA_ROOT}/${TASK_NAME}/_5ch"
TRAIN_FILE="${TASK_DIR}/${TASK_NAME}_chr1_5ch.h5"
VAL_FILE="${TASK_DIR}/${TASK_NAME}_chr2_5ch.h5"
TEST_FILE="${TASK_DIR}/${TASK_NAME}_chr3_5ch.h5"

# 检查数据文件
if [ ! -f "$TRAIN_FILE" ]; then
    echo "错误: 训练文件不存在: $TRAIN_FILE"
    exit 1
fi
if [ ! -f "$VAL_FILE" ]; then
    echo "错误: 验证文件不存在: $VAL_FILE"
    exit 1
fi
if [ ! -f "$TEST_FILE" ]; then
    echo "警告: 测试文件不存在: $TEST_FILE"
    TEST_FILE=""
fi

# ============== 设置输出目录 ==============
if [ -z "$OUTPUT_DIR" ]; then
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    OUTPUT_DIR="${OUTPUT_ROOT}/${TASK_NAME}_${TIMESTAMP}"
fi

mkdir -p "$OUTPUT_DIR"

# ============== 打印配置 ==============
echo "=============================================="
echo "MARCP超参数试验"
echo "=============================================="
echo "任务名称: $TASK_NAME"
echo "输出目录: $OUTPUT_DIR"
echo ""
echo "超参数配置:"
echo "  lr:               $LR"
echo "  bsz:              $BSZ"
echo "  epochs:           $EPOCHS"
echo "  weight_decay:     $WEIGHT_DECAY"
echo "  grad_clip:        $GRAD_CLIP"
echo "  max_stopping_step: $MAX_STOPPING_STEP"
echo "  alpha:            $ALPHA"
echo "  beta:             $BETA"
echo "  gamma:            $GAMMA"
echo "  cuda:             $CUDA"
echo "=============================================="

# ============== 构建并运行训练命令 ==============
cd "$MARCP_DIR"

CMD="python -u marcp_train.py \
  --CUDA $CUDA \
  --task $TASK_NAME \
  --output $OUTPUT_DIR \
  -train $TRAIN_FILE \
  -validation $VAL_FILE \
  --epochs $EPOCHS \
  --bsz $BSZ \
  --lr $LR \
  --weight_decay $WEIGHT_DECAY \
  --grad_clip $GRAD_CLIP \
  --max_stopping_step $MAX_STOPPING_STEP \
  --alpha $ALPHA \
  --beta $BETA \
  --gamma $GAMMA \
  --input_channels 5 \
  --encode OneHot \
  --modelFramework Treat \
  --lossM CombinedProfileLoss \
  --signalMax 0 \
  --seqLength 101"

if [ -n "$TEST_FILE" ]; then
    CMD="$CMD -test $TEST_FILE"
fi

# 保存命令到文件
echo "$CMD" > "${OUTPUT_DIR}/train_command.sh"
chmod +x "${OUTPUT_DIR}/train_command.sh"

# 运行训练
echo ""
echo "开始训练..."
echo ""

$CMD 2>&1 | tee "${OUTPUT_DIR}/train.log"

# ============== 输出结果 ==============
echo ""
echo "=============================================="
echo "训练完成"
echo "=============================================="
echo "输出目录: $OUTPUT_DIR"
echo ""

if [ -f "${OUTPUT_DIR}/val_metric.txt" ]; then
    echo "验证集最终指标:"
    tail -1 "${OUTPUT_DIR}/val_metric.txt"
fi

if [ -f "${OUTPUT_DIR}/test_metric.txt" ]; then
    echo ""
    echo "测试集指标:"
    cat "${OUTPUT_DIR}/test_metric.txt"
fi

echo ""
echo "完成!"
