#!/bin/bash
# 消融实验运行示例（需根据实际数据路径修改 -train / -validation / -test 和 --output）
# 在 MARCP 目录下执行： bash ../说明文档和数据记录/运行示例_消融.sh
# 公平性：所有实验共用相同 COMMON 参数，仅通过单一 --ablation_* 区分，保证可比性。

BASE_OUTPUT="./marcp_trained_model/ablation"   # 所有消融结果会落在 BASE_OUTPUT/no_se、BASE_OUTPUT/no_conservative 等
TRAIN_FILE=""   # 请填写训练集 .h5 路径，如 train_5ch.h5
VAL_FILE=""    # 请填写验证集 .h5 路径
TEST_FILE=""   # 可选，测试集 .h5 路径
TASK="Protein_Cell_Accession"   # 如 AGGF1_HepG2_ENCSR543TPH

if [[ -z "$TRAIN_FILE" || -z "$VAL_FILE" ]]; then
  echo "请先编辑本脚本，填写 TRAIN_FILE、VAL_FILE 和 TASK"
  exit 1
fi

# 资源默认 90% CPU/GPU（--resource_fraction 0.9），运行时间与资源会写入各子目录 run_info.txt。不含 BiLSTM 消融。
COMMON="python marcp_train.py \
  -train $TRAIN_FILE -validation $VAL_FILE -test $TEST_FILE \
  --output $BASE_OUTPUT --task $TASK \
  --modelFramework Treat --lossM CombinedProfileLoss \
  --input_channels 5 --bsz 256 --epochs 100 --lr 2e-4 \
  --weight_decay 0.01 --grad_clip 1.0 \
  --max_stopping_step 8 --resource_fraction 0.9"

# 1) 完整模型（基线）
$COMMON

# 2) 消融：无 SE 通道注意力
$COMMON --ablation_no_se

# 3) 消融：无多尺度（单尺度 k=5）
$COMMON --ablation_no_multiscale

# 4) 消融：无 BiGRU（仅残差+LayerNorm）
$COMMON --ablation_no_gru

# 5) 消融：损失函数（用基线 poissonLoss 替代 CombinedProfileLoss）
$COMMON --ablation_no_combined_loss

# 6) 消融：无第 5 维保守性（输入 5ch 但只用前 4 维）
$COMMON --ablation_no_conservative
# 若要用 MSE 作为基线： $COMMON --ablation_no_combined_loss --ablation_baseline_loss mse

echo "所有消融实验启动完毕，请查看 $BASE_OUTPUT 下各子目录的 val_metric.txt 与 test_metrics_summary.txt"
