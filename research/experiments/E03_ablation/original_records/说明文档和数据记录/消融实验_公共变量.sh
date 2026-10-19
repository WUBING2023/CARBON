#!/bin/bash
# 消融实验公共变量：所有 run_ablation_*.sh 均 source 本文件，保证数据、输出、超参一致，消融公平性。
# 使用前请先填写下方 TRAIN_FILE、VAL_FILE、TASK（及可选 TEST_FILE）。

BASE_OUTPUT="./marcp_trained_model/ablation"
TRAIN_FILE=""    # 必填：训练集 .h5 路径
VAL_FILE=""      # 必填：验证集 .h5 路径
TEST_FILE=""     # 可选：测试集 .h5 路径
TASK="Protein_Cell_Accession"   # 如 AGGF1_HepG2_ENCSR543TPH

# 与运行示例_消融.sh 完全一致的公共参数（仅差在是否加 --ablation_*）
COMMON="python marcp_train.py \
  -train $TRAIN_FILE -validation $VAL_FILE -test $TEST_FILE \
  --output $BASE_OUTPUT --task $TASK \
  --modelFramework Treat --lossM CombinedProfileLoss \
  --input_channels 5 --bsz 256 --epochs 100 --lr 2e-4 \
  --weight_decay 0.01 --grad_clip 1.0 \
  --max_stopping_step 8 --resource_fraction 0.9"
