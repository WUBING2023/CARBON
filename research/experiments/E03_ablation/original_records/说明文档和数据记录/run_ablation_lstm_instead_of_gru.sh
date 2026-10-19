#!/bin/bash
# 消融实验 5：BiLSTM 替换 BiGRU。在 MARCP 目录下执行。
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/消融实验_公共变量.sh"
if [[ -z "$TRAIN_FILE" || -z "$VAL_FILE" ]]; then echo "请先编辑 消融实验_公共变量.sh 填写 TRAIN_FILE、VAL_FILE、TASK"; exit 1; fi
$COMMON --ablation_lstm_instead_of_gru
