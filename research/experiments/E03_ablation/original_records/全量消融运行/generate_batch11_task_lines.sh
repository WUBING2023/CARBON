#!/bin/bash
# 生成本批 11 条任务行（与全量消融同一 file_list 规则），输出到 stdout；不依赖 .run_ablation.state
set -e
BASE="/data/CARBON_PROJECT"
DATA_ROOT="${BASE}/数据/h5data"
TRAIN_CHRS="3,4,5,6,7,10,11,12,13,14,17,18,19,20,21,22"
VAL_CHRS="2,9,16"
TEST_CHRS="1,8,15"

file_list() {
  local task=$1
  local chrs=$2
  local sep="" out="" f c
  for c in ${chrs//,/ }; do
    f="${DATA_ROOT}/${task}/_5ch/${task}_chr${c}_5ch.h5"
    [ -f "$f" ] || continue
    out="${out}${sep}${f}"
    sep=","
  done
  echo "$out"
}

ablation_flag() {
  case "$1" in
    no_se) echo "--ablation_no_se" ;;
    no_multiscale) echo "--ablation_no_multiscale" ;;
    no_gru) echo "--ablation_no_gru" ;;
    no_combined_loss) echo "--ablation_no_combined_loss" ;;
    no_conservative) echo "--ablation_no_conservative" ;;
    baseline) echo "" ;;
    *) echo "" ;;
  esac
}

# 与当时 STATE 中本批 11 条一致：任务名|消融简称
SPECS=(
  "AGGF1_HepG2_ENCSR543TPH|no_conservative"
  "DDX52_HepG2_ENCSR970FEW|no_gru"
  "KHSRP_K562_ENCSR438GZQ|no_conservative"
  "HNRNPC_K562_ENCSR249ROI|no_conservative"
  "GRWD1_HepG2_ENCSR893NWB|no_conservative"
  "SERBP1_K562_ENCSR121GQH|no_se"
  "SERBP1_K562_ENCSR121GQH|no_multiscale"
  "SERBP1_K562_ENCSR121GQH|no_gru"
  "SERBP1_K562_ENCSR121GQH|no_combined_loss"
  "SERBP1_K562_ENCSR121GQH|no_conservative"
  "SERBP1_K562_ENCSR121GQH|baseline"
)

idx=0
for spec in "${SPECS[@]}"; do
  idx=$((idx + 1))
  IFS='|' read -r task ab_short <<< "$spec"
  train_5ch=$(file_list "$task" "$TRAIN_CHRS")
  val_5ch=$(file_list "$task" "$VAL_CHRS")
  test_5ch=$(file_list "$task" "$TEST_CHRS")
  if [ -z "$train_5ch" ] || [ -z "$val_5ch" ] || [ -z "$test_5ch" ]; then
    echo "错误: 任务 $task 缺少 5ch 文件，请检查 ${DATA_ROOT}/${task}/_5ch" >&2
    exit 1
  fi
  flag=$(ablation_flag "$ab_short")
  echo "${task}|${train_5ch}|${val_5ch}|${test_5ch}|${ab_short}|${flag}|${idx}"
done
