#!/bin/bash
# 检查 batch11_tasks.txt 中 11 条是否均在「各个具体实验的指标」下存在 test_metrics_summary.txt
# 使用 awk 取第 1、5 列，避免超长行用 read/here-string 时触达环境限制
set +e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MET="${SCRIPT_DIR}/各个具体实验的指标/MARCP"
LIST="${SCRIPT_DIR}/batch11_tasks.txt"
ok=0 miss=0
while IFS=$'\t' read -r task ablation_short; do
  [[ -z "$task" ]] && continue
  f="${MET}/${task}/${ablation_short}/test_metrics_summary.txt"
  if [[ -f "$f" ]]; then
    echo "OK   ${task} / ${ablation_short}"
    ok=$((ok+1))
  else
    echo "MISS ${task} / ${ablation_short}"
    miss=$((miss+1))
  fi
done < <(awk -F'|' 'NF>=5 {print $1 "\t" $5}' "$LIST")
echo "--- 统计: OK=$ok MISSING=$miss (期望 11 条全 OK)"
[[ "$miss" -eq 0 ]] && exit 0 || exit 1
