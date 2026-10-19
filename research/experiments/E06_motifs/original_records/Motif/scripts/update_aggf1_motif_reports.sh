#!/usr/bin/env bash
# 更新 AGGF1 HepG2 chr22：CISBP 单库对比 + 多源 benchmark（MEME 库 / DEWSeq 目录 / 可选 STREME）
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
# 默认可解释性输出（全染色体/Top-K 等新流水线）：outputs/AGGF1_HepG2_ENCSR543TPH
# 旧版按染色体命名的目录已废弃；统一使用实验任务目录 outputs/<TASK>/
TASK="${TASK:-outputs/AGGF1_HepG2_ENCSR543TPH}"

echo "[1/4] 可选：从 CISBP pwms 同步 MEME 示例子集 …"
python3 scripts/sync_motif_library_examples.py || true

echo "[2/4] CISBP 全库扫描对比 …"
python3 scripts/compare_outputs_to_cisbp.py --task-dir "$TASK"

echo "[3/4] 多源 benchmark（CISBP + meme_rna + dewseq + STREME）…"
python3 scripts/build_motif_benchmark.py --task-dir "$TASK"

echo "[4/4] 完成。请打开："
echo "  $TASK/cisbp_compare/comparison_report.html"
echo "  $TASK/motif_benchmark/benchmark_report.html"
echo "  $TASK/motif_benchmark/benchmark_full.json"
