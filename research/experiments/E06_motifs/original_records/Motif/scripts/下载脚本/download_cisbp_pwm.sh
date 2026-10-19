#!/usr/bin/env bash
# 从 CISBP-RNA 官网「Bulk downloads」获取 Build 2.00：
#   - pwms/*.txt（PWM 矩阵）
#   - RBP_Information.txt、RBP_Information_all_motifs.txt（Motif_ID ↔ RBP 等元数据）
#
# 重要：必须在**两次独立请求**中分别下载 PWM 与 RBP_Information。若在同一 POST 中同时勾选
# PWMs 与 RBP_Information，服务器会生成错误的单文件「PWMs」（非 pwms/ 目录），会破坏比对脚本。
#
# 用法：在 Motif 目录下执行
#   bash scripts/download_cisbp_pwm.sh
# 或指定输出目录：
#   bash scripts/download_cisbp_pwm.sh /path/to/CISBP-RNA

set -euo pipefail

BASE="https://cisbp-rna.ccbr.utoronto.ca"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:-${SCRIPT_DIR}/../reference_data/CISBP-RNA}"
OUTDIR="$(cd "$OUTDIR" && pwd)"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

fetch_zip() {
  local label="$1"
  local zip_name="$2"
  shift 2
  echo "[request] ${label} ..."
  local RESP
  RESP="$(curl -sS -X POST "${BASE}/bulk_archive.php" "$@" -F "submit=Download Entire Archive!")"
  local URL
  URL="$(echo "$RESP" | sed -n 's/.*href="\([^"]*entiredata[^"]*\.zip\)".*/\1/p' | head -1)"
  if [[ -z "${URL}" ]]; then
    echo "Could not parse tmp/*.zip link for: ${label}. Open ${BASE}/bulk.php in a browser." >&2
    exit 1
  fi
  echo "[download] ${BASE}/${URL} -> ${zip_name}"
  curl -fsSL "${BASE}/${URL}" -o "${zip_name}"
}

echo "=== [1/2] PWMs only (pwms/*.txt) ==="
fetch_zip "PWMs only" "cisbp-rna_PWMs_entire_2_00.zip" -F "EntireDataset[]=PWMs"
echo "[unzip] cisbp-rna_PWMs_entire_2_00.zip"
unzip -o -q "cisbp-rna_PWMs_entire_2_00.zip"

echo "=== [2/2] RBP_Information（仅解压两个元数据表，避免重复解压 pwms）==="
fetch_zip "RBP_Information only" "cisbp-rna_RBPinfo_2_00.zip" -F "EntireDataset[]=RBP_Information"
echo "[unzip] RBP_Information.txt + RBP_Information_all_motifs.txt from cisbp-rna_RBPinfo_2_00.zip"
unzip -o -q "cisbp-rna_RBPinfo_2_00.zip" RBP_Information.txt RBP_Information_all_motifs.txt

echo "Done."
echo "  PWM:  ${OUTDIR}/pwms/ (tab-separated, columns Pos A C G U)"
echo "  RBP:  ${OUTDIR}/RBP_Information.txt"
echo "        ${OUTDIR}/RBP_Information_all_motifs.txt"
