#!/usr/bin/env bash
# 从 CISBP-RNA 官网「Bulk downloads → By Species」下载 Homo_sapiens 归档（Build 2.00）。
# 返回的 zip 中常含：
#   - pwms_all_motifs/M*.txt（人类 PWM，目录名与全库 pwms/ 不同）
#   - RBP_Information.txt、RBP_Information_all_motifs.txt
#   - 以及损坏的聚合文件「PWMs_all_motifs」（含 unzip 杂讯，勿使用）
#
# 本脚本仅解压 pwms_all_motifs/* 与两个 RBP 表，删除 PWMs_all_motifs，并把目录改名为 pwms/ 以与全库布局一致。
#
# 若需与全库 pwms **完全一致**的人类 Motif 集合（通常多于物种 zip 内文件数），请改用：
#   python3 scripts/filter_cisbp_homo_sapiens.py
#
# 用法：在 Motif 目录下
#   bash scripts/download_cisbp_homo_sapiens.sh
#   bash scripts/download_cisbp_homo_sapiens.sh /path/to/out_parent   # 将创建 .../homo_sapiens_from_bulk/

set -euo pipefail

BASE="https://cisbp-rna.ccbr.utoronto.ca"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARENT="${1:-${SCRIPT_DIR}/../reference_data/CISBP-RNA}"
PARENT="$(cd "$PARENT" && pwd)"
OUT="${PARENT}/homo_sapiens_from_bulk"
ZIP_NAME="cisbp-rna_Homo_sapiens_bulk_2_00.zip"

mkdir -p "$OUT"
cd "$OUT"

echo "[1/4] Requesting Homo_sapiens archive from ${BASE}/bulk_archive.php ..."
RESP="$(curl -sS -X POST "${BASE}/bulk_archive.php" \
  -F "selSpec=Homo_sapiens" \
  -F "Spec[]=PFMs" \
  -F "Spec[]=RBP_Information" \
  -F "submit=Download Species Archive!")"

URL="$(echo "$RESP" | sed -n 's/.*href="\([^"]*tmp\/[^"]*\.zip\)".*/\1/p' | head -1)"
if [[ -z "${URL}" ]]; then
  echo "Could not parse tmp/*.zip link (expected tmp/Homo_sapiens_*.zip). Open ${BASE}/bulk.php manually." >&2
  exit 1
fi

echo "[2/4] Downloading ${BASE}/${URL} -> ${ZIP_NAME}"
curl -fsSL "${BASE}/${URL}" -o "${ZIP_NAME}"

echo "[3/4] Inspecting zip (物种包有时仅含 RBP 表 + 损坏的 PWMs_all_motifs，无 pwms_all_motifs/) ..."
if ! unzip -l "${ZIP_NAME}" | grep -q 'pwms_all_motifs/'; then
  echo "[error] 本 zip 中未包含 pwms_all_motifs/ 目录（服务器可能返回精简包）。无法从 Bulk 解压人类 PWM。" >&2
  echo "        请改用全库下载 + 本地筛选（与全库 pwms 100% 对齐）：" >&2
  echo "          bash scripts/download_cisbp_pwm.sh" >&2
  echo "          python3 scripts/filter_cisbp_homo_sapiens.py" >&2
  unzip -o -q "${ZIP_NAME}" RBP_Information.txt RBP_Information_all_motifs.txt 2>/dev/null || true
  exit 1
fi

echo "       Unzipping RBP 表 + pwms_all_motifs/*（跳过损坏的 PWMs_all_motifs 单文件） ..."
rm -rf pwms pwms_all_motifs RBP_Information.txt RBP_Information_all_motifs.txt PWMs_all_motifs 2>/dev/null || true
unzip -o -q "${ZIP_NAME}" RBP_Information.txt RBP_Information_all_motifs.txt "pwms_all_motifs/*"
rm -f PWMs_all_motifs

if [[ -d pwms_all_motifs ]]; then
  mv pwms_all_motifs pwms
  echo "       Renamed pwms_all_motifs/ -> pwms/"
fi

echo "[4/4] Done."
echo "  Output: ${OUT}/"
echo "  PWM files: $(find pwms -name 'M*_2.00.txt' 2>/dev/null | wc -l) under pwms/"
echo "  若需与全库 Motif 集合一致的人类子集，仍建议运行: python3 scripts/filter_cisbp_homo_sapiens.py"
