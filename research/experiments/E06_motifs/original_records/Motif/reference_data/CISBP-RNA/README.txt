官网 Bulk downloads：https://cisbp-rna.ccbr.utoronto.ca/bulk.php

推荐：在 Motif 根目录执行
  bash scripts/download_cisbp_pwm.sh

会执行两次独立请求并解压到本目录（Build 2.00）：
1) cisbp-rna_PWMs_entire_2_00.zip → pwms/*.txt（PWM 矩阵）
2) cisbp-rna_RBPinfo_2_00.zip → 仅解压 RBP_Information.txt、RBP_Information_all_motifs.txt

注意：不要在同一 POST 中同时勾选 PWMs 与 RBP_Information，服务器会生成错误的单文件「PWMs」；
脚本已拆成两次请求。用 Motif_ID（如 M29490_2.00）在 RBP_Information.txt 中可查 RBP_Name 等列。

仅人类（Homo_sapiens）子集（推荐，与全库 pwms 一致）：
  python3 scripts/filter_cisbp_homo_sapiens.py
输出：reference_data/CISBP-RNA/homo_sapiens/
  - pwms/（默认硬链接到全库，664 个 Motif，以当前表为准）
  - RBP_Information.txt、RBP_Information_all_motifs.txt（仅人类行）

或从官网「By Species → Homo_sapiens」下载（zip 内另有损坏的 PWMs_all_motifs 聚合文件，脚本会跳过）：
  bash scripts/download_cisbp_homo_sapiens.sh
输出：homo_sapiens_from_bulk/

Motif_ID 列中「.」表示该 RBP 无 motif，无对应 pwms 文件；非全空表。

本机已下载记录（请随重新下载更新）：
- 下载日期：2026-03-27
- Database Build：2.00
- PWM 文件数：约 29502（以 pwms/ 内实际数量为准）

详细说明：docs/CISBP-RNA_参考数据与比对.md
