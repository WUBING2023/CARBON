# 外部功能定位验证：实际研究图交付

这组图已经使用真实的完整分析 CSV 生成。结果不支持“CARBON 在六个 RNA 上总体解释性提升”：101nt 主结果的宏平均位置相关性为 A 0.017736、S 0.066557、同权重 A-flat -0.033184。全部六个位点和负结果保留，图没有用于替换现有论文。

`functional_positions_primary_101.pdf/.svg/.png` 是主图：六个 RNA 的 90 个位置分别展示实验/A/S 三行；每行除以本行最大值只比较空间形状，不能拿颜色深浅比较绝对效应强度。下方同一坐标轴展示 A、S、A-flat 的全部六个位置保真度。三次固定置换的完整结果留在分析表，避免主图堆积。

`U2AF2_signed_variants_primary_101.pdf/.svg/.png` 展示全部 270 个真实单碱基突变。灰色格是 WT 碱基身份，另外三种改变均有实际数据。实验使用一个实测 Δlog2 enrichment 色标，A/S 共同使用另一个原生 Δlog1p-count 色标；两种单位不能互换。主 G4 [31,51) 与 polyC [15,18) 是原论文来源注释，坐标为 0-based half-open。RNA 来源 U2AF2 与结合蛋白 G3BP1 已在图注区分。

`*_sensitivity_90` 是同一全部病例的天然 90nt 长度敏感性副本，没有用它替换主结果。模型输入与真实生化 RNA 的恒定接头/T7 背景不同，101nt 又额外包含天然基因组上下文；这些范围差异与单 checkpoint 限制都在 FIGURE_CAPTIONS.md 明示。

四幅输出宽度均为 89.5mm；最小 PDF 实际文本为 6.3pt；使用 NPG 身份色、统一小写面板字母，没有主标题或阴影。实际 PDF 全部用 PyMuPDF 300dpi 渲染并人工目视复查。48 项全基因/全长度/四条件位置相关性独立重建通过；U2AF2 两长度 A/S 的四项 signed variant rho 也独立重建通过。图源 CSV 保留绝对值、归一化分母与每个实际变体，不只保留图片。

可运行脚本为 `plot_functional_validation.py`，默认读取 `03_protocol_review/functional_map_analysis` 的实际 CSV、STATISTICS.json 与 ANALYSIS_COMPLETION.json。缺少完成标记或真实数据时脚本会停止，没有模拟数据分支。FIGURE_QA.json 记录所有输入和导出的 SHA256、物理尺寸、字体、重建检查与源文件路径。

图注完整英文说明与实际统计结果见 `FIGURE_CAPTIONS.md`。来源为 Martyr et al., Nucleic Acids Research (2026), https://doi.org/10.1093/nar/gkag821。六个基因名表示 RNA 背景，实验的结合蛋白统一为重组人全长 G3BP1；突变是六个位点内部的嵌套观察，不能宣称一千多个独立生物学位点。
