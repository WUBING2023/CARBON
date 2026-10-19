# Figure 1 native single-column repair, v32

交付位于本目录，未改图2/4、共享主稿、CURRENT_REVISION或线上。依据当前Methods和v31权威图，重排为原生85 × 170 mm的可编辑矢量图，保留a–f六个角色。原图/PPTX完整留在baseline/。

## 可读性与物理检查

- 旧源画布400.008 × 346.414 mm，放置为85 mm时高73.611 mm；PDF文字实际缩为2.269–3.830 pt。
- 新源/最终MediaBox均为85 × 170 mm，所有普通矢量文字最小5.5 pt；面板字母9 pt。没有机械缩小宽画布。
- 96/300 dpi完整新图及85 mm前后比较已经实际查看。最终普通文字均在MediaBox内；底部文字余量0.94 mm。已修复首轮乘/加运算符及BiGRU内部箭头的z-order遮挡，并复看真实像素；四个关键运算圆内的×/×/+/+均可见。
- motif字母是科学字形，其相对高度沿用原示例；没有把这些不同高度的图形当作普通字体最小值。外围微小数字/刻度已移除；新增的CWM/PPM标签为5.5 pt。

## 科学关系逐项对应

- a：101 nt、序列4通道+对齐PhyloP 1通道→MSME→RCE；RCE分支进入profile和count，softmax与exp(z)经×合成profile。输入仍只有5通道，无额外结构输入。
- b：宽度3/5/7卷积各32通道、concat96；SE globalmean96→24ReLU→96sigmoid，权重向下作用于Unpooled特征；gate输出先1×1fusion96→96，再加独立5→96输入shortcut。没有在残差相加后添加激活。
- c：96→64投影→双向GRU每方向64→concat128；与独立96→128 shortcut相加→LayerNorm→dropout0.2。分支为profile 1×1 128→64→32→1 logits→softmax，以及count globalmean128→32→1 z→exp(z)。
- d：三条loss分支由α、β、γ权重汇合为sum/s。已核对当前正文：s=1，α=1.4194219、β=0.5402181、γ=1.7983396。参数留在Methods/所引总loss公式中，图内不重复8位小数。没有额外独立机制损失或因果主张。
- e：IG窗口→sequence segments→TF-MoDISco→CWM/PPM示例；原PPTX的image4.svg全部94个彩色motif字形的d路径/色值逐项相同，连同两个原baseline组保留。只移除外围黑色轴字和刻度，未生成新motif或结果。
- f：原流程中的Observed/Predicted示意→A/U单碱基替换→重新预测的WT/Mutated对照；序列片段、profile曲线继续作为流程示意，不是新测量、推理或实验。灰区/虚线分别保留motif区和替换位置的原语义。

颜色沿用本篇PALETTE_V3：A绿、C深蓝#19499C、G橙、U/T红；Full红#E64B35、突变cyan、PhyloP灰蓝及原示意实测mint的含义保留。白底、细灰线、Arial与图2–4统一。已查看local style reference；仅借鉴紧凑层级和连续复合结构，未引入外部数据或模块。

## 文件与可迁移重建

- fig1_v32.pdf / png / svg：供主Agent整合的最终图。
- render_fig1_v32.py：原生mm坐标、可编辑文字/模块/箭头；以__file__定位frozen_inputs，不依赖当前工作目录，不会写原稿目录。
- frozen_inputs/PALETTE_V3.json、image4.svg、motif_letters_only.svg：最终渲染所需冻结资源；Methods/loss snapshots用于科学映射追溯。其余原PPTX嵌入资产保留供追溯。
- FIG1_NATIVE_QA.json：全部节点、箭头坐标/意义、文字物理边界、字体、motif几何和loss对应检查。
- FIG1_BEFORE_AFTER_85mm.pdf/png：完整旧/新图都按85 mm宽度实际放置。
- CAPTION_V32.tex：当前正文图注原文，科学内容无需改动；现有末句已注明序列片段/曲线是流程示意，并将BN、hidden activations和head dropout细节指回Methods/S1.4。

已将脚本和3个所需资源复制到portable_probe/并从无关cwd M:/RBP重建；96 dpi像素与最终图完全相同，实际验证不依赖原PPTX或主稿路径。

重建：`python render_fig1_v32.py`。依赖numpy、matplotlib、PyMuPDF；无需PowerPoint、模型文件、训练/推理或联网。verify_fig1_v32.py另用于本机源稿/图注对应验证和85 mm比较；主Agent完成整合编译后再做主文实际放置及独立审查。
