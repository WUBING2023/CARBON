# BEAPR 整图版式迁移试排：有限独立复核

复核日期：2026-10-15。复核者为独立子 Agent；未参与本轮生产绘图，未修改图、数据、正式稿或 Overleaf。仅写入本文件。本结论限于 `reference_match_20261015` 的候选与左右对照交付，不是整篇审稿、审美认证或投稿许可。

## 范围与结论

依次阅读 `REFERENCE_MAPPING.md`、实际完整 BEAPR Figure 5 原图、实际试排整图、绘图/数据映射、最终图注及两页实际对照预览，并检查保存后的 PDF 页面和对象几何。原图身份为 Yang et al., Nature Communications 10, 1338 (2019), DOI 10.1038/s41467-019-09292-w。本轮没有重新联网核验原文/数据库现状或进行新模型推理。

**版式迁移确实发生**：右侧跨全高的命名对象比较列、左上紧凑队列统计、左下重复等位基因案例列均可在完整原图与试排之间逐一辨认。它比旧 v25 更接近这张参考的整体组织，不只是借用配色。与此同时，CARBON 用的是四个面板、两个案例及两项蛋白终点，原图的生物实验及全部局部几何没有被一比一复制。

在当前明确的探索性/示例性范围内，**未发现需要阻止展示此候选的数据替换、符号、共同尺度、漏掉不利病例或对照 PDF 裁切问题**。这不等于已经满足用户的全部审美要求，也不等于自然变异或临床价值已获正面验证。后两者的证据限制相当实质。

## 原图对应与视觉判断

| 区域 | 实际迁移 | 科学边界 |
|---|---|---|
| 左上 a/b | 两个队列块保持紧凑，和下方案例共享左侧区域 | 原图 KD/TCGA 统计已被 CARBON 的 ASB 方向及疾病注释重叠统计替换；不能借原图实验赋予相同含义 |
| 左下 c | 两个同构的等位基因→七臂预测→作者 reporter 列，统一轴域 | 模型结合预测与作者 RNA 输出是不同终点/细胞，不能把上下排列解读为中介链验证 |
| 右侧 d | 25 个命名蛋白的全高配对比较，两个终点共用对象顺序 | 冻结的蛋白级点估计，不是参考图的组织分布；没有制造箱线图、样本分布或逐蛋白区间 |

实际 178×151 mm 候选的层级比 v25 清楚，病例信息也比泛化的装饰结构具体。所有八个基因位点标注均在；两个 MAP2K3 位点以等位基因区分。

仍有明确的视觉折衷：右列底部四个区间文本只有约 **4.5 pt**，小上标约 3.5 pt；这对最终排印偏小。上方队列与下方案例间仍有留白，右列点线与左列条形并非原图全部视觉语法。临床 n=1/2/3 的 0/50/100% 小条容易显得确定性过强，必须随 n 和完整图注阅读。可在后续采用时把区间写入图注/放大版、增加剩余关键标签字号，而不应加入不存在的分布或新装饰来填空。本复核不把“复杂度增加”视为科学或审美进步。

## 数值、单位与选择核验

对 `PLOTTED_VALUES.csv` 的 **1,221 个唯一键**与冻结来源进行独立核对；没有重算模型输出或新增检验。

- **a：**完整 396 个 RBP–cell–site 事件的 Full/Sequence 方向指示保留；287 位点、85 RBP。五臂 RBP 等权准确率及其全部 15 个点估计/双向聚类区间端点精确对应 `strict_test` 来源，而非错用全体或另一种区间。Full-only 43、Sequence-only 61、both correct 168、neither 124 合计396。最终图与图注均明确“四格正确性”，已清除旧 `joint agreement` 歧义。事件计数与 RBP 等权百分比不是同一估计量。
- **b：**全部 8 位点/13 绑定事件、52 个四臂标量与各位点分母保留；Full/Sequence/Flat/Perm 分别为 **7/13、10/13、12/13、9/13**。所有来源偏好均为 genomic REF，always-REF 为 **13/13**；11 个 training、2 个 validation、0 个独立 strict 事件。这些反证在图注及图内可见。小条采用共同 0–100% 域，表示每个位点事件的来源方向一致比例，不是疾病风险、临床预测准确率或富集检验。
- **c：**纳入恰好三个合格同终点 UPF1 reporter 基因后，在 Full/Sequence 两个相反不一致类别中按字母顺序取首个，得到 ACVR1B 与 SLC39A13。这是已写明的 **outcome-aware、事后、training 位点示例**，不是盲选。PDXK 与全部24绑定事件、66作者 reporter 测量保留在冻结来源/既有完整证据副本，没有用有利案例替换总体。
- 两案例全部 **14 个模型/固定输入控制效应、12 个作者原始 reporter 值**对应来源，包含反向及不利结果；七臂共用线性效应域、reporter 共用 0–2.6 域。ACVR1B genomic A>G，来源偏好 A；SLC39A13 genomic C>T、RNA C>U，来源偏好 **U**。RNA 上行与 reporter 第一组均为 source-preferred，预测是 preferred-minus-other，因此第二例 **U−C** 不能误标为 REF−ALT。最终图注明确此方向。
- ACVR1B 原始作者测量来自 Figure5e 表；SLC39A13 来自 **Figure6c** 表。HeLa qRT-PCR 是 mCherry/eYFP 归一化 RNA 表达、构建约700–1000 nt；模型是 K562 的101 nt 输入。没有把测量写成荧光或模型对 RNA 表型的定量预测。SLC39A13 来源 ASB FDR 约0.0586，不应另称 q<0.05 阳性。
- **d：**完整25蛋白、3臂、2终点的冻结值及各蛋白相对 Frequency 增益保留，顺序一致，负值未丢弃。四个 Full−Frequency/Full−Sequence 的总体均值及冻结 protein-bootstrap 95% 区间对应来源。Full−Sequence 两终点区间均跨零；连接线不是置信区间，底部区间不是逐蛋白不确定性。

## 不能由当前图解决的证据问题

1. 396 队列是先验 ASB-positive 且有开发曝光的方向评价，不含可靠阴性来评价 ASB 检出，也不是新的外部盲测。Full 在该队列弱于 Sequence 和常数 REF 规则；整体不能写成自然突变优越性。
2. ClinVar 致病性注释和 TCGA 重叠不构成新的绑定效应真值；13事件小且偏好全部REF，现有病例/对照不能支持保守性临床优势、致病机制或普遍功能解释改善。
3. 两个 reporter 案例提供实际变异与作者实验的对象衔接，但不能把一正一反的事后示例升级为保守性因果解释；Full/Sequence 匹配训练路径未证实的限制依然存在。
4. 本轮真实修复是证据组织与可比展示，不是新增阳性生物发现。若正文需更强自然变异贡献，仍需独立标签、匹配控制及预先冻结的队列和终点；继续美化不能补齐该缺口。

## 左右对照 PDF

保存后的 PDF 确认为 **2页、420×297 mm**。两页左侧均为完整 BEAPR 原图，宽178 mm、高约144.057 mm；第1页右侧为完整新试排178×151 mm，第2页右侧为完整旧v25、178×184 mm。实际 PDF 图像/表单对象宽均为178 mm（误差<0.001 mm），原始比例保持，没有裁切。实际两页 PNG 已查看；左右位置、正文内容、出处和范围说明均有显示。

## 本次绑定的 SHA-256

| 文件 | SHA-256 |
|---|---|
| `REFERENCE_MAPPING.md` | `688e315883ebc4833b5e9422f674b6dfeb45bdbe6c02ddf7601935a82afbd4ca` |
| `REFERENCE_SELECTION.json` | `d9dd044332f499581738b514b1d9d1d50922b51c7210c05e5187b0bd634a12f7` |
| `references/BEAPR_Yang2019_Figure5_complete.png` | `7f5e3ff8398d89ebf658392ff8651b43746e207c01ae49de288c51e5f191c3bb` |
| `layout_trial/render_trial.py` | `c32d97a2efd76f1a46fd1957ea1b7709c089843e6e491d4e67ad7a3cf8f66104` |
| `layout_trial/data.py` | `a2de92fe26074f36eacb4bdd5f15cd41a334ec30b0c12a2a6b2b5cc9b239e0cb` |
| `layout_trial/CARBON_BEAPR_layout_trial.png` | `ca8d0bc18da89874b9724f6d401647cebcb9705a9a572b1ae1894051ec4bdb8d` |
| `layout_trial/CARBON_BEAPR_layout_trial.pdf` | `ea427ad5f95c8d071f262a86b625beb23aa68e0d275e4cdbc9a8a5c1f9bc9234` |
| `layout_trial/CARBON_BEAPR_layout_trial.svg` | `8f93442c4a6547713f18c76791bbfe37218cc3ea6ef60531354eda4181c3531c` |
| `layout_trial/CAPTION.md` | `7ea2ab8f3621c64e516e3777fc507ec4acb83089e178490d60466df0c62874c5` |
| `layout_trial/PLOTTED_VALUES.csv` | `02af15bbdefff8c86d698255cff77d6c95c0c8687315b7b37138dd07404f205e` |
| `layout_trial/SOURCE_MAPPING.json` | `595e523ea39593c7b436a5df8f80810fed4fe2b85b375e125ede4587bd6b3f87` |
| `build_comparison.py` | `c33933350584ab1c96cb325c643590fe476401724d6969ff50b1e83426de9d5c` |
| `delivery/CARBON_Fig5_reference_left_ours_right.pdf` | `26f6b53f960cd3f45a8ae116d5d412fbea87c4d2d54d3e558a713aa12348f6d9` |
| `delivery/comparison_page_1.png` | `01db6d45a54aca00de234f1ba4025a0a0e3b8124714b5b7bafb68f72af3ff322` |
| `delivery/comparison_page_2.png` | `a30872026aba667129d7c157b93a397bcf3c06f71991e3a11d96f0108077811c` |

数据一致性和对照交付可供本轮审阅；最终图的审美选择仍由用户判断。正式主文、补充材料放置及整篇图文对应不在本次复核范围内。
