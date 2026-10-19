# Figure 2/4 visual repair, v32

仅在本目录制作可编辑图稿；没有改主稿、历史原图、原数据或线上项目。按主Agent已完成的离线preflight沿用当前PaperSpine，不重启工作流。

## 实际改动

- 图2：85 × 185 mm → **85 × 168 mm**，高度减少17 mm（9.2%）。a/b原始尺寸及内容保留，收紧汇总表行距和区域间距，c–f图框由25.5 mm方形改为22.5 mm方形（字体、数值范围、点数、极值选择不变），g热图高度24 → 21.5 mm，h轨道高度/间距小幅压缩。完整a–h仍为一幅复合图。
- 图2b：纵向小提琴、箱线、颜色、alpha、带宽、顶点、轴范围、字体、图例、物理图框均保留。小提琴顶点JSON与v9逐项相同。固定物理裁图的300 dpi平均像素差仅0.014/255，差异集中于PDF变换舍入的抗锯齿边缘；不声称PDF字节或所有像素相同。
- 图4：保持 **85 × 139 mm**。所有定量热图的蓝/红端点统一混入24%白色（新端点约#6b7da5/#ec7665），中点和所有数值→归一化函数保持原版；输入灰阶同样减轻。观测/模型均调用同一函数和同一色尺。没有按模型单独减色、重标定或改变gamma。实测大面积填色alpha为0.55并保留深蓝边线；Full预测线保留NPG红#E64B35，实测身份保留深蓝#19499C。序列碱基颜色、预测/PhyloP线、数值及位置不变。

## 科学内容与映射核对

图2仍包含222任务、150 RBPs、666环形PCC值、每幅散点图222对数据、222 × 12差值矩阵和六条冻结轨道。环形共有power gamma=0.65、PCC[-0.1,0.5]及pooled median白点0.285；g仍为原始差值±0.12，195个饱和值保留。六指标paired median/wins、所有不利比较、极值任务与原轨道数值保留。

图4保留冻结的1,851对、659基因、865 cell–gene行（K562 412、HepG2 453），无筛行、重排或补值。4个865 × 22数组均与归档字节一致；已从冻结pair conditional数组按原pair_rows重建equal-pair gene均值，四臂最大误差<1e-14。原始观测22 nt支持内条件化；预测先三seed均值再对两支持条件化；最后等权pair→gene平均。输入仍是相对于101 nt全窗总input的独立灰色量。

共同mass色尺仍为gamma=0.5，锚点[0, 0.03956027423918412, 0.2]，汇总76,120值、701个高端饱和值。原始support-1 q/p仍以共同[0,0.5,1]线性色尺，两个散点面板保留全部865点及identity线。Raw atlas没有改成input-adjusted主终点：input-adjusted MSE与28.73%结果仍在原S20/正文协议中。冻结adjusted error数组、1,851对映射、Full误差更大的280行（对NoGRU）/275行（对Local+pos）及OCRL失败例全部保留。ZMYM4的主图a为单pair，c中标注行是两pair平均；没有混同。

## 风格参考与实际像素检查

实际查看了v31图1–4及local style reference。以现有图1/3为本论文的配色/字体语法：白底、Arial、细灰轴线、Full红、实测深蓝及固定碱基色。local style reference仅提供紧凑层级和连续复合图结构的辅助参考；未转移任何数据、结论或生物实体，本轮也没有为模仿其结构增删证据。

已检查新图96/300 dpi完整像素及每图85 mm前后对比。首版出现的底部Position标签裁切和c/d字母靠近cell legend的问题已在最终源中修复；最终文本均位于MediaBox内。最小普通标签5.5 pt，图2底部文字余量2.44 mm，图4为0.78 mm。图2散点轴字没有因紧凑化缩小，g每行在600 dpi嵌入图中仍有约2.3像素。弱信号仍可见，共有浅化保持颜色相对关系；近共同白点值本来就应浅。

`BEFORE_AFTER_85mm.pdf`每页完整显示旧图/新图，两个版本均实际放置85 mm宽，适合打印/100%检查。`fig2_v32_85mm_96dpi.png`与`fig4_v32_85mm_96dpi.png`是96 dpi屏幕参考，不能替代主稿放置验收。主稿整合与独立审查由主Agent完成，尚未声称已投稿或已线上更新。

## 文件与重建

- `fig2.pdf/png/svg`、`fig4.pdf/png/svg`：供整合的新图。
- `render_fig2_v32.py`、`render_fig4_v32.py`：可编辑最终绘图脚本，只读取本目录frozen_inputs并写本目录（图2可指定destination）。
- `frozen_inputs/`：全部所需冻结数据、行映射、选择记录和配色。
- `verify_and_compare_v32.py`：核对、85 mm比较PDF/PNG及QA输出；以本机原归档路径核对输入哈希。
- `QA_V32.json`、`FIG2_DISPLAY_V32.json`、`FIG4_DISPLAY_V32.json`：核对详情。
- `CAPTIONS_UNCHANGED.tex`：原图2/4 caption原文；本轮只有几何/ink变化，不需要更改科学图注。
- `source_baseline/`和`baseline/`：保留旧源/旧PDF用于可追溯比较，历史脚本作为.reference.txt，不能当新版入口执行。

在本目录运行：`python render_fig2_v32.py --preview`、`python render_fig4_v32.py`。依赖Python、numpy、matplotlib、PyMuPDF；无需pandas、训练或推理。
