# 独立公开来源与构建资格审查

2026-10-03。本审查直接读取原文、原始图、补充 PDF/XLSX 和 GEO counts/frequencies；没有运行模型，也没有读取本轮新模型输出。写入范围仅为本目录。已阅读 `03_protocol_review/REPORT.md`、`ENDPOINT_SPEC.json` 和 `01_public_binding/CHRNA1_CASE.json`。这次判定来自新来源核验，不是沿用上一轮的自查结论。

**rG4 面板可作为新的外部、逐构建直接 G3BP1FL 结合扰动参考，主端点应锁定 biotin/+guanine 条件。** 六条天然 90nt 均有完整饱和 SNV 覆盖，连续弱效应全部保留。它不是已核准训练未见的确认集，也不是逐碱基接触、唯一 G4 结构核心或精确亲和力的金标准。CHRNA1 可作保持上游 PPT 序列的背景改变导致复合物 PTB 招募降低的一个定性病例。

## 1. rG4：应使用哪一个实验端点

来源为 NAR 2026 `10.1093/nar/gkag821`，GEO `GSE296171`，本地 `02_public_context/sources/rG4_2026.xml` 与补充表。原 Fig5 实际图及补充 FigS9 的实际 PDF 第12页已经观看；本目录保留了 PDF 原页渲染，未重绘原作者实验图。

| 角色 | 精确公开列/条件 | 独立核验结果 |
|---|---|---|
| 主连续扰动端点 | `Mutation deltaWT Scores` 的第12列（0-based 11）：`DeltaWT Score G3BP1 FL (biotin) (+ Guanine)` | 所有1889设计行均等于 ALT 减同基因 WT，最大差 `9.96e-15` |
| 对应绝对结合分数 B | `All Sequences Processed` 的第11列（0-based 10）：`G3BP1 FL (biotin) log2(R) (+ Guanine)` | 所有6000条均可从频率表重现，最大差 `5.11e-15` |
| 原始配对重复 | `G3BP1_FL_WT_PrepA`、`G3BP1_FL_WT_PrepB`，共同 `Input_WT` | 两次独立样本制备；不是三个 RT-stop 重复 |
| 不可代替主终点 | G3BP1FL 7dG、G3BP1 RRM/LCD、FMRP、RT-stop/G4 stability/PDS | 不同蛋白、构建、RNA化学或测量问题；应分别报告 |
| SDS公开 ΔWT | 第14列（0-based 13，末尾有空格）：`DeltaWT Score G3BP1 FL (SDS) (+ Guanine) ` | 与相应 absolute SDS B 的 ALT−WT 有内部不一致；最大差4.43，不能当主结果的重复确认 |

主分数精确计算为：

```
B(s) = 0.5 * [log2(Freq_G3BP1_FL_WT_PrepA(s)/Freq_Input_WT(s))
              + log2(Freq_G3BP1_FL_WT_PrepB(s)/Freq_Input_WT(s))]
ΔB(s) = B(s) - B(same-gene exact WT)
```

这是两个 log-enrichment 的均值，不能替换为先平均频率再取 log。正 ΔB 表示更强富集/结合，负值表示更弱。原图 Fig5B 的结合池与 input 频率比支持该方向。Methods 的主条件是纯化全长人 G3BP1（Q13283-1）200nM、RNA pool 500nM、常规 guanine RNA 在150mM KCl中折叠，使用 SBP/biotin 洗脱；后续 LiCl RT 旨在减少结构引起的读取缺失。这是此条件下的单浓度纯化蛋白结合富集，不是 Kd。

Domains 的 RRM 340–415/LCD 416–466 是1µM条件；鼠 FMRP LCD 478–614 为200nM。相同序列的7dG改变化学及相关折叠条件，现有 DNA/RNA字母+Phi 输入不能表示这些差异。补充 FigS9 点的红蓝颜色来自 RT-stop ΔWT 的 G4稳定/不稳定分类，不能将颜色当 G3BP1FL 增强/降低标签。

SDS 不一致是公开 processed table 的内部核验问题。本审查不能确定它源于 WT baseline、列转换或其他处理；不据此否定实验，不修补其数值，也不将它与 biotin 主端点混合。

## 2. 全覆盖、去重和精度

WT 由补充表1中六条明确90nt天然序列精确匹配，不能依 `1mut`/`AAAmut` 名字排除：有五条 WT 的 ID 仍含 `1mut`，ARHGEF12 的 WT 含 `AAAmut`。主表1889逻辑设计行对应1833个不同实测90nt序列，其中6 WT、1827非WT；56行为设计别名重复，不能给同一构建重复投票。

| 位点 | 独特非WT构建 | SNV | 额外多碱基AAA | 两重复效应符号相反 | 非WT构建中位数 \|ΔB\| |
|---|---:|---:|---:|---:|---:|
| VEGFA | 346 | 270 | 76 | 54 | 0.14814 |
| U2AF2 | 340 | 270 | 70 | 49 | 0.07798 |
| ARHGEF12 | 331 | 270 | 61 | 203 | 0.17243 |
| EWSR1 | 270 | 270 | 0 | 15 | 0.28102 |
| HNRNPA2B1 | 270 | 270 | 0 | 41 | 0.13346 |
| PABPN1 | 270 | 270 | 0 | 10 | 0.65065 |

每条RNA的90个位置、每个位置另三种碱基均存在，缺失0；该覆盖按真实序列差分核验，而非按命名或作者注释推定。只有前三条有额外 AAA 扰动。若 `AAAmut` 序列实际上仅改变一个碱基、或完全等于WT，按真实差分归为 SNV 或 WT/no-op，保持逻辑设计别名记录但不增加构建n。全部多碱基项确认不超过相邻3nt且改变后为A。

`RG4_RAW_ELIGIBILITY.csv` 保存所有1833行、精确 WT/ALT、实际改变位置、别名数、主 ΔB、两重复的原始 count 配对效应和 read counts。原始计数采用公开+1伪计数重算：

```
δ_r = log2((ALT_bound_r+1)/(WT_bound_r+1))
      - log2((ALT_input+1)/(WT_input+1))
```

配对差中池总量归一化因子消去。公开频率经过舍入，因此 raw-count 重算和精确 published processed ΔB 的小差不被覆盖。主分数来自精确处理表；重复符号冲突作为精度诊断完整保留，不能事后用此筛选成功面板。ARHGEF12 两重复的符号冲突尤其多，方向结论应谨慎。两重复并不足以凭本表创造稳定的等价界限或“中性准确率”；接近零/符号不稳项目不能宣称实测中性。

后续效应排序、方向或定位分析须纳入全部技术合格构建，并在每个基因内汇总后对位点等权。1827 variants、90 positions、2重复、seed/window/Phi控制都不是生物n。最多6个来源基因；原作者依据 G4性质选的天然RNA也不是全转录组随机样本。来源论文较新和实验证据外部，不自动证明 CARBON 未见这些坐标/序列。

## 3. native90、完整probe与“core/flank”的合法范围

补充表1的全部1833个实测天然系列构建均为132nt合成probe：

```
5′ common adapter (21nt)  GAGUUCUACAGUCCGACGAUC
native or mutated human sequence (90nt), full-probe indices [21,111)
3′ common adapter (21nt)  UGGAAUUCUCGGGUGUCAAGG
```

因此天然90nt是实测 RNA 的可映射主体，不能叫“90nt裸RNA测量”。固定 natural101 可容纳全部90nt和所有改变，但额外11nt来自基因组延展并替代了原实测的42nt技术尾巴，仍是局部模型范围。天然 reference/Phi 无法合法标注合成adapter；不可发明adapter Phi，不能把输入覆盖全90nt称完整132nt复现。若实际结构/结合受adapter影响，测量与101输入的差异仍是替代解释。

本目录用输入序列中全部 `G{3,}` runs 定义**潜在G-tract**，另区分tract之间的 loops 和 tract hull 以外的区域。U2AF2 的原 Fig5J a–g 七段与独立序列run核对一致，0-based half-open为 `[31,34),[36,40),[42,45),[47,51),[54,57),[66,69),[70,74)`。该mask不利用模型或 binding 标签。

保持这些G段只能叫 potential-G-tract-preserved，不能叫已证明的结构核心保持。原文讨论多种可用G段/拓扑；loop变动也会改变折叠、可及性或蛋白接触。EWSR1/HNRNPA2B1/PABPN1仅有1/2/3段 G>=3 runs，不能强行设定标准四段 rG4 核心。跨天然RNA也不能把不同run数目或不同baseline绑定解释为同一个核心在不同flank的因果比较。

这个面板真正强的是**同一构建骨架的已测结合扰动响应地图**。主可检验量是实测 ΔB 与 count-head WT→ALT 差的方向/位点内排序，以及保留所有实测弱效应的区域效应排序；它不直接提供接触位置、folded-core结构、background恢复或四格交互证据。RT-stop可作为分开的结构辅助标签，不能替换 G3BP1FL 主外部终点。

## 4. CHRNA1：真实图与补充构建核验

原文 Rahman2013 `10.1038/srep02931` 的实际 Fig5c 已观看：lane1 为 beta-globin-MS2 control，lane2 为 iP3Ai-WT-MS2，lane3 为 iP3Ai-mut-MS2。PTB 在 WT lane2 可见，mutant lane3 弱/未检测；hnRNPL、hnRNPLL及 U2AF65/U1 的差异也在对应行。作者明确叙述 WT 的 PTB association 而 mutant 缺失。Methods 是 HeLa nuclear extract 中 MS2-MBP/amylose 分离 EDE 复合物，不是纯化 PTBP1 加RNA测 Kd。

补充 PDF 实际第7–8页引物被独立核对。天然 forward tail `GTTCACACTAGGTCTTTCTCC`；reverse tail `AGGGTGATTACTGACCTCATTCT`，RNA方向反向互补为 `AGAATGAGGTCAGTAATCACCCT`。依据提供的 GRCh38 gene reference 独立定位，129nt insert 唯一匹配负链 `chr2:174757931–174758059`（1-based inclusive）。WT/ALT只有 insert index61（1-based62）G→A，基因组 `chr2:174757998 C→T`，与源病例一致。

101nt anchor 位于天然insert index11开始，改变位于 anchor index50；因此遗漏天然上游11nt、下游17nt，另有未完整重建的vector 5′leader和三个MS2 hairpins。实际测量不是129nt裸RNA，更不是101nt probe。原42nt exon affinity探针属于另一实验，不能用它替换 Fig5c 的 MS2 构建。

37nt上游PPT `TTTCTCCTTTTCTGTGGGTGGACAGGGTGACATGGTA` 位于 insert `[14,51)`、anchor `[3,40)`，在WT/ALT完整保持。机器case的25nt部分core `[3,28)`也保持。改变距37nt片段最后一碱基11nt，作用于下游hnRNPL/LL相关位点 `CACGCCC→CACACCC`。

**判定：core-preserved 的主张在“上游PPT序列保持、辅助结合背景改变而复合物PTB招募下降”这个范围内成立。** Fig5c单独并未建立PTB在PPT的直接RNA接触，PTB/hnRNPL合作机制另由Fig4 PRR/相互作用实验支持。不能写作“PTBP1纯化亲和力因flank改变而下降”，也不能说保留PPT使直接亲和力本应不变。定性 ALT−WT 招募方向为负；没有densitometry、配对误差或等价阈值，不能给精确减幅、零占据、中性标签。一个WT/ALT、一个locus、不同细胞背景与不完整runoff上下文限制泛化。

## 5. 执行与交付状态

- `audit_rg4_raw.py`：可运行源资格/所有SNV覆盖/实际差分/WT/主列/raw重复/processed频率重现/探针范围检查，无模型调用。
- `RG4_RAW_ELIGIBILITY.csv`：全部独特实测序列，无效应筛选、无虚构中性标签。
- `RG4_SOURCE_AUDIT.json`：所有来源SHA256、端点精确列、6000频率及1889ΔWT核验、6位点统计。
- `audit_chrna1.py`、`CHRNA1_SOURCE_AUDIT.json`、`CHRNA1_source_excerpts.txt`：独立引物/参考/构建/core检查与源原文摘录。
- `CHRNA1_supp_actual_page7.png`、`page8.png` 和 `rG4_supp_actual_page1/8/10/12.png`：从真实补充PDF渲染的QA页。

当前没有核模型预测、Phi数值、训练成员或joint OOD；这些不由源资格审查自动成立。后续只读独立数值复算必须明确引用冻结协议/输入/hash和原始模型输出，不能用源核验代替结果审计。

已有证据的问题是：rG4多个位置的真实序列改变能改变纯化全长 G3BP1 的结合富集；CHRNA1保持PPT序列时仍可改变复合物PTB招募。仍待证实的问题是：哪种结构/辅助蛋白造成这些变化、Phi是否帮助未见loci的模型解释、完整probe与局部101nt是否给出相同关系。区分它们可在固定天然位点测完整probe与裸90/101nt的配对结合，并对预定G段保持的loop/flank改变同时测结构和结合；CHRNA1可用纯化PTBP1、加入/去除hnRNPL/LL及同一完整RNA测招募，分开直接亲和力和合作效应。不能凭模型head或Phi重排就宣布找到机制。
