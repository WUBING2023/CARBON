# 公开PTBP1背景机制数据筛选：交付与结论

2026-10-03。仅来源获取、构建与测量审核，本组未做模型推理、训练或新湿实验。web-access已读取并完成依赖检查；访问公开静态原文/补充，不操作用户浏览器或登录资料。先读上一轮SOURCE_LEDGER、CONSTRUCT_SCOPE_REVIEW与FINAL_REPORT；已知旧USP5/ITSN1/TCF3/AXL和ASB总体阴性不能变成新确认集。

## 可实际推进的新证据

最佳新来源是**人CHRNA1 P3A23′G>A**，Rahman等2013，DOI[10.1038/srep02931](https://www.nature.com/articles/srep02931)。上游PTBP1嘧啶富集区保持，较远的外显子G>A改变hnRNPL/LL的招募，Fig5c中从HeLa核提取液组装并MS2纯化的EDE复合体，WT含PTB而ALT的PTB显著弱/不见。它支持“同一个PTB核心因周围辅助序列/合作蛋白变化而改变结合”的真实功能机制。这里的外部标签是复合物中PTB招募方向，**不是纯化PTBP1 Kd、不是剪接PSI预测、不是hnRNPL的测量代替PTB**。

`CHRNA1_CASE.json`可直接交根复核和按既定协议计算。Supplement Methods pp7–8给出了iP3Ai-MS2天然片段的两个PCR边界；二者在下载的GRCh38 CHRNA1基因minus序列各唯一命中。还原自然insert129nt，chr2:174757931–174758059，变位为insert62。人工MS2和vector来源的整个runoff RNA尚未全序列重建，不能称原作者测了101nt RNA。

- GRCh38变位：chr2:174757998 C>T，转录方向minus，对应RNA G>A、P3A外显子23。
- 原自然核心保持：上游25nt intron3片段`TTTCTCCTTTTCTGTGGGTGGACAG`完全相同。嘧啶富集前13nt距离突变35–47nt。
- 101nt模型自然输入中核心完整位于0-based `[3,28)`，突变位于50。其范围chr2:174757948–174758048。
- 与天然位点的匹配仍应由根使用现有本地hg38/链向/Phi独立复核；本组没有访问模型输出。

该来源只有**一个独立locus、一个主要WT/ALT PTB结合方向**，没有方向定量的原始重复数据或CI。HeLa核提取液与CARBON HepG2/K562任务不同，原基因组的训练/验证覆盖本组未核定。短窗口可能遗漏129nt及辅助蛋白的生物作用。可作为新的文献驱动探索性外部方向案例，不能单独证明保守性具有可泛化解释优势。若A/S均正确，或flatten仍正确，不能据“Phi有反应”宣称保守性提高正确性。

## 其余来源的真实价值与缺口

完整六行资格表见`CANDIDATE_LEDGER.csv`，成功下载原件见`sources/`与哈希清单。

1. **Amir-Ahmady2005**，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC1370756/)，Fig1E B/E/F/G最直接地保留两段CUUCUCUCU核心而改上下游嘧啶。E/F单侧改动对Complex1小；G双侧改动减少Complex2且降低Complex1亲和性。WT EMSA的50nM不能搬作所有mutant的Kd；Complex2代表不同化学计量，不能等同全窗口count。这是科学上最贴切的背景面板，但本轮来源Methods尚未证明src RNA为人，也未取得可逐碱基核准的Fig1A完整B/E/F/G序列。全文通过公开原文/作者存档的解析文本阅读，PDF下载202空响应、PMC副本reCAPTCHA，因此**不冒充已保存原图/完整输入**。暂停模型使用，不把human SRC的Phi硬贴上去。

2. **2014 de novo PTBP1 HMM**，[原文](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003442)，纯化人PTBP1对13个69nt **mouse** RNA的EMSA有真正Kd面板。`DENOVO_KD.csv`抄录Fig2A表，未把no binding补成零，probe4仅20%binding不能拟造Kd。Fig2B只公开四个完整probe序列；S3B虽有全部13凝胶与曲线，**没有列全部13序列**。不同RNA且HMM选择，并非核心保持的成对干预。不具备天然人Phi资格，不能作为当前确认样本。Methods还将RNA85°C变性后冰冷；不能把额外结构假说当成已做结构补偿实验。

3. **U1 SL4 Sharma2014**，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4233244/)，Fig7A/B实测HeLa核提取液PTBP1 RNA-pulldown。WT/M10/M10h三个人源24nt构建完整，另有果蝇26nt；见`U1_CONSTRUCTS.json`。M10保留嘧啶并改变茎，PTB相对WT图示约1.8；M10h嘧啶失去bulge呈现，约0.6；果蝇低于0.1。这些仅是图读近似、无原始误差，不能称M10严格中性。天然U1多拷贝/伪基因的唯一坐标与Phi未解析，snRNA机制与mRNA预测任务相异；目前为结构机制参考，不能硬选一个最高Phi拷贝或充作4独立loci。

4. **CHRNA1 Bian2009**，[原文](https://academic.oup.com/hmg/article/18/7/1229/673533)，37nt WT和IVS3-8G>A SPR构建可完全读出并映射同一CHRNA1。作者说PTB无效应但数据未展示，没有等效性区间/Kd/CI，不能作为实测中性正确性标签。可核对上游PTBP1区与天然坐标，且**与2013只能计同一个gene/locus家族**。

5. **ANXA7 Llorian2010**，[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2933513/)，Fig6c/d/e有真实PTB UV-crosslink、EMSA与抗体supershift，突变降低结合。原作者的三个外显子突变故意破坏PTB motif，因此不满足“核心保持、背景改变”。原图取得并看过，但低分辨率完整序列尚未逐碱基核准。可后续发展新的普通结合方向案例，不能归为本次上下文机制确认；KTN1背景插入/删除的主要标签是剪接，不替代直接binding。

## 原图与数据质量复核

实际视觉检查：CHRNA1 Fig5c的PTB条带与WT/ALT列、2014HMM Fig2A/B及S3B全部13凝胶、U1 Fig7比值、ANXA7 Fig6原图。CHRNA1补充PDF实际打开提取15页，pp7–8构建引物核准；fullXML与补充一致。U1图中M10更强，不能用“bound well”歪写等效不变。HMM probe4的20%binding保留未识别Kd状态。

`SOURCE_HASH_MANIFEST.json`将reCAPTCHA HTML及OA404副本标成failed_access_artifact_not_original。成功原件有PDF magic/XML根/实际图片与hash；派生PNG预览和PDF文本单独标derived。可运行脚本为`acquire_sources.py`、`fetch_followup.py`、`fetch_final.py`、`fetch_primary_figures.py`、`map_chrna1.py`、`build_delivery.py`。来源获取受限没有扩大权限/安装/人工联系。第三次获取日志的早期Ensembl长度门槛失败不意味坐标不存在，后续真实40344nt序列已留哈希。

目前可完整落地的是CHRNA1天然129nt WT/ALT的**单点、核心保持、合作招募方向案例**。没有找到可合法立刻用humanPhi完整测试的多locus纯化PTBP1核心保持/补偿面板。本组没有把来源不足包装成新强证据；最终“保守性是否更正确”须看根的全部预定A/S/native/flat结果及失败标准。
