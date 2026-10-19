# 02外部偏好独立审查

审查对象是HNRNPK RNAcompete局部规则候选，而不是重写RBNS主阴性。只读02目录，自己的复算与新公共来源保存在03_functional_rules/cross_audit_02。没有改动执行者文件。

## 实际独立复算

cross_audit_02.py未导入执行者分析代码。从original420矩阵真实16,382个7mer独立重建occurrence-preserving 7→5均值，1024个词全部有47–48次来源贡献，与保存truth和分母一致。逐24个HNRNPK模型条件，从原position_ig以独立convolution重建每个窗口最高5nt归因、winner词与分数；全gene/locus权重为1时的24个top20 recall都精确复现原主点估计。再用直接per-window bincount实现首31个385簇成对重采样，精确复现保存draw；最后独立核差值和percentile区间。39项实际检查通过，结果/哈希在QA.json、INDEPENDENT_POINTS.csv和INDEPENDENT_BOOTSTRAP_SUMMARY.csv。

聚类重采样没有先将六模型ensemble再rank。实际是各fit各算rank/recall，再平均两细胞×三seed；这些条件不是六种独立蛋白。Full−S点差+0.075，95%条件区间[0.0167,0.0917]。Full−Mean+0.0333，区间[−0.0417,0.0333]；Full−Perm+0.0500，区间[−0.0333,0.0417]。后一原点高于bootstrap上端，复算确认不是聚合顺序bug：这个离散top20选择在重采样下明显向下偏（Full−Perm重采样均值约0.0050）。percentile区间只能解释为当前385个原生locus的规则提取稳定性，不能冒称已校准的跨蛋白优越性显著性。三必要Phi对照未全部成立。

## 真实来源与对象

[Ray等Nature2013](https://www.nature.com/articles/nature12311)配套S1 TSV中RNCMPT00026是人HNRNPK，保存构建序列长度实际464aa，三个KH区段44–104、146–209、389–451；这是实际查表，不是从蛋白名推断。当前original420 combined-Z来自[原作者MorrisLab数据仓库](https://github.com/morrislab/RBPbinding)。缺失两个7mer如实不补值；5mer标签是对包含该词的不同7mer Z-score平均，不能改称5mer直接Kd测量。

[ENCODE HNRNPK RBNS ENCSR368NMO](https://www.encodeproject.org/experiments/ENCSR368NMO/)及ENCFF715DIL target已经核对。实际下载并看了该实验附属Protein QC Data图，作者标expected size82kDa；原RBNS实验协议明确随机40nt池、含无蛋白input、多浓度pull-down。QC图及目前元数据没有氨基酸起止/序列，因此**未证明RBNS是某个截短构建**。不能用未经证实的构建差异解释RNAcompete正而RBNS负。能确定的是RNA池设计、读出（RNAcompete归一Z vs RBNS enrichment R）和实验条件不同；构建/测量目标差异仍只是待核解释。

原作者A/B对比文件真含`RNCMPT00026 0.924 0.4 0.54 0.731`，依次为作者Pearson及top10/100/1000重叠。A/B是同一RNAcompete池的两半，不是独立RNA生物研究。我们现有档案没有这两半的原始7mer矩阵；本轮没有对raw A/B重新算CARBON一致性。不能写“我们的结果在A/B独立复核”。静态核查Nature页面重定向、PMC页面captcha、EuropePMC全文XML500都已记录，不绕门获取；这不影响已有真实S1与combined矩阵的结果，限制的是新增A/B确认。

## 简单基线和解释范围

独立重算C-count baseline：按词典顺序断并列，top20 recall=25%；按所有边界并列均匀分配，期望=12.8889%。两种结果都应该留，不能只挑有利的。C-count在全部1024词上与实验标签Spearman=0.81036，明显高于Full IG约0.4018。这说明Full在极高偏好词的提取存在局部候选价值，却尚未胜过简单生化组成对全词偏好的描述。

controlled substitution在同1024宇宙固定16背景，对native log_total_count做central5nt替换；Full全词ρ约0.6482，Mean约0.6716。composition-fixed concordance Full约0.6193，S约0.6240，Perm约0.6286，频率约0.6327。该诊断未建立组成之外的Phi优越性。C-count在同组成配对全为预测tie，因此按同样规则记0.5；不让C-count并列享受不同处理。

IG目标是native log_total_count；序列路径quarter-base→真实序列，Phi固定为各条件值，5nt signed序列IG滑动求和，不是绝对IG或物理接触。执行者已有72条件CPU/CUDA与32/64数值QA，其记录可以作为数值稳定性来源；本审查独立重算的是保存真实IG的规则/统计，没有再发GPU查询。不以可视化好看替代这些门。

## 可发表判断

作为探索性补充或“模型家族学到部分外部序列偏好”的候选，证据有价值，主RBNS阴性不能自动抹去它。若主文使用，必须与RBNS、Mean/Perm、C-count及组成控制共同展示，并称一个HNRNPK assay中的局部top20恢复。当前不能作为“保守性普遍提高解释正确性”“功能机制得到验证”或跨蛋白独立验证的主结论。下一轮最有价值的是实际A/B标签与另一独立HNRNPK偏好来源/匹配构建，而非继续重排同一矩阵的评分。

## 最终绘图与新25蛋白主分析复核

check_02_figure.py独立读取原24模型IG、原作者矩阵与完整统计表，核对250个a/b点与共同坐标、80个c选择格、20个外部Z5/顺序、6个d点和区间及源hash，共620项实际检查全通过，最大差9.77e−15；当前共同navy色标PNG SHA46e838e…，实际85mm PDF已查看，无裁切/碰撞。

check_02_family.py进一步从原作者16,382行TXT重新构建全部assay的5mer真值，从168原response.npz重新读背景均值，核全词ρ/15,240同组成候选配对、每蛋白等权聚合及10,000成对protein bootstrap。新25蛋白主结果及所有necessary contrasts完全复现：Full−frequencyρ+.146187[.038017,.247997]，同组成+.037809[.008987,.067330]。此family证据与局部HN的一蛋白归因证据分开：family支持外部规则能力，Sequence也能恢复，仍未支持Phi独有。

初版严格复核2631项中62微差全部只在已曝光RBFOX2以及含它的all28/exposed3聚合，不在新25。原作者RNCMPT00818 TXT默认pandas浮点解析与Pythonfloat差≤9.02e−17，使ATCCG/GTAGA的5mer汇总发生一个并列变化：默认解析均为−.12273739705577145(rank445.5)，literal分别−.12273739705577143/−.12273739705577145(rank446/445)。最大ρ差3.88e−6；这不是外部实验的真实等值判定。旧严格fail输出保留，已通知执行者以明确literal/round-trip解析契约修复低成本support；不能把最初2631项冒称全PASS。修复后最终读回另存。

执行者保存83项旧parser表/脚本快照后，统一author TXT和derived CSV literal/round-trip契约，没有新推理。2026-10-04 13:55 UTC，独立check_02_family.py从同原始字串Pythonfloat重新运行，2631项全通过，max1.78e−15。新25主值完全不变，已曝光RBFOX2精度并列与支持聚合也精确闭合。旧62fail/原QA/原绘图检查共10文件先存03/cross_audit_02/revisions/parser_before_20261004带hash，不覆盖失败记录。最终FAMILY_QA.json现在指向修复后的真实来源及全PASS；这是明确解析契约的数值修复，不是换终点、抹去阴性或增添生物证据。
