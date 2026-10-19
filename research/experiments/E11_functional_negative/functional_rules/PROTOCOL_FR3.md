# FR3：结合归因能否找到真正影响核保留的关键碱基

2026-10-04，在新核定位标签关联产生前保存。JPX所有binding/IG结果已曝光，本分析是新功能终点、仍单RNA探索，不能称独立RNA确认。FR1/FR2功能全库迁移阴性保持，不在它们上更换评分。

来源EMBOJ2021 DatasetEV1/NucLibC真实四个N/C重复及五个HNRNPK RIP重复；严格同WT、109nt单替代的全部150 SNV、原0-based55..104共50位置，每位置三ALT全部保留。源码来源/WT序列/变异ID用此前真实labels CSV只读核验，再逐变体连接原Excel来源序列；不靠名字假定正确WT。

主功能位置敏感性E_NC[i]=mean_three_ALT abs(mean(N/C_ALT)−mean(N/C_WT))。有方向次级L_NC[i]=mean_three_ALT(mean(N/C_WT)−mean(N/C_ALT))。这分别是突变敏感程度和WT对核保留的有利贡献，不是RNA binding、表达变化或基因KD。对应binding E_RIP/L_RIP独立报告，避免将RIP值改名功能值。

首轮解释方法固定复用已有32step WT count IG、crop4、四核苷酸通道归因绝对和，不重新挑积分目标、crop、seed。全部旧A_native/S/A_flat/permutation0..2以及可用NoGRU/NoSE真实trained模型在两HNRNPK细胞同时报告。SLTM两细胞旧IG全部为固定次级，因为来源功能需要两蛋白作用，不挑更好的蛋白顶替HNRNPK主问题。simple baseline为核心CCTCCC/GCCTCCC覆盖、局部7nt C/CT比例、Phi/absPhi、窗口中心距离；来自旧固定POSITION_MAP或WT重构，与新功能标签无关。全50位置Spearman、top10功能敏感量fraction（边界平均分配ties）、全部49种circular位置错位，对应B库共有25位置重复。B/C共有同一个JPX，不增加RNA n。

必要保守性比较A−S、A−flat、A−mean(permutation rho)，不能只展示胜的对照；NoGRU/NoSE只有旧single checkpoint，训练重复/预算未核齐时只能描述模型组件探测，不当严格机制归因。Bootstrap位置不是跨RNACI；主要结论同时给错位负对照和当前位点实测重复一致性。

第二轮固定补强（无论首轮正负都全部报告）：使用当前matched HNRNPK两细胞四臂×三seed，对WT和全部150ALT产原生count+同定义WT IG；输入从旧精准WT109Phi及ALT sequence重建中央101，ALT保持WT reference Phi，绝不构造ALT专属evolutionary信息。每模型WT与ALT共用本模型相同Phi排列。count的主位置敏感性E_pred=mean_abs(ALT−WT)、signed loss单列，不重跑每位置ALT顺序作为功能主终点。WT IG算法同32step；模型seed是模型条件，不是生物重复。

固定功能/结合分离诊断：全50位置partial rank correlation of I with E_NC controlling E_RIP，至少明确这是同RNA测量误差下的关联而非因果独立。完整位置图将NC与RIP分两行，不能先找模型高分再裁掉bind-only弱功能位置。

如果所有模型输给核心基线、跨读出/复库不稳定，则如实结论功能关键位置的解释优势未建立；既有归因少数局部提升仍须标为探索，不扩展为全模型/跨RNA/保守性普遍优势。
