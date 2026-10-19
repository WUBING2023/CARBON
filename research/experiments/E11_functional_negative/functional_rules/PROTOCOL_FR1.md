# FR1冻结前协议：结合预测能否识别蛋白依赖的核保留

登记：2026-10-04。本任务旧binding预测与14RNA框架已曝光；核定位/KD终点为本波新分析，不称盲测或新实验。分析前不按KD效果、binding强弱、Phi、预测表现筛选RNA或tile。

来源：Lubelsky等，EMBO J 2021，DOI10.15252/embj.2020106357，Dataset EV1作者NucLibB（真实MOESM5 Excel），同109nt AcGFP 3′UTR reporter的Nuc/Cyto、HNRNPK/SLTM RIP与siRNA条件。siNT是非靶向siRNA；HNRNPK、SLTM、SNRNP70敲低后在同一tile库的Nuc/Cyto。三重复列是实验测量，不是独立RNA、训练种子，也不伪称同个reporter分子成对测量。按condition平均再相减；核查作者效应列方向/公式，冲突保留。

## 资格和层级

主范围：旧冻结INPUTS中作者WT=True、NucLibB、109nt精确匹配指定人参考RNA且Phi完整、sequence与本轮Excel构建ID完全一致的全部tiles；模型101nt中央固定crop4:105。原测量物理insert为109nt，转移与crop限制保留。须至少每condition一个finite replicate，主reporter绑定来源不以大效应/显著性筛选。完整三rep条件的固定敏感性另列。去重以gene+109nt序列/坐标实例为键，重复同序列测量同condition等权合并，不增加n；重叠tiles仍按RNA簇汇总。RNA相关计算需>=5唯一tile且终点/预测非恒定；所有资格和undefined保留，不为让平均好看丢弃RNA。

参考匹配不证明原物种。14参考RNA中JPX/PVT1是来源可证人源子集，单列，不用其余大写gene推人源。此次gene-level RNA knockdown读出不与tile reporter混配。

## 主终点与必要对照

主功能依赖D_H = mean(siNT N/C log2) − mean(siHNRNPK N/C log2)，正值为敲低该蛋白后核保留减弱。预测变量是HNRNPK native log(1+count)输出，测试与D_H的RNA内Spearman，各RNA等权宏平均，两checkpoint细胞分别全部报告。不直接把binding输出称定位概率。

所有相同资格比较：A_native、sequence-only、A_flat、固定3个permutation prediction mean；模型seed不作为生物重复。已有single-checkpoint输出先真实重放，明示缺少独立trained Mean/Perm；随后若CPU预算允许，用匹配同预算4trained arms×3seeds HNRNPK/HepG2/K562推断相同冻结输入，作为独立模型控制补强。不得选择最佳细胞或seed。

简单基线固定：C含量、CCC重叠次数、CCCC重叠次数、CCUCCC（DNA CCTCCC）重叠次数、GCCUCCC（GCCTCCC）重叠次数；按同一101nt输入计算；native meanPhi只作外部先验基线。实测RIP为生物正参照，不能叫独立模型基线或理论上限；原始N/C也另读。

主比较A−S、A−flat、A−permutation_mean，以及A相对C-content与motif规则分别报告。RNA簇bootstrap 10,000次，seed20261004；这些RNA同一研究、重叠tiles和共享输入，区间只描述当前库外推而非独立研究总体。主要必要比较三项Bonferroni同时区间；不以p值迫使阳性。

## 预先固定次级诊断

1. D_SLTM、D_SNRNP70与baseline N/C各自独立读，绝不因D_H失败换主列。
2. HNRNPK、SLTM预测对其自身与其他protein KD的相关，区分复合体共同依赖与蛋白专一，不能要求人为严格互斥。
3. 按RNA内rank后，用截距+实测RIP rank残差化预测和D_H，计算partial-rank correlation；需>=8tiles且rank非恒定。它诊断是否仅重述结合量，不把控制噪声RIP的残差解释成新机制。另给baseline N/C控制，因为D含siNT时可能有共同读出/回归到均值。
4. reporter WT功能整体排序与核依赖是一项任务；单RNA突变、双变、核心×背景作为FR2候选登记后另分析，不能把同一RNA多个construct当独立RNA。

所有numeric/results/失败和成本记录保留。不会把来源已有protein KD计为本项目新做湿实验，也不会把A对必要Phi控制未优势包装成保守性独有解释提升。
