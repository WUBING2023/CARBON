# FR1：蛋白结合预测与蛋白依赖核保留

**当前结论为阴性。** 这项真实功能检验完成了1169个独特WT reporter tiles、12个参考RNA，而不是14个RNA；所有资格来自来源和输入匹配，未按敲低效应挑选。实验读出为MCF7中蛋白敲低后同tile reporter的核/胞质log2比变化；CARBON输入为109nt insert中央101nt，用旧已冻结K562/HepG2结合量输出。它检验的是binding模型对定位依赖的迁移，不是直接定位模型。

原主终点D_H = mean(siNT)−mean(siHNRNPK)。RNA等权Spearman宏平均：HepG2 Full −0.1462、Sequence −0.0595；K562 −0.1043、−0.0834。既没有明显正向功能识别，也没有稳定保守性优势。实测HNRNPK RIP与主D_H宏相关也仅0.0359；与原未敲低N/C的相关0.3559，说明“会结合”和“敲低后定位依赖”本来就不是相同量。

作者Effect列逐tile、逐蛋白独立核验3507值，严格为median_j(KD_j−siNT_j)，最大绝对误差4.89e−15。它不是DESeq2、表达量变化或condition均值差。反号后独立保留为来源定义敏感性：HepG2 Full −0.1806、Sequence −0.1223；K562 −0.1054、−0.1272。结果不因选作者定义而转为有力阳性。所有12 RNA、三个蛋白公式/符号均在EFFECT_AXIS_ALL_RNA_FORMULAS.csv；不能只报告冲突大或表现好的RNA。

按文献RNP共同功能关系，仅做一个事后探索：同RNA内HNRNPK与SLTM预测秩的min。D_H的Full宏相关HepG2 −0.0921、K562 0.0089，Sequence 0.0223、0.0003；仍未建立功能解释优势。不再增加这一库的评分公式。

JPX/PVT1有一手来源支持的人源身份，其余来源物种与human reference兼容性分开。1169tiles有重叠，同一研究的12 RNA是汇总/不确定性单位；训练种子、3重复、多个tiles不能扩增为独立RNA。旧flat/permuted仅输入干预，未当成独立trained ablation。Full与Sequence都有GRU，差异不能归因GRU。此项不支持新主文贡献，保留为模型功能迁移边界。

可重放：run_fr1.py、check_fr1_effect_axis.py、finalize_fr1.py。主结果/资格/不确定性在fr1完整CSV，全部来源与历史预测哈希已存；本波新模型查询0，首次CPU重放56.6秒。

来源：[Lubelsky et al., EMBO Journal 2021](https://link.springer.com/article/10.15252/embj.2020106357)，Dataset EV1/NucLibB。来源作者湿实验，不计为本项目新做湿实验。
