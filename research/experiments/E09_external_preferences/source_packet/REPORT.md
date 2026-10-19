# 外部序列偏好：实际结果与能写到什么程度

这轮得到一项可写成独立实验的结果：**CARBON 模型家族在25个此前未探索的真实蛋白上，学到了超过简单词频、不能只用字母组成解释的序列排列规则。** 同时保留一个 HNRNPK RNAcompete 的局部优先级改善。两者都不能改写成“保守性普遍提高功能解释”。没有补造湿实验、没有新训练本目录模型、没有继续挑词长/输出头追逐阳性。

## 新25蛋白的完整确认范围

在作者公开的原始 RNAcompete 420实验集合与我们222任务权重库存之间按真实human protein identity交集，冻结28蛋白/42任务/84原模型权重；其中已查看的HNRNPK、PTBP1、RBFOX2列作曝光探索，新增25蛋白单独统计。每任务固定16个原生背景；中心5位置逐一替换全1024种词，同一背景的WT/全部ALT沿用同一Phi向量。评价外部原7mer Z分数的固定5mer折算，不是直接5mer Kd。所有原实验列、cells均保留，不按外部表现选列或背景。

旧 Full/Sequence 的初始化和训练配对无法证明；因此本范围是家族规则验证，不能把差值当Phi因果增益。frozenFlat/frozenPerm只是完整模型固定权重的输入变换，与前三蛋白真正重训的Mean/Perm不是同一对照。

|新25蛋白等权|Full|Sequence|Frequency|FrozenFlat|FrozenPerm|
|---|---:|---:|---:|---:|---:|
|全1024词 Spearman ρ|0.2574|0.2693|0.1113|0.2335|0.2705|
|相同字母组成的词对排序|56.62%|57.31%|52.84%|56.23%|56.36%|

同组成排序中 Full/S 各有23/25蛋白的点值高于50%，这是点估计计数，不是23个各自显著的检验。四个通用 A/C/G/U 数量规则全部保留；没有给每蛋白事后选最有利的字母。固定Phi-only在当前替换设计下对全部词恒定，ρ不能计算、同组成排序为0.5。

|冻结比较（蛋白配对bootstrap10,000）|均值差|95%区间|正向蛋白|
|---|---:|---:|---:|
|Full−Frequency，全词ρ|0.1462|[0.0380,0.2480]|18/25|
|Sequence−Frequency，全词ρ|0.1581|[0.0733,0.2371]|20/25|
|Full−Frequency，同组成排序|3.78pp|[0.90,6.73]pp|19/25|
|Sequence−Frequency，同组成排序|4.47pp|[1.28,7.51]pp|20/25|
|Full−Sequence，全词ρ|−0.0119|[−0.0665,0.0387]|13/25|
|Full−Sequence，同组成排序|−0.69pp|[−1.91,0.51]pp|9/25|

独立单位是可用的真实蛋白，不是成千上万的词对。区间描述当前可用蛋白范围的稳定性，不能宣称随机抽样全体RBP总体。28全范围、曝光3和新25分开留表；所有单蛋白、单背景与原外部实验列完整保留。

## 三蛋白匹配模型与局部 HNRNPK

前三蛋白使用真正matched六任务×三seed×Full/S/Mean/Perm共72模型，原生101nt输入、同窗口hash抽样及原训练Phi变换完全一致。32步DNA-channel IG沿统一0.25背景积分，Phi固定，count输出头始终未换；每窗取最高5mer IG支持窗口，聚合固定词优先级。主要IG实验在RBNS的HNRNPK/RBFOX与PTBP1 RNAcompete下没有普遍的Full优势，完整阴性见 `results/ALL_METRICS.csv`。受控1024词×16背景分析也不是全阳性：PTBP1全词ρ Full .6528>S .5650，却同组成排序Full .5976<S .6216；不能拿全词相关提升冒充识别词内部排列。

新取得HNRNPK RNAcompete原始标签后，以已冻结IG输出复核，Full top20 recall29.17%，S21.67%，训练Mean25.83%，训练Perm24.17%。2048原生窗口恢复原HDF坐标并逐序列核对，用GENCODE v29 gene body与overlap合并为385簇，配对抽样共同作用于两cell/所有seed；2000次条件bootstrap Full−S +7.50pp，范围[1.67,9.17]pp；Full−Mean +3.33pp，[−4.17,3.33]；Full−Perm +5.00pp，[−3.33,4.17]。这是同一蛋白当前输入范围的条件稳定性，不是新的蛋白复制。

Top-k是非平滑量，原始Full−Perm点5.00pp在中央95% draw范围上界4.17外；逐fit全1cluster权重严格复现原点值，独立reviewer重做31个draw也一致，未改CI类型使其看起来阳性。真实优势来自少数C-rich词的优先级，其中CCCAU/UUCCC等<=3C词有贡献，>=4C词对Full/S recall贡献相同；完整20真值词包括全0和反向差都留在图c/CSV。

简单C-count whole-wordρ为.810，高于Full IG .402；不能说整套规则超过生化共识。C-count lexicographic top20 recall25%，边界tie公平期望12.89%；Full/S没有该边界tie，数字不变。原lex端点和公平tie敏感性同时保留，未以后者替换原端点。原词频recall10%。RBNS外部真值下同一HNRNPK Full−S −4.17pp，[−5.00,0.83]；不隐去跨assay不一致，不凭现有数据把差异归因蛋白构建。

## GRU 外部诊断与上下文审核

六个新PTBP1 NoGRU权重仅CPU8评现成固定16背景/1024词/两RNAcompete列。全词ρFull .6528、NoGRU .6850，Full更低；同组成Full .5976、NoGRU .5742，差+2.34pp主要来自K562，HepG2接近0，Sequence .6216仍更好。新增第二蛋白HNRNPK原同16背景六NoGRU权重完整诊断：RNAcompete Fullρ.648201<NG.668707，同组成.619269<.623272；原冻结RBNS primary20nM Fullρ.618600>NG.545019，同组成.476879/.475786均低于.5。全部五浓度、两细胞、三seed留表，未平均两种assay原始R与Z或选有利列。体外词偏好与原生两个同核心位点的profile分配是不同任务；不能把一项局部顺序差或一个assay相关阳性当整项上下文优势的机制验证。

对A的PTBP1 native主包做独立数值复算：原1851严格窗口/659gene，40点误差<4.8e−13，10,000 paired gene bootstrap误差<4.2e−17。重新fresh初始化NoGRU并逐shared tensor复制，六run hash与原记录完全一致，共同epoch训练样本顺序hash逐一相同。Full/NoGRU使用相同最大100epoch协议，但实际早停epoch不同。完整审核 `CROSS_AUDIT_01.md` 保留44,698对100,826参数容量、101nt截窗边界/中心先验限制。A原仅101位置bias方案在training前被更强 `NoGRU_position_capacity` 替代，旧方案未启动。新局部位置MLP130→216→128及101bias共100,871参数，距Full45参数，并保持RF7/原SE，六fit已真实完成，并通过独立fresh-init/训练顺序/实现审核：loss、validation、inference确实消费加入bias后的logits。新控制原生48点最大误差4.78e−13、10k配对gene bootstrap差最大4.16e−17；Full .079884对新控制 .112084，MSE下降28.73%，差区间[−.0387045,−.0257192]，排序76.74%对66.31%。两cell/三seed一致。它支持GRU路径的定位判断贡献，削弱明显容量/缺绝对坐标的替代解释；不直接证明归因正确性或真实功能机制。

## 数值、数据和失败记录

72模型各16/32/64步检查已完成：32对64步IG相对L1平均0.221%，最差1.406%；32步主要定义不改。CPU32对原CUDA32最大相对L1约6.63e−7。每模型输入/变换、1024词轴、全部WT/ALT和checkpoint SHA都保留。

RNAcompete原矩阵实际16,382个7mer，缺GAAGAGC/GCTCTTC，未补0；所有1024五词仍各有47–48次真实occurrence。主折算occurrence权重已冻结；固定unique-containing敏感性使四已曝光RNAcompete列top20均20/20一致，完整排名ρ>=.999985。作者reported A/B Pearson只是metadata；没有A/B原始分列，不能假称重新验证作者重复实验。

初次wide终端最后一行误导致“168全部完成”的工作状态记录，真实逐manifest核查只有160条。PCBP1两任务合格source窗为632/916，小于原频率样本1024但足够已冻结16背景；此前assert failure未产生预测，也未看标签筛掉PCBP1。先写source-quantity addendum，再只补该8条件，既有160原结果不动，失败日志保留；最终42任务×4条件=168真实摘要/响应、checkpoint和input逐核通过。

原生Phi-word均值不能给未出现词填0；所有模型受控主端点仍全1024，nativePhi association辅助仅在每task共同可观测词域重算，missing词清单完整保留。自动审批曾拒绝潜在占GPU的QA动作；实际发现他人作业后改CPU8并完成，未争卡。NoGRU外部第一次import与审查输出路径错误均已记录修复，错误没有更改旧材料或产生伪结果。

## 可供论文采用的准确主线

本节可以回答“模型学到的序列规则是否在独立体外测量中成立，并超越简单频率/字母组成”。配图与英文结果已产出，原始参考RBPNet Fig4作同宽比较。HNRNPK局部只作为明确后验、assay-specific优先级例子，不能拿来替换新25蛋白范围的总体。论文可以说存在经外部测量支持的序列排列规则；目前不能说保守性普遍提高功能解释，也不能说所有inferential规则比简单canonical规则更好。


原始文本精度闭合：独立family审查发现62项微差只来自已曝光RBFOX2的默认解析人工tie，新增25主结果通过。全局改用round_trip字串忠实契约，保留修复前快照和原失败记录；不做新推理，详PARSER_ADDENDUM.md与最终PARSER_QA.json。
