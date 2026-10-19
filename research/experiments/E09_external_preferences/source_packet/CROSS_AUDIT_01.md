# 独立复核 01_context：真实复算与机制边界

2026-10-04。只读 A 的数据、代码、原生输出及六个新 NoGRU 权重；审查代码与结果仅写在 02。

## 已实际独立复算

从 chr1/X/Y 原生 101 nt 的 `y/control` 与 Full/Sequence/Mean/Perm/NoGRU 的各 seed 原始 profile 数组重新求 softmax 概率总和、两个位点各 11 nt 质量及其份额，而不是读取 A 的最终平均值充当独立验证。完整严格 input-adjusted 集合为 1,851 个窗口、659 个真实 gene ID。先在 task×gene 内平均，再让两细胞等权；没有把窗口或 seed 当生物重复。

独立复算 20 个模型读出 × MSE/排序共 40 个点值，与已存表最大差 4.78×10⁻¹³；用已存采样状态重新生成全部 10,000 次相同 gene ID 的配对重采样，Full−NoGRU bootstrap 最大绝对差 4.16×10⁻¹⁷，95% 区间 [−0.0420575,−0.0285768]。三 seed 和两细胞全保留，Full MSE 优势同向；完整文件 `results/CROSS_AUDIT_01_*`。

原生 profile 是位置概率，两位点 share = mass₁/(mass₁+mass₂)，不会再乘 count；三 seed 先在 float64 平均 profile，再求 ensemble share。实测残差 = signalᵢ − total_signal × inputᵢ/total_input；仅保留残差均非负且和正的预先规定集合，不把负残差截到零。排序仅在实测非平局的行上评，预测平局计 0.5。固定同一坐标、同一位置支持范围；NoGRU 未换 RNA/标签或截窗。

六个本地新权重 SHA256 与各自 complete.json 一致，训练日志 <=100 epochs，best checkpoint 确为该 run 最高验证 PCC；均使用同一 split manifest SHA256。独立 CPU2 从 fresh matched Full 重新初始化每一seed并逐 shared tensor复制，六 run 的 NoGRU初始state hash 与已训练记录完全一致，Full 初始state也与原Full run完全一致。两架构共同训练到的全部 epoch 的样本次序SHA逐一相同（各run共同epoch数60/49/50/38/48/44）。NoGRU实际epoch92/81/91/100/77/44，Full60/49/50/38/48/50；因此准确表述是“相同训练协议与最大预算”，不是实际相同epoch数或计算量。全部真实记录见 `results/CROSS_AUDIT_01_INITIALIZATION.json`。

代码沿用 matched trainer 的 AdamW、损失、batch、order seed、patience、调度与验证选择，仅替换现有 noGRU 架构和共享初始化，末尾历史 heldout evaluation 延后。模型输入继续使用 aligned Phi。NoGRU 44,698 参数，Full 100,826 参数，SE 全局门和 residual/head 保留；所以主结果验证“GRU 路径及其容量有助于当前预测”，不单独证明递归机制。

## 最关键的中心/边界先验混杂

`remote_context.py` 从旧 `_5ch` 和原 raw HDF 读入 101 nt；没有重新建立无偏窗口。现有代码只核 DNA 与坐标一致，原 `peak` 列仍来自旧数据，不能仅因长度101或中心50就断言这是 binding summit 的精确中心。

实际SSH读到原生产源码：`/data/CARBON_PROJECT/数据处理/achieve/preprocessing/RBP_peak_pad.py` 92–195实现了Random与Half两种补边；Random按允许左右transcript范围随机leftPadding，Half才按长度的一半；长peak补到windowSize的整数倍后分窗。`/data/CARBON_PROJECT/数据处理/achieve/ERBP_eCLIP_dataProcess_pipeline_runtype.0_100.py` 86附近默认 `window_Size:101`、`padding_Mode:Random`。这能证明该生产链存在非固定中心方式，但目前没有完整HDF配置链，不能断言当前所有窗口实际使用Random或据此消除位置偏倚；也不能把name_peak当精确summit。第一次尝试的BED2H5路径未找到，未用它补充未验证结论。

Full GRU 的截窗边界可提供绝对位置先验；局部卷积+共享 SE 的 NoGRU 没有同等坐标信息。因此 native Full−NoGRU 效应不能直接称为远端生物序列机制，增强位置/容量对照仍必要。

本轮同严格集合最粗 window-center 排序规则为 52.41%，input-mass 51.56%，低于 Full，但这只排除简单二元规则，不能排除可学习的位置校准。A 已另行冻结并计算 Gaussian priors/相近中心距离诊断。

2026-10-04 13:06 UTC按root实际协调更新：原“仅+101 learnable position bias”的六run方案**在新arm任何training前被替代，未启动，不计作已训练控制**。更强新arm `NoGRU_position_capacity`：局部h128连接固定(t_norm,t_norm²)，经逐位MLP130→216→128 GELU residual（末层初始化0），接原LN/dropout/heads，额外101个profile位置bias初始化0；100,871参数，对Full100,826差45（0.045%），保留RF7与全局共享SE。它用相同freshshared初始化、split与最大100epoch协议，同时针对灵活绝对位置与容量差混杂。由A执行完整六PTBP1 fits，本目录不训练。六fit现已真实完成，实际训练1639.17秒，原生CPU推理5.10秒。独立CPU2重建六run：shared tensors和fresh初值hash完全一致；末层0和101bias0使初始profile/count与原NoGRU误差均0；所有MLP卷积核1、无GRU/LSTM、明确t和t²坐标。实际epoch54/74/62/64/73/94，best45/65/53/55/64/85；同split与共同epoch样本顺序全部吻合。参数差45已核，不再把旧44,698参数差当新控制的混杂。

逐实现核查并数值干预bias：父forward仅返回原始logits/count、没有缓存概率；新bias确加到返回logits，真实训练criterion、validation softmax、inference softmax均消费更新后的logits。人为bias修改造成非零概率改变和非零训练loss梯度，已训练bias也非零。这些是数值实现QA，不是生物功能实验。

新控制原生数组与旧四臂/原NoGRU严格同坐标。重新复算24模型读出×MSE/排序共48点，与A表最大误差4.78×10⁻¹³；按已存状态重新生成全部10,000次paired gene bootstrap，最大差4.16×10⁻¹⁷。Full MSE .079884，近容量位置NoGRU .112084，下降28.73%；Full排序76.74%，新控制66.31%，两cell/三seedMSE和排序都同向。Full−新控制MSE区间[−.0387045,−.0257192]。完整输出为CROSS_AUDIT_01_PC_*及POSITION_CAPACITY_MODEL.json。

该附加控制表明优势不能仅用明显参数数量不足或缺少灵活绝对坐标来解释；没有证明不同架构优化难度完全相同，也没有直接测归因正确性、远端生物机制或功能变化。可准确称GRU路径有助于当前定位判断，不能直接称“已证明远端机制解释”。

## 额外外部诊断

使用相同固定 16 背景、全1024词、两 PTBP1 RNAcompete 列、三seed/两背景的新 NoGRU 体外规则检验已独立完成：Full 全词 ρ=.65280，NoGRU=.68497；同组成顺序 Full=.59762，NoGRU=.57418，其中 HepG2 接近持平、K562 有差。Sequence=.62156 仍高于 Full。这些完整阴性/局部阳性保留，不能把平均2.34pp变成 native 分配优势的普遍分子机制证明。

## 第二蛋白 HNRNPK 的实际独立确认

同样从chr1/X/Y原生profile及y/input读回，保留1418严格pair/527gene；20模型读出×MSE/排序共40个点与原表差<4.86×10⁻¹³，三seed/两cell全部保留。Full/NoGRU MSE为.080893/.106821；从已保存的每次arm means独立计算Full−NG的95%差区间[−.0317861,−.0199898]，与savedCI严格一致。没有为了重复审核而再次生成10,000个相同draw；这项范围审核和PTBP主bootstrap实际重生成是不同check。

全部六HNR权重fresh shared初值和原Full init hash、split manifest、共同epoch sampleorder均独立CPU2复现；NoGRU实际epoch84/76/88/100/88/100，Full37/45/37/49/43/45，同协议/max budget而非相同实际epoch。外部固定word诊断：RNAcompete Fullρ.648201<NG.668707，同组成.619269<.623272；RBNS原primary20nM Fullρ.618600>NG.545019，但同组成.476879/.475786均低于.5。所有RBNS五浓度与RNAcompete分别保留，不能挑RBNS相关阳性替代RNAcompete陰性，也不能因此宣称native分配优势已获得一般体外机制验证。

结论：两蛋白主算式、全部必要训练臂、权重身份、输入坐标、单位、gene 等权和所述成对 bootstrap 均通过独立实际数值复核。科学层面支持当前架构路径的定位判断贡献；PTBP1附加近容量、灵活位置控制仍未追平Full，但不把这一能力测试改写为归因ground-truth或生物功能机制的直接验证。HNRNPK仅完成原NoGRU跨蛋白确认，未把PTBP附加控制冒充HNR也已完成。

最终contrast闭合：按实际`comparison`字段区分Full−control与position_capacity−NoGRU，全部36个Bonf1/2/3/15/45/90区间从已独立重建的10k draws再次核对。Bonferroni两主比较区间Full−Local+pos[−.0395639,−.0247504]，数值完全一致；首次脚本使用了消息内简写而触发assert，修为实际CSV字段后通过，未混淆两种差值。

最终已训练支路身份读回：本地CPU2直接读取六best权重，private五tensor×三norm共90项与A学习QA精确0差；初零MLP末层和profilebias六run均已变成非零。没有再forward或训练。控制源码hash905f265e…与初始审核相同，未漂移。
