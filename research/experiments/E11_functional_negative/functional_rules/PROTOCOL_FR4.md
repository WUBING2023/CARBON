# FR4：PTBP1结合位置能否识别真实敲低剪接响应方向

2026-10-04 11:58 UTC登记；主Agent明确批准在三项NC候选之后追加这一项更匹配的功能路线。45分钟有界执行，截止本片第一次登记后不迟于12:43UTC。已查看源标签数量与作者作用规律，未连接新模型分数；这是操作定义冻结的探索分析，不装完全盲测。旧eCLIP序列可能已见，但本项目未用本组processed PTBP1 KD-splicing endpoint。

来源：[Van Nostrand et al., Nature2020](https://www.nature.com/articles/s41586-020-2077-3) 的公开ENCODE batch-corrected rMATS文件集ENCSR870OLK。K562 ENCFF061DFK来自PTBP1-LV11 ENCSR527IVX，HepG2 ENCFF044JDS来自PTBP1-BGHLV12 ENCSR064DXG。文件metadata明确Sample1=KD/KO、Sample2=non-target-control，包含2+2报告重复。label ΔPSI=KD−control，两重复平均inclusion水平差；非gene表达logFC、非binding强度。文件已真实下载并SHA保存，hg19源坐标。

来源范围：K562 40998 SE rows/39053两组PSI全有限，chr1 4496有限；FDR≤.05且|ΔPSI|≥.05为chr1 50事件/41gene（KD inclusion增加10，减少40）。HepG2 chr1该响应阈值0，明示不可做响应方向复核，不强行改阈值。HepG2全部可测事件作为零响应边界保留；全基因组来源范围不是测试集扩大。

坐标与资格：只用匹配训练完全holdout chr1；旧hg19所有SE及flanking exon六个边界，用UCSC hg19ToHg38 chain逐碱基唯一映射，要求同chr1、同转换方向，目标区间不跨gaps/不改变101nt长度；RNA strand随chain方向变化，所有gene/event身份保存。全部来源行保留资格/排除，不按模型分数、Phi或大效应筛。每条件每重复IJC+SJC至少10，四PSI有限；响应定义如上，非响应操作对照FDR≥.5且|ΔPSI|≤.02（它不证明严格无功能）。中间范围仍保留不强分组。去重复按gene+三exon坐标+strand；多转录本/重叠事件以gene簇汇总，事件不是独立RNA/实验重复。

固定模型输入：cassette exon两侧紧邻、纯intronic101nt，RNA方向上游在3′acceptor前101nt；下游在5′donor后101nt。要求各侧真实intronic长度≥101，避开跨邻接exon。+strand上游[start−101,start)、下游[end,end+101)；−strand对应[exonEnd,exonEnd+101)与[exonStart−101,exonStart)，reverse complement为RNA5′→3′。不因结果移动junction offset、窗口或输出head。native hg38Phi完整、onehot序列完整。**将binding输出与功能读出分别命名。**

主终点：在全部合格K562 responsive events，预测位置指数log_count(upstream)−log_count(downstream)与signed ΔPSI的Spearman。预期正向：上游PTBP1结合更偏向WT抑制exon，KD后inclusion增大；这是事先固定的位置响应规则，不假设所有PTBP1位点只会抑制。HepG2没有合格responsive时终点NA，不借全基因组替代。另全有限事件的关联保留作次级，不能因它更大而换主终点。

模型/必要对照：现成matched PTBP1_K562与HepG2，aligned/sequence_only/mean/permuted，各seed11/23/47；primaryensemble=三seedlog-count平均。三必要对比Full−S/Mean/Perm，同事件资格、同genes bootstrap10000（seed20261004），Bonf3同时CI；种子不作生物重复，cells不是湿实验重复。新可比NoGRU若在登记后可用，作为明确追加组件诊断，真实预算/源头齐才报告，不混旧Full权重。

固定简单基线：两窗口的CT比例、CTCT重叠次数、C比例、native meanPhi各自上游−下游，全部基线用相同事件读出。输入RNAseq coverage与control PSI/|PSI−.5|作为测量动态范围变量。实测eCLIP y/control可用时单独给同坐标实测位置指数；仅覆盖原始h5窗口范围不能当任意位点无结合。binding阳性操作条件不使用模型预测选样，不称未覆盖位点unbound。

次级：response vs nonresponse分类用两侧合计binding（mean log_count），事件AUROC；按gene簇CI、control PSI与RNAseq coverage五分位完整分层和控制rank。若真正同gene有两类事件，可给固定同gene对比；没有不造配对。无eCLIP覆盖/密度支持时只能称功能响应区域分类，不能声称bound/nonresponsive已匹配。不会在这些诊断上新增十种评分。

正反方向、所有有限资格和阴性结果均保留。任何成功最多是当前PTBP1、K562、chr1、该真实KD定义下的局部功能位置规则支持；非普遍Phi、非物理接触、非网络直接输出剪接概率。若主结果失败，完整保留边界，不再本轮换window/head/threshold。
