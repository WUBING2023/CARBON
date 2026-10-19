# FR2：未用RNA中核保留序列区域的跨RNA功能检验

登记2026-10-04（实际文件冻结时刻见fr2/SOURCE_FREEZE.json与INPUT_FREEZE.json），在新CARBON功能预测结果产生前。FR1阴性已知。优先改换真实大库，不用相同失败库重复评分。

来源：Lubelsky & Ulitsky, Nature 2018, DOI10.1038/nature25757，原Supplementary ZIP STable1/2/4和真实全文XML。NucLibA将109nt自然RNA tiles分别放入AcGFP的3′或5′UTR，在MCF7实测三个重复的核/胞质log2比。它是功能定位结果，不是bound/unbound、不直接提供HNRNPK特异KD标签。论文整体叙述5511设计/5153阈值片段；实际STable2/4可用行数按表逐行记录，不强制重写为叙述数。

资格：STable1中RefSeq或hg19自定义human lncRNA设计（前34组），准确ID/109nt序列连接STable4；原STable4已按作者覆盖过滤。本轮保留全部来源资格与排除原因；不按高N/C、SIRLOIN或模型表现筛选。用GencodeV29人源RNA拼接序列精确109nt匹配；仅唯一coordinate map且完整native PhyloP的输入进入Full必要比较。EVF2用已知别名DLX6-AS1，所有映射原名/别名记录；不得将鼠/鱼/蜥蜴MALAT1或鼠3′UTR强映成人源Phi。输出来源34组与reference兼容性单列。

主新RNA范围：排除此前深入探索过JPX、PVT1、MALAT1三个RNA，保留其它全部合格RNA；另全部human设计范围作为明确次级敏感性。未用RNA指本项目功能端未用，不称独立研究或完全盲测：论文机制与STable标签已有公开。必须保存各RNA所有tile，不因结果丢RNA；至少5唯一tile且非恒定才计算相关。重叠tile不当独立RNA。

预测：HNRNPK K562/HepG2四matched训练臂aligned、sequence_only、mean、permuted，各seed11/23/47；101nt固定中央crop4:105，无训练/阈值调优。直接用原生log_total_count（训练目标log1p观测总信号），ensemble为三seed的log输出均值，不把它称定位概率。固定排列用matched_common.permutations(task+'|FR2_NucLibA',0,n)，输入顺序在INPUT_ROWS冻结，各seed一致。三个seed算各自与ensemble全部报告，不选最高seed/细胞。所有model/source sha保存且strict load。

主终点：每RNA内预测结合量对mean(3′UTR N/C log2 replicates)的Spearman，RNA等权宏平均。在给定研究中检验能否识别自然RNA内真正能驱动核保留的区域；不把相关称因果定位机制。

必要对照：Sequence/Mean/Perm训练臂；C-content、CCC、CCCC、CCTCCC、GCCTCCC次数与meanPhi。对相同定义RNA配对比较Full−各对照；简单motif恒定产生undefined时清楚列分母，不能跨分母比较均值。RNA簇bootstrap10000次、seed20261004、Full必要三比较Bonf3同时区间；只外推当前RNA库，不假定跨研究稳健。

固定次级：5′UTR定位同样全部读出作为同109nt序列的placement复核；3′/5′RNA内观测相关可读测量可迁移性，但两位置不是独立RNA。3′UTR WCE/Input表达量与3′NC分离，不能把表达减少当核保留。控制C-content/CCC以及表达量rank的partial-rank仅诊断；如新RNA内主结果失败，不换为事后挑出的单RNA。

即使Full/S等模型家族均较简单基线好，最多说明模型家族捕捉额外序列上下文；只有Full对必要Phi controls同方向且不确定性支持，才考虑保守性额外功能价值。Full/S都有GRU，不能特指GRU贡献。
