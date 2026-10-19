# FR4：真实PTBP1敲低剪接检验的完整结果

这是在前三项核定位候选之后明示批准的第四方向，实测功能终点与模型训练的eCLIP读出不同。原作者ENCODE batch-corrected rMATS，Sample1是KD/KO，Sample2是matched non-target control；已逐文件和实验metadata核target/PTBP1-LV11及K562对照ENCSR129RWD。ΔPSI是真KD−control inclusion变化，不是表达量或结合强度。标签看过后登记操作阈值，因此是探索而非完全盲测。

主范围chr1、固定次范围chr8/15/X/Y，均为matched训练省略染色体，追加全集在连接模型分数前冻结。源hg19 six exon half-open边界与全部101nt intronic flank唯一连续映射至hg38；负链RNA顺序下序列和Phi同时reverse。每RNA-seq条件每重复IJC+SJC≥10；纯intronic长度≥101；完整序列/Phi。完整资格来源、排除原因和原行ID均保存。未找到任意PTBP1 eCLIP BW，不能把这些事件称已结合/未结合。

主问题是预先固定的位置规则：上游−下游binding logcount能否对应signed KD ΔPSI。在K562 chr1最终43响应事件/36gene（9 inclusion、34 exclusion），Fullρ−0.0640，gene-cluster95%CI[−0.4053,0.2738]。S−0.0049，Mean+0.0717，Perm−0.00385。Full与三个必要对照差的Bonf3区间均跨0。NoGRU追加训练组件诊断ρ−0.0842，Full−NoGRU+0.0202[−0.2062,0.2629]，未显示功能机制优势。次级其他holdouts27响应/21gene（5/22），Fullρ−0.0608，S+0.3117但区间跨0。HepG2 chr1响应0，其他范围仅2，不构成可靠方向复核。

类别不均衡明确呈现：K562主恒预测多数方向准确率34/43=79.07%，Full实际55.81%、balanced accuracy51.63%、MCC0.0269。裸准确率没有被用来美化模型。事件点相关是event-weighted，整gene重采样给依赖区间，不误称点估计gene等权。

预先登记的另一个次级问题是**检测哪些事件响应**，使用两侧平均三seed native logcount。主1989事件/608gene（43 responsive、1946 operational nonresponsive）FullAUROC0.5725[0.4706,0.6676]，S0.4882；差+0.0843，Bonf3[0.0148,0.1641]。vsPerm差+0.0564[0.0042,0.1063]是局部差异，必须保留；vsMean差+0.0295同时区间跨0，meanPhi简单基线AUROC0.5851又高于Full。其他holdouts1543事件/542gene FullAUROC0.6063[0.4572,0.7820]，必要差均跨0；HepG2次级仅2正事件，FullAUROC0.2274。它没有确认泛化功能识别，不能替代主位置规则阴性。nonresponsive只是有覆盖/当前阈值未显著变化，不证明严格无功能。

全部4matched arms×3seed native forward233,712次、CPU4 121.70s；NoGRU6fit额外58,428次/9.66s。NoGRU与Full同训练预算/划分但容量44,698 vs100,826，有不同实际early stopping，不解释为纯GRU机制。坐标抽取与模型推理使用各自已有系统库/原torch venv；一次pysam依赖修复发生在模型查询之前，无源/终点变更。简单CT/C比例用准确整数count差除101，避免浮点微差破坏应有并列；模型主数值未变。

资格、全部实测标签、模型输出、必要对照与CI、简单基线、每seed、PSI/coverage固定五分位及partial-rank均在fr4目录。没有改offset/head/响应阈值、没有挑阳性染色体。结论：当前可保留一条局部响应检测差异作为探索边界，**不足以成为功能解释优势节，也不支持保守性独有机制**。
