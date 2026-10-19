# FR4真实可执行性审计

已成功获取作者processed rMATS与原metadata、hg19→hg38 UCSC chain227698bytes，不需下载FASTQ或重新处理RNAseq。K562公开SE表40998行、39053有限PSI；chr1 4496有限、50响应/41genes。HepG2chr1响应0完整保留。两个原实验的target和si/shRNA身份为PTBP1，2+2重复，文件明确KD−control符号。

可用matched两PTBP1任务四臂×三seed，CPU4无冲突；标准101nt上下游intron窗口预计<=108k推理<2min，坐标/测量资格约5min。hg19源测量不直接与hg38输入硬接；全部六exon边界和101nt逐碱基唯一UCSC链映射、strand与intron几何是先过的门。主候选为K562chr1，严格holdout；其有用功能样本最多50事件/41gene，覆盖和intron长度过滤后按实际数报告。

当前共享h5有实测Y/Yctl，但不是完整任意坐标BW入口；不在没有匹配实测binding时称bound nonresponsive。原生binding模型不是剪接预测器，仅检验其位置信号与实际功能方向能否一致。主响应方向与所有事件/简单baseline固定在PROTOCOL_FR4。原三NC阴性仍成立。

来源：[ENCFF061DFK K562](https://www.encodeproject.org/files/ENCFF061DFK/)；[ENCFF044JDS HepG2](https://www.encodeproject.org/files/ENCFF044JDS/)；[正式原文](https://www.nature.com/articles/s41586-020-2077-3)。所有抓取/log/scope/hash在ptbp_splicing_feasibility。
