# FR4固定次级范围与实现门

2026-10-04 12:10 UTC，模型尚未连接事件标签。主Agent允许固定扩展所有其他训练省略染色体chr8/15/X/Y；故本次输入固定chr1/8/15/X/Y全集，chr1仍为主，其他四染色体合并为明确次级确认范围，各细胞完整报告，不挑染色体。总截止仍12:43UTC。原始标签来源相同，不能当第二项独立研究或增加湿实验重复数。

采用gene簇bootstrap，点估计仍event-weighted。方向二分类仅主ρ的辅助诊断，记录恒多数baseline准确率、balanced accuracy和MCC，不能替代主终点。全部范围使用相同coverage、唯一映射和纯intronic101nt资格。

响应01_context/CROSS_AUDIT_03：half-open exon端点映射start与end−1并恢复exclusive end；窗口全部101碱基唯一、连续、同向映射；−strand序列和Phi同时反向；permutation用外部原始row+side身份的固定SHA seed而非冒用HDF row；主分数是三seed平均native log_total_count，上下游相减，非binding/splicing probability。缺少PTBP1任意坐标BigWig，不称任何事件已结合/未结合。是否可由HDF真实观测覆盖，另标支持范围，不改变事件资格。
