# FR3：能对应部分结合关键位置，未定位核保留功能必要位置

同JPX全部50位置、150个单碱基替换，实际四个N/C重复测突变−WT功能变化；RIP对应结合变化独立保留。N/C位置敏感性的同库重复ρ约0.745–0.863，B/C共有25位置ρ0.9177，说明功能标签本身相当可靠。两库共有同一个JPX，nRNA=1，旧binding终点与IG已曝光；新功能关联是探索，不是盲测确认。

首轮保存IG源实际是旧算法32step启动后统一升级的**最终64step**ATTRIBUTIONS，原协议手写“已有32step”不准确；本报告纠正事实，没有修改任何旧文件/归因/排名。重要性按实际代码abs(sum四序列通道IG)，不是sum(abs每通道)。旧HNRNPK Full HG/K的功能敏感性ρ分别−0.2696/−0.0244，结合敏感性ρ+0.3724/+0.2566；核心覆盖基线功能ρ+0.1945。旧NoGRU/NoSE训练条件未严格核同预算，不借它们声称组件机制。

第二轮matched四臂×三seed×两细胞全部真实执行：3624个count forward、792个IG梯度点、24个baseline forward，CPU4 36.87s。WT和全部ALT共用同referencePhi，各模型WT/ALT共用固定Phi排列；不是ALT专属进化信息。平均三seed的native logcount得到ensemble；IG取平均signedIG后绝对值，与平均scalar的IG对应，全部单seed同时保留。

全50位置，matched Full WTIG的功能ρ为HG−0.1326/K−0.1016；同一归因对RIP结合敏感性ρ为HG+0.1808/K+0.4080。count突变敏感性对功能ρ为HG−0.1368/K−0.1189。三个必要trained controls并没有一致劣于Full。即使binding有对应，功能核心优势仍未建立；不能将+0.408的RIP关联改名定位功能正确性。

新的matchedIG固定32step；最大绝对completeness残差0.002124，最大相对0.1810，后者发生在S/K562/seed47很小的输出变化。这是应披露的数值精度限制，不把32step图冒称精确积分；本轮未根据功能结果挑步数、crop或评分。已有count与可靠实测功能阴性共同支持“未建立优势”，不能据此证明所有可能的解释方法永远无效。

fr3/MATCHED_METRICS_ALL.csv、MATCHED_POSITION_MAP_ALL.csv、MATCHED_CONTRASTS_ALL.csv及MATCHED_ANALYSIS_QA.json保存全部条件/分母、错位对照与限制。这个工作有价值的结论是**binding-only解释与binding+function解释需要不同验证**；不宜写功能优势主节。
