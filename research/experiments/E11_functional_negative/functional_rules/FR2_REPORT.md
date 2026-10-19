# FR2：新的31RNA功能库未显示核定位解释优势

真实执行的是Lubelsky与Ulitsky Nature2018的human NucLibA设计与实测N/C reporter读出，排除之前已用JPX/PVT1/MALAT1后，固定31个新参考RNA、3113个合格109nt tile。不是把同RNA增加tile当独立RNA扩大；设计物来自34个human reference group，严格唯一gencode/hg38映射、完整nativePhi的全来源3905tile先固定。原文设计5511/qualified5153叙述，与实际用于这次分析的STable4 4884行不混。

两HNRNPK细胞的四trained arms×三seed真实93,720次native count forward，CPU4 47.254s；中央101nt裁剪固定4。主实测功能读出是3′UTR放置N/C log2，各RNA内Spearman再RNA等权宏平均。HepG2 Fullρ−0.1457，K562−0.0990；Full−S分别+0.00115与+0.01646，必要三对照同时区间全部跨0。5′放置、WCE/Input表达、组成/表达控制和全部seed结果已保留，不换它们当主终点。

这说明当前HN binding输出不能普遍预测这组不同机制RNA的局部核保留。它没有证明模型不能识别HN结合；核定位还需要合作蛋白、RNA背景及细胞过程，且HN输出不是直接训练的定位概率。该方向不宜作为模型功能优势节。

实测表、严格映射及结果均在fr2；PROTOCOL_FR2.md和可重放脚本保留。没有新figure或主文改动。
