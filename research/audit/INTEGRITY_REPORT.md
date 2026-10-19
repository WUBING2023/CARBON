# 来源与一致性核验报告

范围：2026-10-17本地/已授权server15文件盘点、原始小型记录复制、权重字节身份和关键数值复算。不是独立同行评审，也不是科研不端法证认证。没有修改任何原始指标或挑选新的阳性子集。

## 本轮实际通过的检查

| 检查 | 实际结果 | 证据 |
|---|---|---|
| 服务器材料复制 | 24,010文件，零传输哈希差异 | REMOTE_COPY_QA.json |
| 本地代码历史 | 8,185来源条目，3,867个不同内容；0静态凭据格式替换 | LOCAL_CODE_ARCHIVE_QA.json |
| 权重身份 | 2,281文件、3,572,595,106字节；读取前后稳定 | CHECKPOINT_HASH_QA.json |
| 核心代码语法 | 当前CARBON和历史MARCP共52个Python文件可解析 | CORE_SYNTAX_CHECK.csv |
| seed0/1/2旧源码 | 依据历史manifest恢复27/27个哈希完全一致的文件 | S21_SOURCE_RESTORATION.csv |
| 消融 | 1,110配对行、222任务、13,320个原指标核对；最大绝对误差9.71e-17 | NUMERICAL_RECHECK.json |
| 原生分配 | 重算1,851对target/share/loss，659gene、865cell-gene再按细胞等权；MSE与冻结表一致 | NUMERICAL_RECHECK.json |
| 外部偏好 | 180个task/arm/assay组合，1,024词Spearman及15,240同组成词对；最大指标误差0 | NUMERICAL_RECHECK.json |
| 疾病子集 | 全13事件/8变异，Full7/13、Sequence10/13不利结果保留 | NUMERICAL_RECHECK.json |
| ClinVar完整快照 | 找回200,437,060字节原VCF，SHA与保存manifest一致 | E10的SOURCE_MAP与ClinVar manifest |
| 当前图表 | 1–5都有编辑源/入口；2–5源PDF与当前PDF哈希完全一致，图1保留PPTX及原导出 | FIGURE_CODE_MAP.csv |

## 可从证据得到的结果

去除保守性平均Spearman变化为-0.023449478910，219/222任务下降，支持同模型预测价值。原生分配Full MSE=0.079884017654，Local+pos=0.112084039317，误差低约28.73%。外部全词Fullρ=0.257443、Sequenceρ=0.269308、frequencyρ=0.111256；不把Sequence更好的结果隐去。

以上是本轮独立编写的下游算术复算；没有重新从BAM处理信号、没有重训模型、没有重新从7mer原作者表生成全部5mer标签，也没有重跑置信区间。既有独立审查按其原范围保留，不能充当本轮全部工程的独立认证。

## 仍不能由本轮检查证明

不能保证所有原始测量绝无造假；不能宣称三模型已公平统一评价；不能用文件存在证明训练从未看过标签；不能把模型输入依赖变成保守性功能解释优势。原始测序/参考数据仅做位置和体积盘点，没有对全部564GB文件做内容哈希。详见KNOWN_GAPS.md。
