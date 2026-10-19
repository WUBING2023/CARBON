from pathlib import Path
import json,csv,datetime
R=Path(__file__).resolve().parent
(R/'REPORT.md').write_text('''目前得到一个可以保留为论文独立实验的较强结果：**recurrent 的位置耦合处理帮助模型把同一窗口内两处相同短词的真实信号分配得更准确。** 这不是已经证明保守性通道独有优势。

主问题固定为：原生 101 nt 窗口内有两处 PTBP1 CTCT，短词身份无法区分它们，但实测 eCLIP 支持不同；模型是否把两处含端点 11 位邻域的条件份额分得更接近实测？主群体严格扣除 size-matched input，负残差不裁零；短词不等于完整 PWM 亲和力，份额不是结合概率、物理接触、IG 归因或下游功能测量。

PTBP1 两细胞有 1,851 个任务内去重位点对，唯一映射到 659 个真实 GENCODE 基因。Full MSE 为 0.079884，原 NoGRU 为 0.115116。因原对照容量不同，我们新增必要控制 Local+pos：保留局部感受野、全局共享 SE、原 normalization/heads，加入逐位置 130→216→128 残差支路，输入原特征和 t,t²，并有 101 个可学习 profile bias。它允许局部序列×位置交互，但没有 recurrent 跨位置传播。参数 100,871 对 Full 100,826，只差 45（0.0446%）。六 fit 同匹配初始化、全部原训练/验证数据、100 最大 epoch、patience 9、验证选模型规则；实际停止 epoch 可不同。

**Full 相对 Local+pos 的 MSE 下降 28.7%：0.079884 对 0.112084。** Full−Local+pos = −0.032200，主两项比较 Bonferroni-2 配对 joint-gene 聚类区间 [−0.039564, −0.024750]。更高信号位置判断正确率为 76.74% 对 66.31%。两个细胞和三个 seed 均同向；seed 不是生物 RNA 重复。Local+pos 本身比原 NoGRU 改善 −0.003032，调整区间 [−0.005049, −0.001038]，新增控制不是未学习或失效支路。所有六 checkpoint 的原零初始化 MLP 末层和位置 bias 均非零。

支持与限制都保留：原始信号 5,067 对/988 基因的差值为 −0.037844；历史曝光排除 535 对/260 基因为 −0.027153，两者 Bonferroni-90 区间仍低于零。中心距离差≤5 nt 的事后诊断 195 对/150 基因为 −0.027823，Bonferroni-3 区间 [−0.047997, −0.007633]；等距 23 基因和差≤2 nt 的 80 基因均跨零。天然完整相同 17 nt+Phi 子集为零，没有缩词凑例。HNRNPK 原 NoGRU 六 fit 支持 1,418 对/527 基因：Full 0.080893 对 0.106821，下降 24.27%，差 −0.025929，Bonferroni-36 区间 [−0.035348, −0.016202]；正确率 73.95% 对 64.92%。HNRNPK 没有 near-capacity 位置控制，不能借用 PTBP1 的控制完成容量/位置隔离。

保守性必要对照仍没有确立优势。chr1 六任务总体没有稳健 Full−Sequence/Mean/Perm 改善；冻结规则的 chrX/Y 163 对/87 基因有约 8.2% 点趋势，但三必要同时区间均跨零；旧 chr8/15 支持同样未过全部必要对照。它们全部保留，不能改称保守性功能解释优势。

最终主图为 figure/CONTEXT_FINAL.png/pdf/svg（85×125 mm）。a 共享坐标连接实测/Input/Full/NoGRU/Local+pos/Phi，b 用实际数组复算同邻域份额，c 显示全部 1,553 个 task-gene 点真实全范围。LEPR/OCRL 是新控制结果前锁定的原决策转换组中位后验说明，未按新控制曲线换例；原 hash 选定 FOXD2-AS1/KRBOX4 反例保留。最终图注是 FINAL_CAPTION.md，旧 figure/CAPTION.md 只配保留的 precontrol 图。英文正文是 RESULTS.md，并同步到 nogru_matched/RESULTS.md；历史版保留为 RESULTS_PRECONTROL.md。

独立 B 从原数组复算 48 个点和各 10,000 个 bootstrap draws，通过初始化/训练次序/100,871 参数/位置交互/bias 梯度与推理审查，36 项最终区间读回一致；六 checkpoint 的 90 项 private tensor norm 与本地 QA 精确一致。C 的 scope 审查与 root 实际最终 PNG 第二轮视觉复查均通过。机器 QA 验证 PDF 实际 85×125 mm、全部点在视域内、所有图示份额误差≤10⁻¹⁰。

可以明确写：Full 在该 native allocation 问题上优于近同容量、显式建模位置的逐位局部控制。不能写已发现具体远端功能规则、保守性独有解释机制或生化/下游功能已验证。chr1/X/Y 未参与 matched 训练和验证，但研究域有旧分析曝光；历史曝光排除不等于整个研究完全盲测。

实际 batch wall：PTBP1 原 NoGRU 22.85 分钟，HNRNPK 原 NoGRU 21.15 分钟，PTBP1 Local+pos 27.32 分钟。GPU2 每 run 即时 UUID/显存 guard，13:34:10 UTC 已释放，无新训练计划。RBFOX2 K562 完整三 seed 最坏约 35 分钟，在最后批结束后无法在 14:00 前闭合，未启动；HepG2 最坏约 64 分钟也未启动。原因仅预算、不是效果；不称完整三蛋白确认。旧论文、图和权重均只读。本目录为可审查的实验/图文候选包，未替换正式论文。
''',encoding='utf-8')
(R/'REPLAY.md').write_text('''# Replay and source map

Final scientific outputs are RESULTS.md, FINAL_CAPTION.md and REPORT.md. Final image/source/QA are figure/CONTEXT_FINAL.{png,pdf,svg}, FINAL_FIGURE_SOURCE.json, FINAL_FIGURE_QA.json and FINAL_PLOTTED_GENES.csv. figure/CAPTION.md belongs only to CONTEXT_CANDIDATE, not final. Original fixed counterexamples remain in figure/FIXED_COUNTEREXAMPLES.*. The final figure was rendered by render_final.py without any new model query.

The new local-position control is implemented in nogru_position_capacity/control_model.py. PROTOCOL.md was frozen before its pilot or outcomes; remote/TRAINING_MANIFEST.json records all six complete runs. remote/INITIAL_QA.json verifies fresh shared initialization and zero initial output differences. CONTROL_LEARNING_QA.json verifies learned private tensors. remote/INFERENCE_MANIFEST.json maps each native profile and checkpoint identity. A legacy training field named arm=aligned denotes the five-channel input transformation, not the model class; actual control class, architecture, parameter count and checkpoint path are independently audited and preserved.

Native inputs, measured signal, input-control signal, coordinates, original source-row indices, all four trained Phi arms and all seeds are stored under remote/data and phase2_ambiguity/remote/data. Newly trained NoGRU outputs are in nogru_matched/remote/data and nogru_cross_protein/remote/data. Local-position outputs are in nogru_position_capacity/remote/data. Each directory contains source/checkpoint/init hashes. No locally derived sequence material was uploaded to third-party services.

For a reversible replay, copy this entire 01_context directory to a new writable directory before running scripts: paths are resolved relative to each script, and analysis overwrites its tables. With Python and numpy/pandas/matplotlib/torch/pypdf installed, use:

    python nogru_position_capacity/analyze_context.py
    python nogru_cross_protein/analyze_context.py
    python render_final.py
    python finalize_delivery.py

Use the recorded environments and numpy version for exact bootstrap agreement. The local-position analysis runs 10,000 paired joint-gene bootstrap replicates using PCG64 seed 20261008; profile ensembles are float64 means before shares. ANALYSIS_INTERVAL_NOTE.md defines the main two-contrast family and supporting Bonferroni 3/15/45/90 families. All contrasts, methods, cell/seed/range, raw, strict and historical-exposure rows are in full CSVs; bootstrap arrays are retained. Never conflate comparison=full_minus_control with comparison=position_capacity_minus_nogru when filtering control=nogru.

No retraining is needed to reproduce final evidence or image. Authorized server15 train/inference replay wrappers are retained only for provenance. They refer to the existing matched source, native HDFs and GENCODE mappings, use original original full budget, and may launch expensive jobs; do not run them merely to inspect this delivery. A new training replay must use a new remote directory and a passing instantaneous GPU UUID guard, without interrupting other work.

Earlier unsuccessful conservation hypotheses remain at the root CSVs, phase2_ambiguity and test_8_15_support. They are not substituted by the positive component comparison. WINDOW_PROVENANCE_AUDIT.md records why native windows cannot all be called fixed peak-centered windows; RF_AUDIT.md establishes single-position RF 7 nt and union support RF 17 nt, while shared global SE remains.

Independent cross-checks reside in sibling 02_external_preferences/results/CROSS_AUDIT_01_PC_* and POSITION_CAPACITY_MODEL.json, and 05_delivery scope review. They belong to other workers and are read-only here. Final SHA256 identity inventory is DELIVERY_MANIFEST.json; it excludes itself, transport tar archives and generated proof thumbnails.
''',encoding='utf-8')
attempts=[
 {'id':'C1','hypothesis':'Phi distinguishes repeated core context','status':'completed, necessary controls inconclusive','fits':0,'inference_wall_s':219.2,'query_windows':8040,'checkpoint_queries_per_window':12,'notes':'six tasks; all native Phi arms, three seeds; 112392 chr1 input rows audited'},
 {'id':'C2','hypothesis':'frozen Sequence ambiguity on chrX/Y','status':'completed, trend but simultaneous intervals cross zero','fits':0,'inference_wall_s':51.57,'query_windows':2694,'checkpoint_queries_per_window':12,'notes':'primary PTBP1; six-task native inference; three protein identities not six independent proteins'},
 {'id':'C2_support','hypothesis':'same frozen ambiguity rule in old chr8/15','status':'completed, exposed support only, necessary controls inconclusive','fits':0,'inference_wall_s':0,'query_windows':3133,'checkpoint_queries_per_window':0,'notes':'reused cached native profiles; mapping transport 6.31s; local statistics 81.68s'},
]
for folder,task in [('nogru_matched','PTBP1 original NoGRU'),('nogru_cross_protein','HNRNPK original NoGRU'),('nogru_position_capacity','PTBP1 near-capacity local-position control')]:
 m=json.loads((R/folder/'remote/TRAINING_MANIFEST.json').read_text()); im=json.loads((R/folder/'remote/INFERENCE_MANIFEST.json').read_text()); windows=sum(v['windows'] for v in im['tasks']);secs=sum(v.get('seconds',0) for v in m['runs'])
 attempts.append({'id':'C3_'+folder,'hypothesis':task,'status':'all six complete, independent audited','fits':len(m['runs']),'inference_wall_s':im.get('seconds',im.get('elapsed_s','see native manifest')),'query_windows':windows,'checkpoint_queries_per_window':3,'notes':f'six full runs summed wall {secs:.3f}s; all native chr1/X/Y outputs; same maximum/early-stop budget, realized epochs differ'})
attempts.extend([{'id':'C3_position_diagnostic','hypothesis':'fixed Gaussian15/25 and center-distance strata','status':'completed; diagnostics do not replace primary','fits':0,'inference_wall_s':0,'query_windows':0,'checkpoint_queries_per_window':0,'notes':'post hoc registered; both sigma values retained, no fitted share mapping'}, {'id':'C3_RBFOX2_K562','hypothesis':'cross-protein original NoGRU confirmation','status':'not launched, complete three fits could not close by 14:00 budget gate','fits':0,'inference_wall_s':0,'query_windows':0,'checkpoint_queries_per_window':0,'notes':'training-only one-epoch cost pilot 6.688s; worst three complete 100epoch fits about 35min'}, {'id':'C3_RBFOX2_HepG2','hypothesis':'cross-protein original NoGRU confirmation','status':'not launched, exceeds remaining fixed budget','fits':0,'inference_wall_s':0,'query_windows':0,'checkpoint_queries_per_window':0,'notes':'training-only one-epoch cost pilot 12.771s; worst three complete fits about 64min; no heldout result queried'}, {'id':'C3_pure_position_bias','hypothesis':'position-only necessary control proposal','status':'superseded before execution by near-capacity local-position control','fits':0,'inference_wall_s':0,'query_windows':0,'checkpoint_queries_per_window':0,'notes':'no outcome chose final control structure; never run as a separate arm'}])
with (R/'ATTEMPTS.csv').open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=attempts[0].keys());w.writeheader();w.writerows(attempts)
status=(R/'STATUS.md').read_text(encoding='utf-8')
(R/'STATUS.md').write_text('''# Final 01_context status

2026-10-04 14:10 UTC. All authorized batches and native inference in this slice are complete; no further GPU experiments are planned. Six PTBP1 NoGRU, six HNRNPK NoGRU and six PTBP1 near-capacity local-position fits completed with the full matched budget. GPU2 released at 13:34:10 UTC. RBFOX2 remaining confirmation fits were not launched for the frozen completion budget, not for observed effects.

Final PTBP1 primary: 1,851 pairs / 659 true genes; Full MSE .079884 versus Local+pos .112084, delta −.032200, main Bonferroni-2 CI [−.039564,−.024750]. All cells/seeds, raw and exposure exclusion agree. HNRNPK supports Full versus original NoGRU but has no near-capacity control. Prior Phi-specific necessary-control failures remain. Independent B numerical/model/branch checks and C/root actual visual checks passed.

Final Figure 4 candidate: figure/CONTEXT_FINAL.png/pdf/svg, 85×125 mm, a/b/c linked at native coordinates, all 1,553 task-gene points visible. RESULTS.md and FINAL_CAPTION.md match this image. REPORT/REPLAY/ATTEMPTS/QA/source are current; SHA inventory follows source freeze. Old figure/CAPTION.md and REPORT_PRECONTROL are preserved historical outputs and must not be paired with final.

---

'''+status,encoding='utf-8')
print('Final REPORT/REPLAY/ATTEMPTS/STATUS written.')
