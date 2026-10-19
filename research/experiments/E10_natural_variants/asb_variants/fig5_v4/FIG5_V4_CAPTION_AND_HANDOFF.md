# Compact Fig. 5 v4: caption and integration handoff

One native-vector figure, 178 × 142 mm. The v3 figure was 178 × 190 mm; height is reduced by 25.3%. Top: three examples. Bottom: two complete cohorts. No new inference, training, experiments, endpoint changes, example selection or scientific statistics.

## Panel mapping required in manuscript

| v3 | v4 | Evidence object |
|---|---|---|
| a | a | QKI three neighboring edits towards CUAAC |
| b | b | HNRNPA1 UGAGG→UAGGG |
| d | c | Natural C>T locus and all ten RBP events |
| c | d | All 25 external proteins and 1152 paired responses |
| e | e | All 396 natural events and five-arm accuracy |

Update explicit Results references and the caption accordingly: natural example is now 5c, external cohort 5d, and gold event columns in 5e identify 5c. Methods references to 5e remain valid. The root owns these TeX edits; no manuscript prose was edited by this worker.

## Proposed concise English caption

**Figure 5. Sequence-order preferences and natural allele-specific binding.** **a,b**, QKI neighboring edits towards CUAAC and HNRNPA1 UGAGG→UAGGG, independent RNAcompete ΔZ5, and all 32 paired model-background responses per edit (16 HepG2 circles, 16 K562 triangles; ticks denote means). These computational backgrounds are not biological replicates. The PhyloP heatstrip shows the unchanged reference input. **c**, Coverage-first C>T example at GRCh38 chr8:1,844,440 (ARHGEF10/KBTBD11-OT1, K562): aligned RNA alleles, identical reference PhyloP and all ten BEAPR-labelled RBP events for Full, Sequence, flat-Φ and permuted-Φ inputs. Responses are predicted REF-minus-ALT log counts (C−U); positive values favor REF, the observed preferred allele at this locus. **d**, All 25 proteins: the frozen external-first, frequency-opposed, same-composition edit, RNAcompete-derived ΔZ5, native-frequency change, primary all-1024-word Spearman correlation (Full/Sequence/Frequency), and all 1152 paired responses from 36 tasks. H/K denote HepG2/K562; gray denotes an unavailable task. Positive responses favor the word with higher external Z5. The primary Full-minus-Frequency contrast and secondary within-composition concordance contrast appear below; all endpoints remain in Data S23. **e**, All 396 previously exposed chromosome-8/15 events (287 variants, 85 RBPs, 119 RNA-region clusters), ordered by reference-center PhyloP and variant identity independently of prediction outcome. Shared columns link PhyloP, genomic REF/ALT, BEAPR-preferred allele, genomic region and four responses; gold columns locate c. Here responses are predicted observed-preferred-minus-other log counts: positive agrees with BEAPR, irrespective of REF/ALT identity. c–e share the same asinh color mapping, white at zero and original-value range −1.8 to +1.8 (linear width 0.02), without model-specific normalization. PhyloP is identical between alleles within each arm. Frozen RBP-equal quartile and five-arm aggregate accuracies are retained. Overall intervals use two-way RNA-region/RBP multipliers; Full-minus-Sequence is −6.73 percentage points (95% interval −15.97 to +1.64). Flat/permuted Φ are fixed-checkpoint controls; REF always predicts reference preference. These exploratory ASB-positive data establish neither a conservation-specific advantage, ASB detection ability nor pathogenicity. S1.18/S1.20; Data S23/S25.

## 中文逐面板解释

- **a，QKI：** 三个相邻五碱基序列换成 CUAAC，碱基组成不变。中间是独立 RNAcompete 偏好提高的量，右侧是同一个序列背景中 Full 和 Sequence 的预测变化；各点不是新的生物实验。Φ 在替换前后固定。
- **b，HNRNPA1：** UGAGG 换成 UAGGG，外部偏好提高 5.05，但原训练序列中的频率下降 101。两组模型点显示这种“外部更偏好、训练中更少见”的排列能否得到相应预测响应。下方定义 Full、Sequence 和固定 Φ，帮助解码。
- **c，真实天然变异：** 展示一个 C>T 位点的 RNA REF/ALT、相同的参考 Φ，以及全部十个 RBP 的预测。BEAPR 实测都较偏好 C；正值表示模型更偏向 C，负值表示更偏向 U，不能把所有模型响应解释为成功。
- **d，外部完整蛋白队列：** 每行一个蛋白，顺序与原冻结队列一致。左侧的五碱基替换由外部 Z5 和频率选定，右侧保留全部 Full/Sequence 响应；灰色为不存在的任务。彩色三列是全部 1024 词的主终点，底部另保留同组成次终点，二者不能混用。Full 相对 Frequency 有正差异，但并未证明相对 Sequence 的整体优势。
- **e，天然完整事件队列：** 每列一个冻结事件，Φ、REF/ALT、BEAPR 偏好、基因组区段和四个模型共享同一列。按 Φ 和位点身份排序，未按模型收益排序；金框对应 c 的十个事件。这里正值是偏向“BEAPR 实测较偏好的等位基因”，不能统称偏向 REF。下方所有五臂和四分位数保留：Full 52.04%、Sequence 58.77%，差值区间跨零，不能据此宣称 Φ 改善天然变异准确率。

## Frozen scientific values and statistical units

| Item | Frozen value / method |
|---|---|
| External scope | 25 proteins; all 1024 five-mers; 36 available tasks; 16 retained contexts per task; 1152 paired Full/Sequence response cells |
| Primary endpoint | Per-protein all-word Spearman; Full-minus-Frequency mean +0.1461870204, 95% CI [+0.0380167321, +0.2479966561] |
| Secondary endpoint | Within-composition concordance; Full-minus-Frequency +3.7808836395 pp, 95% CI [+0.8987046150, +6.7330172791] |
| External Full-minus-Sequence | Primary −0.0118649541, CI [−0.0665064316, +0.0387161155]; secondary −0.6886701662 pp, CI [−1.9115452756, +0.5054363517] |
| External intervals | 10,000 paired protein bootstrap draws, seed 20261004; no resampling rerun |
| Natural scope | All 396 eligible events / 287 variants / 85 RBPs / 119 RNA-region clusters |
| Natural quartiles | Q1/Q2/Q3/Q4 have 116/101/86/93 events; unique-locus frozen boundaries retain tied Φ values |
| All five overall accuracies | Full 52.0355866221%; Sequence 58.7697258919%; Flat Φ 53.4560185630%; Permuted Φ 52.3647182748%; REF rule 64.7674648044% |
| Full-minus-Sequence | −6.7341392698 pp; two-way CI [−15.9702373697, +1.6411083572] |
| Natural aggregate intervals | 2,000 two-way RNA-region/RBP multiplier draws, seed 20261015; source unchanged |
| Natural quartile intervals | Frozen means shown; original 5,000 RBP-bootstrap intervals, seed 20261012, remain in Data S25 / supplement |
| Shared response display | Same asinh mapping in c/d/e, linear width .02, range [−1.8,+1.8], white zero; not model-specific normalization |
| Shared raw Φ display | TwoSlopeNorm [−7,0,+10]; raw reference PhyloP, unchanged between alleles |

The natural cohort is ASB-positive, previously exposed and in the ENCODE assay family. Genomic region/1000G/ClinVar annotations are context, not binding-effect or pathogenicity truth. The selected natural example is the frozen maximum-RBP-coverage case, not a model-agreement selection. a/b examples remain previously exposed illustrations.

## Exact field bindings

| Panel / object | Source file(s) | Fields / rendering |
|---|---|---|
| a/b | fig5_v3/fig5.pdf and fig5.svg | Literal accepted v3 vector crops; uniform scale 57.5/85 = 0.6764705882 on both dimensions; source axis ticks, shapes, dots, connectors and colors unchanged |
| a reference Φ | b1_reference/INTERVENTION_PACKET.npz | Original QKI_native_conservation and central positions 49–53 already present in the accepted v3 inset; entire inset scaled with a |
| c alleles / Φ | source_inputs/EVENTS_INCLUDED.json.gz | context.reference_sequence_rna_direction_dna, minor_sequence_dna, phi_raw; RNA T→U labels only; positions 46–56, central Φ51 −3.5009999275 |
| c all ten responses | FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv | RBP, event_id, ref_preferred; A_native/S/A_flat/A_permmean_ref_minus_alt |
| d five-mer / ΔZ5 / frequency | fig5_v3/frozen_inputs/ALL_25_SELECTED_EDITS.csv | Original protein order; lower_assay_word, higher_assay_word, external_delta_Z5, native_frequency_delta |
| d contexts | Same directory SOURCE_PACKET.npz, ALL_36_TASKS.csv, ALL_1152_PAIRED_MODEL_RESPONSES.csv | Per-protein *_cells, *_full_paired_delta, *_sequence_paired_delta; H then K, 16 contexts each available task; unavailable rows gray |
| d primary / secondary | manuscript/source_data/S23_external_order_preferences/FROZEN_PRIMARY_METRICS.csv and FROZEN_PRIMARY_CONTRASTS.csv | spearman_all5mer; within_composition_concordance; frozen mean_delta and protein_paired_bootstrap_low/high |
| e full atlas | STRICT396_WITH_ANNOTATIONS.csv | event_id, variant_key, RBP, cell, phi_center, ref, alt, ref_preferred, NCBI_annotation_region; four *_major_minus_minor response fields |
| e display ordering | fig5_v4/EVENT_ATLAS_COLUMNS.csv | Byte-identical to v3 display table; phi_center, numeric chromosome, pos1, ref, alt, cell, RBP, event_id |
| e quartile means | QUARTILE_ASSIGNMENTS_strict_test.csv and DIRECTION_BY_CONSERVATION.csv | Original event assignment; strict_test arm/quartile accuracy |
| e overall / intervals | DIRECTION_SUMMARY.csv and DIRECTION_TWO_WAY_SENSITIVITY.csv | strict_test accuracy; two-way ci_low, ci_high including A_native minus S |

## References actually inspected before drawing

Complete original [PrismNet Fig. 6](https://www.nature.com/articles/s41422-021-00476-y/figures/6), [Fig. 7](https://www.nature.com/articles/s41422-021-00476-y/figures/7), [Scrambler Fig. 3](https://www.nature.com/articles/s42256-021-00428-6/figures/3), and the complete accepted CARBON v3 figure were inspected before the first v4 render. Exact local original images remain in fig5_v3/references; their publisher binding and hashes remain in REFERENCE_PROVENANCE.json.

Fig. 6 supplied the compact three-column top row with cohort evidence below. Fig. 7 supplied literal allele identity linked to coordinate-matched readouts; Scrambler supplied explicit edit identity linked to paired responses. The supported transfer is arrangement and matched identities. PrismNet disease enrichment, experimental structure/splicing validation and Scrambler APA function are not transferred to CARBON. REFERENCE_CURRENT_V4_SAME_WIDTH.pdf provides the complete references and v3/v4 at the same 178-mm width.

## Actual-size and preservation QA

Full 178-mm PDF preview, independent Cairo-rendered SVG preview and 300-dpi top/d/e detail crops were actually viewed. Five-mer letters remain 5.8 pt with 1.4-mm monospaced center spacing; all five letters are distinct. H/K are 4 pt. Original a/b axis labels were uniformly scaled from 6.5 to 4.397 pt; the smallest Φ colorbar labels scale to 3.518 pt. This reduction is the user's authorized whole-panel scaling, not a claim of unchanged physical font size.

At 300 dpi, a/b have zero changed pixels versus independently placed original v3 vector crops outside the panel-letter masks (318985 and 185692 compared pixels). All 17 native heatmaps / 5690 color cells were matched to the frozen source values and shared normalization. All 396 event columns, not just counts, are byte-identical to the previous display table. All five confidence intervals fit their 40–80% axis. No scientific raster content is embedded in the manuscript PDF.

The added b decoder and reset a/b panel letters are the only new overlays on their crops. No source a/b data, axis ticks, point positions or graphical encodings were redrawn.

Deliverables: fig5.pdf, fig5.svg, fig5.png, fig5_178mm_preview.png, independently rendered fig5_svg_178mm_preview.png, render_fig5_v4.py, qa_fig5_v4.py, FIG5_V4_SOURCE_MAPPING.json and FIG5_V4_QA.json. All input/export SHA256 values are recorded in the source mapping; final PDF hash is in QA. Only the final fig5.pdf is installed in manuscript/figures.
