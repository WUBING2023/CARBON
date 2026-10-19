# Candidate Wide: availability-defined external model-family audit

Frozen 2026-10-04 before evaluating new protein labels. This expands Candidate B, not its endpoint or vocabulary. HNRNPK/PTBP1/RBFOX2 are explicitly exposed proteins and are reported separately from the other 25. The original author matrix is public; availability filtering has been completed without reading new model–label agreement.

## Population and source

The full availability-defined intersection is 28 human RNA-binding proteins, 42 CARBON training-background tasks. WIDE_MANIFEST.json freezes every gene, matching RNAcompete experiment ID and Full/Sequence checkpoint path/hash. All local checkpoint hashes match the prior frozen research inventory. Include all same-human-gene author-measured experiment columns, not only a favorable replicate, high-Z experiment or successful model. The author's set-A/set-B quality measures and protein-construct description are reported; poor-quality assay columns remain identifiable. No inferred JPLE profile is used as measurement.

## Fixed input and models

Select 1,024 qualified finite unambiguous one-hot 101-nt native windows per task from raw chr8/15 five-channel files using Candidate A's exact hash(task,original chromosome,original source_row) rule. Use the first 16 as controlled contexts, preserving the same 1,024 five-base substitutions at positions48–52 and scalar delta native log_total_count. No RNA/k-mer length, model output head, context number or seed is optimized on agreement. Cell backgrounds are model contexts, not independent protein replicates. If a task lacks required source files or qualified windows, record that eligibility failure rather than selecting a performance-based replacement.

The old Full and Sequence models are not proven to have matched initializations or a matched random training path. Therefore this is an external model-family rule validation, not a causal estimate of the conservation module. Additional Full frozen-flat and frozen-permuted inference hold its weights unchanged; these are input-dependence controls, never mislabeled separately trained Mean/Perm. Full flat uses native101-nt mean Phi; frozen permutation uses the original task/chromosome/source_row deterministic PCG64 mapping, fixed across all substituted words. Matched trained Mean/Perm for the original three proteins remain a separate evidence layer.

## Two endpoints and universal controls

Primary: complete1,024-word Spearman between mean predicted delta count and true RNAcompete preferences. Secondary: within-exact-mononucleotide-composition pairwise concordance, with true exact ties excluded and predicted ties receiving0.5. Same definitions as Candidate B. Source7→5 mapping is the original occurrence-preserving mean, with unavailable7-mers absent rather than zero-imputed.

For every protein independently report all four universal nucleotide-count baselines A-count/C-count/G-count/T-count, plus native-word frequency and Phi association. No best letter is chosen per protein after seeing the labels. The four base-count predictors all tie within composition and receive0.5; that is expected and useful. Model-independent baselines have no model seeds. No PWM is fitted to evaluation labels and relabeled independent.

Each measured experiment is reported, then experiments and cell backgrounds are averaged within protein with equal weights. The 28 true proteins are the generalization unit; the new25 and exposed3 are separate groups. Biological-population summaries do not count1,024 words, within-composition pairs, sixteen backgrounds, cells or old model seeds as independent n. Retain full model/assay/word tables, all negative comparisons and exact protein-level differences. Report exact counts of proteins where Full exceeds Sequence and both frozen controls; do not call those counts a Phi causal effect.

## Budget and execution

Pilot first lexicographic available task CPEB4_K562_ENCSR046JHH on CPU8. If the projected full four-condition cost is at most90minutes, run the entire frozen population. If cost is longer, freeze a deterministic gene-name priority before any model agreement, truncate by time and report the realized availability prefix. GPU access requires resource-coordinator approval and a process-free guard; do not compete with other tasks. No training. Deadline2026-10-04 14:47:07UTC includes reproducibility/figure/review; no extension after context compaction.

A broad match validates learned intrinsic sequence rules. It cannot establish cellular conservation causality or downstream regulatory function. A composition-only gain is separated from arrangement-aware concordance and is not presented as improved motif syntax.
