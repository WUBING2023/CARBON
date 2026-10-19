# Candidate B: controlled sequence response and motif order

Frozen before Candidate B predictions: 2026-10-04. The Candidate A external labels and results have now been seen. This is a new exploratory diagnostic, not preregistered confirmation. Candidate A's top-20 IG recall and all negative controls remain unchanged and are not replaced.

## Hypothesis and fixed design

Native IG winner rankings combine model sensitivity with native occurrence abundance. A controlled perturbation asks a different question: holding the surrounding RNA sequence and conservation context constant, does the model distinguish measured intrinsic preferences between concrete five-base rules? A useful biological rule must also recognize base order rather than only count cytosines or pyrimidines.

For each of the six matched tasks use the first 16 hash-selected, qualified native windows already frozen in Candidate A. Replace positions 48–52 (zero-based inclusive) by every one of the 1,024 A/C/G/T five-mers. Keep all other sequence positions unchanged. For each matched model (three seeds; aligned, sequence_only, mean, permuted), compute delta native log_total_count relative to that background's original sequence. Average deltas equally over the same 16 contexts. No context is selected by binding, model correctness, IG or external label. This is an in-silico sequence intervention evaluated against real external preferences, not a new measured synthetic-RNA experiment.

Conservation is the native WT background and is fixed while every five-mer is substituted. Mean and Perm are the corresponding separately trained models with their original training-defined conservation transforms; fixed original source_row determines Perm. It is not recomputed per five-mer or chosen by model agreement. The sequence-only model receives four channels. No new output head, k-mer length or context number will be selected by the results.

## Exactly two endpoints

Primary: Spearman correlation of mean predicted delta log_total_count against the complete same-protein 5-mer external preference vector. All 1,024 five-mers are evaluated. RBNS concentration and PTBP1's two RNAcompete experiments follow Candidate A's already frozen source rules.

Secondary: within-composition pairwise concordance. Group the 1,024 words by exact counts of A/C/G/T. Enumerate every unordered pair within every group with at least two words. An experimental pair is eligible when its two measured values are not exactly tied (no tolerance or effect-size selection). A model receives 1 for a correct direction, 0 for an incorrect direction, and 0.5 for an exact prediction tie. Every eligible pair has equal weight; denominators and model tie counts are exported. Both words have identical mononucleotide composition, so C-count/Pyrimidine-count cannot discriminate them and score 0.5. This tests arrangement sensitivity, not binding-site position accuracy. Pair/word counts are not independent biological replicates and do not justify tiny population P values.

Report per-model-fit, background and protein values, paired Full–Sequence/Mean/Perm differences and all source assay columns. Average seed/background only within protein, and PTBP1 experiments without selecting the more favorable experiment. Show the known-consensus and frequency baselines from Candidate A on the same two endpoints; native Phi priority is a context-dependent sequence association baseline and is retained as such, not a genuine artificial sequence-dependent Phi measurement.

Optional conditional uncertainty is a paired bootstrap over the 16 common sequence contexts (1,000 fixed-seed draws). It quantifies context sampling stability within these proteins, not cross-protein generalization. With three proteins, preserve exact values and do not turn contexts, seeds, k-mers, pairs or technical assay conditions into biological replication. No new training or additional measured labels.

## Stop/interpretation rules

Full exceeding Sequence alone supports a family comparison, not aligned conservation specificity. Aligned conservation specificity additionally requires exceeding both matched Mean and Perm on the fixed endpoint. Correct ordering among identical-composition words supports an arrangement-aware learned rule; it does not imply functional RNA regulation or causality in cells. A weak or negative result does not trigger a search over heads, k-mer sizes, favorable proteins or background subsets. This is Candidate B of at most three bounded hypotheses.
