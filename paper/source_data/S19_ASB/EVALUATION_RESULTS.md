# Completed exploratory ASB direction experiment — 2026-09-24

**The frozen primary endpoint did not establish robust allelic-direction accuracy above 50%, nor an improvement from aligned conservation.** These results must be retained as a limitation of extrapolating window-profile performance to variant-direction prediction. The positive signed-score secondary endpoint must not replace the unsuccessful prespecified primary endpoint.

All 396 candidates passed the input gate: reference central nucleotide belongs to the reported allele pair,101 canonical nucleotides, finite PhyloP. Independent reconstruction of five positive- and five negative-strand original training-chr1 H 5 windows matched raw sequence and one-hot encoding exactly; conservation maximum error was 0 in all ten windows. The negative-strand sequence is reverse-complemented, with reversed PhyloP. No source allele labels were relabeled. All 112 checkpoint hashes matched the frozen inventory. Same-input numerical count/profile checks gave zero error. Inference used existing five-channel CARBON models on CPU with two threads and nice 10; it completed all 396×5×2 allele-condition predictions with no exclusion or nonfinite outputs.

| Condition | Event-weighted agreement | Equal-RBP agreement | 95% RBP-cluster interval |
|---|---:|---:|---:|
| Aligned | 53.28% | 52.04% | 44.51–59.36% |
| Window mean | 53.79% | 53.46% | 45.69–60.95% |
| Permutation 17 | 53.54% | 51.41% | 44.02–58.66% |
| Permutation 43 | 54.17% | 54.89% | 47.23–62.23% |
| Permutation 97 | 53.28% | 49.75% | 42.15–57.17% |

Aligned predictions agreed with the supplied favored allele in 211 events and disagreed in 185; no aligned ties. The permutation 43 condition had one numerical tie, scored 0.5 by the predefined 1e-6 rule. The primary aggregation averaged events within task, tasks within RBP, then 85 RBPs equally. Intervals used 5,000 paired RBP-cluster bootstrap replicates.

Aligned minus window-mean agreement was−1.42 percentage points (95% interval−5.23 to+2.18). Aligned minus the three permutations was+0.62 [−4.98,+6.54],−2.85 [−8.86,+2.94] and+2.29 [−2.36,+7.59] percentage points. Every paired interval includes zero. Selecting the most favorable permutation comparison would therefore be inappropriate.

All prespecified sensitivities gave comparable primary estimates with intervals spanning 50%: total named-allele coverage>=20 (395 events,52.33%); each replicate>=10 (310,51.90%); raw direction concordant in both replicates (387,52.33%); no other published same-task ASB within 300 bp (323,52.83%); joint coverage/concordance/proximity (316,52.96%). None creates independent balanced controls or identifies an isolated causal variant.

The aligned raw-event median signed count-head contrast was 0.002094. The equal-RBP mean signed contrast was 0.01058 [0.00089,0.02179]. This secondary amplitude result coexists with the weak primary direction result and does not demonstrate superior direction accuracy, calibration or causal mutation effects.

The model field `log_total_count` is an implementation name: training targets log(1 + observed total count). Accordingly, all scalar effects here refer to count-head output differences. This terminology correction changes no calculation; the frozen protocol and execution hashes are retained unchanged.

## Figure caption

**Exploratory allelic-direction prediction on published ASB events.** (a) Major-minus-minor count-head outputs were evaluated for 396 events representing 287 genomic positions, 112 RBP–context tasks and 85 RBPs. Major/minor labels retain the direction supplied by BEAPR after bias normalization. Circles show task-then-RBP equally weighted direction agreement, intervals show the 2.5th and 97.5th percentiles of 5,000 RBP-cluster bootstrap replicates, and short ticks show event-weighted estimates. Scores within ±1e-6 were assigned half agreement. (b) Paired aligned-minus-control changes in the equally weighted direction agreement, using the same bootstrap draws. No measured non-ASB controls were available, and no AUROC was computed. The source assays belong to the same ENCODE eCLIP family as model training; chr8/15 exclusion limits direct locus leakage but does not constitute an independent assay replication.

## Manuscript-ready Methods text

We performed an exploratory positive-only allelic-direction evaluation using the ASB events released by Yang et al. (2019). Source hg19 coordinates were converted to hg38 using the UCSC chain, requiring a unique mapping, valid allele symbols, a matched RBP–context model, and a position on held-out chromosomes 8 or 15. Each 101-nt reference window was centered on the variant and required canonical sequence, finite PhyloP values, and a central reference base belonging to the reported allele pair. Sequence and conservation were oriented according to the RNA strand. Input reconstruction was checked against five positive- and five negative-strand original training windows before allele inference. The reported bias-normalized major/minor direction was preserved. Frozen CARBON models scored both alleles under aligned conservation, a window-mean track and three within-window permutations; the alleles shared each conservation track. The primary score was the count-head output for the major allele minus that for the minor allele. The head was trained against log(1 + observed total count); its output was not treated as a calibrated log expected count. Direction agreement was 1 for scores > 1e-6, 0 for scores < −1e-6 and 0.5 otherwise. Task-level agreement was averaged within RBP and then across RBPs. We used 5,000 fixed-seed RBP-cluster bootstrap draws for 95% percentile intervals and paired condition differences. Coverage, raw-replicate concordance and published-site proximity were evaluated as prespecified sensitivities. This selected ASB cohort has no experimentally measured balanced-binding negatives and shares ENCODE assay families with the training data.

## Manuscript-ready Results text

The additional ASB analysis exposed a boundary of the current model. Across 396 events from 85 RBPs, aligned-conservation predictions recovered the reported favored allele in 53.3% of events; equally weighting RBP groups yielded 52.0% agreement (95% RBP-cluster interval 44.5–59.4%). Aligned conservation did not improve direction agreement over window-mean conservation or any of the three permutation controls with a paired interval excluding zero. Coverage, replicate-concordance and proximity sensitivities gave similar results. Thus, the window-profile improvements and sequence–conservation sensitivity analyses do not yet establish reliable allelic-direction prediction.

## Files and remaining interpretation limits

`evaluation/event_predictions.csv` contains every event, source label, model hash, allele output and condition contrast. `predictions.npz` preserves all profile vectors and count-head outputs (saved under the implementation key `log_total_count`); `inputs.npz` preserves exactly the screened model inputs. `summary.json`, `task_statistics.csv` and `rbp_statistics.csv` contain primary and sensitivity statistics. The runner and frozen protocol are retained locally and in the isolated server directory. `run_config.json` binds code/input/protocol/checkpoint inventory hashes.

The primary interval is an RBP-cluster uncertainty estimate;287 genomic positions recur across RBPs, so it does not model every form of cross-RBP locus dependence. Source ENCSR accession identity could not be certified from the ASB table, which lists symbols/contexts only. Sparse task denominators preclude strong individual-RBP conclusions. No model retraining, threshold optimization, biological validation or main-manuscript replacement occurred in this subtask.
