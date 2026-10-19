# Supplementary data — CARBON V13, 21 September 2026

S1–S10 retain the approved V12 CSV/TSV contents without numerical changes.
Their figure correspondence is listed in the supplementary PDF, section 4.

S10 is the retained 219-row mapping for Figure 3a. Task IDs and ENCODE accessions identify
the leaves; crop positions are zero-based, end-exclusive and always nine bases.
The displayed pattern is positive pattern 0. Best q values and match counts
refer only to retained comparisons in S6; empty fields mean no retained
comparison, not a tested negative. Columns containing `any_pattern` use all
retained patterns for that task and therefore have different denominators.

Annotation columns preserve the recovered source codes: `Yes` (functional
match) or a domain code such as `RRM`, `ZnF`, `KH` means a matched rule; `Not`
is displayed as **No** in the figure and means no matching cached annotation.
`Unknown` denotes unavailable annotation. No/Not does not demonstrate biological
absence. Operational GO/Pfam rules and the narrow ZnF list are specified in
Supplementary Section S1.7.

The compact main-paper circle retains the V11 layout; S10 supplies the exact leaf identities.


## V13 additions: S11–S14

These tables accompany Supplementary Sections S1.8 and S1.11. The new frozen-model
experiment uses 706,642 windows across 222 tasks / 150 RBP identities, from
chromosomes 8 and 15 only. These chromosomes were previously evaluated; this is
restricted re-evaluation, not a newly independent validation set. `cell_line=other`
in S11/S13/S14 identifies the two adrenal-gland tasks.

### S11 — conservation perturbation (222 task rows)

- `task`, `rbp`, `cell_line`: task identifier, protein identity and context.
- `n_windows`, `n_valid_profile_windows`, `n_chr8`, `n_chr15`: selected and valid
  sample counts. All selected windows had nonconstant observed profiles.
- `primary_median_paired_pcc_drop` and `drop_permutation_mean` are identical:
  median over windows of [true-input within-window PCC minus the mean of the
  three permutation PCC values for that same window]. This is not a difference
  of unpaired medians or a pooled task PCC.
- Other `drop_*` columns: task median of paired true-minus-condition within-window
  PCC. Positive values denote loss under that condition. `no_conservation` uses
  the separately trained four-channel model. `shift_*` denotes a circular shift
  in nucleotides; `perm_17/43/97` identifies permutation seeds.
- `Q1_primary_drop` through `Q4_primary_drop`: the same paired endpoint within
  original whole-window mean-PhyloP quartiles. Stable sorting and four array
  splits define these bins; they are not the peak-core bins used in Figure 4.

Figure 5 summaries first average task effects within RBP and then take the mean
over 150 RBPs. Paired percentile intervals use 5,000 RBP resamples, seed 20260921.
This bootstrap does not estimate variability over training seeds.

### S12 — RBP conservation robustness (4,800 rows)

One row per RBP × stratification × metric × quartile (150 × 2 × 4 × 4).
`n_tasks` is the number of contributing tasks for that protein;
`stratification` is peak-core or whole-window PhyloP; `bin_id=0,1,2,3` means Q1–Q4.
`mean_task_gain` averages the S8 task/bin gains within each RBP.
`mean_task_endpoint_gain` is the within-RBP mean Q4-minus-Q1 gain, repeated across
the four bin rows for convenience; do not count these repeated entries as
independent endpoint measurements. Correlation metrics in S8 are median paired
within-window gains; AP metrics are paired differences of pooled within-bin AP.
Figure 4e uses the MEDIAN over 150 RBP endpoint effects, unlike the MEAN used for
Figure 5. Its percentile intervals use 5,000 paired RBP resamples, seed 20260921.

### S13 — condition scores (3,330 rows)

One row per task × condition (222 × 15). `pooled_pcc` and `pooled_spearman` are
correlations of flattened observed/predicted positions on the selected windows;
`signal_ap` is average precision with observed signal above the task-wide median
as positive. `signal_positive_fraction` reports that positive fraction.
`median_window_pcc` is a separate within-window summary.
`median_paired_pcc_drop` has the true-minus-condition orientation defined for S11;
the true row is zero by construction. The three permutation conditions remain
separate rows. Their score average is not a score of averaged predictions.

### S14 — inverse-count sensitivity (222 rows; metric DELTAS)

IMPORTANT: `pooled_pcc`, `pooled_spearman` and `signal_ap` in this file are changes,
not absolute scores. Each equals CORRECTED MINUS ORIGINAL, where the original
output is p·exp(z) and the corrected output is p·max(exp(z)−1,0) from the same
forward pass. Thus its sign convention differs from the true-minus-control drops
in S11 and S13. Positive means improvement after numerical correction.
`fraction_predicted_total_below_one` is the fraction with exp(z)<1 (zero for all
tasks here). `median_relative_total_change` is the median removed fraction of the
original predicted total: median[min(exp(z),1)/exp(z)], a nonnegative magnitude.
The original benchmark S1 is not replaced by these sensitivity results.

## Re-analysis files

`analysis/summarize_perturbation.py` rebuilds S11/S13/S14 and the Figure 5 summary
from `experiments/outputs/*/summary.json`. `scripts/build_fig4_v13.py` builds S12
from unchanged S8. Per-window endpoints and selected row IDs are retained in the
separate experiment archive. No table establishes biological causality.
