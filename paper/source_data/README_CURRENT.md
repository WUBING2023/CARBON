# Supplementary Data S1–S53: current repository navigation

S1–S27 retain their existing downloadable primary files. Data S40 now includes the task-level table and a standalone reproduction command. S28–S53 are indexed below using the current supplement; rows without public files are explicitly marked as awaiting deposition. The complete material inventory, original manuscript layout, primary hashes and scientific roles are recorded in [SUPPLEMENT_DATA_BINDINGS.csv](../SUPPLEMENT_DATA_BINDINGS.csv).

The original declared manuscript paths remain in declared_path and declared_primary_path. Large identical assets are stored once; all_paths_json can contain declared aliases that materialize_replay.py restores into a separate scratch copy.

The upstream [manuscript-layout README](README_MANUSCRIPT_LAYOUT.md) is preserved byte-for-byte for provenance. Its packet file counts describe the manuscript delivery, not the physical Git layout; no manuscript packet counts are presented as Git counts here.

| Dataset | Direct primary file | Scientific role |
|---|---|---|
| S1 | [Supplementary_Data_S1_benchmark_metrics.csv](S1_S14_original/Supplementary_Data_S1_benchmark_metrics.csv) | Designated manuscript Supplementary Data S1 table. |
| S2 | [Supplementary_Data_S2_ablation_paired_deltas.csv](S1_S14_original/Supplementary_Data_S2_ablation_paired_deltas.csv) | Designated manuscript Supplementary Data S2 table. |
| S3 | [Supplementary_Data_S3_comparator_win_loss_matrix.csv](S1_S14_original/Supplementary_Data_S3_comparator_win_loss_matrix.csv) | Designated manuscript Supplementary Data S3 table. |
| S4 | [Supplementary_Data_S4_cross_cell_complete_score_set.tsv](S1_S14_original/Supplementary_Data_S4_cross_cell_complete_score_set.tsv) | Designated manuscript Supplementary Data S4 table. |
| S5 | [Supplementary_Data_S5_cross_cell_overlap_join_128.tsv](S1_S14_original/Supplementary_Data_S5_cross_cell_overlap_join_128.tsv) | Designated manuscript Supplementary Data S5 table. |
| S6 | [Supplementary_Data_S6_motif_comparisons.csv](S1_S14_original/Supplementary_Data_S6_motif_comparisons.csv) | Designated manuscript Supplementary Data S6 table. |
| S7 | [Supplementary_Data_S7_conservation_task_summary.csv](S1_S14_original/Supplementary_Data_S7_conservation_task_summary.csv) | Designated manuscript Supplementary Data S7 table. |
| S8 | [Supplementary_Data_S8_conservation_bin_summary.csv](S1_S14_original/Supplementary_Data_S8_conservation_bin_summary.csv) | Designated manuscript Supplementary Data S8 table. |
| S9 | [Supplementary_Data_S9_conservation_peak_centred_profile.csv](S1_S14_original/Supplementary_Data_S9_conservation_peak_centred_profile.csv) | Designated manuscript Supplementary Data S9 table. |
| S10 | [Supplementary_Data_S10_motif_leaf_mapping.csv](S1_S14_original/Supplementary_Data_S10_motif_leaf_mapping.csv) | Designated manuscript Supplementary Data S10 table. |
| S11 | [Supplementary_Data_S11_conservation_perturbation.csv](S1_S14_original/Supplementary_Data_S11_conservation_perturbation.csv) | Designated manuscript Supplementary Data S11 table. |
| S12 | [Supplementary_Data_S12_RBP_robustness.csv](S1_S14_original/Supplementary_Data_S12_RBP_robustness.csv) | Designated manuscript Supplementary Data S12 table. |
| S13 | [Supplementary_Data_S13_perturbation_condition_metrics.csv](S1_S14_original/Supplementary_Data_S13_perturbation_condition_metrics.csv) | Designated manuscript Supplementary Data S13 table. |
| S14 | [Supplementary_Data_S14_count_inverse_sensitivity.csv](S1_S14_original/Supplementary_Data_S14_count_inverse_sensitivity.csv) | Designated manuscript Supplementary Data S14 table. |
| S15 | [conditional_curves.csv](S15_conditional_conservation/conditional_curves.csv) | Conditional conservation curves across the recorded signal/conservation strata; paired contrasts remain in all_paths_json. |
| S16 | [fixed_examples_27nt_per_instance.csv](S16_motif_instances/fixed_examples_27nt_per_instance.csv) | Position-level records for all fixed motif examples; motif and population summaries remain in all_paths_json. |
| S17 | [factorial_rbps.csv](S17_factorial/factorial_rbps.csv) | RBP-level factorial perturbation effects, with task/window records retained as companions. |
| S18 | [HNRNPK_block_mass_loss.csv](S18_fixed_cases/HNRNPK_block_mass_loss.csv) | HNRNPK fixed-case block mass-loss profile; RBFOX2 and all conditional/profile companions remain in all_paths_json. |
| S19 | [event_predictions.csv](S19_ASB/event_predictions.csv) | Current reference-flag-audited 396-event ASB prediction table; archived_pre_audit is historical evidence only. |
| S20 | [benchmark.csv](S20_benchmark_annotations/V5_DISPLAY_UPDATE/benchmark.csv) | Original complete task-level benchmark values; earlier display metadata do not govern the current figure palette. |
| S21 | [seed_paired_validity.csv](S21_seed_sensitivity/seed_paired_validity.csv) | Paired seed/control validity and common-mask PCC differences; training and aggregate seed summaries remain companions. |
| S22 | [PAIR_RESULTS_ALL.csv](../../research/experiments/E08_native_allocation/source_packet/nogru_position_capacity/PAIR_RESULTS_ALL.csv) | Complete pair-level native allocation endpoints for the position/capacity control; all arms/gene summaries remain companions. |
| S23 | [WIDE_ALL_METRICS.csv](../../research/experiments/E09_external_preferences/source_packet/results/WIDE_ALL_METRICS.csv) | Complete frozen task/arm/assay primary Spearman and exact-composition endpoints; exploratory adjacent-direction tables are not primary. |
| S24 | [ATtRACT_DATA_S24_V3.zip](S24_motif_reference_audit/ATtRACT_DATA_S24_V3.zip) | Full ATtRACT audit archive containing query/reference matrices, all comparisons and selection records; MEME license notices added, all pre-existing members unchanged. |
| S25 | [COHORT_ALL7262.csv](../../research/experiments/E10_natural_variants/asb_variants/COHORT_ALL7262.csv) | Complete 7,262-event natural-variant cohort; strict subsets, annotations and sensitivities remain companions. |
| S26 | [DISEASE_ALL13_EVENT_RESULTS.csv](S26_disease_annotations/DISEASE_ALL13_EVENT_RESULTS.csv) | Complete 13-event disease-annotation result table, including unfavorable directions and development exposure. |
| S27 | [ALL66_REPORTER_BIOLOGICAL_MEASUREMENTS.csv](S27_reporter_validation/ALL66_REPORTER_BIOLOGICAL_MEASUREMENTS.csv) | All 66 original reporter biological measurements, including unsuccessful loci and external assay units. |


## November manuscript additions (S28–S53)

| Dataset | Contents in the current supplement | Public files / status |
|---|---|---|
| S28 | RBP-weighted ablations, BiGRU/BiLSTM and configuration search. | Retained by the authors; public deposit pending |
| S29 | Native allocation support-radius intersection and outcomes. | Retained by the authors; public deposit pending |
| S30 | Boundary-aware evaluation code and verification. | Retained by the authors; public deposit pending |
| S31 | Motif discovery commands, input provenance and export bindings. | Retained by the authors; public deposit pending |
| S32 | Matched trajectories, selector comparison and fitted-weight bindings. | Retained by the authors; public deposit pending |
| S33 | Matched native conservation controls and gene-level contrasts. | Retained by the authors; public deposit pending |
| S34 | Matched predicted/measured accessibility fits and outputs. | Retained by the authors; public deposit pending |
| S35 | Eight input representations on shared measured windows. | Retained by the authors; public deposit pending |
| S36 | Original-protocol context fits on complete and measured cohorts. | Retained by the authors; public deposit pending |
| S37 | All 18 architectural transfer weights and target predictions. | Retained by the authors; public deposit pending |
| S38 | Matched Sequence transfer on the full native PTBP1 cohort. | Retained by the authors; public deposit pending |
| S39 | Same-cell HepG2 case selection, all candidate scores, fixed K562 profiles and inference-seed sensitivity. | Retained by the authors; public deposit pending |
| S40 | Common-target transfer, complete outputs and case authentication. | [Task metrics, definitions and reproduction](S40_common_target_transfer/README.md) |
| S41 | Preferred-core direction, 7,404 IG maps and fixed-weight controls. | Retained by the authors; public deposit pending |
| S42 | All 14,808 window correlations, nine paired contrasts and native displays. The companion regional record gives the three-part same-cell decomposition and replay. | Retained by the authors; public deposit pending |
| S43 | The 536-variant screen, 136 fits, verified constructs and public context. | Retained by the authors; public deposit pending |
| S44 | Position rankings, all budgets and sequence-matched controls. | Retained by the authors; public deposit pending |
| S45 | Exact edit shortlists, all 47 sites and 94 parents, and RNA outcomes. | Retained by the authors; public deposit pending |
| S46 | Matched-parent ordering and CAG-start-preserving candidate pools. | Retained by the authors; public deposit pending |
| S47 | Frozen MaPSy cohort, complete model profiles and measured splicing. | Retained by the authors; public deposit pending |
| S48 | Signed reporter responses, sequence controls and complete intervals. | Retained by the authors; public deposit pending |
| S49 | Public U2AF2 knockdown and native eCLIP joins. | Retained by the authors; public deposit pending |
| S50 | Physical deletion constructs, paired RNA activity and SNV proxies. | Retained by the authors; public deposit pending |
| S51 | Fixed-ranking selection budgets and original five-choice intervals. | Retained by the authors; public deposit pending |
| S52 | Three-model substitution comparison, rankings and sensitivity analyses. | Retained by the authors; public deposit pending |
| S53 | ClinVar/ENCODE input identities, exact cohort, profile responses, consequences, clinical terms and fixed-case records. | Retained by the authors; public deposit pending |

The S40 package reproduces the principal PCC point estimates. The pending rows are an availability index, not download links; publication of those full packets and the full model-weight archive is still outstanding. The manuscript and supplement links on the repository homepage identify the current version.
