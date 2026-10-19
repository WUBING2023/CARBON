# From code to experiment, result and figure

Each experiment README states its evidence boundary and INCLUDED_FILES.csv enumerates the actual files. research/catalog/SOURCE_MANIFEST.csv retains original source SHA256, current SHA256 and any path-only transformation. EXTERNAL_ASSETS.csv describes omitted originals. paper/FIGURE_BINDINGS.csv enumerates every used main/supplement figure, its output hash, dataset group and verified or unresolved plotting entry.

| Experiment | Code/records and frozen result entry | Manuscript mapping |
|---|---|---|
| E01 data/model | CARBON/, carbon_train.py, preprocessing/; research/historical_model/; archived checkpoint/run catalogs | Figure 1; Methods; S1.2–S1.6 |
| E02 benchmark | E02_benchmark/paper_data/Supplementary_Data_S1_benchmark_metrics.csv; canonical weight availability catalog | Figure 2; historical comparator scores; Supplementary benchmark views |
| E03 ablation | ABLATION_RUN_INDEX.csv → each real test_metrics_summary.txt/config; Supplementary_Data_S2_ablation_paired_deltas.csv; results/ | 222 tasks × baseline and five ablations; Figure S1; Data S2 |
| E04 diagnostics | seed_sensitivity source/manifest and full recorded diagnostics; optimization and BiLSTM/BiGRU records | S1.5/S1.16; Data S21; limited multi-seed checks |
| E05 cross cell | Complete-score and 128-overlap tables; original plotting script | Figure S2; Data S4/S5 |
| E06 motifs | frozen_atlas_inputs queries/Tomtom outputs; reference_provenance; ATtRACT derived matches; current Figure 3 inputs | Figure 3; Figure S21; Data S6/S10; 219-task atlas |
| E07 conservation/ASB | Complete perturbation/ASB frozen metrics and own analysis code | Conservation/ASB supplement; Figure 5 diagnostic context |
| E08 allocation | source_packet/nogru_position_capacity/PAIR_RESULTS_ALL.csv → target/share/loss → METHODS.csv, gene/cell results; own model/control/analysis code | Figure 4; allocation supplement; Data S22 |
| E09 external preferences | source_packet/scripts/evaluate_wide.py; complete WIDE_ALL_5MER_PREDICTIONS.csv + WIDE_EXTERNAL_5MER_SCORES.csv + WIDE_ALL_METRICS.csv | Figure 5a/b; all 25 primary proteins, Full/Sequence/frequency; Data S23 |
| E10 natural variants | disease_annotations/DISEASE_ALL13_EVENT_RESULTS.csv; reporter full measurement/selection tables; complete ASB event data via canonical alias | Figure 5c–e; S26/S27/S31/S32; all favorable and failed directions |
| E11 functional negative | functional_rules/fr1–fr4 full result cohorts and reports; six-RNA functional-map outcomes | Negative/mixed functional exploration, including directions not highlighted in main figures |

Current figure renderer entry points and their small frozen inputs are under paper/editable_figures/. Updated Figure 2/4 share byte-identical packets and Figure 5 reuses full event tables; restore duplicate paths with materialize_replay.py before rendering in scratch. Source data, run metrics, manuscript figures and new verification reports remain distinct artifacts.

Plot entry points marked unresolved in FIGURE_BINDINGS.csv are provenance gaps. A final rendered PDF alone is not evidence that its whole historical upstream workflow was reproduced in this package. The audit checks only the numerical paths listed in its report.
