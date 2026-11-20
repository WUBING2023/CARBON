# Figure 2c reproduction

From the repository root, run:

```shell
python scripts/verify_figure2c.py paper/source_data/S40_common_target_transfer/cross_model_task_metrics.csv
```

Requires Python 3 and its standard library only. Expected output is in `expected_results.json`. This reproduces point estimates from frozen task-level centered sums; it does not rerun inference or reproduce bootstrap confidence intervals.

Each RBP has two transfer directions. Compute task PCC from `sum_pq / sqrt(sum_pp * sum_qq)`, average directions within each RBP, then average equally over RBPs. CARBON, iDeepB and RBPNet rows must have matching task IDs, cohort hashes and window counts. There are 68 RBPs, 136 directed tasks and 1,421,243 target windows.

Data source: Supplementary Data S40, frozen cross_model_task_metrics.csv. SHA256: 63c601b4a0c7eb4a5735ffc3eedbf62913509f44067acfb53d7c35eaab78880d.

## Files

- [Task-level table](cross_model_task_metrics.csv) (408 rows: 136 tasks × 3 models).
- [Expected output](expected_results.json).
- [Calculation script](../../../scripts/verify_figure2c.py).

## Core columns

- `rbp`, `source_cell`, `target_cell`: protein and transfer direction.
- `source_task`, `target_task`: experiment identities.
- `model`: CARBON, iDeepB or RBPNet.
- `cohort_sha256`: identifier of the matched target-window cohort.
- `n_windows`: number of target windows in the directed task.
- `sum_pp`, `sum_qq`, `sum_pq`: centered sums of squares and cross-products for predicted and observed unit-mass profiles, used to recompute task PCC.
- `pcc_unit_mass_task_pooled`: stored task-level PCC, checked against the centered sums.

Companion metrics and output-quality counts are retained unchanged in the table. The reproduction script uses the columns above.
