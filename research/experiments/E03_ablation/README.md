# E03_ablation

222任务×5去除项，另含完整baseline：模型组件预测价值

Paper: **Results 5.2; Data S2** (historical mapping; current binding: paper/MANUSCRIPT_BINDING.json).

Selected original analysis code, frozen results and run metadata are listed in INCLUDED_FILES.csv. Omitted inputs/weights/intermediates are in research/catalog/EXTERNAL_ASSETS.csv. SOURCE_MANIFEST.csv records original and packaged identities.

## Evidence boundary

任务不是生物学重复；noGRU容量改变；no_multiscale容量也改变；reconstruction-only存在指标取舍。保守性预测价值不等于功能解释正确性。

Historical launchers preserve provenance. Account paths are redacted to /data/CARBON_PROJECT or /data/USER_HOME. Set real paths in a new scratch copy before running them. The checked portable entry points are scripts/verify_frozen_results.py and the AARS CPU demo. Exact historical source snapshots, fresh training and end-to-end biological replication are not certified.
