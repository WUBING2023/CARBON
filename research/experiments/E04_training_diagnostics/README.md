# E04_training_diagnostics

保留搜索过程、有限多seed和BiGRU/BiLSTM比较

Paper: **S1.5; S1.16; Data S21** (historical mapping; current binding: paper/MANUSCRIPT_BINDING.json).

Selected original analysis code, frozen results and run metadata are listed in INCLUDED_FILES.csv. Omitted inputs/weights/intermediates are in research/catalog/EXTERNAL_ASSETS.csv. SOURCE_MANIFEST.csv records original and packaged identities.

## Evidence boundary

早期搜索用过chr1；主222任务并非全部多seed。不得把小范围种子诊断当全量重复训练。

Historical launchers preserve provenance. Account paths are redacted to /data/CARBON_PROJECT or /data/USER_HOME. Set real paths in a new scratch copy before running them. The checked portable entry points are scripts/verify_frozen_results.py and the AARS CPU demo. Exact historical source snapshots, fresh training and end-to-end biological replication are not certified.
