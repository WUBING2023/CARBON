# E07_conservation_asb

冻结模型输入依赖、定位、条件扰动与ASB-positive诊断

Paper: **S1.8/S1.11–S1.14; S11–S20; Figure5d** (historical mapping; current binding: paper/MANUSCRIPT_BINDING.json).

Selected original analysis code, frozen results and run metadata are listed in INCLUDED_FILES.csv. Omitted inputs/weights/intermediates are in research/catalog/EXTERNAL_ASSETS.csv. SOURCE_MANIFEST.csv records original and packaged identities.

## Evidence boundary

固定Phi变换不等于另训模型；396 ASB-positive事件不能计算变异检测AUROC；参考Phi不是ALT进化分数。

Historical launchers preserve provenance. Account paths are redacted to /data/CARBON_PROJECT or /data/USER_HOME. Set real paths in a new scratch copy before running them. The checked portable entry points are scripts/verify_frozen_results.py and the AARS CPU demo. Exact historical source snapshots, fresh training and end-to-end biological replication are not certified.
