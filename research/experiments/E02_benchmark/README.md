# E02_benchmark

CARBON、iDeepB、RBPNet历史指标与任务对应

Paper: **Figure 2; S1/S3/S20** (historical mapping; current binding: paper/MANUSCRIPT_BINDING.json).

Selected original analysis code, frozen results and run metadata are listed in INCLUDED_FILES.csv. Omitted inputs/weights/intermediates are in research/catalog/EXTERNAL_ASSETS.csv. SOURCE_MANIFEST.csv records original and packaged identities.

## Evidence boundary

评分规则、ties及输出尺度存在差别。此归档不把旧分数包装成公平SOTA比较。缺失iDeepB权重按本次实查清单记录。

Historical launchers preserve provenance. Account paths are redacted to /data/CARBON_PROJECT or /data/USER_HOME. Set real paths in a new scratch copy before running them. The checked portable entry points are scripts/verify_frozen_results.py and the AARS CPU demo. Exact historical source snapshots, fresh training and end-to-end biological replication are not certified.

Current manuscript clarification: both 0/0 and 4/12 configurations exist in the saved lineage. Do not assign a single threshold contract to all 222 benchmark tasks without run-level evidence.
