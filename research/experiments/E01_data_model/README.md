# E01_data_model

101-nt、4+1通道、profile/count双输出；原始预处理与训练实现

Paper: **Figure 1; Methods; S1.2–S1.6** (historical mapping; current binding: paper/MANUSCRIPT_BINDING.json).

Selected original analysis code, frozen results and run metadata are listed in INCLUDED_FILES.csv. Omitted inputs/weights/intermediates are in research/catalog/EXTERNAL_ASSETS.csv. SOURCE_MANIFEST.csv records original and packaged identities.

## Evidence boundary

示例README里的chr1/2/3分割是历史示例，不能替代正式16/3/3染色体方案；零填Phi不是测得保守性。历史运行使用MARCP命名空间，必须按对应manifest加载。

Historical launchers preserve provenance. Account paths are redacted to /data/CARBON_PROJECT or /data/USER_HOME. Set real paths in a new scratch copy before running them. The checked portable entry points are scripts/verify_frozen_results.py and the AARS CPU demo. Exact historical source snapshots, fresh training and end-to-end biological replication are not certified.
