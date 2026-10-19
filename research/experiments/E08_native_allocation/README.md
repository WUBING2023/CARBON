# E08_native_allocation

PTBP1: Full、NoGRU与近同容量Local+pos，三seed

Paper: **Figure 4; S1.17; Data S22** (historical mapping; current binding: paper/MANUSCRIPT_BINDING.json).

Selected original analysis code, frozen results and run metadata are listed in INCLUDED_FILES.csv. Omitted inputs/weights/intermediates are in research/catalog/EXTERNAL_ASSETS.csv. SOURCE_MANIFEST.csv records original and packaged identities.

## Evidence boundary

主1851对/659gene/865cell-gene；先gene内再细胞等权。所有对照保留Phi和global SE，结果是架构分配证据，不是保守性独有机制证据。

Historical launchers preserve provenance. Account paths are redacted to /data/CARBON_PROJECT or /data/USER_HOME. Set real paths in a new scratch copy before running them. The checked portable entry points are scripts/verify_frozen_results.py and the AARS CPU demo. Exact historical source snapshots, fresh training and end-to-end biological replication are not certified.
