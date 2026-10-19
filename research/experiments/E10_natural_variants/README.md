# E10_natural_variants

BEAPR ASB、ClinVar/TCGA/人群注释及原作者reporter测量

Paper: **Figure5c/d/e; S26/S27/S31/S32** (historical mapping; current binding: paper/MANUSCRIPT_BINDING.json).

Selected original analysis code, frozen results and run metadata are listed in INCLUDED_FILES.csv. Omitted inputs/weights/intermediates are in research/catalog/EXTERNAL_ASSETS.csv. SOURCE_MANIFEST.csv records original and packaged identities.

## Evidence boundary

8个临床/肿瘤位点共13事件，均开发曝光；Full7/13、Sequence10/13。注释不是RBP结合实测；11 reporter变异/66测量与选例规则完整保留。

Historical launchers preserve provenance. Account paths are redacted to /data/CARBON_PROJECT or /data/USER_HOME. Set real paths in a new scratch copy before running them. The checked portable entry points are scripts/verify_frozen_results.py and the AARS CPU demo. Exact historical source snapshots, fresh training and end-to-end biological replication are not certified.
