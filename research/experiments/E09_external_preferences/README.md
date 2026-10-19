# E09_external_preferences

1024个5mer、同组成词序、25蛋白完整主群体

Paper: **Figure5a/b; S1.18; Data S23** (historical mapping; current binding: paper/MANUSCRIPT_BINDING.json).

Selected original analysis code, frozen results and run metadata are listed in INCLUDED_FILES.csv. Omitted inputs/weights/intermediates are in research/catalog/EXTERNAL_ASSETS.csv. SOURCE_MANIFEST.csv records original and packaged identities.

## Evidence boundary

主终点全词Spearman，次终点同组成一致率。Full超过frequency，但不普遍超过Sequence；Full/Sequence训练路径未完全匹配。

Historical launchers preserve provenance. Account paths are redacted to /data/CARBON_PROJECT or /data/USER_HOME. Set real paths in a new scratch copy before running them. The checked portable entry points are scripts/verify_frozen_results.py and the AARS CPU demo. Exact historical source snapshots, fresh training and end-to-end biological replication are not certified.

## Original measured-source matrix

The pinned original RNAcompete archive and author metadata are included at [paper S23 source files](../../../paper/source_data/S23_external_order_preferences/original/sources/). Use python scripts/extract_rnacompete_source.py --verify-only from repository root to check the compressed archive and expanded original-matrix identities. The 115 MB expanded TXT stays outside Git.
