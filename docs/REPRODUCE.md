# Verification and replay

## Checked entry points

1. Install the CPU dependencies. Run python scripts/smoke_model.py for the output contract.
2. Run python scripts/verify_demo_cpu.py. This reads the included checkpoint and test H5, writes to outputs/demo_cpu and compares six saved metrics. Override the destination with --output PATH.
3. Run python scripts/verify_frozen_results.py. This reads frozen tables and writes its numerical report to outputs/frozen_results. It checks all 1,110 paired ablation rows against the 1,332 saved per-run metrics, canonical checkpoint counts, the Figure 4 target/share/loss/gene/cell calculations, all Figure 5 Full/Sequence/frequency RNAcompete 5-mer metrics, and all 13 disease-subset events.
4. Run python scripts/verify_package.py for included file hashes and duplicate targets.

These checks do not train a model, regenerate upstream predictions, reprocess BAM data, rerun bootstrap confidence intervals, or establish a unified fair comparator benchmark. Numerical and hash agreement have a narrower meaning than scientific validity or independent biological replication.

## Restore duplicate paths into scratch

    python scripts/materialize_replay.py --destination ../carbon-replay

The destination must not already exist. This copies the candidate and restores each byte-identical duplicate path from research/catalog/ASSET_ALIASES.csv, checking its SHA256. It never writes to the frozen source repository. Historical scripts may then run in that scratch tree after replacing generic /data/CARBON_PROJECT and /data/USER_HOME anchors and supplying omitted inputs/weights from the external-asset manifests.

Current Figure 2 uses paper/editable_figures/updated_fig2_v33/render_fig2_v33.py. Figure 1/4 source entry points are under paper/editable_figures/updated_fig*_v32; current Figure 5 uses paper/editable_figures/updated_fig5_v34/render_v34.py after materializing its input aliases into scratch; Figure 3 uses accepted_fig3_circlize_v31/render_circos.R. Figure 3 requires R/circlize/ggseqlogo; dependency tarballs and private reference code are not redistributed here. Review docs/EVIDENCE_MAP.md and paper/FIGURE_BINDINGS.csv before replay. Figure scripts can overwrite their own scratch outputs.

Upstream training requires the original experiment-specific namespace, seed, chromosome split, loss, early-stopping and checkpoint-selection rule. A generic default or the small AARS demo must not stand in for those settings. Early 2026-03 per-run source SHA identities were not saved; current archived source associated with launchers is not certified as every historical run's exact source snapshot.

## Original RNAcompete matrix

Run python scripts/extract_rnacompete_source.py --verify-only to verify the included pinned tar and its 115,308,739-byte original measured matrix without placing the expanded file in Git. To extract into scratch, use --destination ../carbon-replay/Zscores_420_origrncmpt.txt. This checks source bytes; it does not rerun the 7-mer-to-5-mer mapping.
