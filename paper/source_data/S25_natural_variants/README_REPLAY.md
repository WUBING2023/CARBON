# Data S25: exposed natural-SNV and conservation analyses

Copy this directory's scripts, CSV/NPZ/JSON files and the complete `source_inputs`, `b1_reference`, `annotation_ncbi`, `annotation_raw`, and `annotation_clinvar` directories together. Paths resolve relative to each script. No original M: drive, server, model weights, GPU, or internet connection is needed to recompute these statistical endpoints. The historical source paths remain in provenance documents. Frozen inputs and original outputs must be retained when moving the package.

The main manuscript addition is `figures/FIG5_ADDITIONAL_d_e_85mm.pdf` (85 × 91 mm), with native SVG and 600-dpi PNG. Keep the existing Figure 5 a/b/c unchanged at 85 × 132 mm; the corresponding byte-exact originals are in `b1_reference/figure.*`. The 85 × 223 mm `FIG5_EXPANDED_85mm.*` is an engineering preview, not a manuscript asset. S22 and S23 are named explicitly in `figures/`.

Use Python 3.11 or newer with NumPy, pandas, SciPy, Matplotlib and PyMuPDF. Execution uses at most four BLAS/OMP threads. An exact package version snapshot is recorded in `DELIVERY_MANIFEST.json`. Arial is used for the delivered figures; on systems without Arial, Matplotlib may substitute an available font. Statistical replay does not depend on fonts.

From the copied directory:

```powershell
python replay_offline.py --metrics-only
python replay_offline.py
```

The first command recomputes scientific tables, cached allele annotations and source audits. The second also regenerates all plots. It reads archived official JSON responses, never invokes the download scripts, and leaves the frozen protocol and acquisition request logs intact. Raw reference conservation and preferred-allele labels are inside `source_inputs/EVENTS_INCLUDED.json.gz`. Original 101-nt allele tensors and predictions are in `INPUTS.npz` and `PREDICTIONS.npz`. `PACKAGED_INPUT_MANIFEST.json` records hashes; the gzipped JSON decompresses to the exact original bytes. No new inference was performed. This package reproduces the analysis of the existing predictions, not their training or GPU inference.

The ordered replay is: `analyze_variants.py`; `verify_weighted_refits.py`; `parse_annotations.py`; `parse_clinvar_cache.py`; `finalize_annotations.py`; `verify_region_sensitivity.py`; `audit_sources.py`; then `render_variants.py`, `render_supplement.py`, `export_continued.py`. Annotation parsing always regenerates the base table before enrichment, avoiding duplicate merge columns. Historical online acquisition scripts are preserved for provenance; do not run them to reproduce the frozen snapshot.

`validate_portable_replay.py` creates an isolated new subdirectory and repeats the full offline replay using only copied package inputs, then compares eleven output tables to the delivered results. `PORTABLE_REPLAY_QA.json` records the completed validation. Its scratch `replay_check_*` directory is not part of the manuscript's Data S25 package.

Scientific scope: all 7,262 ASB-positive events and the historical 396-event/287-locus/85-RBP chromosome subset have already exposed labels and predictions. The latter is not a new blind validation set. Direction is a measured BEAPR preferred-allele label; absolute predicted log-count difference is a computational response magnitude. Clinical significance and population frequency are annotation context only. The source lacks ASB negatives and an independent quantitative binding-effect endpoint. Fixed-input conservation replacements are not retrained models, and old Full/Sequence training is not guaranteed matched. See `INSERTIONS.md` and `OFFICIAL_SOURCES_AND_LIMITS.md` for the complete limits and interval methods.

The local frozen protocol specifies 1e−7 tie credit, whereas the previous S1.14 text used 1e−6. `TIE_THRESHOLD_SUMMARY.csv` verifies that all seven original/derived prediction arms in **this frozen 7,262-event input package** have zero values at either threshold in both complete and strict cohorts. Changing the threshold within this package changes no event, direction credit, or accuracy. The complete (empty) changed-event table is retained. This does not assert that the historical, separately archived S1.14 396-event prediction arrays are identical; their original one-tie report remains valid for those arrays.
