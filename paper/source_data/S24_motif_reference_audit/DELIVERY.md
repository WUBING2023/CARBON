# Figure 3 final delivery — 12 October 2026

The final `fig3.pdf` was copied into `../manuscript/figures/fig3.pdf`. SHA-256: `02910f8c9f57612085b27ab4b3f00656c367f7db846ebe2fbb110c3d7dda42f9`.

The main export is 85 × 199 mm. All ordinary text extracted from the actual PDF is at least 6.5 pt. The final rendered main page was visually opened after the last font increase. Its paired-logo names, motif IDs and p/q values have no clipping or overlap; c retains separated IC/IG ticks. The accepted source evidence is retained: 219 atlas leaves, frozen crop/leaf order, six annotation tracks, stored query/reference matrices and alignments, unchanged p/E/q values, and line/glyph geometry from the old c examples. The central-tree nominal outer radius is 11.5 mm; its center is 3.2 mm left of the canvas midpoint.

## Files for integration

- Main Figure 3: `fig3.pdf` (already installed), `fig3.svg`, `fig3.png`.
- Figure S21: `fig3_supplement_atlas_readable.pdf`, five pages. All 219 labels are present in the full atlas, and each has a complete visible copy in at least one enlarged quadrant. The quadrant label size is approximately 7.84 pt. Boundary fragments are intentional clipped duplicates and do not remove any complete label.
- Data S24 table: `motif_comparisons_all_193.csv`, plus `CISBP_exact_metadata_records.csv`, `paired_logo_matrices.npz`, `MAIN_SELECTION_RULE.json`, `SOURCE_AUDIT.json` and `SOURCE_INPUT_MANIFEST.json`.
- Data S24 independent PDF: `fig3_supplement_all_193_pairs.pdf` (20 pages). The 13-page direct-only inspection booklet is redundant with this complete booklet and need not be appended to the manuscript.
- English integration text: `INSERTIONS.md`. Source-aware audit/display text is labelled S1.19; panel c scaling/provenance still points to S1.7.
- Database URL/version/acquisition record: `CISBP_OFFICIAL_PROVENANCE.md`. Official homepage, bulk and help pages were checked on 12 October 2026. The retained local acquisition record dates the Build 2.00 download to 27 March 2026; this revision did not redownload the database.
- Code and source replay: `README.md`, `draw_fig3.py`, `prepare_motif_sources.py`, `extract_retained_c_vectors.py`, `package_sources.py`, `package_and_qa.py`, cached matrices/geometry and `source_inputs/`. The source copies cover all 78 reference matrices and all 45 task query/Tomtom outputs. There are 177 retained original files, totalling 2,125,243 bytes, each hash-checked.
- Comparison/QA artifacts: `FIGURE_COMPARISONS.pdf`, `FIGURE_QA.json` and `qa_fig3_1.png`.

## Verification and scientific limits

The audit was rerun using only the copied original source files, yielding a bitwise-identical complete CSV. All original S6 p/E/q values match their retained Tomtom rows. Counts remain 193 comparisons / 80 task-pattern pairs / 45 tasks / 31 symbols, with 72 q < 0.05. Direct experimental references comprise 126 comparisons; 45 pass the original threshold. The main six examples include nonsignificant HNRNPA1 and RBFOX2 comparisons. The source audit explicitly separates JPLE-inferred, homology-inferred and database-compiled motifs; it does not call these experimental gold standards.

All pages of the paired-logo booklets and the main/sensitivity PDFs were checked for out-of-page text; none was found. The full atlas has all 219 labels, and the enlarged quadrants jointly contain all 219 complete labels. The main, original, published exemplar, full atlas, selected booklet pages and enlarged sensitivity export were actually viewed. Figure regeneration does not rerun Tomtom, inference or attribution.

Panel c remains a recomposition of accepted drawing geometry. Full per-window checkpoint/IG/substitution provenance remains unresolved, as specified in `INSERTIONS.md`; no mechanistic validation is added. The main atlas retains every leaf but only six numbered representatives; Figure S21 and Data S10 preserve the complete names. No manuscript TeX or other main figure was edited by this worker.
