# Revised Figure 3 and motif source package

The main vector PDF is **85 × 199 mm** and has a verified minimum ordinary-text size of **6.5 pt**. Insert it at `width=85mm`; wider or narrower scaling changes the actual printed type size. The NPG base palette is A navy, C green, G brown and U vermilion. The main circle is shifted 3.2 mm left of the canvas midpoint and its central tree has an explicit nominal outer radius of 11.5 mm.

## Deliverables

- `fig3.pdf`, `fig3.svg`, `fig3.png`: revised main a/b/c. The 219-leaf overview retains the frozen crops, order and annotations. Six numbered representatives link to aligned paired logos. The accepted sensitivity examples are retained as vector geometry.
- `fig3_supplement_atlas_readable.pdf`: Figure S21, five pages comprising the full-name 219-leaf atlas and four enlarged quadrant views. The overview page is 285 mm wide; the quadrants are 210 × 236 mm. Preserve their intended physical dimensions when incorporating them into a larger supplementary document.
- `fig3_supplement_all_193_pairs.pdf`: complete 20-page paired-logo booklet for an independent PDF within Data S24. It includes direct, inferred, compiled and nonsignificant comparisons. It uses full input extents and the retained alignment coordinates.
- `fig3_supplement_direct_pairs.pdf`: 13-page inspection subset of the 126 direct experimental comparisons. It need not be appended alongside the complete booklet, which already contains those comparisons.
- `fig3_supplement_sensitivity.pdf`: enlarged view of the three retained sensitivity examples for source inspection; supplementary numbering is left to the manuscript owner.
- `motif_comparisons_all_193.csv`, `CISBP_exact_metadata_records.csv`, `paired_logo_matrices.npz`, `MAIN_SELECTION_RULE.json`, `SOURCE_AUDIT.json`: complete motif audit, matrices and frozen display rule.
- `INSERTIONS.md`: proposed English Methods, Results, main caption, S1.19, Figure S21 caption and Data S24 description. S1.7 retains the original statistical/provenance discussion and is the main-panel c scaling reference.
- `CISBP_OFFICIAL_PROVENANCE.md`, `source_inputs/`, `SOURCE_INPUT_MANIFEST.json`: official source links, the retained local acquisition record and hashed source files.
- `FIGURE_COMPARISONS.pdf`, `FIGURE_QA.json`, `qa_*.png`: reference/old/new comparisons and checks on the actual exported PDFs. The published exemplar is Sun et al., Cell Research (2021), Figure 5; it informs the paired-logo presentation and circular atlas layout.

## Scientific boundaries

The 193 retained comparisons cover 80 task-pattern pairs, 45 tasks and 31 RBP symbols. Their original p/E/q values are unchanged; 72 pass q < 0.05. Source classification yields 126 direct experimental, 54 JPLE-inferred, 11 database-compiled and two homology-inferred comparisons. The direct experimental subset has 45 passing comparisons. Main examples are fixed named representatives selected lexicographically within each RBP, not chosen by scores. HNRNPA1 and RBFOX2 are nonsignificant examples. The 219 pattern-0 atlas entries are not all experimentally validated.

Panel b applies the stored Tomtom offset/orientation and uses identical nucleotide spacing between each query/reference pair. Its display window is centered on the aligned reference with width max(9, reference width + 2). The full PPMs are preserved, and the significance values refer to the full inputs. HNRNPL genuinely has low information content in the aligned query; the display does not rescale it to appear stronger.

Panel c recomposes the accepted PDF glyphs and line coordinates through linear mappings from the old axis ticks. It does not rerun inference or attribution. eCLIP and wild type retain independent min–max scaling; the mutant uses the wild-type minimum and range. `retained_c_signal_curves_from_pdf.csv` contains extracted drawing geometry, not raw inference arrays. The saved examples still lack full per-window checkpoint, IG target/baseline/step and substitution provenance; the revision does not invent these details or establish a binding mechanism.

## Reproduction

Use the retained matrices/geometry to rebuild all vector outputs without network requests, Tomtom or model inference:

```powershell
& 'D:\Anaconda\envs\rbp\python.exe' -X utf8 .\draw_fig3.py
& 'D:\Anaconda\python.exe' -X utf8 .\package_and_qa.py
```

The two environments are used because the base environment supplies PyMuPDF and the RBP environment supplies compatible h5py/scipy/matplotlib. `draw_fig3.py --main-only` rebuilds the main figure and enlarged sensitivity view. Cached atlas matrices, the frozen S10 mapping/distance matrix and retained c geometry make figure regeneration independent of the external original H5 files.

`prepare_motif_sources.py` rebuilds the exact metadata/PPM audit from the retained source-input copies when available, or from the original bound paths. `extract_retained_c_vectors.py` extracts retained geometry from `fig3_baseline.pdf`. These scripts recover stored results and drawing coordinates; no significance family or model output is recalculated.
