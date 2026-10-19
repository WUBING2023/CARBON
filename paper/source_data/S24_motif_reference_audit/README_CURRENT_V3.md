# Data S24: current v3 entry point

Two separate reference analyses are retained here. Their target collections and
q-value families differ; do not substitute one analysis's q values for the other.

## Current Figure 3b: official ATtRACT comparison

Open `ATtRACT_DATA_S24_V3.zip`. This archive contains the official reference
snapshot, 485 full query matrices, 1,583 target matrices, the complete 767,755
Tomtom comparisons, run provenance, complete candidate tables, post hoc display
selection, and the matrices/coordinates bound to the current paired logos.

The current display retains positive pattern 0 and human, unmutated, same-RBP
references of at most 10 nt. After the complete search, the illustrative subset
requires q < 0.01 and mean query information content at least 0.3 bits in both
matched and displayed intervals. Five best-per-RBP examples prioritize q; the
sixth maximizes displayed query information among the remaining best-per-RBP
candidates. This yields PTBP1, PCBP2, RBM5, HNRNPK, U2AF2 and PUM2. The full
candidate set, including the low-information PCBP1 match, remains available.

These are favorable post hoc illustrations, not a success-rate estimate or a
blind validation set. Logos share the same 0–2-bit height scale. Reference sources
are heterogeneous and are not described as uniformly measured gold standards.
The complete-query q values are unchanged by display filtering.

## Archived CISBP analysis and earlier artwork

All other existing data files in this directory belong to the retained CISBP
analysis or its earlier drawing implementation. In particular,
`MAIN_SELECTION_RULE.json`, `INSERTIONS.md`, `DELIVERY.md`, `FIGURE_QA.json` and
`draw_fig3.py` describe archived selections/artwork, not the current Figure 3b.
The original descriptive README is preserved as `README_ARCHIVED_CISBP.md`.

The archived tables contain all 193 original cognate comparisons, with their own
reference metadata and original statistics. The complete paired-logo booklet and
stored matrices are retained. They are not recomputed or relabeled as ATtRACT
results. The original panel c examples retain their original CISBP statistics.

## Current artwork and named atlas

Current production code lives in `revision_v3/work_fig3/`; source bindings are in
`manuscript/editable_figures/README_V3.md`. The main manuscript uses 85-mm
`fig3_a.pdf` and `fig3_bc.pdf` as a continued figure. All 219 names remain in the
atlas, with separate 178/285-mm vector enlargements and Figure S21. The current
base palette is A green, C deep blue, G orange, U/T red. Numbered 1–6 markers and
the old circle geometry described in the archived README are not current.

For exact commands, input hashes, provenance and selection details, use the
archive README and the current supplementary section S1.19. This package does
not imply that all 219 atlas leaves have significant reference matches and does
not supply missing original IG/inference provenance.
