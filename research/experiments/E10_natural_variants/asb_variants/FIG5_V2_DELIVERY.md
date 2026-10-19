# Fig. 5 v2 delivery

Final manuscript target: `../manuscript/figures/fig5.pdf`. This is one 178 x 183 mm figure; no continued figure is required. It contains a/b/c left and d/e right, each story within an 85 mm column and an 8 mm gutter. No manuscript TeX, baseline, v1, or other agent's figure was edited.

## Artifacts

- `fig5_v2/fig5.pdf`: single native vector manuscript figure.
- `fig5_v2/fig5.svg`: native vector equivalent.
- `fig5_v2/fig5.png`: 600-dpi PDF render.
- `fig5_v2/fig5_178mm_preview.png`: actual-width overview at 150 dpi.
- `fig5_v2/phi_alignment_detail.png`: 300-dpi inset QA crop.
- `render_fig5_v2.py`: deterministic display-only renderer.
- `qa_fig5_v2.py`: geometry, source hash, complete coverage and unchanged-artwork checks.
- `FIG5_V2_SCIENTIFIC_MAPPING.md`: mapping recorded before rendering, revised exact unified caption and interval definitions.
- `fig5_v2/FIG5_V2_SOURCE_MAPPING.json`: exact frozen source hashes and displayed numerical results.
- `fig5_v2/FIG5_V2_QA.json`: machine-readable checks.

## Accepted a/b preservation

Original a/b were not reconstructed. Original `b1_reference/figure.pdf` and `.svg` are embedded at exactly unit scale, clipped only below y=96 mm to remove old c. QKI edits, native sequence, 32 paired backgrounds per example, independent assay endpoints, HNRNPA1 frequency/assay values, axes, letter size and response positions are unchanged. The sole a/b change is an opaque replacement inset at x=35..85 mm, y=5..20.1 mm. It replaces the unclear Phi wiggle with the identical 101 input values as vector heat cells, 49-53 outline, explicit fixed-input text and shared raw PhyloP scale. A 300-dpi raster comparison excludes this declared rectangle and finds **zero changed pixels among 1,032,156 compared pixels**.

Phi uses the NPG navy/white/red endpoints #3C5488 / #FFFFFF / #E64B35, a common raw range -7 to +3 and white at zero. Both plotted input tracks lie inside that range, so no values are clipped. In d, REF and ALT local sequence rows are aligned to the same 11-coordinate Phi strip; a full 101-nt overview retains central 49-53 positions. No variant-driven Phi change is implied.

## Numerical and statistical mapping

The 25 proteins share their original order on primary all-word Spearman and secondary within-composition accuracy axes. New c values are exact cells from `FROZEN_PRIMARY_METRICS.csv`; its accuracy values match the old `INTERVENTION_PACKET.npz` to 1e-10 absolute tolerance. All three arms remain.

Natural example remains chr8:1,844,440 C>T, K562, ARHGEF10/KBTBD11-OT1; all ten alphabetical RBP events are retained. All BEAPR preferred identities are REF C; Full/Sequence stored differences remain unchanged. No clinical/population label is plotted as binding truth.

e retains 396 events, 287 variants, 85 RBPs and 119 RNA-region clusters. Quartile means and 95% intervals come unchanged from `DIRECTION_BY_CONSERVATION.csv` (5,000 paired RBP-bootstrap draws). Aggregate means come from `DIRECTION_SUMMARY.csv`; all five aggregate 95% intervals are from `DIRECTION_TWO_WAY_SENSITIVITY.csv` (2,000 two-way RNA-region/RBP positive-exponential multiplier draws). Full=52.04%, Sequence=58.77%, Flat Phi=53.46%, Permuted Phi=52.36%, REF=64.77%. The paired Full-minus-Sequence mean is -6.734139 pp from `DIRECTION_CONTRASTS.csv`; its **displayed** two-way CI is -15.970237 to +1.641108 pp. The older RBP-only CI remains unmodified in its CSV. Figure and caption distinguish the two interval methods. All quartet and aggregate intervals fit fully within their axes.

The first rendered preview exposed an e-header/legend collision, which was corrected before final rendering. Exposure/implementation sentences were moved out of artwork into the exact caption. Actual-width overview and Phi details were inspected after the final update. Minimum 5.5-pt labels are confined to compact input-coordinate/scale annotations; original a/b text remains 6.5 pt and main new data text is 6.0-7.1 pt. No training, inference, or statistics were newly run.
