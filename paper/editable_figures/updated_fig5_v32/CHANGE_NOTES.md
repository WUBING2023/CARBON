# Figure 5 — v32 visual revision

Accepted scientific baseline: `revision_v31/manuscript/figures/fig5.pdf`, byte-identical to archived `accepted_fig5_v29/fig5_v29.pdf`. All edits are confined to this work directory. No main manuscript, CURRENT_REVISION, original accepted file or online project was changed.

## Changes addressing the feedback

- Matched Figure 3's Arial, black axis/condition labels, plain entity labels, 9 pt bold lowercase panel letters, white canvas and 0.55 pt axis lines. NPG colours keep their original meanings: Full `#E64B35`, Sequence `#4DBBD5`, Frequency `#3C5488`, flat `#00A087`, permuted `#8491B4`.
- Replaced repeated coloured text and heavy panel titles with a shared neutral label hierarchy. Frequency/Sequence names are expanded where space permits. Ordinary figure labels are at least 5.5 pt at the 178 mm export width, matching the scale used by Figure 3's ordinary labels; mathematical subscripts remain smaller.
- Removed panel a's separate large red numeric summaries. The plotted mean/95% interval and every paired protein point remain; exact summaries move to `CAPTION_V32.tex`.
- Retained the QKI lifted $Z_5$ network as the characteristic visual element: 30 nodes, 48 actual edges, four labelled words, split Full/Sequence nodes and the unchanged shared signed response scale. Removed its repeated node/edge count from the artwork; the caption retains it.
- Kept the c/d/e data geometry, shifted the complete lower section upward by 6 mm, and reduced the canvas from 178 × 184 to 178 × 178 mm. Five main panels a–e remain. No separator lines were added.
- Kept both cellular endpoints in c visibly distinct, with seven saved model/control responses and all six published reporter measurements. FDR and the post hoc qualifier are retained in the caption; the training-locus label remains on the figure.

## Evidence and boundaries retained

`render_v32.py` compares every `(panel, kind, id, arm, quantity, value)` with the preserved v29 plotted-value table and asserts exact equality: **2028 / 2028 records**. This is a replot of frozen results, with no training, prediction or new statistical analysis.

- a: all 25 proteins; negative paired differences remain; complete Sequence comparisons remain in the existing supplement. This figure does not claim a general Full-over-Sequence advantage.
- b: 30 QKI words / 48 adjacent exchanges; the original RNAcompete and response values/scales are unchanged.
- c: exact ACVR1B locus, fixed reference conservation, post hoc training exposure, all controls, all measured reporter values and distinct binding/RNA endpoints.
- d: 396 development-exposed ASB events / 85 RBPs / 287 sites; original point coordinates, four quadrant counts (168 / 61 / 43 / 124), five arm estimates and intervals.
- e: all 8 sites / 13 events / 52 effects; Full 7/13, Sequence 10/13, flat 12/13, permuted 9/13, REF rule 13/13. Negative bars are open and hatched. All source labels favour REF. Eleven training and two validation events remain explicit. Disease annotation is exploratory and does not establish clinical or causal validity. Conservation remains the fixed reference input for both alleles.

## Viewed style references and scope

- Current study Figure 3 (`baseline/fig3.pdf`, 85 mm placed width): primary typography, plain black labels, panel lettering, open axes, compact unboxed legends and whitespace. Current Figure 1 (`baseline/fig1.pdf`) was also inspected for the paper's shared colour/line vocabulary, while respecting its different architecture role.
- A local style reference was inspected for consistent lettering, neutral axis labels and restrained figure-wide hierarchy. No scientific statements, results, values or biological entities were transferred. This is a local style correction of accepted CARBON evidence.

The largest remaining difference is justified by evidence density: Figure 5 combines a full external cohort, a lifted word graph, a paired locus illustration, a full ASB cohort and an aligned 13-event table. Its double-column width supports this density. Panel e is still the densest region; its complete unfavourable comparison is retained rather than deleted for appearance.

## Deliverables and verification

- `fig5_v32.pdf` / `.svg` / `.png`: final vector figure and 400 dpi PNG.
- `render_v32.py` + `data.py` + `frozen_inputs/`: active editable producer and task-relative frozen inputs. Run `python render_v32.py` from any working directory. The archived intermediate scripts and `build_revision_source.py` document how the accepted plotting code was restyled; they are not required for a normal replot.
- `PLOTTED_VALUES_V32.csv`, `SOURCE_MAPPING_V32.json`, `INPUT_BINDINGS_V32.json`: numerical equality, cohort/scale bindings and local source hashes.
- `CAPTION_V32.tex`: complete proposed caption; the main Agent must reconcile supplement numbering during integration. It retains a–e mapping and adds the numerical summaries and FDR moved from the artwork.
- `comparison_v29_v32.pdf` / `.png`: complete old/new comparison at the same 178 mm width, with the previous and revised study figure clearly labelled.
- `style_reference_fig3_v32.pdf` / `.png`: current Figure 3 and revised Figure 5 at their actual 85/178 mm widths.

Self-check performed on the actual candidate PNG and on a fresh rasterization of its PDF, then on both complete comparisons. All five panel letters, the four effect columns and every unfavourable bar are visible. There are no duplicate overall legends, off-canvas text, missing points or added dashed separators. PDF page dimensions are 178 × 178 mm; ordinary Arial text has a measured 5.5 pt minimum. Independent review and main-manuscript placement remain the main Agent's acceptance step.
