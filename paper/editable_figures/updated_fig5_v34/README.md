# Current Figure 5 source

Entry: render_v34.py. The current manuscript output is ../../figures/fig5.pdf; fig5_v34.pdf is a byte-identical restorable alias of that canonical PDF. PNG, SVG, caption, plotted values and mapping are supplied here.

All frozen_inputs paths are SHA-checked aliases of existing committed assets. From the repository root run `python scripts/materialize_replay.py --destination ../carbon-v34-replay`, then run `python paper/editable_figures/updated_fig5_v34/render_v34.py` in that scratch tree. The renderer will replace its scratch exports and mapping. It does not rerun model inference, training or source measurements.

PLOTTED_VALUES_V32.csv and PLOTTED_VALUES_V34.csv contain the same 2028 records and bytes. SOURCE_MAPPING_V34.json records input/output hashes and coordinate checks. VISUAL_REVIEW_V34.md is a local rendering review, not user approval. See ../../../docs/FIG5_V34_CHANGES.md for scope and limitations.
