# Third-party notices

Original CARBON software is licensed under [MIT](LICENSE). This does not relicense third-party code, the archived MEME binary, external datasets, source measurement records or papers. Each component retains the terms of its source. Full copies of the code licenses below are included in [licenses/](licenses/); the upstream sources were checked on 2026-10-19.

| Component | Use in CARBON | License and preserved text |
|---|---|---|
| [DeepExplain](https://github.com/marcoancona/DeepExplain) | Attribution inspiration in `CARBON/interpretability/attribution/carbon_attribution.py` and historical MARCP copies | MIT, Copyright (c) 2017 ETH Zurich (Marco Ancona); [full text](licenses/DeepExplain-MIT.txt) |
| [Captum](https://github.com/meta-pytorch/captum) | Adapted DeepLift implementation in `CARBON/interpretability/captum/captum_deeplift_hypothetical.py` and historical copies | BSD 3-Clause, Copyright (c) 2019 PyTorch team; [full text](licenses/Captum-BSD-3-Clause.txt) |
| [TF-MoDISco](https://github.com/kundajelab/tfmodisco) | Adapted sequence visualization in `CARBON/visualization/viz_sequence.py` and historical copies | MIT, Copyright (c) 2018 Kundaje Lab at the cited historical implementation; [full text](licenses/TF-MoDISco-MIT.txt) |
| [MEME Suite / Tomtom 5.5.7](https://meme-suite.org/meme/) | Historical Linux binary inside `paper/source_data/S24_motif_reference_audit/ATtRACT_DATA_S24_V3.zip`, member `tools/meme_bioconda/bin/tomtom` | MEME educational, research and non-profit terms, **not MIT**; [complete same-version COPYING](licenses/MEME-5.5.7-COPYING.txt), also included beside the binary in the ZIP |

## License sources

- DeepExplain: [official LICENSE.md](https://github.com/marcoancona/DeepExplain/blob/master/LICENSE.md).
- Captum: [official LICENSE](https://github.com/meta-pytorch/captum/blob/master/LICENSE). The BSD notice, all three conditions and disclaimer are retained in full.
- TF-MoDISco: [LICENSE at the historical implementation commit](https://github.com/kundajelab/tfmodisco/blob/178db48d9848d4d841ffda8c5d2f0269afa9f78e/LICENSE), associated with the visualization source cited in its header.
- MEME 5.5.7: COPYING from the original same-version source distribution used for the historical analysis; [official release](https://meme-suite.org/meme/meme-software/5.5.7/meme-5.5.7.tar.gz), [official copyright page](https://meme-suite.org/meme/doc/copyright.html). The archived binary is Bioconda `meme 5.5.7 pl5321h1ca524f_3`, SHA256 `998eb82ea2154bd924953fa7f13498b071c21ce16f9c634a8d58ffb9f2e01873`. The package's original `TOMTOM_RUN_PROVENANCE.json` records its source-package URL and run commands.

MEME permits copying and distribution for educational, research and non-profit purposes subject to its complete notice and disclaimer being retained. Incorporation into commercial products or commercial use requires contacting UC San Diego as stated in COPYING. These terms also apply when extracting the archived binary. The license-only ZIP update adds notices while preserving every existing member, including all measurements and the binary. Historical Git objects remain historical; this update does not rewrite them.

## External data and publications

Assay datasets, motif databases, measured-score records, genomes, conservation tracks and published papers retain their original terms and citation requirements. Inclusion of a source record or derived table does not grant a new MIT license to the underlying source material. Follow the manuscript references, per-packet provenance, and [availability guide](docs/DATA_AVAILABILITY.md) to identify sources and obtain unbundled inputs. Author-generated research tables retain their documented provenance; no blanket data license is newly asserted here.
