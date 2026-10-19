# CARBON

**Conservation-aware prediction of base-resolution RNA-binding protein profiles.**

CARBON combines RNA sequence and evolutionary conservation with multi-scale convolutions, channel attention and a bidirectional GRU. It predicts a binding-profile shape and a total signal count for each input window. This research repository contains model and training code, an executable AARS CPU demonstration, frozen experiment results, and the associated manuscript draft.

[Main manuscript](paper/CARBON_main.pdf) · [Supplement](paper/CARBON_supplement.pdf) · [All 27 source-data entries](paper/source_data/README_CURRENT.md) · [Experiment and figure map](docs/EVIDENCE_MAP.md)

## Input and output

The default model consumes a float tensor of shape `(batch, 101, 5)`: four one-hot sequence channels in A/C/G/U order, followed by the aligned PhyloP conservation channel. `T` is encoded as `U`. The HDF5 demonstration stores inputs as `X`, pooled eCLIP signal as `Y`, and input-control signal as `Y_ctl`. Real conservation values require a matching genome assembly and track; zero-filled values are a missing-input convention.

The Treat model returns `logits_profile` of shape `(batch, 101)` and `log_total_count` of shape `(batch, 1)`. The predicted signal is:

```python
prediction = torch.softmax(output["logits_profile"], dim=-1) * torch.exp(output["log_total_count"])
```

See [the model](CARBON/models/carbon_model.py), [5-channel preprocessing](preprocessing/preprocess_5ch.py), and [the smoke check](scripts/smoke_model.py) for the executable contract. The accepted **2026-10-05 Figure 1 is locked**: preserve its exact PDF, editable PowerPoint and 85 mm manuscript placement. [Figure lock](paper/FIGURE_LOCKS.json) · [validation](paper/figure_locks/validate_figure1.py).

## Installation

```bash
git clone https://github.com/WUBING2023/CARBON.git
cd CARBON
conda env create -f environment.portable.yml
conda activate carbon
```

Alternatively, in a Python 3.12 environment:

```bash
python -m pip install -r requirements.cpu.txt
```

[requirements.cpu.txt](requirements.cpu.txt) records versions used for the checked Windows CPU run, including PyTorch 2.4.1. A clean environment installation has not been validated. [environment.yml](environment.yml) is a historical Linux export. Preprocessing and attribution have additional dependencies in [requirements.txt](requirements.txt); some figure workflows require PyMuPDF, R and circlize.

## CPU demonstration

Run from the repository root; these checks require no external data download or GPU:

```bash
python scripts/smoke_model.py
python scripts/verify_demo_cpu.py
python figures/plot_demo_signal_track.py --device cpu --out outputs/demo_signal_track.png
```

The verifier evaluates the included checkpoint on **1,283 AARS test windows** and compares six saved metrics at absolute tolerance `1e-5`. It writes `verification.json` and an evaluation log under `outputs/demo_cpu/`. The plot writes the requested PNG. Its four tracks are CARBON prediction, pooled eCLIP, inverted pooled eCLIP and inverted input control; the inverted tracks are display transforms rather than additional replicates.

## Training and ablations

The included small chr1/chr2/chr3 split supports a short training example:

```bash
bash scripts/run_demo_train.sh outputs/demo_train
python carbon_train.py --help
```

The launcher defaults to three epochs on CPU. Use a new output directory. [Training instructions](docs/TRAINING.md) give the equivalent Python command and the implemented ablation flags. [E01](research/experiments/E01_data_model/README.md) describes preprocessing and model lineage; [E03](research/experiments/E03_ablation/README.md) links the complete 222-task, six-condition archive and recorded configurations.

The demo split is illustrative. Formal experiments require their recorded chromosome splits, seeds, filtering, losses and checkpoint-selection rules. Historical launchers use the MARCP namespace preserved in [research/historical_model](research/historical_model/). A short demo run does not reproduce the full paper.

## Find an experiment, dataset or figure

Follow the [evidence map](docs/EVIDENCE_MAP.md) from a manuscript result to its experiment guide, source data and analysis code. The [S1–S27 index](paper/source_data/README_CURRENT.md) links every primary data entry directly to an included file, including full result cohorts and unfavorable directions. [Figure bindings](paper/FIGURE_BINDINGS.csv) identify the current PDFs and available plotting sources.

| Path | Contents |
|---|---|
| [CARBON/](CARBON/), [carbon_train.py](carbon_train.py), [preprocessing/](preprocessing/) | Model, losses, prediction, training/evaluation and input preparation |
| [data/](data/), [carbon_trained_model/](carbon_trained_model/), [configs/](configs/) | Included AARS HDF5 demonstration and checkpoint |
| [research/experiments/](research/experiments/) | E01–E11 analysis code, frozen results, run records and guides |
| [research/catalog/](research/catalog/) | Per-run ablation index, checkpoint availability, provenance and asset maps |
| [paper/](paper/) | Manuscript TeX/PDF, bibliography, figures, source tables and editable inputs |
| [scripts/](scripts/) | CPU, frozen-result and package verification; scratch replay preparation |
| [docs/](docs/) | Training, experiment navigation, reproduction scope and resource availability |

Verify the included frozen numbers and file bindings:

```bash
python scripts/verify_frozen_results.py
python scripts/verify_package.py
python paper/figure_locks/validate_figure1.py
```

[Reproduction instructions](docs/REPRODUCE.md) explain the checked scope and how to restore duplicate paths into a separate scratch directory. [Known limits](docs/KNOWN_LIMITS.md) document comparator conventions, incomplete training provenance, resource gaps, and negative or mixed functional results.

## Data and weights

The AARS demo checkpoint and test data are included. Paper source tables, selected measured-score records, frozen arrays and analysis code are indexed here. Raw sequencing, reference genomes, large database distributions and full benchmark checkpoint files are external.

A versioned public **full-data/full-weight archive is not currently available**. The canonical 222-task catalog records 221 CARBON, 219 iDeepB and 222 RBPNet checkpoint identities, including explicit missing tasks; these counts describe archived availability, not files bundled here. See [availability](docs/DATA_AVAILABILITY.md), [weight catalog](research/catalog/BENCHMARK_WEIGHT_AVAILABILITY.csv), and [external-asset inventory](research/catalog/EXTERNAL_ASSETS.csv). No pending DOI or private platform link is presented as a public download.

## Citation

Use [CITATION.cff](CITATION.cff) to cite the software, identifying the commit you used. The [associated manuscript](paper/CARBON_main.pdf) is a draft; no journal publication or DOI is claimed here. Cite external assay and data sources described in the manuscript and dataset provenance when using their data.

## License

Original CARBON software is provided under the [MIT license](LICENSE). Adapted code, the archived MEME/Tomtom binary, external datasets and papers retain their own terms; MIT does not relicense them. [Third-party notices](THIRD_PARTY_NOTICES.md) contain the component mapping and full license copies. MEME Suite has educational, research and non-profit use terms, with separate arrangements for commercial use.
