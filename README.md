# CARBON

**Conservation-aware prediction of base-resolution RNA-binding protein profiles.**

CARBON combines RNA sequence and evolutionary conservation with multi-scale convolutions, channel attention and a bidirectional GRU. It predicts a binding-profile shape and a total signal count for each input window. This research repository contains model and training code, an executable AARS example that runs on a standard computer without a GPU, frozen experiment results, and the associated manuscript draft.

**Current manuscript: 20 November 2026.** *CARBON integrates RNA sequence and evolutionary conservation to predict nucleotide-resolution eCLIP profiles*. This update adds the author-reviewed metadata and declarations to the 19 November scientific revision.

[Manuscript (PDF, 16 pages)](paper/revisions/2026-11-20/CARBON_main.pdf) · [Supplement (PDF, 42 pages)](paper/revisions/2026-11-20/CARBON_supplement.pdf) · [Complete LaTeX source](paper/revisions/2026-11-20/CARBON_source.zip) · [Authors and compilation](paper/revisions/2026-11-20/README.md) · [Checksums](paper/revisions/2026-11-20/checksums.json)

[Supplementary Data S1–S53](paper/source_data/README_CURRENT.md) · [Experiment and figure map](docs/EVIDENCE_MAP.md) · [19 November revision](paper/revisions/2026-11-19/) · [October manuscript archive](paper/CARBON_main.pdf)

## Input and output

The default model consumes a float tensor of shape `(batch, 101, 5)`: four one-hot sequence channels in A/C/G/U order, followed by the aligned PhyloP conservation channel. `T` is encoded as `U`. The HDF5 demonstration stores inputs as `X`, pooled eCLIP signal as `Y`, and input-control signal as `Y_ctl`. Real conservation values require a matching genome assembly and track; zero-filled values are a missing-input convention.

The Treat model returns `logits_profile` of shape `(batch, 101)` and `log_total_count` of shape `(batch, 1)`. The predicted signal is:

```python
prediction = torch.softmax(output["logits_profile"], dim=-1) * torch.exp(output["log_total_count"])
```

See [the model](CARBON/models/carbon_model.py), [5-channel preprocessing](preprocessing/preprocess_5ch.py), and [the smoke check](scripts/smoke_model.py) for the executable contract.

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

## Quick start: prediction without a GPU

CPU means the computer’s central processing unit. This small example uses the included AARS checkpoint and data so you can verify installation and prediction on a standard computer without a GPU. Larger training and inference workloads can use GPU hardware; the example is an entry point to CARBON.

Run from the repository root:

```bash
python scripts/smoke_model.py
python scripts/verify_demo_cpu.py
python figures/plot_demo_signal_track.py --device cpu --out outputs/demo_signal_track.png
```

The verifier evaluates the included checkpoint on **1,283 AARS test windows** and compares six saved metrics at absolute tolerance `1e-5`. It writes `verification.json` and an evaluation log under `outputs/demo_cpu/`. The plot writes the requested PNG. Its four tracks are CARBON prediction, pooled eCLIP, inverted pooled eCLIP and inverted input control; the inverted tracks are display transforms rather than additional replicates.

## Reproduce reported results

The [benchmark source table](paper/source_data/S1_S14_original/Supplementary_Data_S1_benchmark_metrics.csv) records the main prediction comparison. The [evidence map](docs/EVIDENCE_MAP.md) links the model comparison, ablations and interpretation analyses to their archived inputs and scripts.

### Figure 2c: common-target transfer

The task-level source table, calculation script and expected results are included. This calculation uses only the Python standard library, with no model download or GPU required:

```bash
python scripts/verify_figure2c.py paper/source_data/S40_common_target_transfer/cross_model_task_metrics.csv
```

The comparison covers 68 RBPs, 136 transfer tasks and 1,421,243 target windows. Averaging the two directions within each RBP and then weighting RBPs equally gives PCC **0.46268 for CARBON, 0.40818 for iDeepB and 0.42047 for RBPNet**. See the [data, column definitions and expected output](paper/source_data/S40_common_target_transfer/README.md). This command recalculates point estimates from task-level sufficient statistics; model inference and confidence-interval resampling are separate analyses.

## Training and ablations

The included small chr1/chr2/chr3 split supports a short training example:

```bash
bash scripts/run_demo_train.sh outputs/demo_train
python carbon_train.py --help
```

The launcher defaults to three epochs on CPU. Use a new output directory. [Training instructions](docs/TRAINING.md) give the equivalent Python command and the implemented ablation flags. [E01](research/experiments/E01_data_model/README.md) describes preprocessing and model lineage; [E03](research/experiments/E03_ablation/README.md) links the complete 222-task, six-condition archive and recorded configurations.

The demo split is illustrative. Formal experiments require their recorded chromosome splits, seeds, filtering, losses and checkpoint-selection rules. Historical launchers use the MARCP namespace preserved in [research/historical_model](research/historical_model/). The full experiment protocols and results are linked in the evidence map.

## Find an experiment, dataset or figure

Follow the [evidence map](docs/EVIDENCE_MAP.md) from a manuscript result to its experiment guide, source data and analysis code. The [S1–S53 index](paper/source_data/README_CURRENT.md) follows the current supplementary material and distinguishes downloadable files from records awaiting public deposition. The current figure PDFs are bundled in the [LaTeX source](paper/revisions/2026-11-20/CARBON_source.zip). The [October figure bindings](paper/FIGURE_BINDINGS.csv) identify the archived plotting sources for that earlier revision.

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

[Reproduction instructions](docs/REPRODUCE.md) explain the checked scope and how to restore duplicate paths into a separate scratch directory. [Known limits](docs/KNOWN_LIMITS.md) document comparator conventions, incomplete training provenance, resource availability and interpretation of the archived analyses.

## Data and weights

The AARS demo checkpoint and test data are included. Paper source tables, selected measured-score records, frozen arrays and analysis code are indexed here. Raw sequencing, reference genomes, large database distributions and full benchmark checkpoint files are external.

A versioned public **full-data/full-weight archive is not currently available**. The canonical 222-task catalog records 221 CARBON, 219 iDeepB and 222 RBPNet checkpoint identities, including explicit missing tasks; these counts describe archived availability, not files bundled here. See [availability](docs/DATA_AVAILABILITY.md), [weight catalog](research/catalog/BENCHMARK_WEIGHT_AVAILABILITY.csv), and [external-asset inventory](research/catalog/EXTERNAL_ASSETS.csv).

## Citation

Use [CITATION.cff](CITATION.cff) to cite the software, identifying the commit you used. The [current manuscript](paper/revisions/2026-11-20/CARBON_main.pdf) and [supplement](paper/revisions/2026-11-20/CARBON_supplement.pdf) describe the November revision; the October draft remains archived. Cite external assay and data sources described in the manuscript and dataset provenance when using their data.

## License

Original CARBON software is provided under the [MIT license](LICENSE). Adapted code, the archived MEME/Tomtom binary, external datasets and papers retain their own terms; MIT does not relicense them. [Third-party notices](THIRD_PARTY_NOTICES.md) contain the component mapping and full license copies. MEME Suite has educational, research and non-profit use terms, with separate arrangements for commercial use.
