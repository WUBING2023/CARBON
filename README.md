# CARBON

CARBON is a conservation-aware deep-learning framework for predicting
base-resolution RNA-binding protein (RBP) binding profiles. It combines
nucleotide one-hot encoding with an evolutionary-conservation channel and
predicts both the profile shape and total count.

This repository is the **public beta code snapshot** accompanying the
manuscript *CARBON: a conservation-aware framework for base-resolution RBP
binding profile prediction*, prepared for submission to **Bioinformatics**.
The manuscript has not yet been peer reviewed.

## Release scope

This beta contains the source code, preprocessing utilities, configuration
examples, figure helpers, and command-line scripts. It intentionally excludes:

- research datasets and derived HDF5 files;
- pretrained weights and checkpoints;
- TensorBoard logs and local caches; and
- archival DOI or Code Ocean records, which remain pending.

The example configuration files document the expected relative data layout.
They do not imply that those data files are bundled in this code-only release.

## Repository layout

| Path | Purpose |
| --- | --- |
| `carbon_train.py` | Training, validation, and evaluation entry point |
| `CARBON/` | Model, loss, operations, utilities, interpretability, and plotting code |
| `preprocessing/` | BED/HDF5 conversion and five-channel preprocessing |
| `configs/` | Example task and path configuration |
| `scripts/` | Training, evaluation, plotting, and regression commands |
| `figures/` | Example plotting and metric helpers |
| `demo/` | Example input manifest; data are not bundled |

## Installation

Python 3.9 is the reference version. Create the portable Conda environment:

```bash
conda env create -f environment.portable.yml
conda activate carbon
```

Alternatively, install the dependencies with pip. Select the PyTorch build
appropriate for your CPU or CUDA platform.

```bash
python -m pip install -r requirements.txt
```

## Data contract

CARBON expects five-channel HDF5 inputs with datasets `X`, `Y`, and `Y_ctl`.
For the default configuration, `X` has shape `(N, 101, 5)`: four nucleotide
channels plus one conservation channel. Paths in `configs/demo.yaml`,
`configs/demo_paths.json`, and `demo/demo_input.json` are examples that become
executable after compatible data are placed at those locations.

To build compatible inputs, see:

```bash
python preprocessing/bed_to_h5.py --help
python preprocessing/preprocess_5ch.py --help
```

## Smoke test without research data

After installing the dependencies, verify model construction and output shapes:

```bash
python scripts/smoke_model.py
```

Expected output:

```text
CARBON smoke test passed: logits_profile=(2, 101), log_total_count=(2, 1)
```

## Training and evaluation

Display all command-line options:

```bash
python carbon_train.py --help
```

Once compatible data are available at the configured paths:

```bash
bash scripts/run_demo_train.sh
bash scripts/run_eval_only.sh <model_params.pth> <output_dir>
bash scripts/run_demo_plot.sh <model_params.pth> <output.png>
```

## Authors and citation

- Wu Bin — College of Artificial Intelligence, Jilin University
- Li Xiangtao — College of Artificial Intelligence, Jilin University
  (corresponding author)

Citation metadata are provided in `CITATION.cff`. Journal identifiers, author
email, ORCID records, data archive URLs, and the final article DOI remain to be
added after confirmation or publication.

## License and third-party notices

CARBON is distributed under the GNU Affero General Public License v3.0 only;
see `LICENSE`. Components adapted from or inspired by permissively licensed
projects are documented in `THIRD_PARTY_NOTICES.md`.

