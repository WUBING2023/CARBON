# CARBON

**CARBON** (Base-resolution RBP binding profile prediction with conservation signals) combines RNA OneHot encoding with evolutionary conservation (e.g. PhyloP) to predict base-resolution RBP binding profiles (profile + total count).

This repository is a complete code release for **GitHub / archival (DOI)**: training, validation, testing, preprocessing, and example plotting scripts, plus a **minimal demo** dataset (`AARS_K562_ENCSR825SVO`, one `_5ch.h5` per chr1/2/3).

**Full datasets:** The complete, detailed datasets used in the paper are archived on **Zenodo**: [URL] (replace `[URL]` with the published Zenodo record before release). This repo ships only the small demo subset for quick sanity checks.

**Figure reproduction:** Full reproduction of main figures, supplementary figures, and associated statistics is primarily provided on **Code Ocean**: [URL] (replace `[URL]` with your capsule URL). This repository offers a minimal local demo and scripts to verify the environment and API; batch plotting paths for publication-scale figures follow the Code Ocean configuration.

## Repository layout (summary)

| Path | Description |
|------|-------------|
| `carbon_train.py` | Single entry for training, validation, **prediction / evaluation** (`--eval_only`) |
| `CARBON/` | Model, losses, utilities, visualization |
| `data/h5data/.../_5ch/*.h5` | Demo 5-channel HDF5 (relative paths) |
| `configs/demo.yaml`, `configs/demo_paths.json` | Demo task id and data/output paths |
| `preprocessing/` | `bed_to_h5.py` (17-col BED→h5), `preprocess_5ch.py` (PhyloP→5ch) |
| `figures/plot_demo_signal_track.py` | **Main demo figure**: first **four tracks** as in paper Figure 2 (CARBON, eCLIP, two inverted controls) |
| `figures/supplementary/plot_train_val_curves.py` | **Supplementary-style plot**: train/val loss and Pearson vs epoch |
| `figures/stats/print_test_metrics.py` | Print `test_metrics_summary.txt` |
| `scripts/run_demo_train.sh` | One-command training on demo data |
| `scripts/run_demo_plot.sh` | One-command demo track figure |
| `scripts/run_eval_only.sh` | Load weights and evaluate on test h5 only |
| `demo/demo_input.json` | Small example input (path list) |
| `exampledata/out_test/` | Example metrics text and small weights (optional) |
| `environment.yml` | Conda env (remove trailing `prefix:` when recreating on another machine) |
| `requirements.txt` | Core pip dependencies |
| `LICENSE` | GPL v3 + academic/non-commercial clause (see file header) |

## Environment

**Conda (recommended)**

```bash
conda env create -f environment.yml
conda activate CARBON
```

If `environment.yml` contains a machine-specific `prefix:` line, delete it before creating the environment elsewhere.

**pip**

```bash
pip install -r requirements.txt
# Install torch/torchvision from https://pytorch.org for your CUDA/CPU setup
```

## Demo data and paths

Demo files:

- `data/h5data/AARS_K562_ENCSR825SVO/_5ch/AARS_K562_ENCSR825SVO_chr1_5ch.h5` (train)
- `data/h5data/AARS_K562_ENCSR825SVO/_5ch/AARS_K562_ENCSR825SVO_chr2_5ch.h5` (validation)
- `data/h5data/AARS_K562_ENCSR825SVO/_5ch/AARS_K562_ENCSR825SVO_chr3_5ch.h5` (test)

Task ID: `AARS_K562_ENCSR825SVO` (matches `task_id` in `configs/demo.yaml`).

When switching data, update `configs/demo.yaml`, `demo/demo_input.json`, and pass new `-train` / `-validation` / `-test` / `--task` on the command line.

## Quick start (minimal demo)

From the repository root:

```bash
# 1) Train (default 3 epochs; set DEMO_EPOCHS=1 for a smoke test)
bash scripts/run_demo_train.sh

# 2) Figure matching Figure 2 first four tracks (needs model_params.pth from step 1)
bash scripts/run_demo_plot.sh
```

Outputs:

- Model and metrics: `carbon_trained_model/AARS_K562_ENCSR825SVO_demo/`
- Figure: `figures/outputs/demo_signal_track_4rows.png`

Without a GPU, training and plotting fall back to CPU (slower).

## Training / validation / evaluation

- **Training** (see `carbon_train.py --help` for all flags):

```bash
python3 -u carbon_train.py \
  --modelFramework Treat \
  --lossM CombinedProfileLoss \
  --input_channels 5 \
  -train <train_5ch.h5> \
  -validation <val_5ch.h5> \
  -test <test_5ch.h5> \
  --task <TASK_ID> \
  --output <output_dir> \
  ... # lr, epochs, etc.: see configs/demo.yaml or the paper
```

- **Evaluation only** (existing `model_params.pth`):

```bash
bash scripts/run_eval_only.sh
# or
python3 carbon_train.py --eval_only --load_model <model_params.pth> -test <test_5ch.h5> \
  --task <TASK_ID> --output <output_dir> ... # same model/data flags as training
```

## Data preprocessing

1. **BED → 9-column h5** (`key=df`): `python3 preprocessing/bed_to_h5.py -i <bed> -o <out.h5>`
2. **Per-chromosome split/merge** can be done upstream (see paper data pipeline).
3. **4-channel table + PhyloP → 5ch HDF5** (`X`, `Y`, `Y_ctl`):  
   `python3 preprocessing/preprocess_5ch.py --input_dir <h5data/TASK> --bw_path <hg38.phyloP100way.bw>`

## Figures and statistics

Full figure and statistics reproduction: **Code Ocean** [URL] (see above).

Zenodo hosts the **full datasets** [URL] (see above).

- **Main demo**: `figures/plot_demo_signal_track.py` — four tracks aligned with Figure 2 (CARBON / eCLIP / inverted tracks); no other comparators required.
- **Supplementary-style curves**: `python3 figures/supplementary/plot_train_val_curves.py --train-dir carbon_trained_model/AARS_K562_ENCSR825SVO_demo`
- **Metrics**: `python3 figures/stats/print_test_metrics.py carbon_trained_model/.../test_metrics_summary.txt`

## Prediction (Python API)

See `CARBON.models.carbon_model.CARBON`: input shape `(batch, 101, 5)`; outputs `logits_profile` and `log_total_count`; profile is `softmax(logits) * exp(log_total_count)`.

## Citation

If you use this code, please cite: **CARBON: Base-resolution RBP binding profile prediction with conservation signals**.

## License

See `LICENSE` (GPL v3, with additional academic/non-commercial terms in the header).
