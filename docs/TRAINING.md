# Training and ablation entry points

Run from the repository root after installing the CPU dependencies. `python carbon_train.py --help` is the authoritative CLI. This command uses the included AARS HDF5 files for a short new training run:

```bash
python carbon_train.py \
  --lr 0.000785 --CUDA -1 --signalMax 0 --seqLength 101 --bsz 64 \
  --encode OneHot --modelFramework Treat --epochs 3 --input_channels 5 \
  -train data/h5data/AARS_K562_ENCSR825SVO/_5ch/AARS_K562_ENCSR825SVO_chr1_5ch.h5 \
  -validation data/h5data/AARS_K562_ENCSR825SVO/_5ch/AARS_K562_ENCSR825SVO_chr2_5ch.h5 \
  -test data/h5data/AARS_K562_ENCSR825SVO/_5ch/AARS_K562_ENCSR825SVO_chr3_5ch.h5 \
  --task AARS_K562_ENCSR825SVO --lossM CombinedProfileLoss \
  --alpha 1.42 --beta 0.54 --gamma 1.80 --weight_decay 0.000114 --grad_clip 2.0 \
  --output outputs/demo_train
```

The equivalent Bash launcher is `bash scripts/run_demo_train.sh outputs/demo_train`. It accepts the output path as its first argument; `DEMO_EPOCHS` changes its epoch count and `CUDA_DEVICE` its device selection. The default `--CUDA -1` selects CPU. GPU use requires a suitable PyTorch/CUDA environment and is outside the CPU validation.

Use a new output path for each condition. Training writes model/checkpoint state, metric tables and run metadata. An existing directory may be used with explicit `--resume` and should not be reused accidentally. Evaluate the supplied checkpoint without training with `python scripts/verify_demo_cpu.py`.

## New ablation runs

Add one of these actual CLI options to the Python command above, and change `--output` to a distinct new directory:

| Option | Implemented change |
|---|---|
| `--ablation_no_se` | Remove SE channel attention |
| `--ablation_no_multiscale` | Use the single k=5 convolution branch |
| `--ablation_no_gru` | Remove the BiGRU sequence mixer |
| `--ablation_lstm_instead_of_gru` | Substitute BiLSTM; mutually exclusive with `--ablation_no_gru` |
| `--ablation_no_combined_loss` | Use `--ablation_baseline_loss` (default `poissonLoss`) |
| `--ablation_no_conservative` | Drop the fifth channel from 5-channel inputs |

These options support new experiments. Their defaults do not reconstruct every historical run. Capacity changes and differences in training history must be considered when interpreting results.

## Historical paper experiments

[E01](../research/experiments/E01_data_model/README.md) describes data/model lineage. [E03](../research/experiments/E03_ablation/README.md) and [ABLATION_RUN_INDEX.csv](../research/catalog/ABLATION_RUN_INDEX.csv) locate the 1,332 archived runs (222 tasks × baseline and five ablations), including saved configurations and metrics. Recover exact upstream data, seeds, chromosome splits, filtering, loss, early stopping and checkpoint-selection rules from the corresponding records.

Historical launchers refer to the [MARCP implementation](../research/historical_model/). Account path prefixes were replaced with `/data/CARBON_PROJECT` and `/data/USER_HOME`; configure real inputs in a separate scratch tree. Omitted HDF5 files, weights and inputs are cataloged in [EXTERNAL_ASSETS.csv](../research/catalog/EXTERNAL_ASSETS.csv). Use [the experiment map](EVIDENCE_MAP.md) and [replay instructions](REPRODUCE.md) to choose a workflow. The small demo does not rerun the complete paper.
