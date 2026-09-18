# Logic equivalence verification

This document explains how to verify that **this refactor did not change the original training/eval logic** (except for two explicit bug fixes), and gives practical verification steps.

## 1. Reference baseline

- **Logic baseline**: the original MARCP tree (e.g. `ablation/MARCP/` in a full checkout; includes `marcp_train.py` and nested `MARCP/MARCP/`).
- The current tree (`CARBON/`, package name `CARBON`, entry `carbon_train.py`) matches that baseline **in logic**; on top of the baseline we reorganized structure, **renamed the model to CARBON**, updated comments, and fixed two bugs.

## 2. Change inventory

### 2.1 Structure / naming / comments only (no logic change)

| Kind | Description | Effect |
|------|-------------|--------|
| Package layout | Removed inner `MARCP/MARCP/`, single package `CARBON/` | Import paths only |
| Imports | `from CARBON.xxx` (or relative) | Reference style only |
| Renames | e.g. `carbon_model.py`, `predict.py` (PEP8) | Filenames and imports |
| Directories | `Attribution/`→`attribution/`, `plot/`→`visualization/` etc. | Paths and imports |
| Visualization | Two `tfmodisco_viz_sequence.py` merged into `visualization/viz_sequence.py` | Single implementation, same behavior |
| Entry point | `carbon_train.py` main under `if __name__ == "__main__":`, `seed_everything(0)` there | No train when imported as module |
| Comments | Chinese removed; English docstrings on modules/main functions | No runtime effect |

None of the above changes: data loading, model definition, loss, optimizer steps, metrics, or `eval_only` vs train branches.

### 2.2 Intentional bug fixes (correctness, same design intent)

| Location | Issue | Fix | Note |
|----------|--------|-----|------|
| **effective_in_channels** | Was only set inside `if ablation_subdir:` but used on main train path | Compute once after ablation config write, before `_test_files_exist`: `effective_in_channels = 4 if getattr(args, "ablation_no_conservative", False) else args.input_channels` | Matches `eff_ch` in `eval_only`; no NameError on train path |
| **writer** | `SummaryWriter` only in `eval_only`, but training loop used `writer` | Before “# Train data”: create `outputDir/tensorboard` and `writer = SummaryWriter(...)` | TensorBoard during training; loss/grad/metrics unchanged |

These are **bug fixes**: the old train path could error or behave undefined; after fix, behavior matches a correct MARCP-style setup.

### 2.3 Other fixes during refactor (no change to correct numerical meaning)

- **utils/u_statistics.py**: Wrong `base_percent(seqlist)` call → use `base_percent_cal` and unpack correctly.
- **loss/loss.py**: `mse_poisson` referenced `self.poisson_loss_test` → use `self.poissonLoss`.
- **models/carbon_model_rnafm.py**: Add `from torch import einsum` where `einsum` is used.

These avoid crashes or wrong outputs under a correct implementation.

## 3. How to check that logic was not unintentionally altered

Recommended: (1)+(2); optionally (3).

### 3.1 Deterministic regression (recommended)

**Idea:** Fixed seed, same data and hyperparameters — run a short training (e.g. 1–2 epochs) with **baseline** and **refactor**, compare key numbers.

1. Run the refactor (example):

```bash
cd /path/to/CARBON
python3 carbon_train.py ... # same data paths, bsz, epochs as baseline
```

2. Run the baseline (needs the original MARCP tree and `marcp_train.py`) with the **same** arguments and data; use a different `--output` (e.g. `exampledata/out_regression_baseline`).

3. **Compare** `train_metric.txt`, `test_metric.txt` / `test_metrics_summary.txt` for the same epoch and metrics (Mean PCC, Spearman, loss). **Exact match** (up to float noise) supports unchanged logic.

Use `scripts/verify_logic_regression.sh` to run the refactor side and write `regression_signature.txt` (see §4).

### 3.2 Module-by-module review (optional)

If you keep baseline MARCP code, compare data flow and formulas (ignore comments/names/imports) for:

- Data: `read_data` / `read_data_5ch`, `_read_data_single`
- Model: `CARBON`, `MultiScaleMotifExtractor`, `ModernBiGRU`, `Conv1dProfileHead`, RNAFM variant `forward`
- Loss: `CombinedProfileLoss`, `loss_fun`
- Training step: optimizer, gradients, Pearson/AUC
- Eval: `calculate_batch_metrics_whole`, `eval_only`, test branch

### 3.3 Unit / small tests (optional)

Same input → same output for:

- `onehot_encode`
- `CombinedProfileLoss.forward`
- `CARBON.forward` with fixed weights
- `read_data_5ch` on one example h5

## 4. Regression script

`scripts/verify_logic_regression.sh`:

- Uses the same fixed seed as `seed_everything(0)` in `carbon_train.py`
- Runs exampledata for 2 epochs (editable)
- Writes `exampledata/out_regression_refactor/regression_signature.txt` from train/test metrics

Run from repo root:

```bash
bash scripts/verify_logic_regression.sh
```

Diff against the baseline’s `regression_signature.txt` with identical args.

## 5. Summary

- **Unchanged:** data flow, model math, loss, optimization, evaluation, `eval_only` vs train split.
- **Changed:** layout, filenames, imports, comments, entry guard, two main bug fixes (`effective_in_channels`, `writer`), three small fixes (u_statistics, loss, einsum).
- **How to confirm:** same-seed regression on key metrics (recommended), plus optional module review or unit tests.
