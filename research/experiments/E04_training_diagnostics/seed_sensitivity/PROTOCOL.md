# Frozen full-data training-seed repeats — 2026-09-24

This protocol is frozen before the four jobs start. It tests training-seed sensitivity for two prespecified representative tasks, not superiority across all 222 tasks. No result-dependent task, epoch or seed selection is allowed. The original sources, checkpoints and datasets remain unchanged.

## Cohort and original training recipe

- HNRNPK_K562_ENCSR268ETU: train 41,990, validation 10,278, test 9,787 windows.
- PTBP1_K562_ENCSR981WKN: train 48,695, validation 12,111, test 14,057 windows.
- Fresh seeds: 1 and 2 for each task; four complete training jobs, not reduced-data demonstrations. Historical seed0 is re-inferred, not retrained.
- Train chromosomes: 3,4,5,6,7,10,11,12,13,14,17,18,19,20,21,22. Validation: 2,9,16. Original test: 1,8,15. Sex and mitochondrial chromosomes are excluded.
- Exact five-channel MARCP/CARBON complete model, 100,826 trainable parameters, dropout 0.2; original source copied from /data/CARBON_PROJECT/MARCP.
- Only scientific-code change: each isolated trainer replaces its one `seed_everything(0)` call with seed1 or seed2. No architecture, loss, loader, schedule, stopping, checkpoint or metric changes.
- AdamW, lr 0.000785, weight decay 0.000114; CombinedProfileLoss alpha=1.42 beta=0.54 gamma=1.80, loss_scale=1; clip norm2; batch256; max100 epochs.
- Preserve original training and validation shuffle=True/drop_last=True. Original early stopping/checkpoint maximizes batch-average validation Peak AP, patience9. The original log incorrectly labels this as val_loss. Validation loss controls lr halving after5 non-improving epochs; min lr1e-7. No cosine/plateau scheduler, pretraining, resume, filtering or log transform.
- Original test loads the best saved checkpoint. Retain all failures, logs and results. Each job has a fresh unique output directory.

## Resources and provenance

CPU only; four processes maximum, four intra-op/BLAS threads per training job, Linux nice10. CUDA_VISIBLE_DEVICES=-1 is set before import. No original launcher, GPU, directory cleanup or unrelated process is touched. The supervisor waits through the original early stopping without a wall-clock truncation being reported as a completed repeat.

Freeze full source hashes, exact argv, HDF5 file hashes and shapes, baseline checkpoint hashes, protocol/runner hashes and Python/library versions before launch. Record PIDs and actual start/end times. Python bytecode writing is disabled.

Original entry SHA256: 125b379f25b3894697e64ce435e5051325969e6e1e4687f3424c3dc3cff221d3.
Original model SHA256: e101164b9bf0d5d8cb4a3471e5cfdf2ac7a6a1cdc4f6b710f4132c3bf2fff55a.
Original loss SHA256: 6f06f377e64f05c33ad3b693cc47957cdd6c29c328039c5f23ea4acc553e3f4e.

## Fixed held-out evaluation

After training, use every window in chr8 and chr15 (HNRNPK 3,327; PTBP1 5,299), independent of training and checkpoint selection. Re-infer historical GPU-trained seed0 and new CPU-trained seeds1/2 with the same CPU code, data and batches. This controls inference hardware but does not remove the historical training-hardware difference. Do not treat the new repeats as a hardware-controlled three-seed retraining study.

Five input conditions: true conservation; its within-window mean; within-window permutations with seeds17,43,97. For each task/chromosome/permutation, a seeded RNG generates a separate permutation of101 positions for each window, and the resulting input is identical across model seeds. Four sequence channels remain unchanged. True targets and input files are immutable. Windows with constant target or predicted profiles have undefined PCC and are excluded pairwise from PCC contrasts with their counts explicitly retained.

Primary descriptive endpoint: mean per-window PCC and true-minus-perturbed per-window PCC difference on pairwise-valid windows. Also export median differences, positive-window counts, pooled six original metrics, count log1p absolute errors, raw per-window condition scores and prediction arrays. No inferential p-values from treating seeds or nested windows as independent biological replicates. Seeds1/2 are sensitivity checks on these two chosen tasks, not additional RBPs.

The original three-chromosome test metrics remain available for historical parity. Primary new analyses use chr8/15 because the historical hyperparameter-search recipe used chr1 for training, chr2 for validation and chr3 for testing.
