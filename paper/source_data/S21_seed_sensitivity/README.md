# S21: limited two-task training-seed sensitivity

This directory is a summary subset, not the full training/replay archive.
It contains the frozen protocol, original training/input manifest, completed
training summary, six task/seed sensitivity results, and independent validation
results (130 checks passed). No checkpoint or per-window prediction NPZ is
included here. The original manifest records original producer paths and hashes.

The complete existing archive is at
`M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/evidence_expansion_r4/training_repeats/`
and its server15 producer directory is
`/data/CARBON_PROJECT/evidence_expansion_r4_20260924/training_repeats/`.
These are local/private paths, not public release links.

Endpoint: task mean of paired within-window PCC differences, averaged over
three permutations; this differs from the main bulk analysis's task median.
Only HNRNPK/PTBP1 K562 have this three-seed check. New seeds were CPU-trained;
historical seed 0 was GPU-trained. This is not benchmark-wide seed validation.
