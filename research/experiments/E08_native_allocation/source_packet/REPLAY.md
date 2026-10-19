# Replay and source map

Final scientific outputs are RESULTS.md, FINAL_CAPTION.md and REPORT.md. Final image/source/QA are figure/CONTEXT_FINAL.{png,pdf,svg}, FINAL_FIGURE_SOURCE.json, FINAL_FIGURE_QA.json and FINAL_PLOTTED_GENES.csv. figure/CAPTION.md belongs only to CONTEXT_CANDIDATE, not final. Original fixed counterexamples remain in figure/FIXED_COUNTEREXAMPLES.*. The final figure was rendered by render_final.py without any new model query.

The new local-position control is implemented in nogru_position_capacity/control_model.py. PROTOCOL.md was frozen before its pilot or outcomes; remote/TRAINING_MANIFEST.json records all six complete runs. remote/INITIAL_QA.json verifies fresh shared initialization and zero initial output differences. CONTROL_LEARNING_QA.json verifies learned private tensors. remote/INFERENCE_MANIFEST.json maps each native profile and checkpoint identity. A legacy training field named arm=aligned denotes the five-channel input transformation, not the model class; actual control class, architecture, parameter count and checkpoint path are independently audited and preserved.

Pair eligibility in the already frozen implementation: native centers use zero-based coordinates, start + floor((core length−1)/2), must lie at 5–95, and are visited in sequence order. Greedy accept the next center at least 12 nt from the previous accepted one and take the first two; supports c−5 through c+5 are therefore separated. This describes existing code, not a changed eligibility rule.

Native inputs, measured signal, input-control signal, coordinates, original source-row indices, all four trained Phi arms and all seeds are stored under remote/data and phase2_ambiguity/remote/data. Newly trained NoGRU outputs are in nogru_matched/remote/data and nogru_cross_protein/remote/data. Local-position outputs are in nogru_position_capacity/remote/data. Each directory contains source/checkpoint/init hashes. No locally derived sequence material was uploaded to third-party services.

For a reversible replay, copy this entire 01_context directory to a new writable directory before running scripts: paths are resolved relative to each script, and analysis overwrites its tables. With Python and numpy/pandas/matplotlib/torch/pypdf installed, use:

    python nogru_position_capacity/analyze_context.py
    python nogru_cross_protein/analyze_context.py
    python render_final.py
    python finalize_delivery.py

Use the recorded environments and numpy version for exact bootstrap agreement. The local-position analysis runs 10,000 paired joint-gene bootstrap replicates using PCG64 seed 20261008; profile ensembles are float64 means before shares. ANALYSIS_INTERVAL_NOTE.md defines the main two-contrast family and supporting Bonferroni 3/15/45/90 families. All contrasts, methods, cell/seed/range, raw, strict and historical-exposure rows are in full CSVs; bootstrap arrays are retained. Never conflate comparison=full_minus_control with comparison=position_capacity_minus_nogru when filtering control=nogru.

No retraining is needed to reproduce final evidence or image. Authorized server15 train/inference replay wrappers are retained only for provenance. They refer to the existing matched source, native HDFs and GENCODE mappings, use original original full budget, and may launch expensive jobs; do not run them merely to inspect this delivery. A new training replay must use a new remote directory and a passing instantaneous GPU UUID guard, without interrupting other work.

Earlier unsuccessful conservation hypotheses remain at the root CSVs, phase2_ambiguity and test_8_15_support. They are not substituted by the positive component comparison. WINDOW_PROVENANCE_AUDIT.md records why native windows cannot all be called fixed peak-centered windows; RF_AUDIT.md establishes single-position RF 7 nt and union support RF 17 nt, while shared global SE remains.

Independent cross-checks reside in sibling 02_external_preferences/results/CROSS_AUDIT_01_PC_* and POSITION_CAPACITY_MODEL.json, and 05_delivery scope review. They belong to other workers and are read-only here. Final SHA256 identity inventory is DELIVERY_MANIFEST.json; it excludes itself, transport tar archives and generated proof thumbnails.
