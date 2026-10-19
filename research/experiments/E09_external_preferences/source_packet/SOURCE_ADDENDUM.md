# Source qualification before model-agreement evaluation

2026-10-04: the authors' official morrislab/RBPbinding archive, pinned to Git tree b67d8fba065b3c010da672a4980c4ed440745710, includes original measured RNAcompete Z-scores for human PTBP1 (ENSG00000011304; RNCMPT00268 and RNCMPT00269). Both matching experiments will be evaluated and averaged, with their individual results retained. This follows the original protocol's conditional same-protein RNAcompete inclusion rule. No PTBP3 measurement is used.

The downloaded 51,769,595-byte archive's Git blob SHA1 is 1bf4ebb9103ea17dfe734e00f26f3e6d29a96508, verified against the author's repository tree. It contains 420 measured-experiment columns and 16,382 unique RNA 7-mers, rather than all 16,384 theoretical 7-mers: GAAGAGC and GCTCTTC are absent after U→T conversion. No missing score is imputed. Each 5-mer is averaged over all actually available overlapping 7-mer occurrences, preserving repeated occurrences. Every 5-mer remains represented. The realized denominator range is exported in EXTERNAL_COLUMNS.csv.

This source qualification was recorded after inspecting the table axes and before successful model-agreement evaluation. It changes no model rank rule, endpoint, output head, seed or window selection. The table's filename explicitly distinguishes original RNAcompete measurements from JPLE-inferred specificity profiles elsewhere in the repository.

Numerical/environment repair: initial GPU run stopped before attribution because deterministic CuBLAS requires CUBLAS_WORKSPACE_CONFIG=:4096:8. The rerun uses that environment. CuDNN GRU training-mode backward is enabled only on the one-layer GRU with zero internal dropout; all other layers remain eval. No parameters change. Primary 32-step IG was retained because 32/64-step maximum-five-mer positions agreed for all 32 QA windows; relative attribution L1 difference was 0.002436.

Local analysis uses the already installed D:/Anaconda/envs/rbp/python.exe (NumPy1.26.4, pandas2.2.3, scipy1.13.1, matplotlib3.10.0), because the base Anaconda has incompatible optional NumPy1 extensions. No package installation was performed.
