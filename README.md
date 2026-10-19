# CARBON

CARBON predicts base-resolution RNA-binding protein profiles from sequence and conservation. This private review candidate combines the existing model/demo code with selected frozen research results and the current manuscript sources.

**Start here:** [S1-S27 data entries](paper/source_data/README_CURRENT.md) · [research evidence map](docs/EVIDENCE_MAP.md) · [replay and verification](docs/REPRODUCE.md) · [known limits](docs/KNOWN_LIMITS.md) · [data and weight availability](docs/DATA_AVAILABILITY.md).

| Path | Purpose |
|---|---|
| CARBON/, carbon_train.py, preprocessing/ | Existing model, losses, training/evaluation and data preparation |
| data/, carbon_trained_model/, configs/ | Small included AARS demonstration and checkpoint |
| research/experiments/E01–E11/ | Original analysis code, complete frozen result cohorts, run records and per-experiment guides |
| research/historical_model/ | MARCP implementation associated with archived historical launchers |
| research/catalog/ | 1,332-run ablation index, canonical weight availability, source hashes and external-asset/duplicate maps |
| paper/ | Current TeX, bibliography, used figure PDFs, source tables and editable figure inputs |
| scripts/verify_demo_cpu.py | Executable CPU check against six saved demo metrics |
| scripts/verify_frozen_results.py | Executable arithmetic checks for archived ablations, Figure 4 and Figure 5 cohorts |

## Install and run the checked paths

Portable Conda setup:

    conda env create -f environment.portable.yml
    conda activate carbon

Or create a Python 3.12 environment and install the CPU dependencies:

    python -m pip install -r requirements.cpu.txt

These dependency versions describe the environment used for validation; a fresh environment install has not been tested. The older environment.yml is a Linux-specific historical lock, not a portable installation recipe. Preprocessing and attribution use additional optional packages listed in requirements.txt; figure replays may also need PyMuPDF, R and circlize.

From the repository root:

    python scripts/smoke_model.py
    python scripts/verify_demo_cpu.py
    python scripts/verify_frozen_results.py
    python figures/plot_demo_signal_track.py --device cpu --out outputs/demo_signal_track.png
    python scripts/verify_package.py

The AARS verifier explicitly passes CPU and task options and evaluates 1,283 included test windows. It compares the six saved metrics with absolute tolerance 1e-5. The four-track plot contains CARBON, pooled eCLIP, inverted pooled eCLIP and inverted input control. The last two tracks are display transforms of Y and Y_ctl, not biological replicates.

The demo chr1/2/3 split is illustrative. Formal archived experiments use their recorded chromosome splits and configurations. Fresh training should use a new output directory, for example:

    bash scripts/run_demo_train.sh outputs/demo_train

## Scope of reproduction

The checked CPU demo and frozen-table calculations do not establish all historical training reproducibility, identical MARCP/refactor behavior, a unified fair three-model benchmark, or independent functional validation. Full/Sequence controls, unfavorable directions, and negative/mixed functional results are retained. Read each experiment's evidence boundary before interpreting its outputs.

A versioned, public full-data/full-weight download is **not available from this package**. Canonical archived checkpoint counts are CARBON 221, iDeepB 219 and RBPNet 222; exact missing tasks and SHA identities are listed. The older local Zenodo package contains a different 222-task set. No pending DOI or private platform URL is represented as a usable public data source. Code Ocean has only a historical demonstration run; complete-paper reproduction there has not been established.

## Provenance and use

The candidate is based on GitHub main 40ca6bfe3b7d601a08aa7aa68abcebd5773f5317, locally checked against 71 downloaded files on 2026-10-19. Source hashes, path-only redactions and current manuscript binding are machine readable. Large identical assets have one canonical copy; materialize a separate scratch tree with scripts/materialize_replay.py before historical figure replays. Packaging preserves supplied scientific data and creates no new experimental observations; file and arithmetic checks do not certify the authenticity of every historical measurement.

The original software and included AARS demo retain their existing licence and third-party notices. External source datasets and papers retain their original terms; this candidate does not apply a new blanket licence to them. Cite the manuscript and external assay/data sources identified in the paper and manifests.
