[![DOI](https://zenodo.org/badge/91751978.svg)](https://doi.org/10.5281/zenodo.17246480)

# Archived Code Collection:

**Associated with the paper:** <br>[A resource of RNA-binding protein motifs across eukaryotes reveals evolutionary dynamics and gene-regulatory function](https://www.nature.com/articles/s41587-025-02733-6)

### This repository is no longer maintained! 

**Please use the actively maintained implementation of JPLE here**: <br>[**https://github.com/morrislab/jple**](https://github.com/morrislab/jple/)

If you want to use JPLE to train your own model for predicting sequence specificity profiles from protein sequence, or to infer RNA sequence specificity for RRM- or KH-domain RNA-binding proteins, please refer to: [https://github.com/morrislab/jple](https://github.com/morrislab/jple/)

All the measured and inferred RNA sequence specificities can be found at: <br> [**cisbp.org**](https://cisbp.org/)

## Repository Contents

This repository serves as an archival collection of scripts and command-line tools used for analyses during the PhD thesis: <br> [Inferring RNA Sequence Specificities from Protein Sequences to Characterize Post-Transcriptional Regulation in Eukaryotes](https://www.proquest.com/openview/d075b03685b79e4426cb11a18d9f19ad/1?pq-origsite=gscholar&cbl=18750&diss=y) 

**Maintenance is limited**, but the repository contains the scripts used to reproduce figures and results presented in the thesis.

#### Recommended setup 

The majority of scripts will run on a virtual environment with python 2.7 adding dependencies listed in [dependencies.txt](https://github.com/LXsasse/RBPbinding/blob/master/dependencies.txt)

Some scripts may require python3 or other dependencies. Please install as required: 

- python3 (agglomerative_clustering.py), scikit-learn==0.23.2)
- Hmmer (http://hmmer.org/)
- conservation_code (Capra JA and Singh M. Predicting functionally important residues from sequence conservation. Bioinformatics, 23(15):1875-82, 2007.
(https://compbio.cs.princeton.edu/conservation/))
- pymol (https://pymol.org/2/) to visualize individual pdbs

#### Reproduce thesis results:

To reconstruct the figures (`fig1-5.sh`), first run the preprocessing scripts in the following order:

1. `rncmpt_data.sh`
2. `fig1.sh`
3. `performance_calc.sh`
4. `interface_importance.sh`
5. `fig2.sh`
6. `jple_reconstruction.sh`
7. `fig3.sh`
8. `cisbp-recstats.sh`
9. `arabidopsis.sh`
10. `fig4.sh`
11. `fig5.sh`

#### License

The repository is licensed under the BSD 3-Clause License. See [LICENSE](https://github.com/lxsasse/RBPbinding/blob/main/LICENSE.md) for details.

#### Citation

If you use this code, please cite the following:

- **Primary paper:**

  > Sasse A., Ray D., Laverty, K.U. et al. **A resource of RNA-binding protein motifs across eukaryotes reveals evolutionary dynamics and gene-regulatory function.**
  > *Nature Biotechnology* (2025). https://doi.org/10.1038/s41587-025-02733-6

- **PhD Thesis:**

  > Sasse A. **Inferring RNA Sequence Specificities from Protein Sequences to Characterize Post-Transcriptional Regulation in Eukaryotes.**
  > PhD Thesis, University of Toronto (2021). [*Proquest link*](https://www.proquest.com/openview/d075b03685b79e4426cb11a18d9f19ad/1?pq-origsite=gscholar&cbl=18750&diss=y)
