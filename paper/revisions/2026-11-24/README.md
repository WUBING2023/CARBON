# CARBON manuscript — 24 November 2026

**Predicting nucleotide-resolution protein–RNA binding signals from RNA sequence and evolutionary conservation using CARBON**

[Main manuscript](CARBON_main.pdf) (16 pages) · [Supplement](CARBON_supplement.pdf) (42 pages) · [Complete LaTeX compilation source](source/) · [Checksums](checksums.json)

This revision uses the author-requested title, the latest five-author metadata, funding and author contributions. It also clarifies the search versus full-scale ablation settings in Supplementary Table S4 and corrects the ACVR1B description to two reversed permutations. Experimental values and figure files are unchanged from the audited source revision.

## Authors

Bin Wu, Yuting Huang, Yuning Yang, Haoran Zhu and Xiangtao Li. Haoran Zhu and Xiangtao Li are the manuscript's corresponding authors. The separately prepared cover letter is signed by Xiangtao Li.

## Author contributions

Bin Wu (Methodology, Investigation, Formal analysis, Writing – original draft), Yuting Huang (Data curation: initial dataset, Methodology), Yuning Yang (Methodology, Writing – review & editing), Haoran Zhu (Methodology, Validation, Writing – review & editing), and Xiangtao Li (Conceptualization, Supervision, Writing – review & editing).

## Funding

This work was supported by the National Natural Science Foundation of China [grant numbers 62472195 and 62076109 to X.L., and 62401131 to Y.Y.], the Natural Science Foundation of Jilin Province [grant number 20260102302JC to X.L.], and the Noncommunicable Chronic Diseases-National Science and Technology Major Project [grant numbers 2026ZD0554000 and 2026ZD0554005 to X.L.].

Conflict of interest: None declared.

## Compile

From the repository root, enter `paper/revisions/2026-11-24/source/` and run:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error CARBON_main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error CARBON_supplement.tex
```

The source directory includes the two manuscripts, their bibliography, class/style and all referenced figure files. It does not contain the full experimental data or model weights. See the [Supplementary Data index](../../source_data/README_CURRENT.md) for data availability. The [20 November revision](../2026-11-20/) remains an archive.
