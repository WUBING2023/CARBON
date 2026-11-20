# CARBON manuscript — 20 November 2026

[Main manuscript](CARBON_main.pdf) (16 pages) · [Supplement](CARBON_supplement.pdf) (42 pages) · [Complete LaTeX source](CARBON_source.zip)

This version retains the scientific text and figures of the 19 November revision and incorporates the author metadata reviewed on 20 November: author affiliations, contributions, ORCID identifiers and the conflict-of-interest declaration. The funding field is left blank in this working manuscript.

## Authors

Bin Wu ([ORCID](https://orcid.org/0009-0002-6711-7230)), Yuting Huang, Haoran Zhu ([ORCID](https://orcid.org/0000-0002-7240-2696)) and Xiangtao Li ([ORCID](https://orcid.org/0000-0002-8716-9823)). Haoran Zhu and Xiangtao Li are corresponding authors.

Bin Wu led the computational analysis and drafted the manuscript. Yuting Huang provided the initial dataset. Haoran Zhu reviewed and revised the manuscript. Xiangtao Li supervised the study.

Conflict of interest: None declared.

## Compile

Extract `CARBON_source.zip` into a new directory. The archive contains both TeX files, the bibliography, journal class and style, and all referenced figure PDFs. With a LaTeX distribution and latexmk installed, run from the extracted directory:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error CARBON_main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error CARBON_supplement.tex
```

The standalone TeX files alongside this README are provided for inspection; compilation requires the companion assets in the ZIP. SHA256 checksums are listed in [checksums.json](checksums.json). The [19 November PDFs](../2026-11-19/) remain available as the prior revision.
