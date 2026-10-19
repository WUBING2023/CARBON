# ClinVar and population annotation: official sources and scope

This branch uses public static files only. The web-access skill was read, and the parent checked Chrome dependencies before network access. No login, external write, contact, upload, patient-level download or new inference occurs in this branch.

## Official source snapshot

- [NCBI GRCh38 VCF directory](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) was accessed on 2026-10-14. It lists the dated `clinvar_20261011.vcf.gz` release posted on 2026-10-12. HTTP HEAD reports **200,437,060 bytes**. The dated file is used instead of the mutable `clinvar.vcf.gz` alias.
- [Dated snapshot](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20261011.vcf.gz) and [official MD5](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20261011.vcf.gz.md5). The source header states `fileDate=2026-10-11` and `reference=GRCh38`. Download headers, full compressed bytes, header text, SHA256/MD5 and exact retained VCF records are stored under `work_disease/clinvar/raw` and its manifest. A network interruption is handled by byte-range continuation; the final digest must match before analysis is declared complete.
- [NCBI file guide](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/) explains variant-level VCV, condition-level RCV, VCF coordinate conventions and coverage. This VCF includes precisely mapped simple alleles under 10 kb; it does not represent every possible ClinVar record or genotype assertion.
- [Classification definitions](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/) and [review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/) were accessed. `CLNSIG`, `ONC` and `SCI` are separate aggregate germline, oncogenicity and somatic clinical-impact classifications. The three associated review fields remain separate. Star counts here apply to germline classification only; a conflict can have one star and is still excluded from P/LP candidates.

## Exact matching and clinical selection

All 4,669 existing GRCh38 chromosome, 1-based position, REF and ALT keys are matched exactly. An absent record is explicitly missing, never benign. Every overlapping row retains Variation ID, unversioned VCV accession, genomic HGVS, Allele ID, disease names and IDs, review status, molecular consequence, allele origin and all three clinical classification types.

Additional candidates are selected **before model inference or inspection of new model effects**: exact A/C/G/T SNV; non-conflicting aggregate germline Pathogenic, Likely pathogenic or Pathogenic/Likely pathogenic; at least one germline review star; splice or UTR molecular consequence; ClinVar gene symbol overlaps an old model-context gene symbol. All qualifying candidates and all B/LB counterparts are retained. A separate file includes every reviewed P/LP SNV inside an existing native 101-nt center ±50 interval. Gene-symbol overlap is an availability screen; transcript coordinates, strand, genomic REF and native input still require verification before inference.

VCV disease-name lists describe conditions represented in the aggregate record. They are not a disease-specific pathogenicity assertion table; such assertions would require RCV/SCV inspection. A clinical P/LP classification is not an RBP-binding measurement, and a large model prediction is not a clinical causal explanation.

## Population background

The prior source-audited 241 exact ALT frequencies from `study_name=1000Genomes` are reused without changing observations. Their original counts, total counts and study version are preserved. [Official IGSR](https://www.internationalgenome.org/data/) identifies the population resource; the actual numerical snapshot remains the archived [NCBI Variation Service](https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/10100210) responses already bound to exact alleles in the v4 package.

The current ClinVar VCF `AF_TGP` header says allele frequencies from TGP. That field is kept separately for current matched records and new clinical candidates. Its population version and denominator are not specified by this VCF and are not invented. Absence is missing rather than zero. ALT frequency is not automatically minor-allele frequency; no biallelic assumption is introduced. Population annotations do not provide a disease or RBP-binding gold label.

## Replay

After the complete dated source passes MD5:

```powershell
D:/Anaconda/python.exe -X utf8 work_disease/clinvar/parse_snapshot.py
D:/Anaconda/python.exe -X utf8 work_disease/population/build_population.py
```

These scripts use the Python standard library, preserve the existing variants and emit new local tables only. No original v4 annotation or prediction file is overwritten.
