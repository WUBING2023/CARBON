# TCGA public somatic annotation: official sources and exact-match results

Scope: authorized local CARBON disease extension, Figure 5 option 2. Source retrieval occurred on 2026-10-14; exact UTC request times are in `work_disease/tcga/requests/*_result.json`. No controlled-access file, authentication token, patient clinical record, model inference, or new experiment was used. The fixed input is all 4,669 existing GRCh38 genomic REF/ALT keys, not a subset chosen by CARBON response.

## Primary sources actually accessed

1. [NCI GDC API Search and Retrieval](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/): documents POST JSON filters, field selection, pagination, and facets. Nineteen batches of up to 250 genomic changes were queried. Every successful response was checked for exact full pagination, no warning, GRCh38, chromosome, 1-based position, REF and ALT. GDC returned six records; there are no unresolved batches or keys.
2. [NCI GDC API Data Analysis](https://docs.gdc.cancer.gov/API/Users_Guide/Data_Analysis/): documents `/ssms`, `/ssm_occurrences`, transcript consequences and `/analysis/mutated_cases_count_by_project`. An occurrence is a mutation in a case; project facets were requested with `size=0`, so the six count queries returned aggregate counts without case records. Final denominator retrieval also used `size=0` and returned no case hits. Case-level observations or aliquot calls were not counted as biological replicates.
3. [GDC Whole Exome and Targeted Sequencing Pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/): explains somatic variant annotation and masking for public release. Open-access masking can omit genuine somatic calls, and exome/targeted coverage is incomplete. Thus an unmatched allele is **not found in the current open GDC index**, rather than proven absent from TCGA or biologically unaffected. This study did not pursue controlled-access data to recover omitted calls.
4. [GDC Aggregated Somatic Mutation](https://docs.gdc.cancer.gov/Encyclopedia/pages/Aggregated_Somatic_Mutation/): describes project-level somatic annotations derived from case-level variants. We retrieved the public API index, not MC3 or a downloaded protected MAF.
5. [Live GDC API status](https://api.gdc.cancer.gov/status): `status=OK`, API tag `8.5.0`, data release **46.0**, release date **2026-08-10**, commit `8f7c2a51ab0084b216ad1b62a3fae8b945439c53`. Saved response: `tcga/raw/status.json`. Live field dictionaries from `/ssms/_mapping` and `/ssm_occurrences/_mapping` are saved with response hashes.

## Fixed matching contract and complete result

The server filter combined `genomic_dna_change in [chrN:g.positionREF>ALT]`, `ncbi_build=GRCh38`, and `occurrence.case.project.program.name=TCGA`. Client-side assertions independently repeated all four allele-key components. No liftover, transcript-strand allele substitution, protein-name match, or nearby-position match was used.

All 4,669 keys were queried, yielding **six exact variants** and **4,663 successful no-match results**. Every existing key has a retained status in `ALL_EXISTING_TCGA_QUERY_STATUS.csv`. All six exact variants are retained; each occurs in **one indexed case in one TCGA project**. These are singleton observations, not recurrent mutations.

| GRCh38 genomic key | Gene | Canonical consequence | TCGA project | Occurrence cases / project cases with SSM |
|---|---|---|---|---|
| chr11:7994518:G:T | EIF3F | splice_donor_variant | TCGA-LUSC | 1 / 504 |
| chr12:103947412:A:T | HSP90B1 | missense_variant | TCGA-UCEC | 1 / 559 |
| chr17:21312245:G:A | MAP2K3 | missense_variant | TCGA-STAD | 1 / 443 |
| chr17:21312240:C:T | MAP2K3 | synonymous_variant | TCGA-BRCA | 1 / 1098 |
| chr19:6686745:C:T | C3 | splice_donor_variant | TCGA-LGG | 1 / 516 |
| chr20:32434485:C:A | ASXL1 | stop_gained | TCGA-BLCA | 1 / 412 |

There are **57 transcript consequence rows**, all preserved, including alternative transcript UTR and noncoding annotations. A variant may be missense in the canonical transcript and UTR/noncoding in another transcript; its alternative annotation must not replace the canonical consequence without explicitly naming that transcript. The case denominators are GDC project cases with indexed SSM data, **not locus-callable sequencing denominators**. Derived fractions describe the current public index and are not epidemiological incidence, germline frequency, or clinical risk estimates. Official cancer project names, not manually assigned disease terms, are in `TCGA_VARIANT_PROJECT_COUNTS.csv`.

## Scientific use and limits

This is an annotation extension of an already explored ASB-positive cohort. TCGA establishes that the exact genomic allele was observed as a public somatic call in the listed tumor cohort. It supplies no independent allele-specific RBP binding measurement, disease-causality label, pathogenicity label, or evidence that a coding consequence acts through RNA binding. The independent binding direction, if shown, remains the existing exposed ASB evidence. CARBON Full/Sequence effects must remain clearly identified as predictions.

Use all six variants, all their existing RBP/cell events, and the full missing-match denominator. Keep splice, missense, synonymous and truncating consequences separate or label them explicitly. With six singleton overlaps, do not report a clinical discrimination AUC, cancer enrichment test, recurrent-driver claim, or conservation-specific explanation advantage. A display example may be selected from these six by a predeclared coverage rule, never by favorable model effect; all other matches remain visible.

## Reproduction and deliverables

- `work_disease/tcga/QUERY_PROTOCOL.json`: input SHA256, filtering and completeness rules.
- `work_disease/tcga/ALL_EXISTING_TCGA_QUERY_STATUS.csv`: all 4,669 keys, including successful no matches.
- `work_disease/tcga/TCGA_EXISTING_VARIANTS.csv`: all six exact variants and canonical/all-transcript consequence summaries.
- `work_disease/tcga/TCGA_ALL_TRANSCRIPT_ANNOTATIONS.csv`: all 57 transcript rows.
- `work_disease/tcga/TCGA_VARIANT_PROJECT_COUNTS.csv`: all six project counts and appropriately labelled SSM availability denominators.
- `work_disease/tcga/TCGA_PROJECT_SSM_DENOMINATORS.csv`: all 33 TCGA project denominator rows.
- `work_disease/tcga/raw/` and `requests/`: official source JSON, exact payloads, HTTP status, headers, timestamp, bytes and SHA256 for each substantive request.
- `work_disease/tcga/scripts/query_existing.py`: existing-cohort query with safe response reuse.
- `work_disease/tcga/scripts/annotate_matches.py --offline`: rebuild annotations from frozen responses without network.
- `work_disease/tcga/scripts/verify_tcga.py`: complete input/result, allele, transcript, count and denominator checks.

Transport note: Python `requests` initially failed with TLS EOF. Windows `curl.exe` succeeded with ordinary HTTPS certificate validation; SSL verification was never disabled. The initial default denominator probe returned empty-source public case UUID hits; they were not used in analysis. The final denominator source is `raw/project_ssm_denominators_aggregate_only.json`, which suppresses all case hits. Ship aggregate CSVs and final aggregate source responses, rather than the unused probe, with any manuscript source-data release.
