# ClinVar exact annotation result, 2026-10-14

The complete dated official GRCh38 snapshot was downloaded (200,437,060 bytes), checked against NCBI's MD5, and streamed through all **4,594,764 VCF records**. The complete existing **4,669 REF/ALT keys** remain in the output, each exactly once. There are **520 exact matches** and **4,149 unavailable matches in this VCF snapshot**.

| Aggregate germline annotation | Existing exact keys |
|---|---:|
| Benign / Likely benign | 453 |
| No germline classification | 52 |
| Conflicting classifications | 7 |
| Uncertain significance | 4 |
| Pathogenic / Likely pathogenic | 3 |
| Other classification | 1 |
| No matching record | 4149 |

The 3 P/LP keys were directly confirmed against current official ESummary canonical SPDI, classification and review status (**3/3**). The genomic-plus allele and RNA-oriented change are different on the ELP1 minus strand.

| Gene | GRCh38 allele | Current VCV | Review | Consequence | Clinical condition boundary |
|---|---|---|---|---|---|
| ASL | chr7:66081825 G>A | VCV000092360.115 | 2 stars, no conflicts | p.Arg12Gln missense | VCV list includes argininosuccinate lyase deficiency; protein mechanism cannot be reassigned to RNA binding. |
| ELP1 | chr9:108880051 C>T | VCV001066986.10 | 2 stars, no conflicts | RNA-oriented c.3460+1G>A, splice donor | VCV list includes familial dysautonomia, medulloblastoma and not provided. These are condition names, not proof that every condition has the aggregate LP classification. |
| ASXL1 | chr20:32434485 C>A | VCV002498883.28 | 2 stars, no conflicts | p.Tyr591Ter nonsense | The submitted condition is only **not provided**. No disease name or cancer-driver status is supplied. |

These keys have 6 old ASB events in total. Their old ASB labels/model responses can be displayed as observed allele preference, with the original development exposure retained. Clinical classifications and disease names are contextual annotations, not binding labels. The clinical labels do not establish an RNA-mediated disease mechanism or a conservation-specific benefit.

Current ESummary replies list no somatic clinical-impact or oncogenicity classification for any of these 3 records. These fields must not be populated from the separate TCGA matches.

The annotation-first availability pools contain 7,668 splice/UTR P/LP candidates in old context-gene symbols, 10,737 corresponding B/LB control candidates, and 364 reviewed P/LP substitutions within old center +/-50 intervals. They are **proposals, not new model experiments**. Parent chose the complete 3 existing P/LP keys for this revision; no additional candidate inference was executed here.

Population output retains the 241 prior exact 1000 Genomes ALT frequencies separately from 428 current matched ClinVar `AF_TGP` fields. Only 48 additional clinical-pool candidates carry that latter field. Missing population annotations remain missing; the VCF header does not specify TGP population release or denominators.

Read `research/CLINVAR_SOURCES.md`, `CLINVAR_SNAPSHOT_MANIFEST.json`, `DIRECT_VERIFICATION_MANIFEST.json`, `ALL3_PLP_DIRECT_VERIFICATION.csv`, and `population/POPULATION_MANIFEST.json` for source identity and replay. Current direct metadata supplement the dated snapshot rather than silently replacing it.
