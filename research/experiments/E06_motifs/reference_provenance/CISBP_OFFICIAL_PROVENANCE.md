# CISBP-RNA reference provenance

Verified against the official site on 12 October 2026:

- [CISBP-RNA homepage](https://cisbp-rna.ccbr.utoronto.ca/): the site identifies the current release as Database Build 2.00 and announces the major 2.0 data release. Its displayed last-update date reads `7-5-2025` (retained as displayed, without assuming a date convention).
- [Official bulk downloads](https://cisbp-rna.ccbr.utoronto.ca/bulk.php): PWM files and RBP annotations are downloadable; the RBP Information file provides the associated metadata. The page identifies Database Build 2.00.
- [Official help](https://cisbp-rna.ccbr.utoronto.ca/help.html): direct motifs represent experimental characterization of an RBP; inferred motifs can be assigned from related proteins. The motif type and study fields identify the assay/source. JPLE reconstructs RNA preference profiles using a joint latent embedding, which supports classifying JPLE entries as inferred profiles rather than direct experimental measurements.

The manuscript may cite these official URLs for the database release and provenance definitions without introducing an unverified bibliographic entry. The exact source records are retained in `CISBP_exact_metadata_records.csv` and the complete human metadata in `source_inputs/reference_data/CISBP-RNA/homo_sapiens/RBP_Information_all_motifs.txt`.

## Local acquisition record

The existing local `Motif/reference_data/CISBP-RNA/README.txt` records a download date of **27 March 2026** and **Database Build 2.00**. This is a retained local acquisition record; the figure worker did not perform a fresh database download. The human metadata files carry 27 March 2026 modification times, which are consistent with that record but are not independently authenticated server timestamps.

The retained acquisition scripts identify the official endpoint `https://cisbp-rna.ccbr.utoronto.ca/bulk_archive.php`. `download_cisbp_pwm.sh` describes separate requests for full-library PWMs and RBP Information, named `cisbp-rna_PWMs_entire_2_00.zip` and `cisbp-rna_RBPinfo_2_00.zip`. `filter_cisbp_homo_sapiens.py` selects `RBP_Species == Homo_sapiens` and obtains PWM files from exact nonempty `Motif_ID` values in the complete metadata. This describes the intended workflow recorded by the scripts; it is not a recovered execution log. Original archive files were not located in the local reference directory.

The present audit uses the existing **homo_sapiens** subset, exact RBP name plus motif ID joins, and all 78 unique Target_ID matrices used by the 193 retained comparisons. The `_2.00` identifiers, acquisition record and source checksums bind the files to the retained Build 2.00 source set. `SOURCE_INPUT_MANIFEST.json` records exact paths, byte lengths and SHA-256 digests, including the acquisition-record files.

## Interpretation used in this revision

`RBP_Status = D` alone is insufficient to call a motif directly experimentally measured. The audit requires a D record and a measured motif type among RNAcompete, HTR-SELEX, RNA-Bind-n-Seq and SEQRS. JPLE is classified as inferred even for D-labelled records; RBPDB-compiled entries and homology-inferred I records are retained separately. All original Tomtom statistics and all nonsignificant comparisons remain in the table. The source classes are annotations of the frozen comparisons, not a rerun against a restricted experimental database.
