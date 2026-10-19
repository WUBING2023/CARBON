# Reference-verification metadata correction

Only the final-event `reference_allele_verified` values were corrected, after
checking the saved 101-nt sequence arrays, hg38 reference base and RNA strand.
All scores, labels, identifiers and other CSV fields are unchanged.

The original bytes of both CSV files are preserved in `archived_pre_audit/`.
Historical `run_config.json` hashes refer to these original files. They have
not been rewritten to pretend the metadata repair was part of the original run.
`reference_flag_repair_manifest.json` records the old/new hashes and the
verification checks. The candidate-stage table remains unchanged.
