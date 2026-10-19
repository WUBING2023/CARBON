Exploratory external functional validation, not a blinded confirmatory experiment. Public experimental directions were known before input freezing. All source-eligible constructs and loci retained, without selection by model success. Primary output is native log_total_count, ALT-WT, not calibrated experimental affinity. A versus sequence-only S and within-A flat/native spatial controls are necessary comparisons; any isolated success does not establish a general conservation advantage. 101nt natural context is primary. Literal 90/32/129nt natural sequence is a separate trained-length extrapolation sensitivity; neither includes artificial assay adapters. Same WT genomic Phi assigned to WT and all ALT. No evolved ALT Phi. One checkpoint per A/S G3BP1/HNRNPC arm; PTBP1 uses 3 genuinely separately trained seeds per cell. Genomic training overlap unknown: all 8 genes are on training chromosomes. No new model training, no labels used for optimization.
{
  "primary_block": "101_0",
  "primary_contrasts": [
    {
      "alt": "GC_to_CC",
      "reference": "WT",
      "expected_sign": 1
    },
    {
      "alt": "CC_to_CG",
      "reference": "GC_to_CC",
      "expected_sign": -1
    }
  ],
  "secondary_contrasts": [
    {
      "alt": "A2577_to_U",
      "reference": "WT",
      "expected_sign": 1
    }
  ],
  "source_measurement": "HNRNPC1 filter binding Kd nM +/- SD, 3 technical replicates; experimental strength = -ln(Kd)",
  "scope": "Native101 and literal32 sensitivity, model trained101. Same-sequence m6A chemical conditions excluded. Natural U5 core remains fixed in all4constructs. No equivalence claim WT vs CG."
}