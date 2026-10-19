# External sequence-preference exploration — frozen 2026-10-04

This is a finite exploratory analysis, not a prospective preregistration. Earlier JPX/PUM2/motif results have been seen. No new external RBNS enrichment table has been scored against the models at this freeze.

## Question and two endpoints

Can sequence rules extracted from CARBON match independent in-vitro protein-binding preferences, and does aligned conservation improve this match relative to sequence-only and matched conservation controls?

Primary: recall of the experimentally strongest 20 unique 5-mers among the 20 highest aggregated model-IG 5-mers. Ties in either top-20 ranking use lexicographic A,C,G,T order; this rule is not adjusted after evaluation. Secondary: Spearman correlation of experimental enrichment with mean 5-mer IG contribution over all native occurrences, using the common observed 5-mer universe for every arm. Both endpoints are preserved; a better secondary does not replace the primary.

## Eligible proteins and external truth

Matched training exists for PTBP1, HNRNPK and RBFOX2 in K562/HepG2, seeds 11/23/47, arms aligned, sequence_only, mean and permuted. HNRNPK ENCSR368NMO and RBFOX2 ENCSR441HLP have released same-protein RBNS 5-mer enrichment tables (ENCFF715DIL and ENCFF002DFE). PTBP1 ENCSR651NGR is revoked with no files. A proper PTBP1 RNAcompete table can be added if independently authenticated; PTBP3 RBNS is not PTBP1 ground truth. Availability, not model performance, determines inclusion. All proteins passing this rule are reported.

For RBNS, use the concentration whose table column has the largest maximum enrichment across all 5-mers, following RBPNet. This choice uses only experimental labels, never model agreement. Report all concentration columns as a sensitivity analysis. File checksums, release status and accession are saved. RNAcompete 7-mer Z scores, if used, map to 5-mer scores by averaging every overlapping occurrence, including repeated occurrences within a 7-mer. Different assays are not interchangeable biological replications.

## Models, input qualification and rule extraction

Use the six previously matched tasks, all three seeds and four separately trained arms. Model source, checkpoints and data are read-only. Same architecture/initialization and training split manifest are checked; sequence-only lacks the conservation channel and has the documented input-layer parameter difference. Mean/Perm are separately trained controls, not merely transformed Full inference.

Use 1,024 native held-out windows per task, selected by smallest SHA256 of task/chromosome/source_row with salt external_preferences_20261004. Qualification: finite 101-nt arrays with unambiguous one-hot DNA at every position. No filtering on measured binding, prediction accuracy, IG strength or experimental k-mer preference. If fewer qualify, use every qualifying row and disclose the count. These windows support model rule extraction; the external binding table is the independent measurement. They are not new RNAs or wet experiments. Each arm sees exactly the same selected sequence windows.

Input arms follow matched_common transformations with the original task/chromosome/source_row permutation, preserving training definitions. Attribute only the four sequence channels conditional on that arm's conservation context. IG uses a neutral sequence baseline (0.25 per base) while holding the conservation vector unchanged at every interpolation point. Fixed scalar output is native log_total_count, midpoint integration with 32 steps, fixed evaluation mode and no dropout. IG attributions are (input minus baseline) times integrated gradient; sum four DNA-channel attributions per position. Numerical checks at 16/32/64 steps are performed on the first 32 hash-selected windows before reading external scores; any numerical repair is documented, never selected by external agreement.

For primary rules, choose the highest summed five-position contribution in each window (earliest offset tie); aggregate this summed IG by the corresponding unique 5-mer. Rank recovered 5-mers by aggregated sum. A negative maximum is retained. If fewer than 20 are recovered, recall still has denominator 20 and no unseen rule is invented. For secondary, average the summed contribution of every sliding 5-mer occurrence. All 97 offsets are counted. In both cases, T becomes U only for display and experimental identifiers are converted back to T.

Simple controls use the same native windows and 5-mer universe: native 5-mer frequency; Phi priority (max-five-position mean Phi, aggregated under the same winner rule for primary, mean occurrence Phi for secondary); and a published-known-consensus composition score specified before external scores: HNRNPK C count, PTBP1 pyrimidine count (C/T), RBFOX2 maximum exact overlap with TGCA(T/G)G consensus fragment. These are transparent weak sequence-prior baselines, not PWMs fitted on the evaluation assay. The RBFOX2 baseline specifically scores the maximum length of a contiguous substring shared with TGCATG; it does not imply a complete six-base motif fits inside a five-base rule.

## Units and comparisons

Three seeds are model fits; K562/HepG2 are training backgrounds; RBNS concentrations are assay conditions within one experiment. The independent protein is the external generalization unit. Report every seed/background, average seeds within background, then backgrounds within protein. Full-minus-Sequence/Mean/Perm differences are paired. With only two or three proteins, show exact protein values and descriptive averages; do not claim population significance from seeds, 1,024 windows, 1,024 k-mers or concentration columns. Window-cluster bootstrap quantifies conditional rule-extraction sampling stability only, and is labeled as such. The full recovered k-mer rankings and all assay columns are delivered.

## Bounded alternatives and scope

If native IG extraction fails numerical checks, fixed 1,024 synthetic 5-mers can be embedded at the central five bases of 16 preselected native contexts and ranked by mean delta log_total_count; this is a separate in-silico response analysis with the same external measurement, not simulated positive ground truth. It will require a dated protocol addendum before scoring. A conservation-stratum scene, if pursued, uses input-defined top/bottom quartiles and is not a new protein replication. At most three scenes and no new training. No new model head or k-mer length chosen because it improves external agreement.

A match to intrinsic in-vitro preferences supports learned sequence-rule correctness. It does not prove conservation is causal, explain RNA localization, or establish organismal function. Conservation-specific benefit requires Full to outperform both Sequence and matched Mean/Perm. Negative or mixed results remain in the report.
