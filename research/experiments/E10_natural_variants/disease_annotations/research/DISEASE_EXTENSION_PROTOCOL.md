# Disease-associated variants: local descriptive protocol

Date: 2026-10-14. User selected Figure 5 option 2 (v12), with ClinVar and TCGA,
no disease/cancer/gene restriction. 1000 Genomes supplies population context.
This is a local analysis plan, not public preregistration or clinical validation.

## Questions and populations

1. Which of the complete existing 4,669 exact GRCh38 REF/ALT variants have current
   ClinVar clinical assertions or open TCGA somatic occurrences? Retain all 7,262
   existing ASB events and their training/validation/chromosome-8/15 provenance.
   Annotation does not change the original binding-direction denominator.
2. How do existing model responses at these annotated sites relate to measured ASB
   preference, clinical annotation, tumor occurrence and reference conservation?
   Germline clinical class and tumor recurrence are contextual annotations, not
   allele-specific RBP-binding endpoints and not evidence of an RNA disease mechanism.
3. If the old ASB population lacks pathogenic ClinVar examples, add an independent,
   annotation-first, small clinical SNV series in already represented transcribed
   genes. Study REF/ALT prediction changes and fixed-input conservation sensitivity,
   without claiming experimentally validated binding or pathogenicity prediction.

## Annotation contract

Use current official ClinVar GRCh38 snapshot and public NCI GDC TCGA aggregates;
save release/date, original replies, commands, hashes and failures. Match chromosome,
one-based coordinate and genomic-plus REF/ALT exactly. Separate germline CLNSIG,
somatic clinical impact and oncogenicity classifications. Preserve review status,
conflicting assertions and missing annotation; missing is not benign, zero frequency
or biological absence. Do not pool TCGA case frequency with population ALT frequency.
No controlled patient data, identifiers, contact or publication are requested.

## Additional clinical-series eligibility and selection, before predictions

Clinical cases: biallelic single-nucleotide substitutions; germline Pathogenic,
Likely pathogenic or the aggregate union; no conflicting classifications; at least
one review star; splice/UTR consequence; gene symbol represented among the previous
ASB gene contexts. Retain the full eligible table. Verify genomic REF, GENCODE v29
gene overlap and RNA strand, and complete 101-nt sequence/PhyloP availability before
inference. A gene-symbol intersection alone is not a valid input.

Use up to 16 clinical cases from distinct genes, sorted by descending review stars,
then SHA-256 of the exact variant key. This is an availability/budget subset, not a
model-performance selection. Seek one B/LB control in the same gene and RNA-region
class, prefer identical REF>ALT class, then nearest coordinate and exact-key tie
break. Keep unmatched clinical cases, failed inputs and the full eligible population
separately. B/LB control matching cannot use prediction magnitude or correctness.
For each gene choose at most one K562 and one HepG2 task from prior same-gene/strand
ASB coverage: descending distinct covered variants, then lexical task ID. The chosen
RBP/cell context does not imply verified occupancy at the new clinical site or that
these cancer-derived cell lines match the disease tissue.

All exact old-cohort TCGA matches are retained. A displayed TCGA case is selected by
maximum distinct aggregate cases, then old RBP coverage, then variant key, before
using prediction direction; all corresponding RBP/cell events must remain visible.

## Predictions, summaries and limits

Reuse existing, hash-verified Full and Sequence checkpoints and implemented class.
For each new REF/ALT pair only the center base changes. Full receives the same raw
reference PhyloP vector for both alleles; Sequence receives sequence alone. Also
evaluate fixed-input flat Phi and three deterministic permutations with the Full
checkpoint. These are sensitivity controls, not retrained ablations. Use eval mode,
record checkpoint hashes and inputs, and validate one previously saved pair before
the new batch. No model retraining, selection by clinical labels, test-label tuning
or pooling of unmatched model scales is permitted.

Report continuous ALT-minus-REF log-total-count responses and, where useful,
profile changes. Clinical-class summaries are descriptive, grouped first by variant
and gene; multiple RBPs/cells are not independent clinical samples. No ASB accuracy
for unlabeled new variants, no pathogenicity AUROC, no cancer-driver or causal
conservation claim. Existing ASB credits retain their frozen definition and exposure
boundaries. Explain the central reference Phi separately from model perturbations.
All original negative Full/Sequence comparisons remain.

## Resources and failure handling

Public data download ceiling 300 MB per source branch; local CPU at most four
threads; remote CPU pilot/batch at most four threads and 20 minutes; new remote
directory only, <=100 MB generated outputs, no training. Inspect server15 current
resources and inputs before execution; do not stop other jobs. GPU use is unnecessary
for this small series. Failed annotation/input/pilot is recorded and remains missing;
no synthetic or hand-filled result substitutes for unavailable evidence.

Changes prompted by source availability are documented before model outcomes.
Figure/caption/Results/Methods must distinguish current database assertions,
measured ASB preference and new exploratory model predictions.
