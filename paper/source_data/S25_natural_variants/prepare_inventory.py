"""Freeze a descriptive extension BEFORE new annotation/association outcomes."""
from pathlib import Path
import csv, json, hashlib, datetime, collections
ROOT = Path(__file__).resolve().parent
SOURCE = Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/expanded_functional_validation_20261003')
P = SOURCE / '06_ASB_expansion'
events = json.loads((P / 'EVENTS_INCLUDED.json').read_text('utf8'))
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
rows=[]
for e in events:
    ref=e['hg38_reference_plus']; major=e['normalized_major_hg38_plus']; minor=e['normalized_minor_hg38_plus']
    assert ref in (major,minor)
    alt=minor if major==ref else major
    k=f"{e['chrom']}:{e['pos0']+1}:{ref}:{alt}"
    rows.append(dict(variant_key=k,chrom=e['chrom'],pos0=e['pos0'],pos1=e['pos0']+1,ref=ref,alt=alt,
      event_id=e['event_id'],RBP=e['RBP'],cell=e['cell'],task=e['task'],split=e['split'],
      RNA_strand=e['strand'],source_ref_is_major=(ref==major),cluster_id=e['cluster_id'],
      phi_center=e['context']['phi_raw'][50],gene_names=';'.join(e['gene_names']),
      gene_ids=';'.join(e['gene_ids']),source_excel_row=e['source_excel_row']))
with (ROOT/'COHORT_ALL7262.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
strict=[r for r in rows if r['split']=='strict_test']
variants={r['variant_key']:r for r in strict}
ordered=sorted(variants, key=lambda k:(int(variants[k]['chrom'][3:]),variants[k]['pos1'],variants[k]['ref'],variants[k]['alt']))
with (ROOT/'ANNOTATION_FROZEN_VARIANTS.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=['variant_key','chrom','pos0','pos1','ref','alt']);w.writeheader()
    w.writerows({a:variants[k][a] for a in w.fieldnames} for k in ordered)
protocol={
  'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
  'status':'local frozen descriptive protocol, not preregistration',
  'question':'How do existing natural ASB direction judgments and predicted REF/ALT perturbations vary with reference conservation and population annotation?',
  'already_exposed':'All historical ASB labels/model outputs already examined; no blind independent validation is claimed.',
  'cohort':{'all_events':len(rows),'strict_events':len(strict),'strict_loci':len(set((r['chrom'],r['pos0']) for r in strict)),
    'annotation_variants':len(ordered),'all_RBPs':len(set(r['RBP'] for r in rows)),'strict_RBPs':len(set(r['RBP'] for r in strict))},
  'eligibility':'All 7262 previously qualified BEAPR ASB-positive events retained; no outcome-dependent filtering. Annotation: ALL unique REF/ALT variants in the prior strict chr8/15 subset.',
  'coordinate_audit':'Existing whole101 unique contiguous GRCh37-to-GRCh38 and reciprocal center mapping retained. Query genome-plus REF/ALT; RNA inputs remain strand oriented.',
  'assay':'Published BEAPR eCLIP-derived significant ASB direction, FDR<0.1; normalized major/minor BASE identity, not numerical normalized binding values. Raw count ratios are not the binding effect size.',
  'primary_descriptive_outputs':['RBP-equal direction accuracy by native center PhyloP quartiles, all and strict separately',
    'Full, Sequence, fixed-input Flat-Phi, fixed-input permutation mean and constant REF preference comparator',
    'center PhyloP association with absolute predicted REF/ALT log-total-count difference for each arm',
    'population frequency and clinical-annotation coverage at allele-matched frozen variants'],
  'analysis':'Direction credits preserved: sign>1e-7 1; abs<=1e-7 .5; else0. Macro: events within task, tasks within RBP, RBPs equal. Quartiles defined from unique loci in each cohort independently, without inspecting outcomes; ties remain in the same bin. Paired RBP resampling 5000 seed20261012. Adjusted association: rank centerPhi or predicted magnitude residualized on task, GC101, GC7, local +/-2 bases, REF>ALT mutation class, strand, log1p raw total depth, and source REF preference as applicable; 2000 two-way positive-exponential gene/locus x RBP multipliers, seed20261013. Refit adjusted rank slopes with weights. These intervals describe sampling sensitivity, not new validation or causal effects.',
  'annotation':'Official Ensembl REST VEP GRCh38 queries; only exact coordinate and REF/ALT matching colocated records. 1000G AF for exact ALT (0..1); minor allele frequency min(AF,1-AF) for biallelic matched records. Missing annotation is not AF0 or benign. Clinical labels cross-check by direct NCBI public ClinVar E-utilities for returned IDs. No clinical label is a measured RBP-binding endpoint.',
  'controls_limits':'The source includes ASB positives only: no ASB detection AUROC/AP. FDR is ascertainment, not effect magnitude. No independent cognate motif-score control is available; local sequence controls are a limited proxy. ENCODE assay overlap and old unmatched Full/Sequence training remain. Flat/permutedPhi are test-time inputs to the Full checkpoint, not retrained models. Same referencePhi for both alleles; no mutant evolution measurement.',
  'TCGA':'Inspect official GDC source suitability. Tumor somatic catalogs do not supply the matched allele-specific RBP-binding endpoint; do not add unlabelled cancer mutations to the binding evaluation. No controlled individual data requested.',
  'resources':{'cpu_max':4,'annotation_loci_max':500,'gpu':'not planned; reuse real existing predictions','bootstrap_max':5000},
  'sources':{str(p):sha(p) for p in [P/'EVENTS_INCLUDED.json',P/'INPUTS.npz',P/'PREDICTIONS.npz',P/'PROTOCOL.json',SOURCE/'07_results/ASB_ALL_EVENT_RESULTS.csv']},
  'deviation_policy':'Any later changes documented before or with analysis; no search for an advantageous subset.'
}
(ROOT/'PROTOCOL_FROZEN.json').write_text(json.dumps(protocol,indent=2,ensure_ascii=False)+'\n','utf8')
(ROOT/'FREEZE_MANIFEST.json').write_text(json.dumps({'frozen_utc':protocol['frozen_utc'],'files':{n:sha(ROOT/n) for n in ['PROTOCOL_FROZEN.json','ANNOTATION_FROZEN_VARIANTS.csv','COHORT_ALL7262.csv','prepare_inventory.py']}},indent=2)+'\n','utf8')
print(json.dumps(protocol['cohort']))
