"""Join current exact public annotations to the complete, exposed ASB archive."""
from pathlib import Path
import csv, json, hashlib, collections

ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
SOURCE=Path(__file__).resolve().parent/'source_inputs/EVENTS_WITH_CONSERVATION_AND_PREDICTIONS.csv'

def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,data):
    assert data
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)

events=rows(SOURCE)
clinical=rows(HERE/'clinvar/ALL4669_EXACT_CLINVAR_ANNOTATIONS.csv')
direct=rows(HERE/'clinvar/ALL3_PLP_DIRECT_VERIFICATION.csv')
tumor=rows(HERE/'tcga/TCGA_EXISTING_VARIANTS.csv')
projects=rows(HERE/'tcga/TCGA_VARIANT_PROJECT_COUNTS.csv')
population=rows(HERE/'population/ALL4669_POPULATION_ANNOTATIONS.csv')
cv={r['variant_key']:r for r in clinical};plp={r['variant_key']:r for r in direct}
tc={r['variant_key']:r for r in tumor};tgpp={r['variant_key']:r for r in population}
assert len(events)==7262 and len(cv)==4669 and len(plp)==3 and len(tc)==6
assert all(r['direct_canonical_spdi_match']=='True' and r['classification_matches_snapshot']=='True' for r in direct)
keys=set(plp)|set(tc)
assert len(keys)==8
included=[r for r in events if r['variant_key'] in keys]
assert len(included)==13
assert all(float(r['ref_preferred'])==1 for r in included)
# Fixed display order by genomic coordinate; no model response sorting.
order=sorted(keys,key=lambda k:(int(k.split(':')[0][3:]),int(k.split(':')[1]),k))
rank={k:i for i,k in enumerate(order)}
included.sort(key=lambda e:(rank[e['variant_key']],e['RBP'],e['cell'],e['task'],e['event_id']))
joined=[]
for e in included:
    k=e['variant_key'];d=plp.get(k,{});c=cv[k];t=tc.get(k,{})
    j=dict(e)
    j.update(database_ClinVar_PLP=int(k in plp),database_TCGA=int(k in tc),
        ClinVar_VCV=d.get('current_VCV_version',''),
        ClinVar_classification=d.get('current_germline_classification',c['CLNSIG']),
        ClinVar_review_stars=d.get('germline_review_stars',c['germline_review_stars']),
        ClinVar_reported_conditions=d.get('current_disease_names',c['CLNDN']),
        condition_assertion_scope='VCV reported conditions; no per-condition pathogenicity claim',
        TCGA_projects=t.get('tcga_projects',''),TCGA_cases=t.get('tcga_case_occurrence_count',''),
        TCGA_canonical_consequence=t.get('canonical_consequences',''),
        full_ALT_minus_REF=-float(e['A_native_ref_minus_alt']),
        sequence_ALT_minus_REF=-float(e['S_ref_minus_alt']),
        observed_ASB_preference='REF' if float(e['ref_preferred'])==1 else 'ALT')
    # Population table is joined intact in a separate table; no missing-as-zero rule.
    joined.append(j)
write(HERE/'DISEASE_ALL13_EVENT_RESULTS.csv',joined)
variant_rows=[]
for k in order:
    es=[e for e in joined if e['variant_key']==k];c=cv[k];d=plp.get(k,{});t=tc.get(k,{})
    phi={float(e['phi_center']) for e in es};assert len(phi)==1
    gene=c['GENEINFO'].split(':')[0] if k in plp else t.get('gene_symbols','')
    z=dict(variant_key=k,gene=gene,ClinVar_PLP=int(k in plp),TCGA=int(k in tc),
      ClinVar_variation_id=d.get('clinvar_variation_id',''),ClinVar_VCV=d.get('current_VCV_version',''),
      ClinVar_classification=d.get('current_germline_classification',''),
      ClinVar_reported_conditions=d.get('current_disease_names',''),
      ClinVar_review_stars=d.get('germline_review_stars',''),
      TCGA_projects=t.get('tcga_projects',''),TCGA_cases=t.get('tcga_case_occurrence_count',''),
      TCGA_canonical_consequence=t.get('canonical_consequences',''),
      ASB_events=len(es),RBPs=';'.join(sorted({e['RBP'] for e in es})),
      splits=';'.join(sorted({e['split'] for e in es})),reference_PhyloP=next(iter(phi)),
      full_direction_matches=sum(float(e['A_native_credit']) for e in es),
      sequence_direction_matches=sum(float(e['S_credit']) for e in es))
    z.update({'population_'+p:v for p,v in tgpp[k].items() if p!='variant_key'})
    variant_rows.append(z)
write(HERE/'DISEASE_ALL8_VARIANT_RESULTS.csv',variant_rows)
case=[e for e in joined if e['variant_key']=='chr9:108880051:C:T']
assert len(case)==2 and {e['RBP'] for e in case}=={'AQR','PRPF8'}
assert all(e['strand']=='-' and e['reference_sequence_rna_dna'][50]=='G' and e['minor_sequence_rna_dna'][50]=='A' for e in case)
write(HERE/'ELP1_CLINICAL_CASE_ALL2.csv',case)
summary=dict(scope='Descriptive exact disease annotation, existing development-exposed ASB predictions; not independent clinical or tumor binding validation',
 all_events=7262,all_variants=4669,ClinVar_exact_variants=520,ClinVar_missing=4149,
 ClinVar_groups=dict(collections.Counter(c['germline_group'] for c in clinical)),
 ClinVar_PLP_variants=3,ClinVar_PLP_events=sum(e['variant_key'] in plp for e in joined),
 TCGA_variants=6,TCGA_events=sum(e['variant_key'] in tc for e in joined),
 union_variants=8,union_events=13,union_RBPs=len({e['RBP'] for e in joined}),
 union_splits=dict(collections.Counter(e['split'] for e in joined)),
 union_full_direction_matches=sum(float(e['A_native_credit']) for e in joined),
 union_sequence_direction_matches=sum(float(e['S_credit']) for e in joined),
 all_TGCA_single_cases=all(int(t['tcga_case_occurrence_count'])==1 for t in tumor),
 strict_chr8_15_disease_variants=0,
 example_selection='Clinical RNA splice consequence first, known reported condition and >=2 review stars; no model response criterion. ELP1 is the sole splice-donor site among the three ClinVar P/LP matches.',
 no_new_training_or_inference=True,
 additional_ClinVar_pool='Availability inventory only; not evaluated because exact old cohort supplies clinical cases.',
 analysis_units='Exact variant for source coverage; event=variant/RBP/cell for descriptive binding direction; no inferential significance test.',
 primary_boundaries=['VCV aggregate classification is not per-condition proof.',
   'TCGA singleton occurrence is not recurrent cancer driver or a locus-callable frequency.',
   'All13 ASB labels/predictions were previously exposed; none belongs to chr8/15 subset.',
   'Reference PhyloP is unchanged across alleles; no mutant evolutionary score.',
   'ClinVar/TCGA classification is not a RBP-binding endpoint; no disease AUROC or RNA-causality claim.'])
(HERE/'DISEASE_SUMMARY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
inputs=[SOURCE,HERE/'clinvar/ALL4669_EXACT_CLINVAR_ANNOTATIONS.csv',HERE/'clinvar/ALL3_PLP_DIRECT_VERIFICATION.csv',
 HERE/'tcga/TCGA_EXISTING_VARIANTS.csv',HERE/'tcga/TCGA_VARIANT_PROJECT_COUNTS.csv',HERE/'population/ALL4669_POPULATION_ANNOTATIONS.csv']
(HERE/'DISEASE_ANALYSIS_MANIFEST.json').write_text(json.dumps({'inputs':[{'path':str(p),'sha256':sha(p)} for p in inputs],
 'script':str(Path(__file__).resolve()),'script_sha256':sha(Path(__file__)),
 'complete_exact_cohort':True,'new_model_run':False,'output_variants':8,'output_events':13},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
