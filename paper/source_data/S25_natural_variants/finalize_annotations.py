"""Keep direct current ClinVar classifications distinct from older RefSNP links."""
import sys
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import pandas as pd,json,collections
R=Path(__file__).resolve().parent
a=pd.read_csv(R/'ANNOTATIONS_ALL287.csv');c=pd.read_csv(R/'CLINVAR_DIRECT_VERIFICATION.csv')
assert c.direct_ClinVar_SPDI_match.all()
agg=c.groupby('variant_key').agg(clinvar_direct_classification=('classification',lambda x:';'.join(sorted(set(x)))),
    clinvar_direct_VCV=('VCV',lambda x:';'.join(sorted(set(x)))),clinvar_direct_review=('review_status',lambda x:';'.join(sorted(set(x)))),clinvar_direct_SPDI_verified=('direct_ClinVar_SPDI_match','all')).reset_index()
a=a.merge(agg,on='variant_key',how='left',validate='one_to_one')
a['NCBI_annotation_region']='unavailable';a['NCBI_consequence_terms']=''
look=json.loads((R/'annotation_ncbi/LOOKUP_ALL_VARIANTS.json').read_text('utf8'))
for item in look:
    if not item['rsids']:continue
    k=item['variant_key'];terms=set();chrom,pos,ref,alt=k.split(':');ac={'chr8':'NC_000008.11','chr15':'NC_000015.10'}[chrom]
    target=dict(seq_id=ac,position=int(pos)-1,deleted_sequence=ref,inserted_sequence=alt)
    for rid in item['rsids']:
        p=R/f'annotation_ncbi/rs{rid}.json'
        if not p.exists():continue
        z=json.loads(p.read_text('utf8'))['primary_snapshot_data'];top=[p for p in z['placements_with_allele'] if p['seq_id']==ac and p['is_ptlp']]
        ids=[i for i,x in enumerate(top[0]['alleles']) if x['allele']['spdi']==target] if len(top)==1 else []
        if len(ids)!=1:continue
        ann=z['allele_annotations'][ids[0]]
        for q in ann.get('assembly_annotation',[]):
            if q['seq_id']!=ac:continue
            for g in q.get('genes',[]):
                terms.update(t['name'] for t in g.get('sequence_ontology',[]))
                for r in g.get('rnas',[]):terms.update(t['name'] for t in r.get('sequence_ontology',[]))
    t=';'.join(sorted(terms));a.loc[a.variant_key==k,'NCBI_consequence_terms']=t
    if any(q in terms for q in ['missense_variant','synonymous_variant','stop_gained','stop_lost','frameshift_variant','coding_sequence_variant']):region='coding'
    elif any('splice' in q for q in terms):region='splice'
    elif any('UTR' in q or 'utr' in q for q in terms):region='UTR'
    elif 'non_coding_transcript_exon_variant' in terms:region='noncoding_exon'
    elif 'intron_variant' in terms:region='intronic'
    elif terms:region='other_annotated'
    else:region='unavailable'
    a.loc[a.variant_key==k,'NCBI_annotation_region']=region
a.to_csv(R/'ANNOTATIONS_ALL287.csv',index=False)
d=pd.read_csv(R/'EVENTS_WITH_CONSERVATION_AND_PREDICTIONS.csv');s=d[d.split=='strict_test'].merge(a,on=['variant_key','chrom','pos0','pos1','ref','alt'],how='left',validate='many_to_one');s.to_csv(R/'STRICT396_WITH_ANNOTATIONS.csv',index=False)
summary=json.loads((R/'ANNOTATION_SUMMARY.json').read_text('utf8'))
summary['ClinVar_direct_exact_SPDI_records']=len(c);summary['ClinVar_current_direct_classifications']=c.classification.value_counts().to_dict()
summary['NCBI_annotation_regions']=a.NCBI_annotation_region.value_counts().to_dict();summary['HTTP404_variants']=8
summary['clinical_snapshot_note']='Direct ClinVar esummary classifications and review statuses are current retrieved records; older RefSNP RCV terms are preserved but not substituted for current VCV classifications. Source conflict/uncertainty is retained.'
(R/'ANNOTATION_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n','utf8')
print(json.dumps(summary,indent=2))
