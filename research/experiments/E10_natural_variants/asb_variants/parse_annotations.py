"""Bind public annotations to exact GRCh38 REF/ALT and preserve all misses."""
import sys
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import json,csv,collections,numpy as np,pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parent;A=R/'annotation_ncbi';V=R/'annotation_raw'
lookups=json.loads((A/'LOOKUP_ALL_VARIANTS.json').read_text('utf8'))
variants=list(csv.DictReader((R/'ANNOTATION_FROZEN_VARIANTS.csv').open(encoding='utf8')))
vep={}
for p in list(V.glob('vep_*.json')):
    if 'request' in p.name:continue
    for z in json.loads(p.read_text('utf8')):
        q=z['input'].split();k=f"chr{q[0]}:{q[1]}:{q[3]}:{q[4]}";vep[k]=z
rows=[];freqrows=[];clinicalrows=[];errors=[]
for v,lookup in zip(variants,lookups):
    assert v['variant_key']==lookup['variant_key']
    target=dict(seq_id=lookup['spdi'].split(':')[0],position=int(v['pos0']),deleted_sequence=v['ref'],inserted_sequence=v['alt'])
    r=dict(v,ncbi_status=lookup['status'],rsids=';'.join('rs'+str(x) for x in lookup['rsids']),ncbi_exact_GRCh38_REF_ALT=False,
      ncbi_alt_af_1000g=None,ncbi_1000g_version=None,ncbi_1000g_alt_count=None,ncbi_1000g_total_count=None,ncbi_1000g_observed_alleles=None,
      clinvar_measure_ids='',clinvar_linked_terms='',clinvar_linked_review_status='',ensembl_success=(v['variant_key'] in vep),ensembl_exact_REF_ALT=False,
      ensembl_alt_af_1000g=None,ensembl_rsids='',ensembl_clinical_significance='',ensembl_most_severe_consequence='')
    matched_freq=[];clin=[]
    for rid in lookup['rsids']:
        p=A/f'rs{rid}.json'
        if not p.exists():continue
        obj=json.loads(p.read_text('utf8'));snapshot=obj.get('primary_snapshot_data',{})
        placements=snapshot.get('placements_with_allele',[]);annotations=snapshot.get('allele_annotations',[])
        top=[p for p in placements if p.get('is_ptlp') and p['seq_id']==target['seq_id']]
        if len(top)!=1:errors.append(dict(variant_key=v['variant_key'],reason='no_unique_primary_GRCh38_placement',rsid=rid));continue
        alleles=top[0]['alleles'];ids=[i for i,z in enumerate(alleles) if z['allele']['spdi']==target]
        if len(ids)!=1 or len(alleles)!=len(annotations):errors.append(dict(variant_key=v['variant_key'],reason='exact_allele_or_annotation_index_mismatch',rsid=rid));continue
        idx=ids[0];r['ncbi_exact_GRCh38_REF_ALT']=True
        # Each frequency observation must correspond to this allele index in a retained placement.
        equivalent=[p['alleles'][idx]['allele']['spdi'] for p in placements if len(p['alleles'])==len(alleles)]
        ann=annotations[idx]
        for f in ann.get('frequency',[]):
            if f['observation'] not in equivalent:
                errors.append(dict(variant_key=v['variant_key'],reason='frequency_observation_not_equivalent_to_exact_allele',rsid=rid,study=f['study_name']));continue
            if f['total_count']<=0:continue
            af=f['allele_count']/f['total_count'];assert 0<=af<=1
            freqrows.append(dict(variant_key=v['variant_key'],rsid=rid,study=f['study_name'],study_version=f['study_version'],allele_count=f['allele_count'],total_count=f['total_count'],ALT_AF=af,observation=json.dumps(f['observation'],sort_keys=True)))
            if f['study_name']=='1000Genomes':matched_freq.append((f['study_version'],af,f['allele_count'],f['total_count'],f['local_row_id'],annotations))
        for c in ann.get('clinical',[]):
            clin.append(c);clinicalrows.append(dict(variant_key=v['variant_key'],rsid=rid,allele_id=c.get('allele_id'),measure_set_id=c.get('measure_set_id'),RCV=c.get('accession_version'),terms=';'.join(c.get('clinical_significances',[])),review_status=c.get('review_status',''),diseases=';'.join(c.get('disease_names',[])),last_evaluated_date=c.get('last_evaluated_date','')))
    if matched_freq:
        best=max(q[0] for q in matched_freq);z=[q for q in matched_freq if q[0]==best]
        if len(set((q[1],q[2],q[3]) for q in z))==1:
            version,af,count,total,rowid,anns=z[0]
            r.update(ncbi_alt_af_1000g=af,ncbi_1000g_version=version,ncbi_1000g_alt_count=count,ncbi_1000g_total_count=total)
            counts=[ff['allele_count'] for aa in anns for ff in aa.get('frequency',[]) if ff['study_name']=='1000Genomes' and ff['study_version']==version and ff['local_row_id']==rowid]
            r['ncbi_1000g_observed_alleles']=sum(z>0 for z in counts)
        else:errors.append(dict(variant_key=v['variant_key'],reason='conflicting_1000g_frequencies'))
    r['clinvar_measure_ids']=';'.join(str(i) for i in sorted(set(c['measure_set_id'] for c in clin)))
    r['clinvar_linked_terms']=';'.join(sorted(set(t for c in clin for t in c.get('clinical_significances',[]))))
    r['clinvar_linked_review_status']=';'.join(sorted(set(c.get('review_status','') for c in clin)))
    if v['variant_key'] in vep:
        z=vep[v['variant_key']];r['ensembl_exact_REF_ALT']=z['assembly_name']=='GRCh38' and z['seq_region_name']==v['chrom'][3:] and z['start']==int(v['pos1']) and z['end']==int(v['pos1']) and z['allele_string']==v['ref']+'/'+v['alt']
        assert r['ensembl_exact_REF_ALT']
        r['ensembl_most_severe_consequence']=z.get('most_severe_consequence','')
        cols=[c for c in z.get('colocated_variants',[]) if c['start']==int(v['pos1']) and c['end']==int(v['pos1']) and v['ref'] in c.get('allele_string','').split('/') and v['alt'] in c.get('allele_string','').split('/') and c.get('strand',1)==1]
        r['ensembl_rsids']=';'.join(sorted(set(c['id'] for c in cols)))
        ff=[c['frequencies'][v['alt']]['af'] for c in cols if 'frequencies' in c and 'af' in c['frequencies'].get(v['alt'],{})]
        if len(set(ff))==1:r['ensembl_alt_af_1000g']=ff[0]
        r['ensembl_clinical_significance']=';'.join(sorted(set(t for c in cols for t in c.get('clin_sig',[]))))
    # No source-priority swapping to favor a result; NCBI exact allele primary, VEP coverage separately.
    rows.append(r)
ann=pd.DataFrame(rows)
ann['pos0']=ann.pos0.astype(int);ann['pos1']=ann.pos1.astype(int)
ann.to_csv(R/'ANNOTATIONS_ALL287.csv',index=False)
pd.DataFrame(freqrows).to_csv(R/'NCBI_ALL_ALLELE_FREQUENCIES.csv',index=False);pd.DataFrame(clinicalrows).to_csv(R/'NCBI_CLINVAR_LINKS.csv',index=False)
(R/'ANNOTATION_PARSE_ERRORS.json').write_text(json.dumps(errors,indent=2)+'\n','utf8')
events=pd.read_csv(R/'EVENTS_WITH_CONSERVATION_AND_PREDICTIONS.csv');strict=events[events.split=='strict_test'].merge(ann,on=['variant_key','chrom','pos0','pos1','ref','alt'],how='left',validate='many_to_one')
strict.to_csv(R/'STRICT396_WITH_ANNOTATIONS.csv',index=False)
unique=strict.drop_duplicates('variant_key');f=unique[unique.ncbi_alt_af_1000g.notna()];freqcorr=float(spearmanr(f.phi_center,f.ncbi_alt_af_1000g).statistic) if len(f)>2 else None
both=ann[ann.ncbi_alt_af_1000g.notna() & ann.ensembl_alt_af_1000g.notna()]
summary=dict(total_variants=len(ann),ncbi_lookup_status=ann.ncbi_status.value_counts().to_dict(),exact_NCBI_GRCh38_REF_ALT=int(ann.ncbi_exact_GRCh38_REF_ALT.sum()),
  NCBI_1000Genomes_AF_variants=int(ann.ncbi_alt_af_1000g.notna().sum()),NCBI_1000Genomes_AF_events=int(strict.ncbi_alt_af_1000g.notna().sum()),
  NCBI_linked_ClinVar_variants=int(sum(ann.clinvar_measure_ids.astype(str).str.len()>0)),
  NCBI_linked_ClinVar_terms=dict(collections.Counter(t for x in ann.clinvar_linked_terms if x for t in x.split(';'))),
  Ensembl_success=int(ann.ensembl_success.sum()),Ensembl_missing_or_failed=int(sum(~ann.ensembl_success)),Ensembl_1000g_AF_variants=int(ann.ensembl_alt_af_1000g.notna().sum()),
  both_sources_AF_variants=len(both),AF_source_difference_max=(float(np.abs(both.ncbi_alt_af_1000g-both.ensembl_alt_af_1000g).max()) if len(both) else None),
  exploratory_unique_variant_phi_vs_ALT_AF_spearman=freqcorr,parse_errors=len(errors),
  limits='Exact NCBI allele annotations are primary; unavailable population AF or ClinVar links remain missing. Linked ClinVar absence is not comprehensive proof of absence in current ClinVar. ClinVar linked terms must be directly checked. ALT AF does not measure binding effect or cancer risk. Complete exposed strict cohort retained.')
(R/'ANNOTATION_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n','utf8')
print(json.dumps(summary,indent=2))
