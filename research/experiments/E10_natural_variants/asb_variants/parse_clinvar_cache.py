"""Replay exact ClinVar allele verification from archived JSON; no network."""
from pathlib import Path
import csv,json
R=Path(__file__).resolve().parent
links=list(csv.DictReader((R/'NCBI_CLINVAR_LINKS.csv').open(encoding='utf8')))
objects={}
for p in sorted((R/'annotation_clinvar').glob('esummary_*.json')):
    objects.update({k:v for k,v in json.loads(p.read_text('utf8'))['result'].items() if k!='uids'})
rows=[]
for key in sorted(set(x['variant_key'] for x in links)):
    chrom,pos1,ref,alt=key.split(':');ac={'chr8':'NC_000008.11','chr15':'NC_000015.10'}[chrom];spdi=f'{ac}:{int(pos1)-1}:{ref}:{alt}'
    for mid in sorted(set(x['measure_set_id'] for x in links if x['variant_key']==key)):
        o=objects.get(mid,{});match=any(v.get('canonical_spdi')==spdi for v in o.get('variation_set',[]));c=o.get('germline_classification',{})
        rows.append(dict(variant_key=key,measure_set_id=mid,query_spdi=spdi,direct_ClinVar_SPDI_match=match,VCV=o.get('accession_version',''),classification=c.get('description',''),review_status=c.get('review_status',''),last_evaluated=c.get('last_evaluated',''),query_status='returned' if o else 'failed_or_not_returned',canonical_spdis=';'.join(v.get('canonical_spdi','') for v in o.get('variation_set',[]))))
assert len(rows)==31 and all(x['direct_ClinVar_SPDI_match'] for x in rows)
with (R/'CLINVAR_DIRECT_VERIFICATION.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('Offline ClinVar SPDI matches: 31/31')
