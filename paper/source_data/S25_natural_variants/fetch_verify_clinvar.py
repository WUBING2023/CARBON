"""Direct ClinVar record confirmation using canonical SPDI; archive stale links."""
from pathlib import Path
import csv,json,urllib.request,datetime,hashlib
R=Path(__file__).resolve().parent;A=R/'annotation_clinvar';A.mkdir(exist_ok=True)
links=list(csv.DictReader((R/'NCBI_CLINVAR_LINKS.csv').open(encoding='utf8')))
ids=sorted(set(x['measure_set_id'] for x in links));objects={};log=[]
for i in range(0,len(ids),50):
    chunk=ids[i:i+50];out=A/f'esummary_{i:03d}.json';url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id='+','.join(chunk)+'&retmode=json'
    try:
        if out.exists():raw=out.read_bytes()
        else:
            with urllib.request.urlopen(url,timeout=35) as f:raw=f.read();headers=dict(f.headers)
            out.write_bytes(raw);log.append(dict(url=url,headers=headers,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sha256=hashlib.sha256(raw).hexdigest()))
        z=json.loads(raw);objects.update({k:v for k,v in z['result'].items() if k!='uids'})
    except Exception as e:log.append(dict(url=url,error=repr(e)))
rows=[]
for key in sorted(set(x['variant_key'] for x in links)):
    chrom,pos1,ref,alt=key.split(':');ac={'chr8':'NC_000008.11','chr15':'NC_000015.10'}[chrom];spdi=f'{ac}:{int(pos1)-1}:{ref}:{alt}'
    for mid in sorted(set(x['measure_set_id'] for x in links if x['variant_key']==key)):
        o=objects.get(mid,{});match=any(v.get('canonical_spdi')==spdi for v in o.get('variation_set',[]))
        classification=o.get('germline_classification',{});rows.append(dict(variant_key=key,measure_set_id=mid,query_spdi=spdi,direct_ClinVar_SPDI_match=match,VCV=o.get('accession_version',''),classification=classification.get('description',''),review_status=classification.get('review_status',''),last_evaluated=classification.get('last_evaluated',''),query_status='returned' if o else 'failed_or_not_returned',canonical_spdis=';'.join(v.get('canonical_spdi','') for v in o.get('variation_set',[]))))
with (R/'CLINVAR_DIRECT_VERIFICATION.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(A/'REQUEST_LOG.json').write_text(json.dumps(log,indent=2)+'\n','utf8')
print('Direct ClinVar matches',sum(x['direct_ClinVar_SPDI_match'] for x in rows),'of',len(rows),'records')
