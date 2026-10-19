"""Bounded fallback after repeat batch HTTP500: smaller requests without HGVS."""
from pathlib import Path
import csv,json,urllib.request,time,datetime,hashlib
R=Path(__file__).resolve().parent;A=R/'annotation_raw'
v=list(csv.DictReader((R/'ANNOTATION_FROZEN_VARIANTS.csv').open(encoding='utf8')))
logs=[]
for i in range(0,len(v),50):
    if (A/f'vep_{i:03d}.json').exists():continue
    for j in range(i,min(i+50,len(v)),10):
        out=A/f'vep_small_{j:03d}.json'
        if out.exists():continue
        p={'variants':[f"{x['chrom'][3:]} {x['pos1']} . {x['ref']} {x['alt']} . . ." for x in v[j:j+10]],'af':1,'af_1kg':1,'check_existing':1,'symbol':1}
        (A/f'request_small_{j:03d}.json').write_text(json.dumps(p,indent=2)+'\n','utf8')
        req=urllib.request.Request('https://rest.ensembl.org/vep/human/region',data=json.dumps(p).encode(),headers={'Content-Type':'application/json','Accept':'application/json'})
        try:
            with urllib.request.urlopen(req,timeout=45) as f:raw=f.read();headers=dict(f.headers)
            out.write_bytes(raw);logs.append(dict(batch=j,n=len(p['variants']),status=200,headers=headers,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sha256=hashlib.sha256(raw).hexdigest()));print('OK',j,len(raw),flush=True)
        except Exception as e:
            logs.append(dict(batch=j,error=repr(e),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
            if hasattr(e,'read'):(A/f'error_small_{j:03d}.txt').write_bytes(e.read())
            print('ERROR',j,repr(e),flush=True)
        (A/'SMALL_REQUEST_LOG.json').write_text(json.dumps(logs,indent=2)+'\n','utf8')
        time.sleep(.35)
