"""Official public annotation, cached raw responses; never fill absent AF with zero."""
from pathlib import Path
import csv,json,urllib.request,urllib.parse,time,datetime,hashlib
R=Path(__file__).resolve().parent
A=R/'annotation_raw';A.mkdir(exist_ok=True)
variants=list(csv.DictReader((R/'ANNOTATION_FROZEN_VARIANTS.csv').open(encoding='utf8')))
BASE='https://rest.ensembl.org'
def request(url,payload=None):
    data=None if payload is None else json.dumps(payload).encode()
    req=urllib.request.Request(url,data=data,headers={'Content-Type':'application/json','Accept':'application/json','User-Agent':'CARBON-public-variant-research/1.0'})
    with urllib.request.urlopen(req,timeout=90) as f:return f.read(),dict(f.headers),f.status
log=[]
for path in ['info/rest','info/software','info/data']:
    name=path.replace('/','_')+'.json'
    try:
        raw,h,status=request(BASE+'/'+path+'?content-type=application/json');(A/name).write_bytes(raw)
        log.append({'url':BASE+'/'+path,'status':status,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    except Exception as e:log.append({'url':BASE+'/'+path,'error':repr(e)})
for i in range(0,len(variants),50):
    chunk=variants[i:i+50];out=A/f'vep_{i:03d}.json'
    lines=[f"{v['chrom'][3:]} {v['pos1']} . {v['ref']} {v['alt']} . . ." for v in chunk]
    payload={'variants':lines,'hgvs':1,'af':1,'af_1kg':1,'check_existing':1,'variant_class':1,'symbol':1,'canonical':1}
    (A/f'request_{i:03d}.json').write_text(json.dumps(payload,indent=2)+'\n','utf8')
    if out.exists():print('CACHE',i,flush=True);continue
    url=BASE+'/vep/human/region'
    for attempt in range(2):
        try:
            raw,h,status=request(url,payload);out.write_bytes(raw)
            log.append({'url':url,'batch_start':i,'n':len(chunk),'status':status,'headers':h,'sha256':hashlib.sha256(raw).hexdigest(),'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
            print('OK',i,len(chunk),len(raw),flush=True);break
        except Exception as e:
            log.append({'url':url,'batch_start':i,'attempt':attempt,'error':repr(e)})
            print('ERROR',i,repr(e),flush=True)
            if attempt==0:time.sleep(3)
    (A/'REQUEST_LOG.json').write_text(json.dumps(log,indent=2)+'\n','utf8')
    time.sleep(.5)
(A/'REQUEST_LOG.json').write_text(json.dumps(log,indent=2)+'\n','utf8')
