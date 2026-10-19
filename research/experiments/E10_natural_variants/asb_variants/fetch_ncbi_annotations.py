"""NCBI Variation API fallback: allele-level SPDI lookup and RefSNP annotation."""
from pathlib import Path
import csv,json,urllib.request,concurrent.futures,threading,time,datetime,hashlib
R=Path(__file__).resolve().parent;A=R/'annotation_ncbi';A.mkdir(exist_ok=True)
variants=list(csv.DictReader((R/'ANNOTATION_FROZEN_VARIANTS.csv').open(encoding='utf8')))
accessions={'chr8':'NC_000008.11','chr15':'NC_000015.10'}
lock=threading.Lock();last=[0.];logs=[]
def get(url,path):
    if path.exists():return json.loads(path.read_text('utf8'))
    with lock:
        delay=max(0,.36-(time.monotonic()-last[0]));time.sleep(delay);last[0]=time.monotonic()
    req=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'CARBON-research-public-annotation/1.0'})
    with urllib.request.urlopen(req,timeout=35) as f:raw=f.read();headers=dict(f.headers)
    obj=json.loads(raw);path.write_bytes(raw)
    with lock:logs.append(dict(url=url,file=path.name,headers=headers,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sha256=hashlib.sha256(raw).hexdigest()))
    return obj
def work(pair):
    index,v=pair
    spdi=f"{accessions[v['chrom']]}:{v['pos0']}:{v['ref']}:{v['alt']}"
    r=dict(v,spdi=spdi,status='pending',rsids=[])
    try:
        obj=get('https://api.ncbi.nlm.nih.gov/variation/v0/spdi/'+spdi+'/rsids',A/f'spdi_{index:03d}.json')
        r['rsids']=obj.get('data',{}).get('rsids',[])
        r['status']='found' if r['rsids'] else 'no_refsnp_record'
        for rid in r['rsids']:
            get('https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/'+str(rid),A/f'rs{rid}.json')
    except Exception as e:r['status']='request_failed';r['error']=repr(e)
    if index%20==0:print(index,r['status'],flush=True)
    return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(work,enumerate(variants)))
(A/'LOOKUP_ALL_VARIANTS.json').write_text(json.dumps(results,indent=2)+'\n','utf8')
(A/'REQUEST_LOG.json').write_text(json.dumps(logs,indent=2)+'\n','utf8')
from collections import Counter
print(dict(Counter(x['status'] for x in results)),flush=True)
