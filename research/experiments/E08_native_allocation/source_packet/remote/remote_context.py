"""New chr1 inference, no fitting or modification of historical files."""
import os
os.environ['PYTHONDONTWRITEBYTECODE']='1'
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
import sys,json,time,hashlib,re,ast,resource
from pathlib import Path
import numpy as np,h5py,pandas as pd,torch
R=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched')
OUT=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context')
OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(R))
from matched_common import make_model,transformed,permutations,sha
M=json.loads((R/'manifest.json').read_text())
torch.set_num_threads(4)
torch.set_num_interop_threads(1)
START=time.time()
def emit(event,**kw):
 rec={'event':event,'unix':time.time(),'elapsed_s':time.time()-START,**kw}
 print(json.dumps(rec,ensure_ascii=False),flush=True)
 with open(OUT/'events.jsonl','a',encoding='utf-8') as f:f.write(json.dumps(rec,ensure_ascii=False)+'\n')
 (OUT/'status.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2),encoding='utf-8')

def attrs(s):
 return dict(p.split('=',1) for p in s.rstrip().split(';') if '=' in p)
GFF=Path('/data/CARBON_PROJECT/数据处理/gencodeV29/gencode.v29.annotation.gff3')
genes=[]
with open(GFF,encoding='utf-8') as f:
 for line in f:
  if line.startswith('#'):continue
  row=line.rstrip().split('\t')
  if len(row)<9:continue
  if row[0]=='chr2' and genes:break
  if row[0]=='chr1' and row[2]=='gene':
   a=attrs(row[8]);genes.append({'id':a.get('gene_id',a.get('ID')),'name':a.get('gene_name'),
      'start':int(row[3])-1,'end':int(row[4]),'strand':row[6],'type':a.get('gene_type')})
assert len(genes)>100,'No usable GENCODE v29 chr1 gene mapping'
emit('gene_annotation',genes=len(genes),source=str(GFF))
source={'matched_manifest_sha256':sha(R/'manifest.json'),'protocol_sha256':sha(OUT/'PROTOCOL.md'),
 'gencode_source':str(GFF),'gene_annotation_chr1_sha256':hashlib.sha256(json.dumps(genes,sort_keys=True).encode()).hexdigest(),
 'source_scripts':{p:sha(R/p) for p in ['matched_train.py','matched_common.py','source/MARCP/models/marcp_model.py']},
 'initialization':'make_model generates seeded random MARCP; shared S weight slices are copied from same random full state. Only own matched resume checkpoint is optionally restored by trainer; no old full-model pretrained weight is loaded.',
 'tasks':[],'cpu_threads':4,'device':'cpu','torch':torch.__version__,'numpy':np.__version__}
(OUT/'genes_chr1.json').write_text(json.dumps(genes,ensure_ascii=False,indent=2),encoding='utf-8')

for task in M['tasks']:
 began=time.time()
 d=Path('/data/CARBON_PROJECT/数据/h5data')/task
 p=d/'_5ch'/(task+'_chr1_5ch.h5');raw=d/(task+'_chr1.h5')
 with h5py.File(p,'r') as h:x=np.asarray(h['X'],dtype=np.float32);y=np.asarray(h['Y'],dtype=np.float32);ctl=np.asarray(h['Y_ctl'],dtype=np.float32)
 df=pd.read_hdf(raw,key='df');df.columns=['chr','chromStart','chromEnd','peak','score','strand','signal','signal_ctl','seq']
 df=df[~df['seq'].str.contains('N',case=False)].dropna()
 assert len(df)==len(x);assert x.shape==(len(df),101,5)
 seq=[''.join('ACGT'[v] for v in row[:,:4].argmax(1)) for row in x]
 assert all(s==str(t).upper().replace('U','T') for s,t in zip(seq,df['seq']))
 word={'PTBP1':'CTCT','HNRNPK':'CCCC','RBFOX2':'TGCATG'}[task.split('_')[0]]
 eligible=[]
 for i,s in enumerate(seq):
  centers=[]
  for k in range(101-len(word)+1):
   c=k+(len(word)-1)//2
   if s[k:k+len(word)]==word and 5<=c<=95 and (not centers or c-centers[-1]>=12):centers.append(c)
  if len(centers)<2:continue
  row=df.iloc[i];c1,c2=centers[:2]
  g=[v for v in genes if v['strand']==str(row['strand']) and v['start']<=int(row['chromStart']) and v['end']>=int(row['chromEnd'])]
  eligible.append({'index':i,'chromosome':str(row['chr']),'start':int(row['chromStart']),'end':int(row['chromEnd']),
   'strand':str(row['strand']),'centers':[c1,c2],'hash':hashlib.sha256((task+'|1|'+str(i)).encode()).hexdigest(),
   'gene_ids':[v['id'] for v in g],'gene_names':[v['name'] for v in g],
   'same7':s[c1-3:c1+4]==s[c2-3:c2+4],'sevenmers':[s[c1-3:c1+4],s[c2-3:c2+4]],
   'phi_support':[float(x[i,c1-5:c1+6,4].mean()),float(x[i,c2-5:c2+6,4].mean())]})
 comps=[]
 for strand in ['+','-']:
  block=sorted([v for v in eligible if v['strand']==strand],key=lambda v:(v['start'],v['end'],v['index']));cur=[];end=-1
  for v in block:
   if cur and v['start']>=end:comps.append(cur);cur=[]
   cur.append(v);end=max(end,v['end']) if len(cur)>1 else v['end']
  if cur:comps.append(cur)
 reps=sorted([min(c,key=lambda v:v['hash']) for c in comps],key=lambda v:v['hash'])[:10000]
 idx=np.array([v['index'] for v in reps],dtype=int)
 z={'x':x[idx],'y':y[idx],'control':ctl[idx],'index':idx,'centers':np.asarray([v['centers'] for v in reps],dtype=int).reshape(-1,2)}
 pm=permutations(task,1,len(x))[idx]
 od=OUT/'data'/task;od.mkdir(parents=True,exist_ok=True)
 (od/'metadata.json').write_text(json.dumps(reps,ensure_ascii=False,indent=2),encoding='utf-8')
 checkpoint_sources=[]
 emit('selected',task=task,source_windows=len(x),sequence_eligible=len(eligible),loci=len(comps),retained=len(reps),annotation_unique=sum(len(v['gene_ids'])==1 for v in reps))
 for seed in M['seeds']:
  for arm in M['arms']:
   ts=time.time();ck=R/'outputs'/task/('seed'+str(seed))/arm/'best.pth'
   model,init=make_model(seed,arm)
   run=json.loads((ck.parent/'run.json').read_text())
   assert init['actual_initial_state_sha256']==run['actual_initial_state_sha256']
   model.load_state_dict(torch.load(ck,map_location='cpu',weights_only=True),strict=True);model.eval()
   xx=transformed(z['x'],pm,arm);pr=[];counts=[]
   with torch.no_grad():
    for st in range(0,len(xx),128):
     pred=model(torch.from_numpy(np.ascontiguousarray(xx[st:st+128])))
     pr.append(torch.softmax(pred['logits_profile'],dim=1).numpy());counts.append(pred['log_total_count'].numpy().ravel())
   prob=np.concatenate(pr) if pr else np.empty((0,101),dtype=np.float32)
   z[arm+'_seed'+str(seed)]=prob;z[arm+'_seed'+str(seed)+'_logcount']=np.concatenate(counts) if counts else np.array([])
   assert np.isfinite(prob).all() and (not len(prob) or np.allclose(prob.sum(1),1,atol=1e-6))
   elapsed=time.time()-ts
   checkpoint_sources.append({'seed':seed,'arm':arm,'path':str(ck),'sha256':sha(ck),'initial_state_sha256':init['actual_initial_state_sha256'],'prediction_seconds':elapsed,'windows':len(xx),'windows_per_s':len(xx)/max(elapsed,1e-9)})
   emit('inferred',task=task,seed=seed,arm=arm,windows=len(xx),seconds=elapsed,windows_per_s=len(xx)/max(elapsed,1e-9))
   del model
 np.savez_compressed(od/'data.npz',**z)
 source['tasks'].append({'task':task,'motif':word,'input_path':str(p),'input_sha256':sha(p),'raw_path':str(raw),'raw_sha256':sha(raw),'source_windows':len(x),'sequence_eligible':len(eligible),'retained_loci':len(reps),'checkpoints':checkpoint_sources,'data_sha256':sha(od/'data.npz'),'metadata_sha256':sha(od/'metadata.json'),'task_seconds':time.time()-began})
 (OUT/'SOURCE_MANIFEST.json').write_text(json.dumps(source,ensure_ascii=False,indent=2),encoding='utf-8')
 emit('task_complete',task=task,seconds=time.time()-began)
 del x,y,ctl,df,z
source['total_elapsed_seconds']=time.time()-START
source['cpu_user_seconds']=resource.getrusage(resource.RUSAGE_SELF).ru_utime
source['cpu_system_seconds']=resource.getrusage(resource.RUSAGE_SELF).ru_stime
(OUT/'SOURCE_MANIFEST.json').write_text(json.dumps(source,ensure_ascii=False,indent=2),encoding='utf-8')
emit('complete',seconds=time.time()-START,cpu_user_seconds=source['cpu_user_seconds'],cpu_system_seconds=source['cpu_system_seconds'])
