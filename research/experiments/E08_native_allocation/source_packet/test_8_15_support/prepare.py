from pathlib import Path
import subprocess,shlex,json,time,hashlib,sys
import numpy as np
O=Path(__file__).resolve().parent;R=O.parent;BASE=R.parent.parent
OLD=BASE/'two_directions_exploration_20261004/direction1/data'
code=r'''
from pathlib import Path
import json,hashlib,time
p=Path('/data/CARBON_PROJECT/数据处理/gencodeV29/gencode.v29.annotation.gff3');out=[];t=time.time()
with open(p) as f:
 for line in f:
  if line.startswith('#'):continue
  r=line.rstrip().split('\t')
  if len(r)<9 or r[0] not in ['chr8','chr15'] or r[2]!='gene':continue
  a=dict(v.split('=',1) for v in r[8].split(';') if '=' in v)
  out.append(dict(chromosome=r[0],start=int(r[3])-1,end=int(r[4]),strand=r[6],id=a.get('gene_id',a.get('ID')),name=a.get('gene_name')))
print(json.dumps(dict(genes=out,path=str(p),seconds=time.time()-t),ensure_ascii=False))
'''
cmd=['C:/Windows/System32/OpenSSH/ssh.exe','-F','/data/USER_HOME/.ssh/config','-o','BatchMode=yes','-o','ConnectTimeout=15','server15',shlex.join(['/data/CARBON_PROJECT/消融实验/全量消融运行/venv/bin/python3','-c',code])]
t=time.perf_counter();p=subprocess.run(cmd,capture_output=True,check=True);ann=json.loads(p.stdout);(O/'ANNOTATION.json').write_text(json.dumps(ann,ensure_ascii=False,indent=2),encoding='utf-8');genes=ann['genes'];manifest={'tasks':[],'old_directory':str(OLD),'scope':'Previously exposed old chr8/15 cached representatives; not a new independent primary'}
for task in ['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']:
 src=OLD/task
 with np.load(src/'data.npz') as f:z={k:f[k] for k in f.files}
 meta=json.loads((src/'metadata.json').read_text(encoding='utf-8'));assert len(meta)==len(z['x'])
 for ch in ['8','15']:
  ii=[i for i,v in enumerate(meta) if v['chromosome']=='chr'+ch];mm=[]
  for i in ii:
   v=meta[i].copy();c1,c2=v['centers'];seq=''.join('ACGT'[k] for k in z['x'][i,:,:4].argmax(1));g=[q for q in genes if q['chromosome']==v['chromosome'] and q['strand']==v['strand'] and q['start']<=v['start'] and q['end']>=v['end']]
   v.update(gene_ids=[q['id'] for q in g],gene_names=[q['name'] for q in g],same7=seq[c1-3:c1+4]==seq[c2-3:c2+4],sevenmers=[seq[c1-3:c1+4],seq[c2-3:c2+4]],phi_support=[float(z['x'][i,c-5:c+6,4].mean()) for c in [c1,c2]])
   mm.append(v)
  od=O/'remote/data'/task/ch;od.mkdir(parents=True,exist_ok=True);np.savez_compressed(od/'data.npz',**{k:v[ii] for k,v in z.items()});(od/'metadata.json').write_text(json.dumps(mm,indent=2),encoding='utf-8')
  manifest['tasks'].append(dict(task=task,chromosome=ch,windows=len(ii),original_array_sha256=hashlib.sha256((src/'data.npz').read_bytes()).hexdigest(),original_metadata_sha256=hashlib.sha256((src/'metadata.json').read_bytes()).hexdigest(),unique_gene_windows=sum(len(v['gene_ids'])==1 for v in mm)))
manifest['seconds']=time.perf_counter()-t;(O/'remote/SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
(O/'HISTORICAL_EXPOSURE.json').write_bytes((R/'HISTORY_REPLAY_ALL_CHROMOSOMES.json').read_bytes())
s=(R/'phase2_ambiguity/analyze_context.py').read_text().replace('20261005','20261009').replace("for p in ['PTBP1','HNRNPK','RBFOX2']","for p in ['PTBP1']").replace('for fam in [1,3,27,54]:','for fam in [1,3,9,27,54]:')
(O/'analyze_context.py').write_text(s,encoding='utf-8');print(json.dumps(manifest,indent=2))
