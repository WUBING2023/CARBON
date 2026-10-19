from pathlib import Path
import csv,json,hashlib,datetime
import numpy as np
R=Path(__file__).resolve().parent;O=R/'fr1'
def read(p):return list(csv.DictReader(p.open(encoding='utf8')))
def rank(x):
 x=np.asarray(x,float);u,ii,n=np.unique(x,return_inverse=True,return_counts=True);e=np.cumsum(n);return ((e+e-n-1)/2+1)[ii]
def rho(x,y):
 if len(x)<5:return None
 a=rank(x);b=rank(y);a-=a.mean();b-=b.mean();d=np.linalg.norm(a)*np.linalg.norm(b);return float(a@b/d) if d>1e-12 else None
def write(name,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (O/name).open('w',encoding='utf8',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
T=read(O/'MODEL_FUNCTION_JOIN.csv');axis=read(O/'EFFECT_AXIS_PER_TILE.csv')
effects={(r['construct_id'],r['protein']):-float(r['author_effect']) for r in axis}
methods=['A_native','S','A_flat','A_permutation_mean'];genes=sorted({r['gene'] for r in T});rows=[]
for task in sorted({r['task'] for r in T}):
 for gene in genes:
  rr=[r for r in T if r['task']==task and r['gene']==gene]
  for p in ['HNRNPK','SLTM','SNRNP70']:
   y=[np.mean([effects[(cid,p)] for cid in json.loads(r['construct_ids'])]) for r in rr]
   for m in methods+['C_content','CCC','CCTCCC','measured_HNRNPK_RIP','measured_SLTM_RIP']:
    rows.append(dict(task=task,gene=gene,endpoint='author_paired_median_D_'+p,method=m,n=len(rr),rho=rho([float(r[m]) for r in rr],y)))
joint=[]
for cell in ['HepG2','K562']:
 for gene in genes:
  a={r['unique_key']:r for r in T if r['task'].startswith('HNRNPK_'+cell+'_') and r['gene']==gene}
  b={r['unique_key']:r for r in T if r['task'].startswith('SLTM_'+cell+'_') and r['gene']==gene};assert set(a)==set(b) and len(a)>=5
  keys=sorted(a)
  for m in methods:
   hn=rank([a[k][m] for k in keys])/len(keys);sl=rank([b[k][m] for k in keys])/len(keys);j=np.minimum(hn,sl)
   for ep in ['D_HNRNPK','D_SLTM','baseline_NC']:
    y=[float(a[k][ep]) for k in keys]
    joint.append(dict(cell=cell,gene=gene,endpoint=ep,method=m,n=len(keys),joint_min_rho=rho(j,y),HN_single_rho=rho(hn,y),SLTM_single_rho=rho(sl,y)))
write('AUTHOR_MEDIAN_PER_RNA.csv',rows);write('JOINT_SUPPORT_PER_RNA.csv',joint)
summary=[]
for data,vfields,group in [(rows,['rho'],['task','endpoint','method']),(joint,['joint_min_rho','HN_single_rho','SLTM_single_rho'],['cell','endpoint','method'])]:
 for key in sorted({tuple(r[k] for k in group) for r in data}):
  rs=[r for r in data if tuple(r[k] for k in group)==key]
  z=dict(zip(group,key));z['RNA_total']=len(rs)
  for v in vfields:
   vv=[r[v] for r in rs if r[v] is not None];z[v+'_RNA_defined']=len(vv);z[v+'_macro']=float(np.mean(vv)) if vv else None
  summary.append(z)
write('SOURCE_AND_JOINT_SUMMARY.csv',summary)
(O/'FINAL_AUX_EXECUTION.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'new_model_queries':0,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [O/'MODEL_FUNCTION_JOIN.csv',O/'EFFECT_AXIS_PER_TILE.csv',R/'PROTOCOL_FR1_SOURCE_ADDENDUM.md',Path(__file__)]},'FR1_frozen':True,'primary_retained':True},indent=2),encoding='utf8')
for r in summary:
 if (r.get('task','').startswith('HNRNPK') and r['endpoint']=='author_paired_median_D_HNRNPK' and r['method'] in methods) or (r.get('cell') and r['endpoint']=='D_HNRNPK' and r['method'] in ['A_native','S']):print(json.dumps(r))
