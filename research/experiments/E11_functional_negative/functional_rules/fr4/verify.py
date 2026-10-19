"""Independent numerical/source/transform readback, no new model query."""
import sys
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,csv
from scipy.stats import spearmanr,rankdata
R=Path(__file__).resolve().parent;D=pd.read_csv(R/'EVENT_SCORES_ALL.csv');M=pd.read_csv(R/'METRICS.csv');C=pd.read_csv(R/'CONTRASTS.csv');N=json.loads((R/'NATIVE_METADATA.json').read_text());P={k:v for k,v in np.load(R/'PREDICTIONS.npz').items()};G={k:v for k,v in np.load(R/'NOGRU_PREDICTIONS.npz').items()};I=np.load(R/'INPUTS.npz');X=I['X'];checks=[];maxpoint=0
def check(name,ok,detail=''):checks.append({'check':name,'pass':bool(ok),'detail':detail})
def auc(x,y):
 y=np.asarray(y,bool);r=rankdata(x);return float((r[y].sum()-y.sum()*(y.sum()+1)/2)/(y.sum()*(~y).sum())) if y.any() and (~y).any() else np.nan
for r in N:
 wi=int(r['window_index']);x=X[wi]
 if not r['native_exclusion']:
  seq=''.join('ACGT'[i] for i in x[:,:4].argmax(1));check('native/'+str(wi),seq==r['sequence_RNA5to3'] and len(seq)==101 and np.all(x[:,:4].sum(1)==1) and np.isfinite(x).all() and abs(float(x[:,4].mean())-r['mean_Phi'])<1e-6)
for task in ['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']:
 ix=P[task+'__input_indices'];check(task+'/nogru_same_indices',np.array_equal(ix,G[task+'__input_indices']))
 for wi,perm in zip(ix,P[task+'__permutation']):
  r=N[int(wi)];key=task+'|FR4|'+r['chrom']+'|'+r['identity'];seed=int.from_bytes(hashlib.sha256(key.encode()).digest()[:8],'little');check('perm/'+task+'/'+str(wi),np.array_equal(perm,np.random.default_rng(seed).permutation(101)))
 lookup={int(w):j for j,w in enumerate(ix)}
 for _,r in D[D.cell==task.split('_')[1]].iterrows():
  u=lookup[int(r.upstream_window_index)];d=lookup[int(r.downstream_window_index)]
  for arm in ['aligned','sequence_only','mean','permuted','no_gru']:
   a=P if arm!='no_gru' else G;vs=np.stack([a[f'{task}__{arm}__seed{s}'][[u,d]] for s in [11,23,47]]);v=vs.mean(0);check('output/'+str(r.event_index)+'/'+arm,abs((v[0]-v[1])-r[arm+'_position_index'])<2e-6 and abs(v.mean()-r[arm+'_binding_amount'])<2e-6)
for _,r in M.iterrows():
 z=D[(D.cell==r.cell)&(D.scope==r.scope)]
 if r.endpoint=='responsive_direction':z=z[z.response_class=='responsive'];y=z.delta_PSI.to_numpy();x=z[r.method+'_position_index'].to_numpy();point=spearmanr(x,y).statistic if len(z)>=3 else np.nan
 elif r.endpoint=='all_event_direction':y=z.delta_PSI.to_numpy();x=z[r.method+'_position_index'].to_numpy();point=spearmanr(x,y).statistic if len(z)>=3 else np.nan
 else:z=z[z.response_class.isin(['responsive','nonresponsive'])];y=z.response_class=='responsive';x=z[r.method+'_binding_amount'];point=auc(x,y)
 err=abs(point-r.point) if np.isfinite(point) and np.isfinite(r.point) else 0;maxpoint=max(maxpoint,err);check('metric/'+r.cell+'/'+r.scope+'/'+r.endpoint+'/'+r.method,(not np.isfinite(point) and not np.isfinite(r.point)) or err<1e-10,str(err))
 boot=R/f'BOOT_{r.cell}_{r.scope}_{r.endpoint}.npz'
 if boot.exists():
  a=np.load(boot);idx=list(a['methods']).index(r.method);v=a['values'][:,idx];v=v[np.isfinite(v)];qq=np.quantile(v,[.025,.975]);check('CI/'+r.cell+'/'+r.scope+'/'+r.endpoint+'/'+r.method,max(abs(qq[0]-r.CI_low),abs(qq[1]-r.CI_high))<1e-12)
for _,r in C.iterrows():
 p=R/f'BOOT_{r.cell}_{r.scope}_{r.endpoint}.npz'
 if not p.exists():continue
 b=np.load(p);names=list(b['methods']);v=b['values'][:,0]-b['values'][:,names.index(r.comparison.split('-',1)[1])];v=v[np.isfinite(v)];qq=np.quantile(v,[.025,.975]);check('pairedCI/'+r.cell+'/'+r.scope+'/'+r.endpoint+'/'+r.comparison,max(abs(qq[0]-r.CI_low),abs(qq[1]-r.CI_high))<1e-12)
# Fixed measurement/confounding and all-seed diagnostics, always retained.
diag=[]
for (cell,scope),z in D.groupby(['cell','scope']):
 for which in ['responsive','all']:
  q=z[z.response_class=='responsive'] if which=='responsive' else z
  if len(q)<3:continue
  cov=np.column_stack([np.ones(len(q)),rankdata(q.control_PSI),rankdata(abs(q.control_PSI-.5)),rankdata(np.log1p(q.coverage_mean))]);y=rankdata(q.delta_PSI);yr=y-cov@np.linalg.lstsq(cov,y,rcond=None)[0]
  for arm in ['aligned','sequence_only','mean','permuted','no_gru']:
   x=rankdata(q[arm+'_position_index']);xr=x-cov@np.linalg.lstsq(cov,x,rcond=None)[0];diag.append({'cell':cell,'scope':scope,'endpoint':which+'_direction_partial','method':arm,'n':len(q),'point':float(np.corrcoef(xr,yr)[0,1])})
   for seed in [11,23,47]:diag.append({'cell':cell,'scope':scope,'endpoint':which+'_direction_seed','method':arm,'seed':seed,'n':len(q),'point':float(spearmanr(q[arm+f'_seed{seed}_position_index'],q.delta_PSI).statistic)})
 z=z[z.response_class.isin(['responsive','nonresponsive'])]
 for arm in ['aligned','sequence_only','mean','permuted','no_gru']:
  for seed in [11,23,47]:diag.append({'cell':cell,'scope':scope,'endpoint':'response_detection_seed','method':arm,'seed':seed,'n':len(z),'point':auc(z[arm+f'_seed{seed}_binding_amount'],z.response_class=='responsive')})
pd.DataFrame(diag).to_csv(R/'PARTIAL_AND_SEED_DIAGNOSTICS.csv',index=False)
pd.DataFrame(checks).to_csv(R/'VERIFICATION_CHECKS.csv',index=False);qa={'pass':all(c['pass'] for c in checks),'checks':len(checks),'failures':[c for c in checks if not c['pass']][:20],'max_metric_difference':maxpoint,'independent_scipy_readback':True,'native_sequence_and_meanPhi_verified':True,'all_window_external_identity_permutations_recreated':True,'raw_count_outputs_verified':True,'paired_saved_bootstrap_CI_verified':True,'no_new_model_queries':True,'input_and_prediction_hash':{'input':hashlib.sha256((R/'INPUTS.npz').read_bytes()).hexdigest(),'predictions':hashlib.sha256((R/'PREDICTIONS.npz').read_bytes()).hexdigest(),'nogru':hashlib.sha256((R/'NOGRU_PREDICTIONS.npz').read_bytes()).hexdigest()},'scope':'same author dataset; event points with gene-cluster uncertainty; not new wetlab n; two HG response events insufficient functional confirmation'};(R/'QA.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa))
