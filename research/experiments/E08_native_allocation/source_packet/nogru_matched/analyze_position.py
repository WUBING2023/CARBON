from pathlib import Path
import os,sys,json,time,hashlib
os.environ['OMP_NUM_THREADS']='4';sys.modules['numexpr']=None;sys.modules['bottleneck']=None
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent;O=R/'position_diagnostic';O.mkdir(exist_ok=True)
d=pd.read_csv(R/'PAIR_RESULTS_ALL.csv');d=d[d.adj_good&d.unique_gene].copy()
tasks=sorted(d.task.unique());rng=np.random.default_rng(20261007);B=10000
d['center_distance_difference']=np.abs(np.abs(d.c1-50)-np.abs(d.c2-50))
gpos=np.arange(101)
for sig in [15,25]:
 g=np.exp(-.5*((gpos-50)/sig)**2)
 m1=np.array([g[int(c)-5:int(c)+6].sum() for c in d.c1]);m2=np.array([g[int(c)-5:int(c)+6].sum() for c in d.c2])
 d[f'gaussian{sig}_share']=m1/(m1+m2);d[f'gaussian{sig}_loss']=(d[f'gaussian{sig}_share']-d.q_adjusted)**2
 d[f'gaussian{sig}_correct']=np.where(d.q_adjusted==.5,np.nan,np.where(d[f'gaussian{sig}_share']==.5,.5,(np.sign(d[f'gaussian{sig}_share']-.5)==np.sign(d.q_adjusted-.5)).astype(float)))
py={}
for t in tasks:
 for ch in ['1','X','Y']:
  src=R.parent/'remote/data'/t if ch=='1' else R.parent/'phase2_ambiguity/remote/data'/t/ch
  if not (src/'data.npz').exists():continue
  with np.load(src/'data.npz') as z:x=z['x']
  meta=json.loads((src/'metadata.json').read_text())
  for i,v in enumerate(meta):
   c1,c2=v['centers'];a=x[i,c1-5:c1+6,[1,3]].sum();b=x[i,c2-5:c2+6,[1,3]].sum()
   py[(t,ch,v['index'])]=np.sign(a-b)
d['rule_py']=[py[(v.task,str(v.chromosome).removeprefix('chr'),v.index)] for v in d.itertuples()]
methods=['aligned_ensemble','nogru_ensemble']+[f'{a}_seed{s}' for s in [11,23,47] for a in ['aligned','nogru']]
records=[];contr=[];counts=[]
for scene in ['all','equal_distance','difference_le2','difference_le5']:
 sub=d.copy() if scene=='all' else d[d.center_distance_difference.le({'equal_distance':0,'difference_le2':2,'difference_le5':5}[scene])].copy()
 m=methods+['gaussian15','gaussian25'];cols=[v+('_loss' if v.startswith('gaussian') else '_adjusted_loss') for v in m]
 tab=sub.groupby(['task','gene_id'])[cols].mean();names=sorted(set(tab.index.get_level_values('gene_id')))
 xx=np.full((len(names),len(tasks),len(m)),np.nan)
 for ti,t in enumerate(tasks):
  if t not in tab.index.get_level_values('task'):continue
  tt=tab.xs(t,level='task')
  for ni,n in enumerate(names):
   if n in tt.index:xx[ni,ti]=tt.loc[n].to_numpy()
 den=np.isfinite(xx[:,:,0]).sum(0);valid=len(names)>0 and np.all(den>0)
 point=np.nanmean(np.nanmean(xx,axis=0),axis=0) if valid else np.full(len(m),np.nan)
 state=rng.bit_generator.state;bs=[]
 if len(names):
  for b in range(0,B,50):
   ix=rng.integers(0,len(names),(min(50,B-b),len(names)));sample=xx[ix];dd=np.isfinite(sample[:,:,:,0]).sum(1)
   mm=(np.nansum(sample,axis=1)/np.maximum(dd[:,:,None],1)).mean(1);mm[np.any(dd==0,axis=1)]=np.nan;bs.append(mm)
  bs=np.concatenate(bs)
 else:bs=np.full((B,len(m)),np.nan)
 np.savez_compressed(O/(scene+'_bootstrap.npz'),means=bs,methods=np.array(m),genes=np.array(names),tasks=np.array(tasks),rng_state_json=json.dumps(state))
 counts.append(dict(scene=scene,pairs=len(sub),genes=len(names),task_gene_units=len(tab),genes_per_task=json.dumps(dict(zip(tasks,map(int,den)))),ties=int((sub.q_adjusted==.5).sum()),valid_bootstraps=int(np.isfinite(bs[:,0]).sum())))
 for j,v in enumerate(m):
  cc=v+('_correct' if v.startswith('gaussian') else '_adjusted_correct');acc=sub.groupby(['task','gene_id'])[cc].mean().groupby('task').mean().mean()
  records.append(dict(scene=scene,method=v,gene_mse=point[j],pooled_mse=sub[cols[j]].mean(),gene_ranking=acc,pooled_ranking=sub[cc].mean()))
 for a,b in [('aligned_ensemble','nogru_ensemble')]+[(f'aligned_seed{s}',f'nogru_seed{s}') for s in [11,23,47]]+[(methods[0],f'gaussian{sig}') for sig in [15,25]]:
  ai=m.index(a);bi=m.index(b);vv=bs[:,ai]-bs[:,bi];vv=vv[np.isfinite(vv)]
  rec=dict(scene=scene,full=a,control=b,delta=point[ai]-point[bi],pairs=len(sub),genes=len(names))
  for fam in [1,2,3]:
   rec[f'bonf{fam}_low']=np.quantile(vv,.05/(2*fam)) if len(vv) else None;rec[f'bonf{fam}_high']=np.quantile(vv,1-.05/(2*fam)) if len(vv) else None
  contr.append(rec)
rules=[]
for rule in ['phi','gc','center','input','py']:
 target=np.sign(d.q_adjusted-.5);s=d['rule_'+rule];nt=d.q_adjusted!=.5;tt=d.loc[nt,['task','gene_id']].copy();tt['correct']=np.where(s[nt]==0,.5,(s[nt]==target[nt]).astype(float))
 for t in tasks+['PTBP1']:
  sub=tt if t=='PTBP1' else tt[tt.task==t]
  rules.append(dict(rule=rule,task=t,non_tied_pairs=len(sub),genes=sub.gene_id.nunique(),gene_ranking=sub.groupby(['task','gene_id']).correct.mean().groupby('task').mean().mean(),pooled_ranking=sub.correct.mean()))
for file,obj in [('COHORT',counts),('METHODS',records),('CONTRASTS',contr),('SIMPLE_RULES',rules)]:pd.DataFrame(obj).to_csv(O/(file+'.csv'),index=False,float_format='%.12g')
d.to_csv(O/'PAIR_RESULTS.csv',index=False,float_format='%.12g')
(O/'MANIFEST.json').write_text(json.dumps(dict(registered_protocol_sha256=hashlib.sha256((R/'PROTOCOL_POSITION_DIAGNOSTIC.md').read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),B=B,seed=20261007,posthoc=True),indent=2))
print(pd.DataFrame(counts).to_string(index=False));print(pd.DataFrame(records).to_string(index=False));print(pd.DataFrame(contr).to_string(index=False));print(pd.DataFrame(rules).to_string(index=False))
