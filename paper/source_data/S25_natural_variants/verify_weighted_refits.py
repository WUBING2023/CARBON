"""Bounded sensitivity: 200 real weighted nuisance refits, no endpoint selection."""
import os
os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['OMP_NUM_THREADS']='4';os.environ['MKL_NUM_THREADS']='4'
import sys
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.stats import rankdata
R=Path(__file__).resolve().parent
d=pd.read_csv(R/'EVENTS_WITH_CONSERVATION_AND_PREDICTIONS.csv');f=d[d.split=='strict_test'].reset_index(drop=True)
cols=[np.ones(len(f)),f.GC101.to_numpy(),f.GC7.to_numpy(),np.log1p(f.raw_total_depth),f.ref_preferred.to_numpy()]
for j in range(5):
    for b in 'ACG':cols.append((f.local5.str[j]==b).to_numpy().astype(float))
cat=pd.get_dummies(f[['task','strand']].astype(str),drop_first=True,dtype=float);mutation=pd.get_dummies(f.ref+'>'+f.alt,drop_first=True,dtype=float)
x=rankdata(f.phi_center);x=(x-x.mean())/x.std();X=np.column_stack([x]+cols+[cat.to_numpy(),mutation.to_numpy()])
names=[];ys=[]
for a in ['A_native','S','A_flat','A_permmean']:
    for outcome in ['abs_delta','credit']:
        y=f[a+'_'+outcome].to_numpy()
        if outcome=='abs_delta':y=rankdata(y);y=(y-y.mean())/y.std()
        names.append((a,outcome));ys.append(y)
Y=np.column_stack(ys);point=np.linalg.lstsq(X,Y,rcond=None)[0][0]
gene=pd.factorize(f.cluster_id)[0];rbp=pd.factorize(f.RBP)[0];rng=np.random.default_rng(20261014);samples=[]
for i in range(200):
    w=rng.exponential(size=gene.max()+1)[gene]*rng.exponential(size=rbp.max()+1)[rbp];sw=np.sqrt(w)
    samples.append(np.linalg.lstsq(X*sw[:,None],Y*sw[:,None],rcond=None)[0][0])
s=np.array(samples);rows=[]
for i,(a,o) in enumerate(names):
    lo,hi=np.quantile(s[:,i],[.025,.975]);rows.append(dict(cohort='strict_test',arm=a,outcome=o,point=float(point[i]),refit_ci_low=float(lo),refit_ci_high=float(hi),replicates=200,seed=20261014,zero_in_interval=bool(lo<=0<=hi)))
pd.DataFrame(rows).to_csv(R/'WEIGHTED_REFIT_SENSITIVITY_200.csv',index=False)
np.savez_compressed(R/'WEIGHTED_REFIT_SAMPLES_200.npz',samples=s,point=point)
print(pd.DataFrame(rows).to_string(index=False))
# Gene x RBP multiplier sensitivity for macro direction contrasts.
out=[]
for cohort in ['strict_test','all']:
    ff=(d if cohort=='all' else d[d.split==cohort]).reset_index(drop=True)
    tasks,ti=np.unique(ff.task,return_inverse=True);rbps,ri=np.unique(ff.RBP,return_inverse=True);genes,gi=np.unique(ff.cluster_id,return_inverse=True)
    taskrbp=np.array([np.where(rbps==ff.loc[ff.task==t,'RBP'].iloc[0])[0][0] for t in tasks]);taskn=np.bincount(taskrbp,minlength=len(rbps))
    arms=['A_native','S','A_flat','A_permmean','ref'];values=[ff.ref_preferred.to_numpy() if a=='ref' else ff[a+'_credit'].to_numpy() for a in arms]
    rng=np.random.default_rng(20261015);boot=[]
    for b in range(2000):
        wg=rng.exponential(size=len(genes))[gi];wr=rng.exponential(size=len(rbps));den=np.bincount(ti,weights=wg,minlength=len(tasks));v=[]
        for y in values:
            ta=np.bincount(ti,weights=wg*y,minlength=len(tasks))/den;ra=np.bincount(taskrbp,weights=ta,minlength=len(rbps))/taskn;v.append((ra*wr).sum()/wr.sum())
        boot.append(v)
    boot=np.asarray(boot)
    for i,a in enumerate(arms):
        lo,hi=np.quantile(boot[:,i],[.025,.975]);out.append(dict(cohort=cohort,arm=a,ci_low=lo,ci_high=hi,kind='macro direction two-way multiplier',replicates=2000))
    for i,a in enumerate(arms[1:],1):
        lo,hi=np.quantile(boot[:,0]-boot[:,i],[.025,.975]);out.append(dict(cohort=cohort,arm='A_native minus '+a,ci_low=lo,ci_high=hi,kind='paired macro contrast two-way multiplier',replicates=2000))
pd.DataFrame(out).to_csv(R/'DIRECTION_TWO_WAY_SENSITIVITY.csv',index=False)
print(pd.DataFrame(out).query("cohort=='strict_test'").to_string(index=False))
