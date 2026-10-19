"""Registered event-weighted functional endpoint; gene-cluster uncertainty."""
import sys
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,time
from scipy.stats import spearmanr,rankdata
R=Path(__file__).resolve().parent;t=time.perf_counter();Q=pd.read_csv(R/'QUALIFICATION_ALL.csv').fillna('');N=pd.DataFrame(json.loads((R/'NATIVE_METADATA.json').read_text()));P=np.load(R/'PREDICTIONS.npz');I=np.load(R/'INPUTS.npz');res=json.loads((R/'RESULT.json').read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(R/'PREDICTIONS.npz')==res['prediction_sha256'];assert sha(R/'INPUTS.npz')==res['input_sha256']
Q['native_exclusion']='not_model_eligible_before_native_input'
for e,z in N.groupby('event_index',sort=False):
 assert len(z)==2 and set(z.side)=={'upstream','downstream'}
 Q.loc[int(e),'native_exclusion']=z.native_exclusion.iloc[0]
Q.to_csv(R/'QUALIFICATION_FINAL.csv',index=False)
E=Q[(Q.exclusion=='') & (Q.native_exclusion=='')].copy();out=[];PN=np.load(R/'NOGRU_PREDICTIONS.npz')
assert sha(R/'NOGRU_PREDICTIONS.npz')==json.loads((R/'NOGRU_RESULT.json').read_text())['prediction_sha256']
for cell,task in [('K562','PTBP1_K562_ENCSR981WKN'),('HepG2','PTBP1_HepG2_ENCSR384KAN')]:
 ix=P[task+'__input_indices'];assert np.all(I['valid'][ix]);lookup={int(w):i for i,w in enumerate(ix)}
 for _,r in E[E.cell==cell].iterrows():
  u=int(float(r.upstream_window_index));d=int(float(r.downstream_window_index));u0=lookup[u];d0=lookup[d];rec=r.to_dict()
  for arm in ['aligned','sequence_only','mean','permuted','no_gru']:
   uv=[];dv=[]
   for seed in [11,23,47]:
    v=(PN if arm=='no_gru' else P)[f'{task}__{arm}__seed{seed}'];uv.append(float(v[u0]));dv.append(float(v[d0]));rec[f'{arm}_seed{seed}_position_index']=float(v[u0]-v[d0]);rec[f'{arm}_seed{seed}_binding_amount']=float((v[u0]+v[d0])/2)
   rec[arm+'_position_index']=float(np.mean(uv)-np.mean(dv));rec[arm+'_binding_amount']=float((np.mean(uv)+np.mean(dv))/2)
  for name in ['CT_fraction','C_fraction','CTCT_count','mean_Phi']:
   uv=float(N.iloc[u][name]);dv=float(N.iloc[d][name])
   if name in ['CT_fraction','C_fraction']:uv=round(uv*101);dv=round(dv*101);rec[name+'_position_index']=(uv-dv)/101;rec[name+'_binding_amount']=(uv+dv)/202
   else:rec[name+'_position_index']=uv-dv;rec[name+'_binding_amount']=(uv+dv)/2
  out.append(rec)
D=pd.DataFrame(out);D.to_csv(R/'EVENT_SCORES_ALL.csv',index=False)
def groups(x):
 o=np.argsort(x,kind='stable');z=x[o];g=np.r_[0,np.flatnonzero(z[1:]!=z[:-1])+1];ends=np.r_[g[1:],len(x)];inv=np.repeat(np.arange(len(g)),ends-g);return o,g,inv
def wrank(x,w):
 o,g,inv=groups(x);ws=w[:,o];m=np.add.reduceat(ws,g,axis=1);rk=np.cumsum(m,axis=1)-m/2+.5;v=np.empty_like(w);v[:,o]=rk[:,inv];return v
def corr(x,y,w):
 rx=wrank(x,w);ry=wrank(y,w);s=w.sum(axis=1);xm=(rx*w).sum(axis=1)/s;ym=(ry*w).sum(axis=1)/s;xc=rx-xm[:,None];yc=ry-ym[:,None];den=np.sqrt((w*xc*xc).sum(axis=1)*(w*yc*yc).sum(axis=1));return np.divide((w*xc*yc).sum(axis=1),den,out=np.full(len(s),np.nan),where=den>0)
def auc(x,y,w):
 o,g,_=groups(x);p=np.add.reduceat((w*y[None,:])[:,o],g,axis=1);n=np.add.reduceat((w*(1-y)[None,:])[:,o],g,axis=1);den=p.sum(axis=1)*n.sum(axis=1);return np.divide((p*(np.cumsum(n,axis=1)-n/2)).sum(axis=1),den,out=np.full(len(w),np.nan),where=den>0)
methods=['aligned','sequence_only','mean','permuted','no_gru','CT_fraction','C_fraction','CTCT_count','mean_Phi'];stats=[];contr=[];diagn=[];drawmeta=[];rng=np.random.default_rng(20261004)
for cell in ['K562','HepG2']:
 for scope in ['chr1_primary','other_holdouts_secondary']:
  allz=D[(D.cell==cell)&(D.scope==scope)].copy()
  for endpoint in ['responsive_direction','all_event_direction','response_detection']:
   if endpoint=='responsive_direction':z=allz[allz.response_class=='responsive'];suffix='_position_index';y=z.delta_PSI.to_numpy(float);fn=corr
   elif endpoint=='all_event_direction':z=allz;suffix='_position_index';y=z.delta_PSI.to_numpy(float);fn=corr
   else:z=allz[allz.response_class.isin(['responsive','nonresponsive'])];suffix='_binding_amount';y=(z.response_class=='responsive').to_numpy(float);fn=auc
   n=len(z);gene=z.gene.to_numpy();gg,gi=np.unique(gene,return_inverse=True);ns=len(gg);xs=[z[m+suffix].to_numpy(float) for m in methods];point=[float(fn(x,y,np.ones((1,n)))[0]) if n and len(np.unique(y))>1 else np.nan for x in xs]
   # Never create a correlation from the two HepG2 secondary events.
   valid=n>=3 and len(np.unique(y))>1
   boots=np.full((10000,len(methods)),np.nan)
   if valid:
    state=rng.bit_generator.state
    for b in range(0,10000,200):
     gw=rng.multinomial(ns,np.full(ns,1/ns),size=min(200,10000-b));w=gw[:,gi].astype(float)
     for j,x in enumerate(xs):boots[b:b+len(w),j]=fn(x,y,w)
    np.savez_compressed(R/f'BOOT_{cell}_{scope}_{endpoint}.npz',values=boots,methods=np.array(methods),gene_ids=gg)
    drawmeta.append({'cell':cell,'scope':scope,'endpoint':endpoint,'n_events':n,'genes':ns,'RNG_state_before':state,'bootstrap':'10000 whole-gene resamples; event-weighted point statistic'})
   for j,m in enumerate(methods):
    a=boots[:,j];a=a[np.isfinite(a)];ci=np.quantile(a,[.025,.975]) if len(a) else [np.nan,np.nan]
    stats.append({'cell':cell,'scope':scope,'endpoint':endpoint,'method':m,'n_events':n,'n_genes':ns,'point':point[j] if valid else np.nan,'CI_low':ci[0],'CI_high':ci[1],'valid_bootstrap':len(a),'source_positive':int(sum(y>0)) if endpoint!='response_detection' else int(y.sum()),'source_negative':int(sum(y<0)) if endpoint!='response_detection' else int((1-y).sum())})
   for k in [1,2,3,4]:
    a=boots[:,0]-boots[:,k];a=a[np.isfinite(a)];ci=np.quantile(a,[.025,.975]) if len(a) else [np.nan,np.nan];adj=np.quantile(a,[.05/6,1-.05/6]) if len(a) else [np.nan,np.nan]
    contr.append({'cell':cell,'scope':scope,'endpoint':endpoint,'comparison':'aligned-'+methods[k],'point':point[0]-point[k] if valid else np.nan,'CI_low':ci[0],'CI_high':ci[1],'Bonf3_low':adj[0] if k<4 else np.nan,'Bonf3_high':adj[1] if k<4 else np.nan,'multiplicity_scope':'three necessary Phi contrasts' if k<4 else 'registered extra component diagnostic; unadjusted CI','valid_bootstrap':len(a),'n_events':n,'n_genes':ns})
   if endpoint=='responsive_direction' and n:
    yt=y>0;majority=max(yt.mean(),1-yt.mean())
    for m,x in zip(methods,xs):
     pred=x>0;tp=np.sum(pred&yt);tn=np.sum(~pred&~yt);fp=np.sum(pred&~yt);fneg=np.sum(~pred&yt);den=np.sqrt(float((tp+fp)*(tp+fneg)*(tn+fp)*(tn+fneg)))
     diagnostic={'cell':cell,'scope':scope,'method':m,'n':n,'positive':int(yt.sum()),'negative':int((~yt).sum()),'accuracy':float(np.mean(pred==yt)),'balanced_accuracy':float(.5*(tp/max(1,tp+fneg)+tn/max(1,tn+fp))),'MCC':float((tp*tn-fp*fneg)/den) if den else np.nan,'constant_majority_accuracy':float(majority),'constant_majority_balanced_accuracy':.5,'constant_majority_MCC':0}
     diagn.append(diagnostic)
  # Complete fixed quantile strata for secondary detection. No selected stratum becomes primary.
  z=allz[allz.response_class.isin(['responsive','nonresponsive'])].copy()
  for cov in ['control_PSI','coverage_mean']:
   z['stratum']=pd.qcut(pd.to_numeric(z[cov]),5,labels=False,duplicates='drop')
   for k,a in z.groupby('stratum'):
    yy=(a.response_class=='responsive').to_numpy(float)
    for m in methods:
     v=float(auc(a[m+'_binding_amount'].to_numpy(),yy,np.ones((1,len(a))))[0]) if len(np.unique(yy))>1 else np.nan
     diagn.append({'cell':cell,'scope':scope,'diagnostic':'fixed_quantile_detection','covariate':cov,'stratum':int(k),'method':m,'n':len(a),'positive':int(yy.sum()),'AUROC':v})
pd.DataFrame(stats).to_csv(R/'METRICS.csv',index=False);pd.DataFrame(contr).to_csv(R/'CONTRASTS.csv',index=False);pd.DataFrame(diagn).to_csv(R/'DIAGNOSTICS.csv',index=False)
den=[]
for (cell,scope),a in D.groupby(['cell','scope']):
 b=a[a.response_class=='responsive'];den.append({'cell':cell,'scope':scope,'events':len(a),'genes':a.gene.nunique(),'responsive':len(b),'responsive_genes':b.gene.nunique(),'inclusion':int((b.delta_PSI>0).sum()),'exclusion':int((b.delta_PSI<0).sum()),'nonresponsive':int((a.response_class=='nonresponsive').sum()),'intermediate':int((a.response_class=='intermediate').sum())})
pd.DataFrame(den).to_csv(R/'DENOMINATORS_FINAL.csv',index=False)
(R/'BOOTSTRAP_PROVENANCE.json').write_text(json.dumps(drawmeta,indent=2))
(R/'ANALYSIS_EXECUTION.json').write_text(json.dumps({'wall_seconds':time.perf_counter()-t,'script_sha256':sha(__file__),'pred_sha256':sha(R/'PREDICTIONS.npz'),'inputs_sha256':sha(R/'INPUTS.npz'),'quantity':'event-weighted Spearman and AUROC with paired gene-cluster resampling','majority_direction_reported':True,'outputs':[{'path':p.name,'sha256':sha(p)} for p in R.glob('*.csv')]},indent=2))
print(pd.DataFrame(stats).query("endpoint=='responsive_direction' and cell=='K562'").to_string(index=False));print(pd.DataFrame(contr).query("endpoint=='responsive_direction' and cell=='K562'").to_string(index=False))
