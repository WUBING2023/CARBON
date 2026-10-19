from pathlib import Path
import os,sys,json,time,hashlib
os.environ['OMP_NUM_THREADS']='4';os.environ['MKL_NUM_THREADS']='4'
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent
ARMS=['aligned','sequence_only','mean','permuted','nogru'];SEEDS=[11,23,47]
METHODS=[a+'_'+s for a in ARMS for s in ['ensemble','seed11','seed23','seed47']]+['fixed_half']
SCENES=['all_core','long_gap','identical_local_input']
RNG=np.random.default_rng(20261006);B=10000
START=time.perf_counter();rows=[]
history=json.loads((ROOT.parent/'HISTORY_REPLAY_ALL_CHROMOSOMES.json').read_text(encoding='utf-8'))['intervals']
spans=sorted(set((v['chromosome'],v['start'],v['end']) for v in history))
old_names=set(json.loads((ROOT.parent/'HISTORY_REPLAY_ALL_CHROMOSOMES.json').read_text(encoding='utf-8'))['gene_names'])
source=json.loads((ROOT/'remote/INFERENCE_MANIFEST.json').read_text(encoding='utf-8'))
for spec in source['tasks']:
 task=spec['task'];d=ROOT.parent/'remote/data'/task if spec['chromosome']=='1' else ROOT.parent/'phase2_ambiguity/remote/data'/task/spec['chromosome']
 with np.load(d/'data.npz') as f:z={k:f[k] for k in f.files}
 meta=json.loads((d/'metadata.json').read_text(encoding='utf-8'))
 with np.load(ROOT/'remote/data'/task/spec['chromosome']/'data.npz') as f:nz={k:f[k] for k in f.files}
 assert np.array_equal(nz['index'],z['index']) and np.array_equal(nz['centers'],z['centers'])
 z.update({k:v for k,v in nz.items() if k.startswith('nogru_')})
 assert len(meta)==len(z['x'])
 profiles={a+'_seed'+str(s):z[a+'_seed'+str(s)].astype(np.float64) for a in ARMS for s in SEEDS}
 profiles.update({a+'_ensemble':np.mean([profiles[a+'_seed'+str(s)] for s in SEEDS],axis=0,dtype=np.float64) for a in ARMS})
 for i,v in enumerate(meta):
  c1,c2=v['centers'];ix1=slice(c1-5,c1+6);ix2=slice(c2-5,c2+6)
  y=z['y'][i].astype(np.float64);ctl=z['control'][i].astype(np.float64)
  Y=y.sum();C=ctl.sum();y1=y[ix1].sum();y2=y[ix2].sum();b1=ctl[ix1].sum();b2=ctl[ix2].sum()
  raw_good=y1+y2>0
  r1=(y1*C-Y*b1)/C if C>0 else np.nan;r2=(y2*C-Y*b2)/C if C>0 else np.nan
  adj_good=Y>0 and C>0 and r1>=0 and r2>=0 and r1+r2>0
  qraw=y1/(y1+y2) if raw_good else np.nan;qadj=r1/(r1+r2) if adj_good else np.nan
  unique_gene=len(v['gene_ids'])==1
  gene=v['gene_ids'][0] if unique_gene else None
  exposed_locus=any(v['chromosome']==ch and v['start']<end and v['end']>start for ch,start,end in spans)
  exposed_gene=any(name in old_names for name in v['gene_names'])
  identical_local=8<=c1<=92 and 8<=c2<=92 and np.array_equal(z['x'][i,c1-8:c1+9,:4],z['x'][i,c2-8:c2+9,:4]) and np.max(np.abs(z['x'][i,c1-8:c1+9,4]-z['x'][i,c2-8:c2+9,4]))<=1e-6
  row={'range':('chr1' if spec['chromosome']=='1' else 'XY'),'long_gap':c2-c1>=31,'identical_local_input':identical_local,'task':task,'rbp':task.split('_')[0],'cell':task.split('_')[1],'index':v['index'],
   'chromosome':v['chromosome'],'start':v['start'],'end':v['end'],'strand':v['strand'],'gene_id':gene,'gene_name':v['gene_names'][0] if unique_gene else None,
   'gene_ids_json':json.dumps(v['gene_ids']),'hash':v['hash'],'c1':c1,'c2':c2,
   'same7':v['same7'],'sevenmer1':v['sevenmers'][0],'sevenmer2':v['sevenmers'][1],
   'phi1':v['phi_support'][0],'phi2':v['phi_support'][1],
   'all_core':True,'identical7':v['same7'],'identical7_phi05':v['same7'] and abs(v['phi_support'][0]-v['phi_support'][1])>=.5,
   'signal_total':Y,'input_total':C,'y1':y1,'y2':y2,'input1':b1,'input2':b2,'residual1':r1,'residual2':r2,
   'raw_good':raw_good,'raw_tie':raw_good and qraw==.5,'adj_good':adj_good,'adj_tie':adj_good and qadj==.5,
   'unique_gene':unique_gene,'historical_locus_exposed':exposed_locus,'historical_gene_exposed':exposed_gene,
   'q_raw':qraw,'q_adjusted':qadj,'rule_phi':np.sign(v['phi_support'][0]-v['phi_support'][1]),
   'rule_center':np.sign(abs(c2-50)-abs(c1-50)),
   'rule_gc':np.sign(z['x'][i,ix1,1:3].sum()-z['x'][i,ix2,1:3].sum()),'rule_input':np.sign(b1-b2)}
  for method in METHODS:
   if method=='fixed_half':p=.5;m1=m2=np.nan
   else:
    profile=profiles[method][i];m1=profile[ix1].sum();m2=profile[ix2].sum();assert m1+m2>0 and np.isfinite(m1+m2);p=m1/(m1+m2)
   row[method+'_share']=p;row[method+'_mass1']=m1;row[method+'_mass2']=m2
   for co,q in [('raw',qraw),('adjusted',qadj)]:
    row[method+'_'+co+'_loss']=(p-q)**2 if np.isfinite(q) else np.nan
    row[method+'_'+co+'_correct']=((float(np.sign(p-.5)==np.sign(q-.5)) if p!=.5 else .5) if np.isfinite(q) and q!=.5 else np.nan)
  rows.append(row)


df=pd.DataFrame(rows)
df['ambiguous_S']=(df['sequence_only_ensemble_share']-.5).abs()<=.1
df['ambiguous_S_phi05']=df['ambiguous_S']&(df['phi1']-df['phi2']).abs().ge(.5)
df['confident_S']=~df['ambiguous_S']
for co in ['raw','adjusted']:
 q=df['q_raw' if co=='raw' else 'q_adjusted'];s=df['sequence_only_ensemble_share'];a=df['aligned_ensemble_share']
 df['shift_'+co+'_toward_term']=2*(a-s)*(q-s);df['shift_'+co+'_squared_penalty']=(a-s)**2
 df['shift_'+co+'_net_improvement']=df['shift_'+co+'_toward_term']-df['shift_'+co+'_squared_penalty']
 assert np.allclose(df['shift_'+co+'_net_improvement'],df['sequence_only_ensemble_'+co+'_loss']-df['aligned_ensemble_'+co+'_loss'],equal_nan=True,atol=1e-14)
df.to_csv(ROOT/'PAIR_RESULTS_ALL.csv',index=False,float_format='%.12g')
cohorts=[];methods=[];contrasts=[];genesout=[];simple=[];transitions=[]
bootdir=ROOT/'bootstrap';bootdir.mkdir(exist_ok=True)
tasks=list(dict.fromkeys(s['task'] for s in source['tasks']))
for scene in SCENES:
 for cohort in ['adjusted','raw','adjusted_unexposed_gene']:
  co='raw' if cohort=='raw' else 'adjusted';good='raw_good' if co=='raw' else 'adj_good'
  mask=df[scene]&df[good]&df['unique_gene']
  if cohort.endswith('unexposed_gene'):mask &= ~df['historical_gene_exposed']&~df['historical_locus_exposed']
  sd=df[mask].copy();cols=[m+'_'+co+'_loss' for m in METHODS]
  g=sd.groupby(['task','gene_id'],sort=True)[cols].mean()
  for group,tt,ranges in [(t,[t],['chr1','XY']) for t in tasks]+[(g,tasks,rs) for g,rs in [('HNRNPK',['chr1','XY']),('HNRNPK_chr1',['chr1']),('HNRNPK_XY',['XY'])]]:
   sub=sd[sd['task'].isin(tt)&sd['range'].isin(ranges)];gg=sub.groupby(['task','gene_id'],sort=True)[cols].mean();names=sorted(set(gg.index.get_level_values('gene_id')))
   xx=np.full((len(names),len(tt),len(METHODS)),np.nan)
   for ti,t in enumerate(tt):
    if t not in gg.index.get_level_values('task'):continue
    tab=gg.xs(t,level='task')
    for ni,n in enumerate(names):
     if n in tab.index:xx[ni,ti]=tab.loc[n].to_numpy()
   den=np.isfinite(xx[:,:,0]).sum(0) if len(xx) else np.zeros(len(tt));estimable=bool(len(xx) and (den>0).all())
   point=np.nanmean(np.nanmean(xx,axis=0),axis=0) if estimable else np.full(len(METHODS),np.nan)
   state=RNG.bit_generator.state;boots=[]
   if len(xx):
    for b in range(0,B,25):
     idx=RNG.integers(0,len(xx),(min(25,B-b),len(xx)));sampled=xx[idx];dd=np.isfinite(sampled[:,:,:,0]).sum(1)
     mm=(np.nansum(sampled,axis=1)/np.maximum(dd[:,:,None],1)).mean(1);mm[np.any(dd==0,axis=1)]=np.nan;boots.append(mm)
    bs=np.concatenate(boots)
   else:bs=np.full((B,len(METHODS)),np.nan)
   np.savez_compressed(bootdir/(scene+'_'+cohort+'_'+group+'.npz'),means=bs,methods=np.array(METHODS),genes=np.array(names),tasks=np.array(tt),rng_state_json=json.dumps(state))
   cohorts.append({'scene':scene,'cohort':cohort,'group':group,'pair_loci':len(sub),'genes':len(names),'task_gene_units':len(gg),'genes_per_task_json':json.dumps(dict(zip(tt,map(int,den)))),'observed_ties':int((sub['q_raw' if co=='raw' else 'q_adjusted']==.5).sum()),'valid_bootstraps':int(np.isfinite(bs[:,0]).sum()),'estimable':estimable})
   for j,m in enumerate(METHODS):
    vv=bs[:,j];vv=vv[np.isfinite(vv)]
    methods.append({'scene':scene,'cohort':cohort,'group':group,'method':m,'gene_equal_mse':point[j],'pair_pooled_mse':sub[cols[j]].mean(),'ci95_low':np.quantile(vv,.025) if len(vv) else None,'ci95_high':np.quantile(vv,.975) if len(vv) else None,'genes':len(names),'loci':len(sub),'gene_equal_ranking':sub.groupby(['task','gene_id'])[m+'_'+co+'_correct'].mean().groupby('task').mean().mean(),'pair_ranking':sub[m+'_'+co+'_correct'].mean()})
   for seed in ['ensemble','seed11','seed23','seed47']:
    ai=METHODS.index('aligned_'+seed)
    for control in ['nogru','sequence_only','mean','permuted']:
     ci=METHODS.index(control+'_'+seed);v=bs[:,ai]-bs[:,ci];v=v[np.isfinite(v)]
     rec={'scene':scene,'cohort':cohort,'group':group,'seed':seed,'control':control,'gene_equal_delta':point[ai]-point[ci],'pair_pooled_delta':(sub['aligned_'+seed+'_'+co+'_loss']-sub[control+'_'+seed+'_'+co+'_loss']).mean(),'genes':len(names),'loci':len(sub),'valid_bootstraps':len(v)}
     for fam in [1,12,36]:
      rec['bonf'+str(fam)+'_low']=np.quantile(v,.05/(2*fam)) if len(v) else None;rec['bonf'+str(fam)+'_high']=np.quantile(v,1-.05/(2*fam)) if len(v) else None
     contrasts.append(rec)
   if len(tt)==1:
    t=tt[0]
    for n in names:
     for j,m in enumerate(METHODS):genesout.append({'scene':scene,'cohort':cohort,'task':t,'gene_id':n,'method':m,'mse':gg.loc[(t,n),cols[j]]})
    q=sub['q_raw' if co=='raw' else 'q_adjusted'];nt=q!=.5;target=np.sign(q-.5)
    for rule in ['phi','gc','center','input']:
     p=sub['rule_'+rule];corr=np.where(p==0,.5,(p==target).astype(float));tab=sub.loc[nt,['gene_id']].copy();tab['correct']=corr[nt]
     simple.append({'scene':scene,'cohort':cohort,'task':t,'rule':rule,'loci':int(nt.sum()),'gene_equal_ranking':tab.groupby('gene_id')['correct'].mean().mean(),'pair_ranking':tab['correct'].mean()})
    sok=sub['sequence_only_ensemble_'+co+'_correct'];aok=sub['aligned_ensemble_'+co+'_correct']
    transitions.append({'scene':scene,'cohort':cohort,'task':t,'non_tied':int(nt.sum()),'S_wrong_A_right':int(((sok==0)&(aok==1)).sum()),'S_right_A_wrong':int(((sok==1)&(aok==0)).sum()),'both_right':int(((sok==1)&(aok==1)).sum()),'both_wrong':int(((sok==0)&(aok==0)).sum())})
for name,records in [('COHORT',cohorts),('METHODS',methods),('CONTRASTS',contrasts),('GENE_RESULTS',genesout),('SIMPLE_RULES',simple),('TRANSITIONS',transitions)]:pd.DataFrame(records).to_csv(ROOT/(name+'.csv'),index=False,float_format='%.12g')
record={'rows':len(df),'seconds':time.perf_counter()-START,'B':B,'rng_seed':20261006,'numpy':np.__version__,'methods':METHODS,'task_order':tasks,'scenes':SCENES,'primary':'HNRNPK combinedchr1XY / allcore / adjusted, Full−NoGRU;95%CI, samegeneIDclusters joint acrosscells','supporting':'3scenes×4controls Bonf12 for combinedrange; Bonf36 conservativeall3ranges. Bothranges also reportedseparately. ModelsseedsnotRNAreplicates.','ensemble':'Float64 average native3seed profiles before ratio','exposure':'15table conservative coordinate/gene-name ledger; not whole study blinding'}
(ROOT/'ANALYSIS_MANIFEST.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record,indent=2));c=pd.DataFrame(contrasts);print(c[(c['group'].str.startswith('HNRNPK'))&(c['cohort']=='adjusted')&(c['seed']=='ensemble')].to_string(index=False));print(pd.DataFrame(cohorts).query("group=='HNRNPK' and cohort=='adjusted'").to_string(index=False))
