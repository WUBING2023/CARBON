from pathlib import Path
import os,sys,json,time,hashlib
os.environ['OMP_NUM_THREADS']='4';os.environ['MKL_NUM_THREADS']='4'
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent
ARMS=['aligned','sequence_only','mean','permuted'];SEEDS=[11,23,47]
METHODS=[a+'_'+s for a in ARMS for s in ['ensemble','seed11','seed23','seed47']]+['fixed_half']
SCENES=['all_core','identical7','identical7_phi05']
RNG=np.random.default_rng(20261004);B=10000
START=time.perf_counter();rows=[]
history=json.loads((ROOT/'HISTORICAL_EXPOSURE.json').read_text(encoding='utf-8'))['intervals']
spans=sorted(set((v['start'],v['end']) for v in history))
old_names=set(v['gene'] for v in history if v.get('gene'))
source=json.loads((ROOT/'remote/SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
for spec in source['tasks']:
 task=spec['task'];d=ROOT/'remote/data'/task
 with np.load(d/'data.npz') as f:z={k:f[k] for k in f.files}
 meta=json.loads((d/'metadata.json').read_text(encoding='utf-8'))
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
  exposed_locus=any(v['start']<end and v['end']>start for start,end in spans)
  exposed_gene=any(name in old_names for name in v['gene_names'])
  row={'task':task,'rbp':task.split('_')[0],'cell':task.split('_')[1],'index':v['index'],
   'start':v['start'],'end':v['end'],'strand':v['strand'],'gene_id':gene,'gene_name':v['gene_names'][0] if unique_gene else None,
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
df=pd.DataFrame(rows);df.to_csv(ROOT/'PAIR_RESULTS_ALL.csv',index=False,float_format='%.12g')
summaries=[];contrasts=[];cohorts=[];gene_rows=[];rules=[];macros=[];macro_contrasts=[];transitions=[]
bootdir=ROOT/'bootstrap';bootdir.mkdir(exist_ok=True)
for scene in SCENES:
 for cohort in ['adjusted','raw','adjusted_unexposed_gene']:
  co='raw' if cohort=='raw' else 'adjusted';good='raw_good' if co=='raw' else 'adj_good'
  mask=df[scene]&df[good]&df['unique_gene']
  if cohort.endswith('unexposed_gene'):mask &= ~df['historical_gene_exposed']&~df['historical_locus_exposed']
  sub=df[mask].copy()
  columns=[m+'_'+co+'_loss' for m in METHODS]
  g=sub.groupby(['task','gene_id'],sort=True)[columns].mean()
  assert np.isfinite(g.to_numpy()).all()
  for task in [t['task'] for t in source['tasks']]:
   sd=sub[sub['task']==task];n=len(sd)
   x=g.xs(task,level='task').to_numpy() if task in g.index.get_level_values('task') else np.empty((0,len(METHODS)))
   genes=g.xs(task,level='task').index.tolist() if len(x) else []
   cohorts.append({'scene':scene,'cohort':cohort,'task':task,'sequence_loci':int((df['task']==task).sum()),'scene_loci':int(((df['task']==task)&df[scene]).sum()),'measurable_loci':n,'genes':len(genes),'observed_ties':int(sd[co=='raw' and 'raw_tie' or 'adj_tie'].sum()),'signal_support_median':float((sd['y1']+sd['y2']).median()) if n else None})
   if not len(x):continue
   state=RNG.bit_generator.state;boots=[]
   for b in range(0,B,50):
    idx=RNG.integers(0,len(x),(min(50,B-b),len(x)));boots.append(x[idx].mean(axis=1))
   bs=np.concatenate(boots)
   np.savez_compressed(bootdir/(scene+'_'+cohort+'_'+task+'.npz'),means=bs,methods=np.array(METHODS),genes=np.array(genes),rng_state_json=json.dumps(state))
   for j,m in enumerate(METHODS):
    summaries.append({'scene':scene,'cohort':cohort,'task':task,'method':m,'loci':n,'genes':len(x),'gene_equal_mse':x[:,j].mean(),'pair_pooled_mse':sd[columns[j]].mean(),'ci95_low':np.quantile(bs[:,j],.025),'ci95_high':np.quantile(bs[:,j],.975),'gene_equal_ranking':sd.groupby('gene_id')[m+'_'+co+'_correct'].mean().mean(),'pair_ranking':sd[m+'_'+co+'_correct'].mean()})
   for seed in ['ensemble','seed11','seed23','seed47']:
    ai=METHODS.index('aligned_'+seed)
    for control in ['sequence_only','mean','permuted']:
     ci=METHODS.index(control+'_'+seed);delta=bs[:,ai]-bs[:,ci];d=x[:,ai]-x[:,ci]
     contrasts.append({'scene':scene,'cohort':cohort,'task':task,'seed':seed,'control':control,'genes':len(x),'loci':n,'gene_equal_delta':d.mean(),'ci95_low':np.quantile(delta,.025),'ci95_high':np.quantile(delta,.975),'bonf18_low':np.quantile(delta,.05/36),'bonf18_high':np.quantile(delta,1-.05/36),'fraction_genes_improved':np.mean(d<0),'fraction_genes_harmed':np.mean(d>0),'pair_pooled_delta':(sd['aligned_'+seed+'_'+co+'_loss']-sd[control+'_'+seed+'_'+co+'_loss']).mean()})
   q=sd['q_'+co if co=='raw' else 'q_adjusted'];not_tie=q!=.5;target=np.sign(q-.5)
   for rule in ['phi','center','gc','input']:
    valid=not_tie;pred=sd['rule_'+rule];corr=np.where(pred==0,.5,(pred==target).astype(float));tab=sd.loc[valid,['gene_id']].copy();tab['correct']=corr[valid]
    rules.append({'scene':scene,'cohort':cohort,'task':task,'rule':rule,'non_tied_loci':int(valid.sum()),'genes':tab['gene_id'].nunique(),'pair_ranking':tab['correct'].mean(),'gene_equal_ranking':tab.groupby('gene_id')['correct'].mean().mean(),'prediction_ties':int((pred[valid]==0).sum())})
   for gi,name in enumerate(genes):
    for j,m in enumerate(METHODS):gene_rows.append({'scene':scene,'cohort':cohort,'task':task,'gene_id':name,'method':m,'mse':x[gi,j]})
   s_ok=sd['sequence_only_ensemble_'+co+'_correct'];a_ok=sd['aligned_ensemble_'+co+'_correct']
   transitions.append({'scene':scene,'cohort':cohort,'task':task,'non_tied':int(not_tie.sum()),'S_wrong_A_right':int(((s_ok==0)&(a_ok==1)).sum()),'S_right_A_wrong':int(((s_ok==1)&(a_ok==0)).sum()),'both_right':int(((s_ok==1)&(a_ok==1)).sum()),'both_wrong':int(((s_ok==0)&(a_ok==0)).sum())})
  tasks=[t['task'] for t in source['tasks']];names=sorted(set(g.index.get_level_values('gene_id')))
  if not names:continue
  xx=np.full((len(names),len(tasks),len(METHODS)),np.nan)
  for ti,t in enumerate(tasks):
   if t not in g.index.get_level_values('task'):continue
   tab=g.xs(t,level='task')
   for ni,name in enumerate(names):
    if name in tab.index:xx[ni,ti]=tab.loc[name].to_numpy()
  point=np.nanmean(np.nanmean(xx,axis=0),axis=0) if all(np.isfinite(xx[:,:,0]).sum(0)>0) else np.full(len(METHODS),np.nan)
  state=RNG.bit_generator.state;boot=[]
  for b in range(0,B,25):
   idx=RNG.integers(0,len(xx),(min(25,B-b),len(xx)))
   sampled=xx[idx];valid=np.isfinite(sampled[:,:,:,0]);den=valid.sum(1)
   means=np.nansum(sampled,axis=1)/np.maximum(den[:,:,None],1)
   mm=means.mean(1);mm[np.any(den==0,axis=1)]=np.nan;boot.append(mm)
  bs=np.concatenate(boot);np.savez_compressed(bootdir/(scene+'_'+cohort+'_macro.npz'),means=bs,methods=np.array(METHODS),genes=np.array(names),tasks=np.array(tasks),rng_state_json=json.dumps(state))
  for j,m in enumerate(METHODS):
   macros.append({'scene':scene,'cohort':cohort,'method':m,'unique_genes':len(names),'six_task_macro_gene_mse':point[j],'pair_pooled_mse':sub[columns[j]].mean(),'valid_bootstraps':int(np.isfinite(bs[:,j]).sum())})
  for seed in ['ensemble','seed11','seed23','seed47']:
   ai=METHODS.index('aligned_'+seed)
   for control in ['sequence_only','mean','permuted']:
    ci=METHODS.index(control+'_'+seed);delta=bs[:,ai]-bs[:,ci];valid=np.isfinite(delta);d=delta[valid]
    macro_contrasts.append({'scene':scene,'cohort':cohort,'seed':seed,'control':control,'gene_equal_macro_delta':point[ai]-point[ci],'pair_pooled_delta':(sub['aligned_'+seed+'_'+co+'_loss']-sub[control+'_'+seed+'_'+co+'_loss']).mean(),'ci95_low':np.quantile(d,.025) if len(d) else None,'ci95_high':np.quantile(d,.975) if len(d) else None,'bonf3_low':np.quantile(d,.05/6) if len(d) else None,'bonf3_high':np.quantile(d,1-.05/6) if len(d) else None,'bonf9_low':np.quantile(d,.05/18) if len(d) else None,'bonf9_high':np.quantile(d,1-.05/18) if len(d) else None,'valid_bootstraps':len(d),'unique_genes':len(names)})
for name,data in [('COHORT',cohorts),('METHODS',summaries),('CONTRASTS',contrasts),('GENE_RESULTS',gene_rows),('SIMPLE_RULES',rules),('MACRO',macros),('MACRO_CONTRASTS',macro_contrasts),('TRANSITIONS',transitions)]:pd.DataFrame(data).to_csv(ROOT/(name+'.csv'),index=False,float_format='%.12g')
record={'rows':len(df),'elapsed_s':time.perf_counter()-START,'bootstrap_B':B,'rng_seed':20261004,'numpy':np.__version__,'methods':METHODS,'scene_order':SCENES,'cohort_order':['adjusted','raw','adjusted_unexposed_gene'],'task_order':[t['task'] for t in source['tasks']],'ensemble':'Average three native probability profiles in float64 before taking support shares','cluster':'Same genomic gene_id is resampled jointly across task-specific gene mean losses; each task averages its represented genes; tasks are then equally averaged. No repeated-seed biological replication.'}
(ROOT/'ANALYSIS_MANIFEST.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record,indent=2));print(pd.DataFrame(cohorts).to_string(index=False));print(pd.DataFrame(macro_contrasts).query("seed=='ensemble'").to_string(index=False))
