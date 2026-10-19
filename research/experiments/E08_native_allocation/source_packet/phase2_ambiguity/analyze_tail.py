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
  for group,tt in [(t,[t]) for t in tasks]+[(p,[t for t in tasks if t.startswith(p+'_')]) for p in ['PTBP1','HNRNPK','RBFOX2']]:
   sub=sd[sd['task'].isin(tt)];gg=g[g.index.get_level_values('task').isin(tt)];names=sorted(set(gg.index.get_level_values('gene_id')))
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
    for control in ['sequence_only','mean','permuted']:
     ci=METHODS.index(control+'_'+seed);v=bs[:,ai]-bs[:,ci];v=v[np.isfinite(v)]
     rec={'scene':scene,'cohort':cohort,'group':group,'seed':seed,'control':control,'gene_equal_delta':point[ai]-point[ci],'pair_pooled_delta':(sub['aligned_'+seed+'_'+co+'_loss']-sub[control+'_'+seed+'_'+co+'_loss']).mean(),'genes':len(names),'loci':len(sub),'valid_bootstraps':len(v)}
     for fam in [1,3,27,54]:
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
record={'rows':len(df),'seconds':time.perf_counter()-START,'B':B,'rng_seed':20261005,'numpy':np.__version__,'methods':METHODS,'task_order':tasks,'scenes':SCENES,'primary':'PTBP1 / ambiguous_S / adjusted; protein-level same geneID clusters jointly across cells;3controls','supporting':'protein×scene×control27; task×scene×control54. confident_S complement descriptive only, not confirmatory family.','ensemble':'Float64 average native3seed profiles before ratio','exposure':'15table conservative coordinate/gene-name ledger; not whole study blinding'}
(ROOT/'ANALYSIS_MANIFEST.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record,indent=2));c=pd.DataFrame(contrasts);print(c[(c['group']=='PTBP1')&(c['cohort']=='adjusted')&(c['seed']=='ensemble')].to_string(index=False));print(pd.DataFrame(cohorts).query("group=='PTBP1' and cohort=='adjusted'").to_string(index=False))
