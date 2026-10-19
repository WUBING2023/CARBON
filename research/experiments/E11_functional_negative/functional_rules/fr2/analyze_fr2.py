from pathlib import Path
import csv,json,hashlib,datetime,time
import numpy as np
R=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(p.open(encoding='utf8')))
def num(v):
 try:v=float(v);return v if np.isfinite(v) else None
 except (ValueError,TypeError):return None
def ranks(v):
 v=np.asarray(v,float);_,ii,n=np.unique(v,return_inverse=True,return_counts=True);end=np.cumsum(n);return ((end+end-n-1)/2)[ii]
def corr(a,b):
 a=np.asarray(a,float);b=np.asarray(b,float);a=a-a.mean();b=b-b.mean();d=np.linalg.norm(a)*np.linalg.norm(b);return float(a@b/d) if d>1e-12 else None
def rho(a,b):return corr(ranks(a),ranks(b)) if len(a)>=5 else None
def partial(a,b,cs):
 x=ranks(a);y=ranks(b);z=np.column_stack([np.ones(len(x))]+[ranks(c) for c in cs]);return corr(x-z@np.linalg.lstsq(z,x,rcond=None)[0],y-z@np.linalg.lstsq(z,y,rcond=None)[0])
def write(n,rs):
 with (R/n).open('w',encoding='utf8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rs for k in r)));w.writeheader();w.writerows(rs)
start=time.perf_counter();inp=np.load(R/'INPUTS.npz',allow_pickle=False);pred=np.load(R/'PREDICTIONS.npz',allow_pickle=False);rr=read(R/'INPUT_ROWS.csv');res=json.loads((R/'RESULT.json').read_text('utf8'))
assert res['inputs_sha256']==sha(R/'INPUTS.npz') and res['predictions_sha256']==sha(R/'PREDICTIONS.npz');assert inp['construct_ids'].tolist()==[r['construct_id'] for r in rr];assert len({r['unique_key'] for r in rr})==len(rr)
freeze=json.loads((R/'INPUT_FREEZE.json').read_text('utf8'))
for path,h in freeze['source_files_sha256'].items():assert sha(Path(path))==h,path
tasks=sorted({m['task'] for m in res['matched_models']});arms=['aligned','sequence_only','mean','permuted'];base=['C_content','CCC','CCCC','CCTCCC','GCCTCCC','mean_Phi'];methods=arms+base
allrows=[];per=[];summary=[];contrasts=[];partialrows=[];place=[]
for task in tasks:
 scores={}
 for arm in arms:
  scores[arm]=np.mean([pred[task+'__'+arm+'__seed'+str(s)] for s in [11,23,47]],axis=0)
  for seed in [11,23,47]:scores[arm+'_seed'+str(seed)]=pred[task+'__'+arm+'__seed'+str(seed)]
 scores.update({b:np.array([float(r[b]) for r in rr]) for b in base})
 for i,r in enumerate(rr):allrows.append({**r,'task':task,**{m:float(v[i]) for m,v in scores.items()}})
 for cohort in ['new_RNAs_primary','all_human_source_sensitivity']:
  ids=[i for i,r in enumerate(rr) if cohort=='all_human_source_sensitivity' or r['new_RNA_primary']=='True'];genes=sorted({rr[i]['official_gene'] for i in ids})
  for ep in ['3_NC','5_NC','3_expression']:
   records=[]
   for g in genes:
    ix=[i for i in ids if rr[i]['official_gene']==g and num(rr[i][ep]) is not None];y=[float(rr[i][ep]) for i in ix];z={'task':task,'cohort':cohort,'endpoint':ep,'gene':g,'n_unique_tiles':len(ix),'RNA_eligible_min5':len(ix)>=5}
    for m,v in scores.items():z[m+'__rho']=rho(v[ix],y)
    per.append(z);records.append(z)
    if ep=='3_NC':
     cc=[i for i in ix if num(rr[i]['3_expression']) is not None]
     pp={'task':task,'cohort':cohort,'gene':g,'n':len(cc),'endpoint':ep}
     for m in arms:
      pp[m+'__partial_C_CCC']=partial(scores[m][cc],[float(rr[i][ep]) for i in cc],[[float(rr[i][b]) for i in cc] for b in ['C_content','CCC']]) if len(cc)>=8 else None
      pp[m+'__partial_expression']=partial(scores[m][cc],[float(rr[i][ep]) for i in cc],[[float(rr[i]['3_expression']) for i in cc]]) if len(cc)>=8 else None
     partialrows.append(pp)
   for m in scores:
    vv=[z[m+'__rho'] for z in records if z[m+'__rho'] is not None];a=np.asarray(vv);rng=np.random.default_rng(20261004);boot=a[rng.integers(0,len(a),(10000,len(a)))].mean(1) if len(a) else []
    summary.append({'task':task,'cohort':cohort,'endpoint':ep,'method':m,'RNA_total':len(records),'RNA_defined':len(vv),'undefined_genes':json.dumps([z['gene'] for z in records if z[m+'__rho'] is None]),'macro_rho':float(a.mean()) if len(a) else None,'CI95_low':float(np.quantile(boot,.025)) if len(a) else None,'CI95_high':float(np.quantile(boot,.975)) if len(a) else None})
   for m in ['sequence_only','mean','permuted']+base:
    pairs=[z for z in records if z['aligned__rho'] is not None and z[m+'__rho'] is not None];d=np.array([z['aligned__rho']-z[m+'__rho'] for z in pairs]);rng=np.random.default_rng(20261004);boot=d[rng.integers(0,len(d),(10000,len(d)))].mean(1) if len(d) else []
    contrasts.append({'task':task,'cohort':cohort,'endpoint':ep,'comparison':'aligned_minus_'+m,'RNA_paired':len(d),'paired_genes':json.dumps([z['gene'] for z in pairs]),'delta_rho':float(d.mean()) if len(d) else None,'CI95_low':float(np.quantile(boot,.025)) if len(d) else None,'CI95_high':float(np.quantile(boot,.975)) if len(d) else None,'Bonf3_low':float(np.quantile(boot,.05/6)) if len(d) else None,'Bonf3_high':float(np.quantile(boot,1-.05/6)) if len(d) else None,'Bonf6_two_cells_low':float(np.quantile(boot,.05/12)) if len(d) else None,'Bonf6_two_cells_high':float(np.quantile(boot,1-.05/12)) if len(d) else None,'Full_win_RNAs':int((d>0).sum()),'Full_loss_RNAs':int((d<0).sum())})
  for g in genes:
   ix=[i for i in ids if rr[i]['official_gene']==g and num(rr[i]['3_NC']) is not None and num(rr[i]['5_NC']) is not None];place.append({'task':task,'cohort':cohort,'gene':g,'n':len(ix),'measured_same_insert_3prime_vs_5prime_rho':rho([float(rr[i]['3_NC']) for i in ix],[float(rr[i]['5_NC']) for i in ix])})
write('MODEL_FUNCTION_JOIN.csv',allrows);write('PER_RNA_METRICS.csv',per);write('SUMMARY.csv',summary);write('CONTRASTS.csv',contrasts);write('PARTIAL_RANK_DIAGNOSTICS.csv',partialrows);write('PLACEMENT_REPLICATION.csv',place)
facts={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'analysis_wall_seconds':time.perf_counter()-start,'new_RNAs':freeze['RNA_new'],'primary_tiles':freeze['qualified_new_RNA_rows'],'primary_RNA_count':len(freeze['RNA_new']),'all_tiles':len(rr),'prediction_query_count':res['model_queries'],'experimental_unit':'RNA; overlapping tiles clustered; within one study','source_species_gate':'human lncRNA source design AND unique exact109human reference','functional_endpoint':'MCF7 reporter N/C; not binding or endogenous protein-dependent localization','score':'HNRNPK native predicted log-count; not N/C probability','cells_both_reported':tasks,'primary':'newRNAs 3primeNC ensemble','sha256':{str(p):sha(p) for p in [R/'INPUTS.npz',R/'PREDICTIONS.npz',R/'RESULT.json',R/'INPUT_FREEZE.json',R.parent/'PROTOCOL_FR2.md',Path(__file__)]}}
(R/'ANALYSIS.json').write_text(json.dumps(facts,indent=2),encoding='utf8')
for z in summary:
 if z['cohort']=='new_RNAs_primary' and z['endpoint']=='3_NC' and z['method'] in methods:print(json.dumps(z))
for z in contrasts:
 if z['cohort']=='new_RNAs_primary' and z['endpoint']=='3_NC' and z['comparison'] in ['aligned_minus_sequence_only','aligned_minus_mean','aligned_minus_permuted']:print(json.dumps(z))
