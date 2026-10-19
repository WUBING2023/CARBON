from pathlib import Path
import csv,json,hashlib,numpy as np
from scipy.stats import spearmanr,rankdata
R=Path(__file__).resolve().parent;P=np.load(R/'PREDICTIONS_MATCHED.npz');meta=json.loads((R/'INPUT_METADATA.json').read_text());labels=list(csv.DictReader((R/'LABELS_FUNCTION_ALL.csv').open()));res=json.loads((R/'RESULT_MATCHED.json').read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();assert sha(R/'PREDICTIONS_MATCHED.npz')==res['predictions_sha256'];assert sha(R/'INPUTS_MATCHED.npz')==res['input_sha256']
def rho(x,y):return float(spearmanr(x,y).statistic)
def partial(x,y,c):
 a=rankdata(x);b=rankdata(y);z=np.column_stack([np.ones(len(c)),rankdata(c)]);return float(np.corrcoef(a-z@np.linalg.lstsq(z,a,rcond=None)[0],b-z@np.linalg.lstsq(z,b,rcond=None)[0])[0,1])
def top_share(s,e):
 cut=np.sort(s)[-10];g=s>cut;ti=s==cut;w=g.astype(float);w[ti]=(10-g.sum())/ti.sum();return float(w@e/e.sum())
metrics=[];points=[];contr=[]
for lib in ['NucLibC','NucLibB']:
 lab=sorted([r for r in labels if r['library']==lib],key=lambda r:(int(r['position0']),r['ALT_base']));positions=np.array(sorted({int(r['position0']) for r in lab}));delta=np.array([float(r['NC_delta_repmean']) for r in lab]).reshape(len(positions),3);e=abs(delta).mean(1);l=-delta.mean(1);b=abs(np.array([float(r['HNRNPK_RIP_delta_repmean']) for r in lab]).reshape(len(positions),3)).mean(1)
 for task in ['HNRNPK_HepG2_ENCSR828ZID','HNRNPK_K562_ENCSR268ETU']:
  mm={}
  for arm in ['aligned','sequence_only','mean','permuted']:
   counts=[P[f'{task}__{arm}__seed{s}__logcounts'] for s in [11,23,47]];signed=[P[f'{task}__{arm}__seed{s}__IG_signed'] for s in [11,23,47]]
   for seed,co,ig in [(s,c,v) for s,c,v in zip([11,23,47],counts,signed)]+[('ensemble',np.mean(counts,axis=0),np.mean(signed,axis=0))]:
    allpos=np.array([int(q['position0']) for q in meta[1:]]);pred_delta=co[1:]-co[0];ss=np.array([abs(pred_delta[allpos==p]).mean() for p in positions]);imp=abs(ig[positions-4])
    for method,v in [('WT_IG',imp),('ALT_count_sensitivity',ss)]:
     r=rho(v,e);shift=[rho(np.roll(v,k),e) for k in range(1,len(v))];row={'library':lib,'task':task,'arm':arm,'seed':seed,'method':method,'n_positions':len(v),'RNA_n':1,'rho_function_E_NC':r,'rho_function_L_NC':rho(v,l),'rho_binding_E_RIP':rho(v,b),'partial_function_E_controlling_binding':partial(v,e,b),'top10_function_E_share':top_share(v,e),'circular_shift_p':(1+sum(z>=r for z in shift))/len(v)};metrics.append(row)
     if seed=='ensemble':mm[(arm,method)]=row
     for pos,value,ev,lv,bv in zip(positions,v,e,l,b):points.append({'library':lib,'task':task,'arm':arm,'seed':seed,'method':method,'position0':int(pos),'importance':float(value),'function_E_NC':float(ev),'function_L_NC':float(lv),'binding_E_RIP':float(bv)})
  for method in ['WT_IG','ALT_count_sensitivity']:
   for ctrl in ['sequence_only','mean','permuted']:
    a=mm[('aligned',method)];c=mm[(ctrl,method)];contr.append({'library':lib,'task':task,'method':method,'comparison':'aligned-'+ctrl,'delta_rho_function':a['rho_function_E_NC']-c['rho_function_E_NC'],'delta_rho_binding':a['rho_binding_E_RIP']-c['rho_binding_E_RIP'],'RNA_n':1})
def save(p,rs):
 with (R/p).open('w',encoding='utf8',newline='') as f:w=csv.DictWriter(f,list(dict.fromkeys(k for r in rs for k in r)));w.writeheader();w.writerows(rs)
save('MATCHED_METRICS_ALL.csv',metrics);save('MATCHED_POSITION_MAP_ALL.csv',points);save('MATCHED_CONTRASTS_ALL.csv',contr)
qa={'n_conditions':len(metrics),'seed_conditions_all_retained':True,'source_hashes_unchanged':True,'n_RNA':1,'ensemble_definition':'mean log counts and mean signed nucleotideIG acrossseed; ensemble importance abs(mean signedIG), equivalent to IG of average scalar, not mean absoluteIG','IG_steps':32,'IG_max_absolute_completeness_error':max(q['absolute_error'] for q in res['IG_completeness']),'IG_max_relative_error':max(q['relative_error'] for q in res['IG_completeness']),'IG_conditions_relative_error_gt_002':[q for q in res['IG_completeness'] if q['relative_error']>.02],'source_old_IG_fact':'reused old ATTRIBUTIONS were final64-step upgrade, not 32-step; source old prepare_freeze/run_ig and QUERY_LOG.json establish this. Original protocol wording corrected in report, no score rerank','script_sha256':sha(__file__)}
(R/'MATCHED_ANALYSIS_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf8')
for q in metrics:
 if q['library']=='NucLibC' and q['seed']=='ensemble':print(json.dumps(q))
