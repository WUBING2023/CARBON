from pathlib import Path
import json,sys,shutil
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
import pandas as pd,numpy as np
R=Path(__file__).resolve().parent;O=R/'figure';S=R/'nogru_matched'
for e in ['png','pdf','svg']:
 shutil.copy2(O/('CONTEXT_CANDIDATE.'+e),O/('FIXED_COUNTEREXAMPLES.'+e))
shutil.copy2(O/'CAPTION.md',O/'FIXED_COUNTEREXAMPLES_CAPTION.md');shutil.copy2(O/'FIGURE_SOURCE.json',O/'FIXED_COUNTEREXAMPLES_SOURCE.json')
d=pd.read_csv(S/'PAIR_RESULTS_ALL.csv');d=d[d.adj_good&d.unique_gene].copy();d['improvement']=d.nogru_ensemble_adjusted_loss-d.aligned_ensemble_adjusted_loss
rec=[]
for label,mask in [('corrected',(d.aligned_ensemble_adjusted_correct==1)&(d.nogru_ensemble_adjusted_correct==0)),('worsened',(d.aligned_ensemble_adjusted_correct==0)&(d.nogru_ensemble_adjusted_correct==1))]:
 sub=d[mask].copy();med=sub.improvement.median();sub['median_distance']=(sub.improvement-med).abs();v=sub.sort_values(['median_distance','hash']).iloc[0].to_dict()
 rec.append({'group':label,'eligible_group_pairs':len(sub),'group_median_improvement':med,'selected':v})
(O/'ILLUSTRATION_SELECTION.json').write_text(json.dumps({'selection':'Post hoc two decision-transition groups, nearest median improvement and hash tie breaker; no extremal or support selection','examples':rec},indent=2),encoding='utf-8');print(json.dumps(rec,indent=2))
