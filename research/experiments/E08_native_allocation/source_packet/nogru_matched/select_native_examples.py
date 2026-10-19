from pathlib import Path
import json,sys,datetime
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
import numpy as np,pandas as pd
O=Path(__file__).resolve().parent;P=O.parent
rows=[]
for range_,root in [('chr1',P),('XY',P/'phase2_ambiguity')]:
 f=pd.read_csv(root/'PAIR_RESULTS_ALL.csv')
 f=f[f['rbp']=='PTBP1'].copy();f['range']=range_
 if range_=='chr1':f['chromosome']='chr1'
 rows.append(f)
f=pd.concat(rows,ignore_index=True)
same=[]
for task in f['task'].unique():
 for chrom in ['1','X','Y']:
  d=P/'remote/data'/task if chrom=='1' else P/'phase2_ambiguity/remote/data'/task/chrom
  with np.load(d/'data.npz') as z:x=z['x'];indices=z['index'];cs=z['centers']
  for i,(idx,(c1,c2)) in enumerate(zip(indices,cs)):
   good=8<=c1<=92 and 8<=c2<=92 and np.array_equal(x[i,c1-8:c1+9,:4],x[i,c2-8:c2+9,:4]) and np.max(np.abs(x[i,c1-8:c1+9,4]-x[i,c2-8:c2+9,4]))<=1e-6
   if good:same.append({'task':task,'chromosome':'chr'+chrom,'index':int(idx)})
examples=[]
for cell in ['K562','HepG2']:
 tab=f[(f['cell']==cell)&f['unique_gene']&f['adj_good']&((f['c2']-f['c1'])>=31)].sort_values('hash')
 if len(tab):
  v=tab.iloc[0];examples.append({k:(v[k].item() if hasattr(v[k],'item') else v[k]) for k in ['task','cell','range','chromosome','index','start','end','strand','gene_id','gene_name','hash','c1','c2','y1','y2','input1','input2','q_adjusted','q_raw']})
j={'frozen_before_new_nogru_predictions':datetime.datetime.now(datetime.timezone.utc).isoformat(),'rule':'First lexicographic inputidentityhash within original measurable strict-adjusted gene-qualified longgap>=31 population in eachcell; eligibility uses existing measurement viability, no Full/NoGRU success or new model outputs. Not called pureinput-only selection. All population stats unchanged.','examples':examples,'native_identical17_phi_all_candidate_count':len(same),'native_identical17_phi_candidates':same}
(O/'EXAMPLE_SELECTION.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(j,ensure_ascii=False,indent=2))
