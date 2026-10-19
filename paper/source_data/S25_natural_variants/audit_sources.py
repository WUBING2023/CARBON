"""Source integrity, allele-input equivalence, numerical tie and B1 direction audit."""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='4';sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import json,gzip,hashlib,numpy as np,pandas as pd
R=Path(__file__).resolve().parent;P=R/'source_inputs'
raw=gzip.decompress((P/'EVENTS_INCLUDED.json.gz').read_bytes());events=json.loads(raw)
manifest=json.loads((R/'PACKAGED_INPUT_MANIFEST.json').read_text('utf8'))
assert hashlib.sha256(raw).hexdigest()==manifest['events_source_sha256']
for name,entry in manifest['files'].items():
    assert hashlib.sha256((R/Path(name)).read_bytes()).hexdigest()==entry['sha256'],name
z=np.load(P/'INPUTS.npz',allow_pickle=False)
shape={k:list(z[k].shape) for k in z.files}
# INPUTS contains the exact historical tensor; the paired rows differ only at the SNV base.
X=z['X_101'] if 'X_101' in z.files else z['101'] if '101' in z.files else None
if X is None:
    X=next(z[k] for k in z.files if z[k].ndim==3 and 101 in z[k].shape)
pair_diff=[]
for e in events:
    a,b=X[e['major_input_row']],X[e['minor_input_row']]
    delta=a!=b
    positions=np.flatnonzero(np.any(delta,axis=1 if a.shape[0]==101 else 0))
    assert np.array_equal(positions,[50]),e['event_id']
    channel_delta=np.any(delta,axis=0 if a.shape[0]==101 else 1)
    assert not channel_delta[-1],e['event_id']
d=pd.read_csv(R/'EVENTS_WITH_CONSERVATION_AND_PREDICTIONS.csv');rows=[];summary=[]
for cohort in ['all','strict_test']:
    f=d if cohort=='all' else d[d.split==cohort]
    for arm in ['A_native','S','A_flat','A_perm0','A_perm1','A_perm2','A_permmean']:
        vals=f[arm+'_major_minus_minor'].to_numpy();changed=(np.abs(vals)>1e-7)&(np.abs(vals)<=1e-6)
        summary.append(dict(cohort=cohort,arm=arm,new_threshold=1e-7,old_threshold=1e-6,old_ties=int(sum(np.abs(vals)<=1e-6)),new_ties=int(sum(np.abs(vals)<=1e-7)),changed_events=int(sum(changed))))
        for idx in np.flatnonzero(changed):
            e=f.iloc[idx];v=vals[idx];rows.append(dict(cohort=cohort,arm=arm,event_id=e.event_id,RBP=e.RBP,task=e.task,variant_key=e.variant_key,major_minus_minor=float(v),old_credit=.5,new_credit=float(v>0)))
pd.DataFrame(summary).to_csv(R/'TIE_THRESHOLD_SUMMARY.csv',index=False)
pd.DataFrame(rows,columns=['cohort','arm','event_id','RBP','task','variant_key','major_minus_minor','old_credit','new_credit']).to_csv(R/'TIE_THRESHOLD_CHANGED_EVENTS.csv',index=False)
b=np.load(R/'b1_reference/INTERVENTION_PACKET.npz',allow_pickle=True)
qki={k:b[k].tolist() for k in b.files if 'QKI' in k}
(R/'SOURCE_QA_REPORT.json').write_text(json.dumps(dict(events=len(events),compressed_source_hash_verified=True,packaged_files_hash_verified=True,input_array_shapes=shape,all_pairs_differ_only_at_center=True,all_pairs_have_identical_Phi=True,tie_threshold_summary=summary,QKI_source_arrays=qki,QKI_direction='neighbor to CUAAC; assay contrast CUAAC Z minus neighbor Z'),indent=2)+'\n','utf8')
print(pd.DataFrame(summary).to_string(index=False));print('Source integrity and all paired SNP/Phi input audits passed.')
