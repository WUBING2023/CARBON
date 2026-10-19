import os,sys,json,time,base64,io,hashlib,datetime
from pathlib import Path
import numpy as np,torch
R=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/03_functional_rules/fr2')
OLD=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched')
sys.path.insert(0,str(OLD));from matched_common import make_model,transformed,permutations,sha
torch.set_num_threads(4);torch.set_num_interop_threads(1)
assert not (R/'RESULT.json').exists(),'Refuse overwriting results'
inp=np.load(R/'INPUTS.npz',allow_pickle=False);x=inp['X'];assert x.shape[1:]==(101,5) and np.isfinite(x).all()
tasks=['HNRNPK_HepG2_ENCSR828ZID','HNRNPK_K562_ENCSR268ETU'];arms=['aligned','sequence_only','mean','permuted'];seeds=[11,23,47]
start=time.time();out={};sources=[]
for task in tasks:
 pm=permutations(task+'|FR2_NucLibA',0,len(x))
 out[task+'__permutation']=pm
 for seed in seeds:
  for arm in arms:
   ck=OLD/'outputs'/task/('seed'+str(seed))/arm/'best.pth';model,init=make_model(seed,arm);run=json.loads((ck.parent/'run.json').read_text());assert init['actual_initial_state_sha256']==run['actual_initial_state_sha256'];model.load_state_dict(torch.load(ck,map_location='cpu',weights_only=True),strict=True);model.eval()
   xx=transformed(x,pm,arm);counts=[]
   with torch.no_grad():
    for st in range(0,len(xx),128):
     y=model(torch.from_numpy(np.ascontiguousarray(xx[st:st+128])));counts.append(y['log_total_count'].numpy().ravel())
   val=np.concatenate(counts);assert val.shape==(len(x),) and np.isfinite(val).all();out[task+'__'+arm+'__seed'+str(seed)]=val
   sources.append({'task':task,'seed':seed,'arm':arm,'checkpoint':str(ck),'sha256':sha(ck),'run_json_sha256':sha(ck.parent/'run.json'),'initialization':init})
   print(json.dumps({'task':task,'seed':seed,'arm':arm,'elapsed_seconds':time.time()-start}),file=sys.stderr,flush=True)
   del model
buf=io.BytesIO();np.savez_compressed(buf,**out);data=buf.getvalue()
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.time()-start,'model_queries':len(x)*24,'n_rows':len(x),'device':'CPU','threads':4,'strict_load':True,'matched_models':sources,'source_scripts_sha256':{str(p):sha(p) for p in [OLD/'matched_common.py',OLD/'matched_train.py',OLD/'source/MARCP/models/marcp_model.py',R/'remote_predict_fr2.py']},'inputs_sha256':sha(R/'INPUTS.npz'),'protocol_sha256':sha(R/'PROTOCOL_FR2.md'),'predictions_sha256':hashlib.sha256(data).hexdigest(),'output':'native log_total_count trained against log1p observed total; ensemble mean of log output, not nuclear localization probability','permutation':'row-key task|FR2_NucLibA, chromosome0, frozen qualified input row; same each seed'}
(R/'PREDICTIONS.npz').write_bytes(data);(R/'RESULT.json').write_text(json.dumps(result,indent=2))
print(json.dumps({'result':result,'arrays_base64':base64.b64encode(data).decode()}))
