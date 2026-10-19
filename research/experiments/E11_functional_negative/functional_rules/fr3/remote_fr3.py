import sys,time,io,json,base64,datetime,hashlib
from pathlib import Path
import numpy as np,torch
R=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/03_functional_rules/fr3');M=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched')
sys.path.insert(0,str(M));from matched_common import make_model,transformed,permutations,sha
torch.set_num_threads(4);torch.set_num_interop_threads(1);assert not (R/'RESULT_MATCHED.json').exists()
x=np.load(R/'INPUTS_MATCHED.npz')['X'];assert x.shape==(151,101,5) and np.all(x[:,:,4]==x[0,:,4]);start=time.time();out={};sources=[];checks=[]
for task in ['HNRNPK_HepG2_ENCSR828ZID','HNRNPK_K562_ENCSR268ETU']:
 pm=np.repeat(permutations(task+'|FR3_JPX_WT',0,1),151,axis=0)
 out[task+'__permutation']=pm[0]
 for seed in [11,23,47]:
  for arm in ['aligned','sequence_only','mean','permuted']:
   ck=M/'outputs'/task/('seed'+str(seed))/arm/'best.pth';model,init=make_model(seed,arm);run=json.loads((ck.parent/'run.json').read_text());assert init['actual_initial_state_sha256']==run['actual_initial_state_sha256'];model.load_state_dict(torch.load(ck,map_location='cpu',weights_only=True),strict=True);model.eval();xx=transformed(x,pm,arm);key=task+'__'+arm+'__seed'+str(seed)
   with torch.no_grad():counts=model(torch.from_numpy(np.ascontiguousarray(xx)))['log_total_count'].numpy().ravel()
   out[key+'__logcounts']=counts
   wt=torch.from_numpy(np.ascontiguousarray(xx[:1]));base=wt.clone();base[:,:,:4]=.25;integral=torch.zeros_like(wt)
   for k in range(33):
    z=(base+(k/32)*(wt-base)).detach().requires_grad_(True);f=model(z)['log_total_count'].ravel();grad=torch.autograd.grad(f.sum(),z)[0];integral+=(.5 if k in [0,32] else 1)*grad.detach()
   ig=((wt-base)*integral/32).numpy();signed=ig[:,:,:4].sum(2).ravel();out[key+'__IG_signed']=signed;out[key+'__IG_importance']=abs(signed)
   with torch.no_grad():bcount=float(model(base)['log_total_count'].item())
   d=float(counts[0]-bcount);err=float(ig.sum()-d);checks.append({'task':task,'arm':arm,'seed':seed,'IG_steps':32,'f_WT_minus_f_base':d,'sum_IG':float(ig.sum()),'absolute_error':abs(err),'relative_error':abs(err)/max(abs(d),1e-5)})
   sources.append({'task':task,'arm':arm,'seed':seed,'checkpoint':str(ck),'checkpoint_sha256':sha(ck),'run_sha256':sha(ck.parent/'run.json'),'initialization':init});print(json.dumps({'task':task,'arm':arm,'seed':seed,'wall_seconds':time.time()-start,'IG_relative_error':checks[-1]['relative_error']}),file=sys.stderr,flush=True);del model
buf=io.BytesIO();np.savez_compressed(buf,**out);data=buf.getvalue();assert all(np.isfinite(v).all() for v in out.values())
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.time()-start,'count_forward_queries':151*24,'IG_gradient_evaluations':33*24,'baseline_forward_queries':24,'cpu_threads':4,'device':'CPU','strict_load':True,'models':sources,'IG_completeness':checks,'input_sha256':sha(R/'INPUTS_MATCHED.npz'),'protocol_sha256':sha(R/'PROTOCOL_FR3.md'),'script_sha256':sha(R/'remote_fr3.py'),'matched_common_sha256':sha(M/'matched_common.py'),'predictions_sha256':hashlib.sha256(data).hexdigest(),'WT_ALT_same_condition_Phi':True,'permutation_definition':'permutations(task|FR3_JPX_WT,0,1), repeated151; all seeds same permutation','IG_scalar':'native log_total_count; sequence path .25 to WT; Phi kept current-condition fixed; abs(sum nucleotideIG)'}
(R/'PREDICTIONS_MATCHED.npz').write_bytes(data);(R/'RESULT_MATCHED.json').write_text(json.dumps(result,indent=2));print(json.dumps({'result':result,'arrays_base64':base64.b64encode(data).decode()}))
