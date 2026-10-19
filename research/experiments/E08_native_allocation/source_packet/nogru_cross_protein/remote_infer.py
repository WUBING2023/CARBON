import os
os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
import sys,json,time,hashlib,resource
import numpy as np,torch
R=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched');P=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context');O=P/'nogru_cross_protein'
sys.path.insert(0,str(R));from matched_common import make_model,sha,state_hash
from MARCP.models.marcp_model import MARCP
torch.set_num_threads(4);torch.set_num_interop_threads(1)
T=time.time();out={'tasks':[],'device':'CPU4','started_unix':T,'source_sha256':sha(O/'remote_infer.py'),'protocol_sha256':sha(O/'PROTOCOL.md')}
for task in ['HNRNPK_K562_ENCSR268ETU','HNRNPK_HepG2_ENCSR828ZID']:
 for chrom in ['1','X','Y']:
  src=P/'data'/task if chrom=='1' else P/'phase2_ambiguity/data'/task/chrom
  zfile=src/'data.npz';metafile=src/'metadata.json'
  with np.load(zfile) as f:x=f['x'];index=f['index'];centers=f['centers']
  z={'index':index,'centers':centers};checkpoints=[]
  for seed in [11,23,47]:
   ck=O/'outputs'/task/('seed'+str(seed))/'best.pth';run=json.loads((ck.parent/'run.json').read_text());complete=json.loads((ck.parent/'complete.json').read_text());assert complete['checkpoint_sha256']==sha(ck)
   full,init=make_model(seed,'aligned');model=MARCP(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True);model.load_state_dict({k:full.state_dict()[k].clone() for k in model.state_dict()},strict=True);del full
   assert state_hash(model.state_dict())==run['actual_initial_state_sha256'];model.load_state_dict(torch.load(ck,map_location='cpu',weights_only=True),strict=True);model.eval();p=[];cnt=[];t=time.time()
   with torch.no_grad():
    for st in range(0,len(x),128):
     pred=model(torch.from_numpy(np.ascontiguousarray(x[st:st+128])));p.append(torch.softmax(pred['logits_profile'],dim=1).numpy());cnt.append(pred['log_total_count'].numpy().ravel())
   prob=np.concatenate(p) if p else np.empty((0,101),np.float32);count=np.concatenate(cnt) if cnt else np.empty(0,np.float32)
   assert np.isfinite(prob).all() and (not len(prob) or np.allclose(prob.sum(1),1,atol=1e-6))
   z['nogru_seed'+str(seed)]=prob;z['nogru_seed'+str(seed)+'_logcount']=count
   checkpoints.append({'seed':seed,'checkpoint_path':str(ck),'checkpoint_sha256':sha(ck),'run_sha256':sha(ck.parent/'run.json'),'complete_sha256':sha(ck.parent/'complete.json'),'seconds':time.time()-t,'best_epoch':complete['best_epoch'],'parameters':run['parameters'],'initial_state_sha256':run['actual_initial_state_sha256']})
   del model
  od=O/'data'/task/chrom;od.mkdir(parents=True,exist_ok=True);np.savez_compressed(od/'data.npz',**z)
  rec={'task':task,'chromosome':chrom,'windows':len(x),'native_data_path':str(zfile),'native_data_sha256':sha(zfile),'native_metadata_path':str(metafile),'native_metadata_sha256':sha(metafile),'native_input_sha256':hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest(),'checkpoints':checkpoints,'output_data_sha256':sha(od/'data.npz')};out['tasks'].append(rec);(O/'INFERENCE_MANIFEST.json').write_text(json.dumps(out,indent=2));print(json.dumps(rec),flush=True)
out.update(seconds=time.time()-T,cpu_user_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime,cpu_system_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_stime);(O/'INFERENCE_MANIFEST.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
