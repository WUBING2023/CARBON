from pathlib import Path
import sys,json,time,io,base64,hashlib
import numpy as np,torch
R=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/03_functional_rules/fr4');N=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_matched');sys.path.insert(0,'/data/CARBON_PROJECT/conservation_upgrade_20260929_matched/source');from MARCP.models.marcp_model import MARCP
torch.set_num_threads(4);torch.set_num_interop_threads(1);sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();start=time.time();x=np.load(R/'INPUTS.npz')['X'];pred=np.load(R/'PREDICTIONS.npz');out={};sources=[]
assert not (R/'NOGRU_RESULT.json').exists()
for task in ['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']:
 ix=pred[task+'__input_indices'];out[task+'__input_indices']=ix
 for seed in [11,23,47]:
  ck=N/'outputs'/task/f'seed{seed}'/'best.pth';m=MARCP(in_channels=5,dropOut=0.2,control=False,ablation_no_gru=True);m.load_state_dict(torch.load(ck,map_location='cpu',weights_only=True),strict=True);m.eval();vals=[]
  with torch.no_grad():
   for st in range(0,len(ix),128):vals.append(m(torch.from_numpy(np.ascontiguousarray(x[ix[st:st+128]])))['log_total_count'].numpy().ravel())
  val=np.concatenate(vals);assert len(val)==len(ix) and np.isfinite(val).all();out[f'{task}__no_gru__seed{seed}']=val;sources.append({'task':task,'seed':seed,'parameters':sum(p.numel() for p in m.parameters()),'checkpoint':str(ck),'checkpoint_sha256':sha(ck),'run_json_sha256':sha(ck.parent/'run.json'),'complete_json_sha256':sha(ck.parent/'complete.json'),'run_metadata':json.loads((ck.parent/'run.json').read_text())});del m
buf=io.BytesIO();np.savez_compressed(buf,**out);data=buf.getvalue();(R/'NOGRU_PREDICTIONS.npz').write_bytes(data);res={'wall_seconds':time.time()-start,'model_queries':sum(len(pred[t+'__input_indices']) for t in ['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN'])*3,'strict_load':True,'device':'CPU','threads':4,'input_sha256':sha(R/'INPUTS.npz'),'prediction_sha256':sha(R/'NOGRU_PREDICTIONS.npz'),'models':sources,'scientific_scope':'extra registered trained-component diagnostic, same budget/split but capacity differs Full100826/NoGRU44698; not pure mechanistic attribution'};(R/'NOGRU_RESULT.json').write_text(json.dumps(res,indent=2));print(json.dumps({'files':{p:base64.b64encode((R/p).read_bytes()).decode() for p in ['NOGRU_PREDICTIONS.npz','NOGRU_RESULT.json']},'wall_seconds':res['wall_seconds'],'model_queries':res['model_queries']}))
