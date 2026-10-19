import os
os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['CUDA_VISIBLE_DEVICES']='2';os.environ['CUBLAS_WORKSPACE_CONFIG']=':4096:8'
from pathlib import Path
import json,time,subprocess,sys,hashlib
import numpy as np,torch
R=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched');O=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_pilot');O.mkdir(parents=True,exist_ok=True)
uuid=subprocess.check_output(['nvidia-smi','-i','2','--query-gpu=uuid','--format=csv,noheader'],text=True).strip()
apps=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid','--format=csv,noheader'],text=True).strip()
guard='\n'.join(v for v in apps.splitlines() if v.split(',')[0].strip()==uuid)
mem=subprocess.check_output(['nvidia-smi','-i','2','--query-gpu=memory.used','--format=csv,noheader,nounits'],text=True).strip()
if guard or int(mem)>10:raise RuntimeError('GPU2 resource guard rejected: live apps or nontrivial memory; no waiting/contending')
sys.path.insert(0,str(R));from matched_common import make_model,set_seed,state_hash,load_split,CombinedProfileLoss,sha
from matched_train import validate
from MARCP.models.marcp_model import MARCP
torch.set_num_threads(4);torch.set_num_interop_threads(1)
m=json.loads((R/'manifest.json').read_text());result={'resource_guard':{'memory_MiB':int(mem),'apps':guard},'tasks':[],'device':'GPU2','protocol_sha256':sha(O/'PROTOCOL.md'),'source_sha256':sha(O/'remote_pilot.py')}
for task in ['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']:
 seed=11;t=time.time();full,initial=make_model(seed,'aligned');model=MARCP(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True)
 shared={k:full.state_dict()[k].clone() for k in model.state_dict()};model.load_state_dict(shared,strict=True);del full
 params=sum(p.numel() for p in model.parameters());h=state_hash(model.state_dict());model.to('cuda:0')
 tr=load_split(m,task,'train','aligned',verify=True);va=load_split(m,task,'validation','aligned',verify=True)
 tx=torch.as_tensor(tr['x'],device='cuda:0');ty=torch.as_tensor(tr['y'],device='cuda:0');vx=torch.as_tensor(va['x'],device='cuda:0');vy=torch.as_tensor(va['y'],device='cuda:0');del tr,va
 torch.cuda.synchronize();load=time.time()-t
 opt=torch.optim.AdamW(model.parameters(),lr=.000785,weight_decay=.000114);criterion=CombinedProfileLoss(alpha=1.42,beta=.54,gamma=1.80,loss_scale=1)
 gen=torch.Generator(device='cpu').manual_seed(seed+10000);set_seed(seed+20000);bsz=256;order=torch.randperm(len(tx),generator=gen);order_hash=hashlib.sha256(order.numpy().tobytes()).hexdigest();n=(len(tx)//bsz)*bsz;order=order[:n].to('cuda:0')
 model.train();t=time.time()
 for start in range(0,n,bsz):
  ix=order[start:start+bsz];pred=model(tx[ix]);loss=criterion(pred['logits_profile'],pred['log_total_count'],ty[ix])['loss'];assert torch.isfinite(loss)
  opt.zero_grad(set_to_none=True);loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),2,error_if_nonfinite=True);opt.step()
 v=validate(model,vx,vy,criterion,bsz);torch.cuda.synchronize();secs=time.time()-t
 t=time.time();torch.save({k:v.detach().cpu() for k,v in model.state_dict().items()},O/(task+'_one_epoch_COST_ONLY.pth'));io=time.time()-t
 previous=[]
 for s in [11,23,47]:
  p=R/'outputs'/task/('seed'+str(s))/'aligned';ck=json.loads((p/'complete.json').read_text());epochs=[json.loads(l) for l in (p/'epochs.jsonl').read_text().splitlines()]
  previous.append({'seed':s,'epochs':len(epochs),'best_epoch':ck['best_epoch'],'training_seconds':ck['training_seconds'],'median_epoch_seconds':float(np.median([x['seconds'] for x in epochs]))})
 rec={'task':task,'parameters':params,'full_parameters':initial['parameters'],'common_full_initial_state_sha256':initial['actual_initial_state_sha256'],'actual_nogru_initial_state_sha256':h,'load_seconds':load,'epoch_seconds':secs,'checkpoint_io_seconds':io,'training_examples':len(tx),'epoch_training_examples':n,'validation_examples':len(vx),'validation':v,'sampler_order_sha256':order_hash,'old_full_runs':previous}
 result['tasks'].append(rec);(O/'PILOT.json').write_text(json.dumps(result,indent=2));print(json.dumps(rec),flush=True)
 del tx,ty,vx,vy,model,opt;torch.cuda.empty_cache()
result['upper6runs100epochs_seconds']=3*100*sum(t['epoch_seconds'] for t in result['tasks'])+3*sum(t['load_seconds']+100*t['checkpoint_io_seconds'] for t in result['tasks']);result['feasible_under90min']=result['upper6runs100epochs_seconds']<=5400
(O/'PILOT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
