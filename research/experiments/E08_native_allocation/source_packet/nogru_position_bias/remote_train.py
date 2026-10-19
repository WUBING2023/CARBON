import os
os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['CUDA_VISIBLE_DEVICES']='2';os.environ['CUBLAS_WORKSPACE_CONFIG']=':4096:8'
from pathlib import Path
import subprocess,json,time,sys,hashlib
R=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched');O=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_position_bias');O.mkdir(parents=True,exist_ok=True)
source=(R/'matched_train.py').read_text()
override='''
def make_model(seed, arm):
    assert arm=='aligned'
    from MARCP.models.marcp_model import MARCP
    class NoGRUPosition(MARCP):
        def __init__(self):
            super().__init__(in_channels=5,dropOut=0.2,control=False,ablation_no_gru=True)
            self.position_profile_bias=torch.nn.Parameter(torch.zeros(101))
        def forward(self,x):
            result=super().forward(x)
            result['logits_profile']=result['logits_profile']+self.position_profile_bias.unsqueeze(0)
            return result
    full,original=_matched_make_model(seed,'aligned')
    model=NoGRUPosition()
    shared={k:(full.state_dict()[k].clone() if k in full.state_dict() else torch.zeros_like(v)) for k,v in model.state_dict().items()}
    model.load_state_dict(shared,strict=True)
    return model, {'common_five_channel_state_sha256':original['common_five_channel_state_sha256'],'actual_initial_state_sha256':state_hash(model.state_dict()),'parameters':sum(p.numel() for p in model.parameters()),'full_parameters':original['parameters'],'architecture':'native noGRU+101 zero-init learned absolute-position profile logits bias; count/SE/heads retained'}
'''
source=source.replace('    make_model, set_seed, state_hash,','    make_model as _matched_make_model, set_seed, state_hash,')
source=source.replace('\ndef save_torch(',override+'\ndef save_torch(')
source=source.replace("    evaluation(model,manifest,args.task,args.arm,outdir,device)","    # Heldout evaluation deferred; no test labels used to select the checkpoint.")
source="import sys\nsys.path.insert(0,"+repr(str(R))+')\n'+source
(O/'train_nogru.py').write_text(source)
manifest={'original_matched_train_sha256':hashlib.sha256((R/'matched_train.py').read_bytes()).hexdigest(),'new_train_source_sha256':hashlib.sha256(source.encode()).hexdigest(),'protocol_sha256':hashlib.sha256((O/'PROTOCOL.md').read_bytes()).hexdigest(),'started_unix':time.time(),'runs':[],'gpu':2,'max_batch_seconds':2100}
(O/'TRAINING_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
for task in ['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']:
 for seed in [11,23,47]:
  if time.time()-manifest['started_unix']>2100:raise RuntimeError('Registered35min cap exhausted; no newrun')
  uuid=subprocess.check_output(['nvidia-smi','-i','2','--query-gpu=uuid','--format=csv,noheader'],text=True).strip()
  apps=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid','--format=csv,noheader'],text=True).strip()
  guard=[v for v in apps.splitlines() if v.split(',')[0].strip()==uuid]
  mem=int(subprocess.check_output(['nvidia-smi','-i','2','--query-gpu=memory.used','--format=csv,noheader,nounits'],text=True).strip())
  if guard or mem>10:raise RuntimeError('GPU2 guard rejected newrun; no waiting or contention: '+str(guard)+' memory='+str(mem))
  out=O/'outputs'/task/('seed'+str(seed));out.mkdir(parents=True,exist_ok=True);t=time.time()
  cmd=[sys.executable,str(O/'train_nogru.py'),'--task',task,'--seed',str(seed),'--arm','aligned','--output',str(out)]
  with open(out/'console.log','w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
  rec={'task':task,'seed':seed,'exit_code':p.returncode,'seconds':time.time()-t,'guard':{'uuid':uuid,'apps':guard,'memory_MiB':mem}}
  if (out/'complete.json').exists():rec['complete']=json.loads((out/'complete.json').read_text())
  manifest['runs'].append(rec);(O/'TRAINING_MANIFEST.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(rec),flush=True)
  if p.returncode:raise RuntimeError('NoGRU run failed; see new console.log; no effect-dependent continuation')
manifest['elapsed_seconds']=time.time()-manifest['started_unix'];(O/'TRAINING_MANIFEST.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2),flush=True)
