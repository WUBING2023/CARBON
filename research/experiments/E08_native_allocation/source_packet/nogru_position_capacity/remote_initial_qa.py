import os
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['PYTHONDONTWRITEBYTECODE']='1'
from pathlib import Path
import sys,json,time
import torch
R=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched');O=Path(__file__).resolve().parent
sys.path.insert(0,str(R));from matched_common import make_model,state_hash,sha
from MARCP.models.marcp_model import MARCP
from control_model import NoGRUPositionCapacity
torch.set_num_threads(4);out={'device':'CPU4','protocol_sha256':sha(O/'PROTOCOL.md'),'control_source_sha256':sha(O/'control_model.py'),'seeds':[]}
for seed in [11,23,47]:
 full,ini=make_model(seed,'aligned');model=NoGRUPositionCapacity();model.load_state_dict({k:(full.state_dict()[k].clone() if k in full.state_dict() else v.clone()) for k,v in model.state_dict().items()},strict=True);del full
 base=MARCP(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True);base.load_state_dict({k:model.state_dict()[k].clone() for k in base.state_dict()},strict=True);base.eval();model.eval()
 x=torch.linspace(-1,1,16*101*5).reshape(16,101,5)
 with torch.no_grad():aa=base(x);bb=model(x)
 error={k:float((aa[k]-bb[k]).abs().max()) for k in aa};assert all(v==0 for v in error.values())
 assert sum(p.numel() for p in model.parameters())==100871
 assert model.mixer.position_mlp[0].weight.abs().sum()>0 and model.mixer.position_mlp[2].weight.abs().sum()==0 and model.position_profile_bias.abs().sum()==0
 out['seeds'].append({'seed':seed,'parameters':100871,'initial_state_sha256':state_hash(model.state_dict()),'common_full_initial_state_sha256':ini['actual_initial_state_sha256'],'exact_output_max_error':error,'conv_kernels':{n:list(m.kernel_size) for n,m in model.named_modules() if isinstance(m,torch.nn.Conv1d)}})
(O/'INITIAL_QA.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
