import os
os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
import sys,json,time
import numpy as np,torch
R=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched');O=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_matched')
sys.path.insert(0,str(R/'source'));from MARCP.models.marcp_model import MARCP
torch.set_num_threads(4);torch.set_num_interop_threads(1)
out={'definition':'Numerical QA only, artificialcopies withinnativeinputs; no newbiology/functionalmeasurement. GlobalSEsamewithinwindow.', 'tasks':[]}
for task in ['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']:
 p=O.parent/'data'/task/'data.npz'
 with np.load(p) as f:x=f['x'][:16].copy()
 x[:,67:84,:]=x[:,17:34,:];c1,c2=25,75
 assert np.array_equal(x[:,c1-8:c1+9],x[:,c2-8:c2+9])
 for seed in [11,23,47]:
  ck=O/'outputs'/task/('seed'+str(seed))/'best.pth'
  if not ck.exists():continue
  model=MARCP(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True);model.load_state_dict(torch.load(ck,map_location='cpu',weights_only=True),strict=True);model.eval()
  with torch.no_grad():logits=model(torch.from_numpy(x))['logits_profile'].numpy()
  err=float(np.max(np.abs(logits[:,c1-5:c1+6]-logits[:,c2-5:c2+6])));assert err<1e-5
  out['tasks'].append({'task':task,'seed':seed,'artificial_QA_inputs':len(x),'max_absolute_corresponding11logit_difference':err,'conv1d_kernels':{n:list(m.kernel_size) for n,m in model.named_modules() if isinstance(m,torch.nn.Conv1d)},'LayerNorm_shapes':{n:list(m.normalized_shape) for n,m in model.named_modules() if isinstance(m,torch.nn.LayerNorm)}})
(O/'RF_QA.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
