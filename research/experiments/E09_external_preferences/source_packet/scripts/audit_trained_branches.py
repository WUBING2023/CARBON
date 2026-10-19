"""Independent CPU readback of all six fitted position/capacity branches."""
from pathlib import Path
import hashlib,json
import torch
torch.set_num_threads(2)
R=Path(__file__).resolve().parents[1];A=R.parent/'01_context'/'nogru_position_capacity'
saved=json.loads((A/'CONTROL_LEARNING_QA.json').read_text());checks=[]
for row in saved['checkpoints']:
    p=A/'remote'/'outputs'/row['task']/row['seed']/'best.pth'
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row['checkpoint_sha256']
    state=torch.load(p,map_location='cpu',weights_only=True)
    for name,reported in row['learned_tensors'].items():
        tensor=state[name].float()
        observed={'l1':float(tensor.abs().sum()),'l2':float(tensor.norm()),'max_abs':float(tensor.abs().max())}
        for metric,value in observed.items():
            assert abs(value-reported[metric])<1e-6
            checks.append(dict(task=row['task'],seed=row['seed'],tensor=name,metric=metric,independent=value,reported=reported[metric],abs_difference=abs(value-reported[metric])))
        if name in ['position_profile_bias','mixer.position_mlp.2.weight','mixer.position_mlp.2.bias']:
            assert observed['l1']>0
out=dict(device='local CPU2',no_forward_no_training=True,six_checkpoint_hashes_verified=True,
         all_zero_initialized_private_tensors_trained_nonzero=True,n_independent_norm_checks=len(checks),
         max_abs_error=max(v['abs_difference'] for v in checks),checks=checks,
         control_source_sha256=hashlib.sha256((A/'control_model.py').read_bytes()).hexdigest(),
         audit_is_implementation_validity_not_biological_mechanism=True)
(R/'results'/'CROSS_AUDIT_01_PC_TRAINED_BRANCHES.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='checks'},indent=2))
