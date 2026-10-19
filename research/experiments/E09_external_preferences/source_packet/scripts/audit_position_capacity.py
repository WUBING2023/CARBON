"""Independent CPU structural, fresh-state, logit-consumption and run audit."""
import os
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['PYTHONDONTWRITEBYTECODE']='1'
from pathlib import Path
import sys,json
import numpy as np,torch
M=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched');P=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_position_capacity');R=P.parent.parent/'02_external_preferences'
sys.path.insert(0,str(M));sys.path.insert(0,str(P))
from matched_common import make_model,state_hash,sha,CombinedProfileLoss
from MARCP.models.marcp_model import MARCP
from control_model import NoGRUPositionCapacity
torch.set_num_threads(2);records=[]
for task in ['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN']:
    native=np.load(M/'test_data'/task/'observed.npz');idx=np.flatnonzero(native['y'].sum(1)>0)[:4];xx=torch.from_numpy(native['x'][idx].astype(np.float32));yy=torch.from_numpy(native['y'][idx].astype(np.float32))
    for seed in [11,23,47]:
        full,init=make_model(seed,'aligned');model=NoGRUPositionCapacity();model.load_state_dict({k:full.state_dict()[k].clone() if k in full.state_dict() else v.clone() for k,v in model.state_dict().items()},strict=True)
        d=P/'outputs'/task/f'seed{seed}';run=json.loads((d/'run.json').read_text());complete=json.loads((d/'complete.json').read_text());fresh=state_hash(model.state_dict());assert fresh==run['actual_initial_state_sha256']
        old=MARCP(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True);old.load_state_dict({k:full.state_dict()[k].clone() for k in old.state_dict()},strict=True);old.eval();model.eval()
        assert sum(p.numel() for p in model.parameters())==100871 and sum(p.numel() for p in full.parameters())==100826
        assert not any(isinstance(m,(torch.nn.GRU,torch.nn.LSTM)) for m in model.modules());assert all(m.kernel_size==(1,) for m in model.mixer.position_mlp.modules() if isinstance(m,torch.nn.Conv1d))
        assert np.allclose(model.position_coordinates[0,0].numpy(),(np.arange(101)-50)/50) and np.allclose(model.position_coordinates[0,1].numpy(),((np.arange(101)-50)/50)**2)
        with torch.no_grad():a=model(xx);b=old(xx)
        assert set(a)=={'logits_profile','log_total_count'};errors={k:float((a[k]-b[k]).abs().max()) for k in a};assert max(errors.values())==0
        bias=torch.linspace(-1,1,101)
        with torch.no_grad():model.position_profile_bias.copy_(bias);modified=model(xx)
        delta=modified['logits_profile']-a['logits_profile'];biaserr=float((delta-bias[None]).abs().max());assert biaserr<3e-7
        criterion=CombinedProfileLoss(alpha=1.42,beta=.54,gamma=1.80,loss_scale=1.0);model.zero_grad();now=model(xx);loss=criterion(now['logits_profile'],now['log_total_count'],yy)['loss'];loss.backward();grad=float(model.position_profile_bias.grad.abs().sum());assert grad>0
        probchange=float((torch.softmax(modified['logits_profile'],1)-torch.softmax(a['logits_profile'],1)).abs().max());assert probchange>0
        model.load_state_dict(torch.load(d/'best.pth',map_location='cpu',weights_only=True),strict=True);assert sha(d/'best.pth')==complete['checkpoint_sha256'];trainedbias=float(model.position_profile_bias.detach().abs().sum());assert trainedbias>0
        logs=[json.loads(line) for line in (d/'epochs.jsonl').read_text().splitlines()];flogs=[json.loads(line) for line in (M/'outputs'/task/f'seed{seed}'/'aligned'/'epochs.jsonl').read_text().splitlines()];shared=min(len(logs),len(flogs));assert all(logs[i]['train_order_sha256']==flogs[i]['train_order_sha256'] for i in range(shared));assert complete['best_epoch']==max(logs,key=lambda v:v['validation_pcc'])['epoch'] and len(logs)<=100
        records.append(dict(task=task,seed=seed,parameters=100871,full_parameters=100826,shared_tensor_init_exact=True,fresh_init_hash_matches=True,initial_oldNoGRU_output_errors=errors,has_no_RNN=True,MLP_pointwise=True,explicit_t_t2_coordinates=True,additive_bias_logit_max_error=biaserr,bias_effect_consumed_in_softmax=probchange,bias_nonzero_training_loss_gradient_L1=grad,trained_bias_L1=trainedbias,all_real_common_epoch_orders_match=True,n_common_epoch_orders=shared,epochs=len(logs),best_epoch=complete['best_epoch'],checkpoint_hash_verified=True,validation_checkpoint_selection_verified=True))
        del full,model,old
trainer=(P/'train_nogru.py').read_text();inference=(P/'remote_infer.py').read_text()
assert "criterion(pred['logits_profile'],pred['log_total_count'],yy)" in trainer
assert "torch.softmax(out['logits_profile'],dim=1)" in trainer
assert "torch.softmax(pred['logits_profile'],dim=1)" in inference
j=dict(device='CPU2',runs=records,loss_consumes_updated_logit=True,validation_consumes_updated_logit=True,inference_consumes_updated_logit=True,parent_forward_has_no_cached_probability=True,source_sha256=sha(__file__),control_source_sha256=sha(P/'control_model.py'),trainer_sha256=sha(P/'train_nogru.py'),inference_sha256=sha(P/'remote_infer.py'),synthetic_bias_is_numerical_QA_not_biological_evidence=True)
(R/'CROSS_AUDIT_01_POSITION_CAPACITY_MODEL.json').write_text(json.dumps(j,indent=2));print(json.dumps(j,indent=2))
