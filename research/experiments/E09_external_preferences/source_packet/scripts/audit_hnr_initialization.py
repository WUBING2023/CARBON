"""Independent CPU audit; all existing runs/weights/source are read-only."""
import os
os.environ['CUDA_VISIBLE_DEVICES']='';os.environ['PYTHONDONTWRITEBYTECODE']='1'
from pathlib import Path
import sys,json,hashlib
import torch
M=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched');A=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_cross_protein');O=A.parent.parent/'02_external_preferences'
sys.path.insert(0,str(M));from matched_common import make_model,state_hash,sha
from MARCP.models.marcp_model import MARCP
torch.set_num_threads(2);rows=[]
for task in ['HNRNPK_K562_ENCSR268ETU','HNRNPK_HepG2_ENCSR828ZID']:
    for seed in [11,23,47]:
        f,init=make_model(seed,'aligned');fullrun=json.loads((M/'outputs'/task/f'seed{seed}'/'aligned'/'run.json').read_text());ngr=json.loads((A/'outputs'/task/f'seed{seed}'/'run.json').read_text())
        assert init['actual_initial_state_sha256']==fullrun['actual_initial_state_sha256']
        model=MARCP(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True)
        assert all(k in f.state_dict() and f.state_dict()[k].shape==v.shape for k,v in model.state_dict().items())
        model.load_state_dict({k:f.state_dict()[k].clone() for k in model.state_dict()},strict=True)
        assert state_hash(model.state_dict())==ngr['actual_initial_state_sha256'] and fullrun['manifest_sha256']==ngr['manifest_sha256']
        nglog=[json.loads(v) for v in (A/'outputs'/task/f'seed{seed}'/'epochs.jsonl').read_text().splitlines()]
        flog=[json.loads(v) for v in (M/'outputs'/task/f'seed{seed}'/'aligned'/'epochs.jsonl').read_text().splitlines()]
        shared=min(len(nglog),len(flog));orders=[nglog[i]['train_order_sha256']==flog[i]['train_order_sha256'] for i in range(shared)];assert all(orders)
        ck=A/'outputs'/task/f'seed{seed}'/'best.pth';model.load_state_dict(torch.load(ck,map_location='cpu',weights_only=True),strict=True)
        rows.append(dict(task=task,seed=seed,full_init_matches_original=True,all_nogru_shared_tensors_equal_initial_full=True,nogru_initial_hash_matches=True,same_data_manifest=True,shared_epoch_orders_equal=True,n_shared_epoch_orders=shared,nogru_epochs=len(nglog),full_epochs=len(flog),full_parameters=init['parameters'],nogru_parameters=sum(p.numel() for p in model.parameters()),nogru_final_state_hash=state_hash(model.state_dict()),nogru_checkpoint_sha256=sha(ck)))
        del f,model
j=dict(rows=rows,device='CPU2',no_existing_mutation=True,source_sha256=sha(__file__),training_source_original_sha256=sha(M/'matched_train.py'),training_source_nogru_sha256=sha(A/'train_nogru.py'))
(O/'CROSS_AUDIT_01_HNR_INITIALIZATION.json').write_text(json.dumps(j,indent=2));print(json.dumps(j,indent=2))
