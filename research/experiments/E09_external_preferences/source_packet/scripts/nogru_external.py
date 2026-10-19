"""Fixed PTBP1 controlled-word diagnostic; CPU8, no external labels read."""
import os
os.environ['CUDA_VISIBLE_DEVICES']=''
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import sys,json,itertools,time
from pathlib import Path
import numpy as np,torch
sys.path.insert(0,str(Path(__file__).resolve().parent))
sys.path.insert(0,str(Path('/data/CARBON_PROJECT/supervised_expansion_20261004/02_external_preferences/scripts')))
from extract_rules import MATCHED,ROOT,select,sha,save_json
sys.path.insert(0,str(MATCHED/'source'))
from MARCP.models.marcp_model import MARCP
WEIGHTS=ROOT.parent/'01_context'/'nogru_matched'/'outputs'
torch.set_num_threads(8);torch.set_num_interop_threads(1)
for task in ['PTBP1_HepG2_ENCSR384KAN','PTBP1_K562_ENCSR981WKN']:
    x,pm,seq,meta=select(task);x=x[:16];seq=seq[:16]
    kmers=[''.join(k) for k in itertools.product('ACGT',repeat=5)]
    onehot=np.eye(4,dtype=np.float32)[np.array([['ACGT'.index(c) for c in k] for k in kmers])]
    for seed in [11,23,47]:
        start=time.time();out=ROOT/'nogru_external_outputs'/task/f'seed{seed}';out.mkdir(parents=True,exist_ok=True)
        if (out/'summary.json').exists():print('SKIP',task,seed,flush=True);continue
        ck=WEIGHTS/task/f'seed{seed}'/'best.pth';complete=json.loads((ck.parent/'complete.json').read_text());assert complete['checkpoint_sha256']==sha(ck)
        model=MARCP(in_channels=5,dropOut=.2,control=False,ablation_no_gru=True);model.load_state_dict(torch.load(ck,map_location='cpu',weights_only=True),strict=True);model.eval()
        with torch.no_grad():
            wt=model(torch.from_numpy(x))['log_total_count'].ravel().numpy();pred=[]
            for row in x:
                z=np.repeat(row[None],1024,axis=0);z[:,48:53,:4]=onehot
                pred.append(np.concatenate([model(torch.from_numpy(z[i:i+256]))['log_total_count'].ravel().numpy() for i in range(0,1024,256)]))
        pred=np.stack(pred);np.savez_compressed(out/'responses.npz',kmers=np.array(kmers),delta_log_count=pred-wt[:,None],predicted_log_count=pred,wt_log_count=wt,sequence=seq,native_conservation=x[:,:,4],conditional_conservation=x[:,:,4],chromosome=np.array(meta['chromosome'][:16]),source_row=np.array(meta['source_row'][:16]))
        info=dict(task=task,seed=seed,condition='nogru',device='CPU8',parameters=sum(p.numel() for p in model.parameters()),checkpoint_sha256=sha(ck),input_sha256=meta['input_sha256'],script_sha256=sha(__file__),n_backgrounds=16,n_words=1024,elapsed_seconds=time.time()-start)
        save_json(out/'summary.json',info);print('DONE',json.dumps(info),flush=True);del model
