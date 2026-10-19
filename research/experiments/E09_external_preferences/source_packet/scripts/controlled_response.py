"""All five-mer substitutions in 16 fixed native contexts; no external labels read."""
import sys, os, json, itertools, time, argparse
from pathlib import Path
import numpy as np
import torch
sys.path.insert(0,str(Path(__file__).resolve().parent))
from extract_rules import MATCHED,ROOT,TASKS,ARMS,SEEDS,select,sha,save_json
sys.path.insert(0,str(MATCHED))
from matched_common import make_model,transformed

def run(task,seed,arm):
    out=ROOT/'controlled_outputs'/task/f'seed{seed}'/arm;out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists():print('SKIP',task,seed,arm,flush=True);return
    start=time.time();x,pm,seq,meta=select(task);x=x[:16];pm=pm[:16];seq=seq[:16]
    kmers=[''.join(k) for k in itertools.product('ACGT',repeat=5)]
    nk=np.asarray([['ACGT'.index(c) for c in k] for k in kmers]);onehot=np.eye(4,dtype=np.float32)[nk]
    model,init=make_model(seed,arm);ckpt=MATCHED/'outputs'/task/f'seed{seed}'/arm/'best.pth';model.load_state_dict(torch.load(ckpt,map_location='cpu',weights_only=True),strict=True);model=model.to('cuda').eval()
    xx=transformed(x,pm,arm)
    wt=[];pred=[]
    with torch.no_grad():
        wt=model(torch.from_numpy(xx).to('cuda'))['log_total_count'].ravel().cpu().numpy()
        for row in xx:
            z=np.repeat(row[None],1024,axis=0);z[:,48:53,:4]=onehot
            yy=[]
            for lo in range(0,1024,256):yy.append(model(torch.from_numpy(z[lo:lo+256]).to('cuda'))['log_total_count'].ravel().cpu().numpy())
            pred.append(np.concatenate(yy))
    predicted=np.stack(pred);delta=predicted-wt[:,None]
    np.savez_compressed(out/'responses.npz',kmers=np.array(kmers),delta_log_count=delta,predicted_log_count=predicted,wt_log_count=wt,sequence=seq,native_conservation=x[:,:,4],conditional_conservation=xx[:,:,4] if arm!='sequence_only' else np.empty((0,101)),chromosome=np.array(meta['chromosome'][:16]),source_row=np.array(meta['source_row'][:16]))
    meta.update(seed=seed,arm=arm,task=task,backgrounds=16,kmers=1024,replaced_zero_based_positions=[48,52],initialization=init,checkpoint_sha256=sha(ckpt),script_sha256=sha(__file__),elapsed_seconds=time.time()-start)
    save_json(out/'summary.json',meta);print('DONE',task,seed,arm,round(time.time()-start,3),flush=True)
    del model;torch.cuda.empty_cache()

if __name__=='__main__':
    torch.set_num_threads(8)
    for task in TASKS:
        for seed in SEEDS:
            for arm in ARMS:run(task,seed,arm)
