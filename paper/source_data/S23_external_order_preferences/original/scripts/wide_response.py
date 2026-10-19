import os,sys,time,json,itertools,hashlib,argparse
from pathlib import Path
import numpy as np
import h5py,torch
MATCHED=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched');ROOT=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/02_external_preferences')
sys.path.insert(0,str(MATCHED));from matched_common import make_model,row_key,transformed
ARMS=['full','sequence','frozen_flat','frozen_perm']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,j):p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(j,indent=2),encoding='utf-8')

def selected(task):
    out=ROOT/'wide_inputs'/f'{task}.npz';out.parent.mkdir(exist_ok=True)
    if out.exists():
        a=np.load(out);return a['x'],a['pm'],a['sequence'],a['chromosome'],a['source_row']
    xx=[];cc=[];rr=[];files=[]
    for chrom in [8,15]:
        p=Path('/data/CARBON_PROJECT/数据/h5data')/task/'_5ch'/f'{task}_chr{chrom}_5ch.h5'
        with h5py.File(p) as f:x=f['X'][:].astype(np.float32)
        assert x.shape[1:]==(101,5);xx.append(x);cc.extend([chrom]*len(x));rr.extend(range(len(x)));files.append(dict(path=str(p),bytes=p.stat().st_size,n_windows=len(x),sha256=sha(p)))
    x=np.concatenate(xx);c=np.array(cc);r=np.array(rr);good=np.isfinite(x).all((1,2))&np.isin(x[:,:,:4],[0,1]).all((1,2))&(x[:,:,:4].sum(2)==1).all(1)
    ids=sorted(np.flatnonzero(good),key=lambda i:hashlib.sha256(f'external_preferences_20261004|{task}|{int(c[i])}|{int(r[i])}'.encode()).digest())[:1024]
    ids=np.array(ids);assert len(ids)>=16,'Fewer than16 qualified native contexts'
    chosen=x[ids];cs=c[ids];rs=r[ids]
    pm=np.array([np.random.default_rng(int.from_bytes(row_key(task,int(ci),int(ri))[:8],'big')).permutation(101) for ci,ri in zip(cs,rs)],dtype=np.uint8)
    seq=np.asarray([''.join('ACGT'[j] for j in ss) for ss in chosen[:,:,:4].argmax(2)])
    np.savez_compressed(out,x=chosen,pm=pm,sequence=seq,chromosome=cs,source_row=rs)
    save(out.with_suffix('.json'),dict(task=task,input_files=files,n_total=len(x),n_eligible=int(good.sum()),n_selected=len(ids),selection='external_preferences_20261004 hash; min1024/all qualified; fixed16 controlled contexts',cache_sha256=sha(out)))
    return chosen,pm,seq,cs,rs

def run(spec,arm,device):
    task=spec['task'];out=ROOT/'wide_outputs'/task/arm;out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists():print('SKIP',task,arm,flush=True);return
    start=time.time();x,pm,seq,chrom,row=selected(task);x=x[:16];pm=pm[:16]
    cp=spec['sequence_checkpoint_path' if arm=='sequence' else 'full_checkpoint_path'];expected=spec['sequence_checkpoint_sha256' if arm=='sequence' else 'full_checkpoint_sha256'];assert sha(cp)==expected
    model,init=make_model(11,'sequence_only' if arm=='sequence' else 'aligned');model.load_state_dict(torch.load(cp,map_location='cpu',weights_only=True),strict=True);model=model.to(device).eval()
    mapping={'full':'aligned','sequence':'sequence_only','frozen_flat':'mean','frozen_perm':'permuted'};xx=transformed(x,pm,mapping[arm]);words=[''.join(v) for v in itertools.product('ACGT',repeat=5)];onehot=np.eye(4,dtype=np.float32)[np.array([['ACGT'.index(b) for b in k] for k in words])]
    pred=[]
    with torch.no_grad():
        wt=model(torch.from_numpy(xx).to(device))['log_total_count'].ravel().cpu().numpy()
        for context in xx:
            z=np.repeat(context[None],1024,0);z[:,48:53,:4]=onehot;yy=[]
            for lo in range(0,1024,256):yy.append(model(torch.from_numpy(z[lo:lo+256]).to(device))['log_total_count'].ravel().cpu().numpy())
            pred.append(np.concatenate(yy))
    pred=np.stack(pred);np.savez_compressed(out/'responses.npz',kmers=np.array(words),delta_log_count=pred-wt[:,None],wt_log_count=wt,sequence=seq[:16],native_conservation=x[:,:,4],conditional_conservation=xx[:,:,4] if arm!='sequence' else np.empty((0,101)),chromosome=chrom[:16],source_row=row[:16])
    result=dict(task=task,protein=spec['protein'],arm=arm,n_contexts=16,n_kmers=1024,device=device,threads=torch.get_num_threads(),elapsed_seconds=time.time()-start,checkpoint_sha256=expected,checkpoint_path=cp,model_source_sha256=sha(MATCHED/'source/MARCP/models/marcp_model.py'),input_cache_sha256=sha(ROOT/'wide_inputs'/f'{task}.npz'),scalar='log_total_count',strict_load=True,trained_matching_not_verified=True,initialization_not_used=init)
    save(out/'summary.json',result);print('DONE',task,arm,result['elapsed_seconds'],flush=True)
    del model

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--task');ap.add_argument('--arm');ap.add_argument('--device',default='cpu');args=ap.parse_args();torch.set_num_threads(8)
    specs=json.loads((ROOT/'WIDE_MANIFEST.json').read_text())
    for spec in specs:
        if args.task and spec['task']!=args.task:continue
        for arm in ([args.arm] if args.arm else ARMS):
            try:run(spec,arm,args.device)
            except Exception as exc:
                save(ROOT/'wide_errors'/f'{spec["task"]}_{arm}.json',dict(task=spec['task'],arm=arm,error=repr(exc)));print('ERROR',spec['task'],arm,repr(exc),flush=True)
