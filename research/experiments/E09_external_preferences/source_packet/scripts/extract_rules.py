"""Bounded matched-model IG extraction. External preference tables are not read here."""
import os, sys, json, hashlib, time, argparse
from pathlib import Path
import numpy as np
import torch

MATCHED=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched')
ROOT=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/02_external_preferences')
sys.path.insert(0,str(MATCHED))
from matched_common import make_model, transformed, row_key

TASKS=['HNRNPK_HepG2_ENCSR828ZID','HNRNPK_K562_ENCSR268ETU','PTBP1_HepG2_ENCSR384KAN','PTBP1_K562_ENCSR981WKN','RBFOX2_HepG2_ENCSR987FTF','RBFOX2_K562_ENCSR756CKJ']
ARMS=['aligned','sequence_only','mean','permuted']
SEEDS=[11,23,47]

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save_json(p,j):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(j,indent=2,allow_nan=False),encoding='utf-8')

def select(task):
    path=MATCHED/'test_data'/task/'observed.npz'
    d=np.load(path);x=d['x'];chrom=d['chromosome'];rows=d['source_row']
    good=np.isfinite(x).all((1,2))&np.isin(x[:,:,:4],[0,1]).all((1,2))&(x[:,:,:4].sum(2)==1).all(1)
    ids=np.flatnonzero(good)
    ids=sorted(ids,key=lambda i:hashlib.sha256(f'external_preferences_20261004|{task}|{int(chrom[i])}|{int(rows[i])}'.encode()).digest())[:1024]
    ids=np.array(ids,dtype=np.int64)
    pm=np.stack([np.random.default_rng(int.from_bytes(row_key(task,int(c),int(r))[:8],'big')).permutation(101) for c,r in zip(chrom[ids],rows[ids])]).astype(np.uint8)
    chosen=x[ids].astype(np.float32)
    seq=np.asarray([''.join('ACGT'[int(v)] for v in row) for row in chosen[:,:,:4].argmax(2)])
    return chosen,pm,seq,dict(input_sha256=sha(path),n_total=len(x),n_eligible=int(good.sum()),n_selected=len(ids),chromosome=chrom[ids].tolist(),source_row=rows[ids].tolist(),original_rows=ids.tolist(),input_shape=list(x.shape),selected_sequence_sha256=hashlib.sha256('\n'.join(seq).encode()).hexdigest())

def ig(model,x,nsteps=32,batch=16):
    baseline=x.clone();baseline[:,:,:4]=0.25
    delta=x-baseline;res=torch.zeros_like(x)
    sums=[]
    for lo in range(0,len(x),batch):
        xx=x[lo:lo+batch];base=baseline[lo:lo+batch];dd=delta[lo:lo+batch]
        grad=torch.zeros_like(xx)
        for ast in range(0,nsteps,8):
            aa=(torch.arange(ast,min(ast+8,nsteps),device=x.device,dtype=x.dtype)+.5)/nsteps
            z=(base[None]+aa[:,None,None,None]*dd[None]).reshape(-1,101,x.shape[2]).detach().requires_grad_(True)
            out=model(z)['log_total_count'].sum()
            g=torch.autograd.grad(out,z)[0].reshape(len(aa),len(xx),101,x.shape[2]).sum(0)
            grad+=g/nsteps
        rr=dd*grad;res[lo:lo+batch]=rr
        with torch.no_grad():
            truth=model(xx)['log_total_count'].ravel()-model(base)['log_total_count'].ravel()
        sums.append(torch.stack([rr.sum((1,2)),truth],1).detach().cpu().numpy())
    return res[:,:,:4].sum(2).detach().cpu().numpy(),np.concatenate(sums)

def aggregate(seq,values):
    win=np.lib.stride_tricks.sliding_window_view(values,5,axis=1).sum(2)
    kmers=np.asarray([[s[j:j+5] for j in range(97)] for s in seq])
    keys=sorted(set(kmers.ravel()));loc={s:i for i,s in enumerate(keys)}
    total=np.zeros(len(keys));cnt=np.zeros(len(keys));winner_sum=np.zeros(len(keys));winner_n=np.zeros(len(keys))
    for i in range(len(seq)):
        w=int(np.argmax(win[i]));kw=loc[kmers[i,w]];winner_sum[kw]+=float(win[i,w]);winner_n[kw]+=1
        for k,v in zip(kmers[i],win[i]):
            ii=loc[k];total[ii]+=float(v);cnt[ii]+=1
    return dict(kmers=keys,mean_occurrence=(total/cnt).tolist(),occurrence_n=cnt.tolist(),winner_sum=winner_sum.tolist(),winner_n=winner_n.tolist())

def run(task,seed,arm,args):
    started=time.time();x,pm,seq,meta=select(task)
    out=ROOT/'outputs'/task/f'seed{seed}'/arm;out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists() and not args.overwrite and not args.qa:
        print('SKIP',task,seed,arm,flush=True);return
    model,initial=make_model(seed,arm);checkpoint=MATCHED/'outputs'/task/f'seed{seed}'/arm/'best.pth'
    model.load_state_dict(torch.load(checkpoint,map_location='cpu',weights_only=True),strict=True)
    model=model.to(args.device).eval()
    # cuDNN RNN needs training-mode backward; this single-layer GRU has dropout=0.
    # Every BatchNorm/Dropout outside GRU stays in evaluation mode.
    model.mixer.gru.train(True)
    xx=torch.from_numpy(transformed(x,pm,arm)).to(args.device)
    if args.qa:
        qa={}
        for st in [16,32,64]:
            vals,cpl=ig(model,xx[:32],st,args.batch)
            np.savez_compressed(out/f'qa_steps{st}.npz',position_ig=vals,completeness=cpl)
            qa[str(st)]=dict(mean_abs_error=float(np.abs(cpl[:,0]-cpl[:,1]).mean()),max_abs_error=float(np.abs(cpl[:,0]-cpl[:,1]).max()),mean_abs_output_difference=float(np.abs(cpl[:,1]).mean()))
        q32=np.load(out/'qa_steps32.npz')['position_ig'];q64=np.load(out/'qa_steps64.npz')['position_ig']
        qa['32_64']=dict(relative_l1=float(np.abs(q32-q64).sum()/max(np.abs(q64).sum(),1e-12)),winner_agreement=float((np.argmax(np.lib.stride_tricks.sliding_window_view(q32,5,axis=1).sum(2),axis=1)==np.argmax(np.lib.stride_tricks.sliding_window_view(q64,5,axis=1).sum(2),axis=1)).mean()))
        save_json(out/'NUMERICAL_QA.json',qa);print('QA',task,seed,arm,qa,flush=True)
        if args.qa_only:return
    vals,cpl=ig(model,xx,args.steps,args.batch)
    np.savez_compressed(out/'attributions.npz',position_ig=vals,completeness=cpl,sequence=seq,conservation=x[:,:,4],chromosome=np.array(meta['chromosome']),source_row=np.array(meta['source_row']))
    rules=aggregate(seq,vals);save_json(out/'rules.json',rules)
    meta.update(initialization=initial,checkpoint_sha256=sha(checkpoint),arm=arm,task=task,seed=seed,steps=args.steps,batch=args.batch,scalar='log_total_count',integration='midpoint',model_eval_except_zero_dropout_GRU=True,elapsed_seconds=time.time()-started,mean_abs_completeness_error=float(np.abs(cpl[:,0]-cpl[:,1]).mean()),mean_abs_output_difference=float(np.abs(cpl[:,1]).mean()),script_sha256=sha(__file__))
    save_json(out/'summary.json',meta);print('DONE',task,seed,arm,round(time.time()-started,2),meta['mean_abs_completeness_error'],flush=True)
    del model,xx;torch.cuda.empty_cache()

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--task');ap.add_argument('--seed',type=int);ap.add_argument('--arm');ap.add_argument('--steps',type=int,default=32);ap.add_argument('--batch',type=int,default=16);ap.add_argument('--device',default='cuda');ap.add_argument('--qa',action='store_true');ap.add_argument('--qa-only',action='store_true');ap.add_argument('--overwrite',action='store_true');args=ap.parse_args()
    torch.set_num_threads(8)
    save_json(ROOT/('QA_ENVIRONMENT.json' if args.qa_only else 'RUN_ENVIRONMENT.json'),dict(start_utc=time.strftime('%Y-%m-%d %H:%M:%S',time.gmtime()),torch=torch.__version__,cuda_visible_devices=os.environ.get('CUDA_VISIBLE_DEVICES'),device=args.device,cuda_name=torch.cuda.get_device_name(0) if args.device=='cuda' else None,script_sha256=sha(__file__)))
    for task in ([args.task] if args.task else TASKS):
        for seed in ([args.seed] if args.seed else SEEDS):
            for arm in ([args.arm] if args.arm else ARMS):run(task,seed,arm,args)
