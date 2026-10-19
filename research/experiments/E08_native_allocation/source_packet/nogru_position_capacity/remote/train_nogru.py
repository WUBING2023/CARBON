import sys
sys.path.insert(0,'/data/CARBON_PROJECT/conservation_upgrade_20260929_matched')
sys.path.insert(0,'/data/CARBON_PROJECT/supervised_expansion_20261004/01_context/nogru_position_capacity')
"""Matched retraining; no test-dependent selection. Supports exact epoch resume."""
import argparse, csv, fcntl, json, os, platform, sys, time
from pathlib import Path
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
import numpy as np
import torch
from matched_common import (ROOT, ARMS, sha, write_json, verify_source, load_split,
    make_model as _matched_make_model, set_seed, state_hash, CombinedProfileLoss, pcc_tensor,
    transformed, metrics, covariates, row_key)

def make_model(seed, arm):
    assert arm=='aligned'
    from control_model import NoGRUPositionCapacity
    full,original=_matched_make_model(seed,'aligned')
    model=NoGRUPositionCapacity()
    shared={k:(full.state_dict()[k].clone() if k in full.state_dict() else v.clone()) for k,v in model.state_dict().items()}
    model.load_state_dict(shared,strict=True)
    return model, {'common_five_channel_state_sha256':original['common_five_channel_state_sha256'],'actual_initial_state_sha256':state_hash(model.state_dict()),'parameters':sum(p.numel() for p in model.parameters()),'full_parameters':original['parameters'],'architecture':'noGRU pointwise capacity+position control; 130→216→128 zero-final MLP +101 profile bias'}

def save_torch(path,value):
    tmp=Path(str(path)+'.tmp'); torch.save(value,tmp); tmp.replace(path)

def reconcile_epoch_log(path,committed_epoch):
    if not path.exists(): return
    text=path.read_text(encoding='utf-8'); rows={}
    for line in text.splitlines():
        try: row=json.loads(line)
        except json.JSONDecodeError: continue
        if row['epoch']<=committed_epoch: rows[row['epoch']]=row
    corrected=''.join(json.dumps(rows[k])+'\n' for k in sorted(rows))
    if corrected!=text:
        path.with_name('epochs_uncommitted_%d.jsonl'%time.time_ns()).write_text(text,encoding='utf-8')
        tmp=path.with_suffix('.tmp'); tmp.write_text(corrected,encoding='utf-8'); tmp.replace(path)

@torch.no_grad()
def validate(model,x,y,criterion,bsz):
    model.eval(); loss=0.; score=0.; nscore=0; invalid=0
    for start in range(0,len(x),bsz):
        out=model(x[start:start+bsz]); yy=y[start:start+bsz]
        ll=criterion(out['logits_profile'],out['log_total_count'],yy)['loss']
        assert torch.isfinite(ll).item(),'Nonfinite validation loss'
        loss+=float(ll)*len(yy)
        r,good,observed=pcc_tensor(torch.softmax(out['logits_profile'],dim=1),yy)
        score+=float(torch.where(good,r,torch.zeros_like(r))[observed].sum())
        nscore+=int(observed.sum()); invalid+=int((observed&~good).sum())
    assert nscore>0
    return dict(loss=loss/len(x),pcc=score/nscore,target_valid=nscore,predicted_constant=invalid)

@torch.no_grad()
def predict(model,x,device,bsz=512):
    model.eval(); ps=[]; zs=[]; logs=[]
    for s in range(0,len(x),bsz):
        out=model(torch.as_tensor(np.ascontiguousarray(x[s:s+bsz]),device=device))
        ps.append(torch.softmax(out['logits_profile'],dim=1).cpu().numpy())
        logs.append(torch.log_softmax(out['logits_profile'],dim=1).cpu().numpy())
        zs.append(out['log_total_count'].cpu().numpy().ravel())
    return np.concatenate(ps),np.concatenate(zs),np.concatenate(logs)

def evaluation(model,manifest,task,arm,outdir,device):
    data=load_split(manifest,task,'test',verify=True)
    inputs=['sequence_only'] if arm=='sequence_only' else ['aligned','mean','permuted']
    if arm=='aligned': inputs += ['shift_%d'%d for d in [-25,-10,-5,-1,1,5,10,25]]
    summaries={}
    for inp in inputs:
        p,z,lp=predict(model,transformed(data['x'],data['permutation'],inp),device)
        assert np.isfinite(p).all() and np.isfinite(z).all()
        assert np.allclose(p.sum(1),1,atol=1e-6)
        mt=metrics(p,z,data['y'],logprob=lp)
        np.savez_compressed(outdir/('test_'+inp+'.npz'),profile=p,log_profile=lp,log_count=z,chromosome=data['chromosome'],source_row=data['source_row'],**mt)
        safe_mean=lambda v: float(np.nanmean(v)) if np.isfinite(v).any() else None
        summaries[inp]={'n':len(p),'target_valid':int(mt['target_valid'].sum()),'pcc_valid':int(mt['valid'].sum()),'predicted_constant':int((mt['target_valid']&~mt['valid']).sum()),'mean_pcc':safe_mean(mt['pcc']),'median_pcc':float(np.nanmedian(mt['pcc'])) if mt['valid'].any() else None,'mean_pcc_zero':safe_mean(mt['pcc_zero']),'mean_jsd':safe_mean(mt['jsd']),'mean_log_count_abs_error':safe_mean(mt['log_count_abs_error'])}
    write_json(outdir/'evaluation.json',summaries)
    return summaries

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--task',required=True); ap.add_argument('--seed',type=int,required=True)
    ap.add_argument('--arm',choices=ARMS,required=True); ap.add_argument('--output',required=True)
    args=ap.parse_args(); outdir=Path(args.output); outdir.mkdir(parents=True,exist_ok=True)
    run_lock=open(outdir/'run.lock','a')
    fcntl.flock(run_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    write_json(outdir/'process.json',{'pid':os.getpid(),'parent_pid':os.getppid(),'started_unix':time.time(),'argv':sys.argv})
    manifest=json.loads((ROOT/'manifest.json').read_text()); verify_source(manifest)
    assert args.seed in manifest['seeds'] and args.task in manifest['tasks']
    if (outdir/'complete.json').exists(): raise RuntimeError('Already completed; refuse overwrite')
    torch.set_num_threads(4)
    assert torch.cuda.is_available(),'Production experiment requires CUDA'
    device=torch.device('cuda:0'); started=time.time()
    model,initial=make_model(args.seed,args.arm); model.to(device)
    info={'task':args.task,'seed':args.seed,'arm':args.arm,'manifest_sha256':sha(ROOT/'manifest.json'),**initial,'python':sys.version,'torch':torch.__version__,'numpy':np.__version__,'host':platform.node(),'gpu':torch.cuda.get_device_name(0),'cuda_visible_devices':os.environ.get('CUDA_VISIBLE_DEVICES'),'started_unix':started}
    if (outdir/'run.json').exists():
        previous=json.loads((outdir/'run.json').read_text())
        for k in ('task','seed','arm','manifest_sha256','actual_initial_state_sha256','gpu','cuda_visible_devices'): assert info[k]==previous[k],k
    else: write_json(outdir/'run.json',info)
    print(json.dumps({'event':'loading','task':args.task,'seed':args.seed,'arm':args.arm}),flush=True)
    tr=load_split(manifest,args.task,'train',args.arm,verify=True)
    va=load_split(manifest,args.task,'validation',args.arm,verify=True)
    tx=torch.as_tensor(tr['x'],device=device); ty=torch.as_tensor(tr['y'],device=device)
    vx=torch.as_tensor(va['x'],device=device); vy=torch.as_tensor(va['y'],device=device)
    del tr,va
    opt=torch.optim.AdamW(model.parameters(),lr=.000785,weight_decay=.000114)
    criterion=CombinedProfileLoss(alpha=1.42,beta=.54,gamma=1.80,loss_scale=1.0)
    generator=torch.Generator(device='cpu').manual_seed(args.seed+10000)
    set_seed(args.seed+20000)
    best_pcc=-float('inf'); best_loss=float('inf'); bad_pcc=bad_loss=0; best_epoch=0; epoch0=0; elapsed=0.; finished=False
    resume=outdir/'resume.pth'
    if resume.exists():
        ck=torch.load(resume,map_location='cpu',weights_only=False)
        model.load_state_dict(ck['model'],strict=True); opt.load_state_dict(ck['optimizer'])
        for state in opt.state.values():
            for k,v in state.items():
                if torch.is_tensor(v): state[k]=v.to(device)
        generator.set_state(ck['sampler']); torch.set_rng_state(ck['rng_cpu']); torch.cuda.set_rng_state(ck['rng_cuda'])
        np.random.set_state(ck['numpy_rng'])
        epoch0=ck['epoch']; best_pcc=ck['best_pcc']; best_loss=ck['best_loss']; best_epoch=ck['best_epoch']
        bad_pcc=ck['bad_pcc']; bad_loss=ck['bad_loss']; elapsed=ck['elapsed']; finished=ck['finished']
        print(json.dumps({'event':'resume','next_epoch':epoch0+1,'training_finished':finished}),flush=True)
    reconcile_epoch_log(outdir/'epochs.jsonl',epoch0)
    bsz=256
    for epoch in range(epoch0+1,101) if not finished else []:
        epoch_start=time.time(); model.train(); losses=torch.zeros((),device=device)
        order=torch.randperm(len(tx),generator=generator); order_hash=sha_bytes(order.numpy().tobytes())
        ntrain=(len(tx)//bsz)*bsz; order=order[:ntrain].to(device)
        for start in range(0,ntrain,bsz):
            idx=order[start:start+bsz]; pred=model(tx[idx]); yy=ty[idx]
            loss=criterion(pred['logits_profile'],pred['log_total_count'],yy)['loss']
            if not torch.isfinite(loss).item(): raise RuntimeError('Nonfinite training loss')
            opt.zero_grad(set_to_none=True); loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(),2.0,error_if_nonfinite=True)
            opt.step(); losses+=loss.detach()*bsz
        v=validate(model,vx,vy,criterion,bsz)
        if v['pcc']>best_pcc:
            best_pcc=v['pcc']; bad_pcc=0; best_epoch=epoch
            save_torch(outdir/'best.pth',{k:val.detach().cpu() for k,val in model.state_dict().items()})
        else: bad_pcc+=1
        if v['loss']<best_loss: best_loss=v['loss']; bad_loss=0
        else:
            bad_loss+=1
            if bad_loss>=5:
                for group in opt.param_groups: group['lr']/=2
                bad_loss=0
        finished=bad_pcc>=9 or opt.param_groups[0]['lr']<1e-7 or epoch==100
        torch.cuda.synchronize(); secs=time.time()-epoch_start; elapsed+=secs
        line={'epoch':epoch,'train_loss':float(losses)/ntrain,'validation_loss':v['loss'],'validation_pcc':v['pcc'],'validation_target_valid':v['target_valid'],'validation_predicted_constant':v['predicted_constant'],'best_epoch':best_epoch,'lr':opt.param_groups[0]['lr'],'seconds':secs,'train_order_sha256':order_hash,'training_finished':finished}
        with open(outdir/'epochs.jsonl','a',encoding='utf-8') as f: f.write(json.dumps(line)+'\n')
        save_torch(resume,dict(model=model.state_dict(),optimizer=opt.state_dict(),epoch=epoch,best_pcc=best_pcc,best_loss=best_loss,best_epoch=best_epoch,bad_pcc=bad_pcc,bad_loss=bad_loss,elapsed=elapsed,finished=finished,sampler=generator.get_state(),rng_cpu=torch.get_rng_state(),rng_cuda=torch.cuda.get_rng_state(),numpy_rng=np.random.get_state()))
        write_json(outdir/'progress.json',line)
        print(json.dumps({'event':'epoch',**line}),flush=True)
        if finished: break
    del tx,ty,vx,vy; torch.cuda.empty_cache()
    model.load_state_dict(torch.load(outdir/'best.pth',map_location=device,weights_only=True),strict=True)
    print(json.dumps({'event':'evaluating','best_epoch':best_epoch}),flush=True)
    # Heldout evaluation deferred; no test labels used to select the checkpoint.
    write_json(outdir/'complete.json',{'task':args.task,'seed':args.seed,'arm':args.arm,'best_epoch':best_epoch,'best_validation_pcc':best_pcc,'training_seconds':elapsed,'checkpoint_sha256':sha(outdir/'best.pth'),'final_state_sha256':state_hash(model.state_dict()),'manifest_sha256':sha(ROOT/'manifest.json'),'completed_unix':time.time()})
    print(json.dumps({'event':'complete','task':args.task,'arm':args.arm,'seed':args.seed}),flush=True)

def sha_bytes(b):
    import hashlib
    return hashlib.sha256(b).hexdigest()

if __name__=='__main__': main()
