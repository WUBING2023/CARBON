"""Shared, result-independent transforms and metrics for matched retraining."""
import hashlib, json, os, random, sys
from pathlib import Path
import numpy as np
import h5py
import torch

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'source'))
from MARCP.models.marcp_model import MARCP
from MARCP.loss.loss import CombinedProfileLoss as OriginalCombinedProfileLoss
ARMS = ('aligned','mean','permuted','sequence_only')

class CombinedProfileLoss(OriginalCombinedProfileLoss):
    """Same CDF-MSE objective, using deterministic triangular matrix products.

    Torch 2.4.1 CUDA cumsum refuses strict determinism. The original source is
    preserved; this wrapper changes only the cumulative-sum implementation.
    """
    def jitter_loss(self,logits_profile,y_true):
        p=torch.softmax(logits_profile,dim=-1)
        q=y_true/(y_true.sum(dim=-1,keepdim=True)+1e-6)
        key=(p.shape[-1],p.device,p.dtype)
        if getattr(self,'_cdf_key',None)!=key:
            self._cdf_upper=torch.triu(torch.ones((p.shape[-1],p.shape[-1]),device=p.device,dtype=p.dtype))
            self._cdf_key=key
        return torch.nn.functional.mse_loss(p@self._cdf_upper,q@self._cdf_upper)

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def write_json(path, obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False),encoding='utf-8')
    tmp.replace(path)

def row_key(task, chromosome, row):
    return hashlib.sha256(('20260929|%s|%s|%s'%(task,chromosome,row)).encode()).digest()

def permutations(task, chromosome, n):
    p=np.empty((n,101),dtype=np.uint8)
    for i in range(n):
        p[i]=np.random.default_rng(int.from_bytes(row_key(task,chromosome,i)[:8],'big')).permutation(101)
    return p

def set_seed(seed):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic=True
    torch.backends.cudnn.benchmark=False
    torch.backends.cuda.matmul.allow_tf32=False
    torch.backends.cudnn.allow_tf32=False
    torch.use_deterministic_algorithms(True)

def state_hash(state):
    h=hashlib.sha256()
    for k,v in sorted(state.items()):
        h.update(k.encode()); h.update(str(tuple(v.shape)).encode())
        h.update(v.detach().cpu().contiguous().numpy().tobytes())
    return h.hexdigest()

def make_model(seed, arm):
    set_seed(seed)
    full=MARCP(in_channels=5,dropOut=0.2,control=False)
    initial={k:v.clone() for k,v in full.state_dict().items()}
    parent_hash=state_hash(initial)
    if arm=='sequence_only':
        model=MARCP(in_channels=4,dropOut=0.2,control=False)
        target=model.state_dict()
        for k,v in target.items():
            src=initial[k]
            if src.shape!=v.shape:
                assert k in ['motif.branches.%d.0.weight'%i for i in range(3)]+['motif.shortcut.weight'],k
                assert src.shape[1]==5 and v.shape[1]==4
                src=src[:,:4,:]
            target[k]=src.clone()
        model.load_state_dict(target,strict=True)
    else: model=full
    return model, {'common_five_channel_state_sha256':parent_hash,'actual_initial_state_sha256':state_hash(model.state_dict()),'parameters':sum(p.numel() for p in model.parameters())}

def verify_source(manifest):
    for p,h in manifest['execution_sha256'].items():
        assert sha(ROOT/p)==h, 'Frozen source changed: '+p
    assert sha(ROOT/'PROTOCOL.md')==manifest['protocol_sha256']

def load_split(manifest, task, split, arm='aligned', verify=False):
    xs=[]; ys=[]; ctls=[]; chrs=[]; rows=[]; ps=[]
    for spec in manifest['tasks'][task]['inputs'][split]:
        p=Path(spec['path'])
        if verify: assert sha(p)==spec['sha256'],str(p)
        with h5py.File(p,'r') as h:
            x=np.asarray(h['X'],dtype=np.float32)
            y=np.asarray(h['Y'],dtype=np.float32)
            ctl=np.asarray(h['Y_ctl'],dtype=np.float32)
        assert x.shape==(spec['n'],101,5) and y.shape==ctl.shape==(spec['n'],101)
        assert np.isfinite(x).all() and np.isfinite(y).all() and np.isfinite(ctl).all()
        assert (y>=0).all() and (ctl>=0).all()
        perm_path=ROOT/spec['permutation_path']
        assert sha(perm_path)==spec['permutation_sha256']
        pm=np.load(perm_path)
        if arm=='mean': x[:,:,4]=x[:,:,4].mean(1,keepdims=True)
        elif arm=='permuted': x[:,:,4]=np.take_along_axis(x[:,:,4],pm,axis=1)
        elif arm=='sequence_only': x=x[:,:,:4]
        xs.append(x); ys.append(y); ctls.append(ctl); ps.append(pm)
        chrs.append(np.full(len(x),spec['chromosome'],dtype=np.int16))
        rows.append(np.arange(len(x),dtype=np.int32))
    return dict(x=np.ascontiguousarray(np.concatenate(xs)),y=np.concatenate(ys),control=np.concatenate(ctls),chromosome=np.concatenate(chrs),source_row=np.concatenate(rows),permutation=np.concatenate(ps))

def transformed(x, perm, arm):
    if arm=='sequence_only': return x[:,:,:4].copy()
    z=x.copy()
    if arm=='mean': z[:,:,4]=z[:,:,4].mean(1,keepdims=True)
    elif arm=='permuted': z[:,:,4]=np.take_along_axis(z[:,:,4],perm,axis=1)
    elif arm.startswith('shift_'): z[:,:,4]=np.roll(z[:,:,4],int(arm.split('_')[1]),axis=1)
    elif arm!='aligned': raise ValueError(arm)
    return z

def pcc_tensor(prob,y):
    a=prob-prob.mean(1,keepdim=True); b=y-y.mean(1,keepdim=True)
    aa=(a*a).sum(1); bb=(b*b).sum(1)
    target_valid=y.amax(1)>y.amin(1)
    good=(prob.amax(1)>prob.amin(1))&target_valid
    r=(a*b).sum(1)/(aa*bb).sqrt().clamp_min(1e-30)
    return r.clamp(-1,1), good, target_valid

def metrics(prob,z,y,logprob=None):
    p=np.asarray(prob,dtype=np.float64); y=np.asarray(y,dtype=np.float64)
    a=p-p.mean(1,keepdims=True); b=y-y.mean(1,keepdims=True)
    den=np.sqrt((a*a).sum(1)*(b*b).sum(1))
    target_valid=np.ptp(y,axis=1)>0
    valid=(np.ptp(p,axis=1)>0)&target_valid&(den>0)
    r=np.full(len(y),np.nan); r[valid]=np.clip((a*b).sum(1)[valid]/den[valid],-1,1)
    total=y.sum(1); nz=total>0
    q=y/np.maximum(total[:,None],1e-30); mix=(p+q)/2
    eps=1e-30
    ce=-(q*(np.log(np.maximum(p,eps)) if logprob is None else np.asarray(logprob,dtype=np.float64))).sum(1)
    jsd=.5*((q*np.log(np.maximum(q,eps)/np.maximum(mix,eps))).sum(1)+(p*np.log(np.maximum(p,eps)/np.maximum(mix,eps))).sum(1))
    ce[~nz]=np.nan; jsd[~nz]=np.nan
    return dict(pcc=r,valid=valid,target_valid=target_valid,pcc_zero=np.where(target_valid,np.nan_to_num(r,nan=0),np.nan),jsd=jsd,cross_entropy=ce,log_count_abs_error=np.abs(np.asarray(z).ravel()-np.log1p(total)))

def covariates(data):
    x=data['x']; y=data['y']; phi=x[:,:,4]; composition=x[:,:,:4].mean(1)
    return dict(mean_phylop=phi.mean(1),sd_phylop=phi.std(1),signal_total=y.sum(1),signal_max=y.max(1),control_total=data['control'].sum(1),gc=composition[:,1]+composition[:,2],sequence_entropy=-(composition*np.log(np.maximum(composition,1e-30))).sum(1))
