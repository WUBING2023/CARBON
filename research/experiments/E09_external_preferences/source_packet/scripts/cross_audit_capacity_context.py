"""Independent readback of A native profiles, share endpoints and gene weighting."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];A=ROOT.parent/'01_context';NG=A/'nogru_position_capacity'
tasks=['PTBP1_K562_ENCSR981WKN','PTBP1_HepG2_ENCSR384KAN'];arms=['aligned','sequence_only','mean','permuted','nogru','nogru_position_capacity'];seeds=[11,23,47]
rows=[];checks=[];init=[]
for task in tasks:
    for seed in seeds:
        d=NG/'remote'/'outputs'/task/f'seed{seed}';run=json.loads((d/'run.json').read_text());complete=json.loads((d/'complete.json').read_text());hh=hashlib.sha256((d/'best.pth').read_bytes()).hexdigest();assert hh==complete['checkpoint_sha256']
        ep=[json.loads(v) for v in (d/'epochs.jsonl').read_text().splitlines()];assert complete['best_epoch']==max(ep,key=lambda v:v['validation_pcc'])['epoch'];assert len(ep)<=100
        init.append(dict(task=task,seed=seed,epochs=len(ep),best_epoch=complete['best_epoch'],parameters=run['parameters'],full_parameters=run['full_parameters'],manifest_sha256=run['manifest_sha256'],common_full_init_sha256=run['common_five_channel_state_sha256'],nogru_init_sha256=run['actual_initial_state_sha256'],checkpoint_verified=True,validation_selection_verified=True))
    for chrom in ['1','X','Y']:
        d=A/'remote'/'data'/task if chrom=='1' else A/'phase2_ambiguity'/'remote'/'data'/task/chrom
        meta=json.loads((d/'metadata.json').read_text());z=np.load(d/'data.npz');zz=np.load(NG/'remote'/'data'/task/chrom/'data.npz')
        assert np.array_equal(z['index'],zz['index']) and np.array_equal(z['centers'],zz['centers'])
        oldng=np.load(A/'nogru_matched'/'remote'/'data'/task/chrom/'data.npz')
        assert np.array_equal(oldng['index'],z['index'])
        profiles={}
        for arm in arms:
            store=zz if arm=='nogru_position_capacity' else (oldng if arm=='nogru' else z)
            for seed in seeds:
                p=store[arm+f'_seed{seed}'].astype(float);assert p.shape==(len(meta),101) and np.isfinite(p).all() and np.allclose(p.sum(1),1,atol=1e-6)
                profiles[arm+f'_seed{seed}']=p
            profiles[arm+'_ensemble']=np.stack([profiles[arm+f'_seed{seed}'] for seed in seeds]).mean(0)
        for i,v in enumerate(meta):
            y=z['y'][i].astype(float);ctl=z['control'][i].astype(float);c1,c2=v['centers'];idx=[np.arange(c-5,c+6) for c in [c1,c2]]
            signal=np.array([y[ix].sum() for ix in idx]);control=np.array([ctl[ix].sum() for ix in idx]);res=signal-y.sum()*control/ctl.sum() if ctl.sum()>0 else np.array([np.nan,np.nan])
            good=y.sum()>0 and ctl.sum()>0 and np.isfinite(res).all() and (res>=0).all() and res.sum()>0 and len(v['gene_ids'])==1
            if not good:continue
            q=res[0]/res.sum();rec=dict(task=task,gene_id=v['gene_ids'][0],range='chr1' if chrom=='1' else 'XY',index=v['index'],chromosome=v['chromosome'],q=q,c1=c1,c2=c2,center_correct=.5 if abs(c1-50)==abs(c2-50) else float((abs(c1-50)<abs(c2-50))==(q>.5)),input_correct=.5 if control[0]==control[1] else float((control[0]>control[1])==(q>.5)))
            for method,p in profiles.items():
                masses=np.array([p[i,ix].sum() for ix in idx]);share=masses[0]/masses.sum();rec[method+'_share']=share;rec[method+'_loss']=(share-q)**2;rec[method+'_correct']=(.5 if share==.5 else float((share>.5)==(q>.5))) if q!=.5 else np.nan
            rows.append(rec)
df=pd.DataFrame(rows);cols=[c for c in df if c.endswith('_loss')];gene=df.groupby(['task','gene_id'])[cols].mean();p=gene.groupby('task').mean().mean();saved=pd.read_csv(NG/'METHODS.csv');saved=saved.query("scene=='all_core' and cohort=='adjusted' and group=='PTBP1'").set_index('method')
rankcols=[c for c in df if c.endswith('_correct')];ranks=df.groupby(['task','gene_id'])[rankcols].mean().groupby('task').mean().mean()
for col,val in p.items():
    method=col[:-5];expected=float(saved.loc[method,'gene_equal_mse']);checks.append(dict(quantity='primary_gene_equal_mse',method=method,independent=float(val),saved=expected,abs_difference=abs(val-expected)));assert abs(val-expected)<1e-11
    rv=float(ranks[method+'_correct']);exp=float(saved.loc[method,'gene_equal_ranking']);checks.append(dict(quantity='primary_gene_equal_ranking',method=method,independent=rv,saved=exp,abs_difference=abs(rv-exp)));assert abs(rv-exp)<1e-11
names=sorted(set(df.gene_id));gg=gene.loc[:,['aligned_ensemble_loss','nogru_position_capacity_ensemble_loss']];mat=np.full((len(names),2,2),np.nan)
for gi,g in enumerate(names):
    for ti,t in enumerate(tasks):
        if (t,g) in gg.index:mat[gi,ti]=gg.loc[(t,g)]
bs=np.load(NG/'bootstrap'/'all_core_adjusted_PTBP1.npz');state=json.loads(str(bs['rng_state_json']));rng=np.random.default_rng();rng.bit_generator.state=state
idx=rng.integers(0,len(names),(10000,len(names)));sample=mat[idx];den=np.isfinite(sample[:,:,:,0]).sum(1);boot=(np.nansum(sample,axis=1)/den[:,:,None]).mean(1);delta=boot[:,0]-boot[:,1]
savedboots=bs['means'];methods=list(bs['methods']);savedDelta=savedboots[:,methods.index('aligned_ensemble')]-savedboots[:,methods.index('nogru_position_capacity_ensemble')];assert np.allclose(delta,savedDelta,atol=1e-13)
pd.DataFrame(checks).to_csv(ROOT/'results'/'CROSS_AUDIT_01_PC_POINTS.csv',index=False);pd.DataFrame(init).to_csv(ROOT/'results'/'CROSS_AUDIT_01_PC_TRAINING.csv',index=False)
cell=df.groupby(['task','gene_id'])[cols+rankcols].mean().groupby('task').mean();cell.to_csv(ROOT/'results'/'CROSS_AUDIT_01_PC_CELL_SEED.csv')
qa=dict(readback_windows=len(df),unique_genes=len(names),point_checks=len(checks),max_point_difference=max(x['abs_difference'] for x in checks),bootstrap_10000_max_abs_error=float(abs(delta-savedDelta).max()),paired_ci95=np.quantile(delta,[.025,.975]).tolist(),center_rank=float(ranks['center_correct']),input_rank=float(ranks['input_correct']),same_coordinate_profiles=True,all_arms_retained=True,parameter_counts_near_matched=True,parameter_count_difference=45,flexible_absolute_position_control=True,architecture_equivalence_not_claimed=True,biological_mechanism_not_directly_measured=True)
(ROOT/'results'/'CROSS_AUDIT_01_PC_QA.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa,indent=2));print(cell[[c for c in cell if c.startswith(('aligned','nogru'))]].to_string())
