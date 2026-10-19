import pathlib,json,itertools,tarfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,rankdata
from evaluate_rules import ROOT,KMERS,ARMS,truth,consensus_score

def pairs():
    groups={}
    for i,k in enumerate(KMERS):groups.setdefault(tuple(k.count(b) for b in 'ACGT'),[]).append(i)
    pair=np.array([p for members in groups.values() for p in itertools.combinations(members,2)])
    assert len(groups)==56
    return pair,groups

def concordance(pred,actual,pp):
    t=actual[pp[:,0]]-actual[pp[:,1]];p=pred[pp[:,0]]-pred[pp[:,1]];good=t!=0
    score=np.where(p[good]==0,.5,(np.sign(p[good])==np.sign(t[good])).astype(float))
    return float(score.mean()),int(good.sum()),int((p[good]==0).sum()),int((~good).sum())

def run():
    folder=ROOT/'results'/'controlled_extracted';folder.mkdir(exist_ok=True)
    if not (folder/'controlled_outputs').exists():
        with tarfile.open(ROOT/'results'/'controlled_outputs.tar.gz') as tf:
            for m in tf:
                assert (folder/m.name).resolve().is_relative_to(folder.resolve());tf.extract(m,folder)
    ts,primary=truth();pp,groups=pairs();meta=dict(composition_groups=56,possible_pairs=len(pp),criteria='same exact A/C/G/T counts; experimental non-tie; prediction tie=0.5',independent_unit='protein',n_proteins=3,native_contexts=16,kmers=1024,external_labels_exposed=True)
    rows=[];scorerows=[];audits=[]
    for task in sorted((folder/'controlled_outputs').iterdir()):
        protein,cell=task.name.split('_')[:2];id0=None
        for seed in [11,23,47]:
            for arm in ARMS:
                path=task/f'seed{seed}'/arm;a=np.load(path/'responses.npz');j=json.loads((path/'summary.json').read_text())
                assert list(a['kmers'])==KMERS and a['delta_log_count'].shape==(16,1024)
                ids=(a['sequence'].tobytes(),a['chromosome'].tobytes(),a['source_row'].tobytes())
                if id0 is None:id0=ids
                assert ids==id0
                pred=a['delta_log_count'].mean(0)
                for k,v in zip(KMERS,pred):scorerows.append(dict(task=task.name,protein=protein,cell=cell,seed=seed,arm=arm,kmer=k,mean_delta_log_count=float(v)))
                for label,tru in ts[protein].items():
                    actual=np.array([tru[k] for k in KMERS]);cc,n,ties,eties=concordance(pred,actual,pp)
                    rows.append(dict(task=task.name,protein=protein,cell=cell,seed=seed,arm=arm,label=label,primary_column=label in primary,spearman_all5mer=float(spearmanr(pred,actual).statistic),within_composition_concordance=cc,n_pairs=n,prediction_ties=ties,experimental_ties=eties,n_contexts=16,n_kmers=1024))
                phi=a['native_conservation'];cp=a['conditional_conservation']
                if arm=='aligned':assert np.array_equal(phi,cp)
                elif arm=='mean':assert np.allclose(cp,phi.mean(1,keepdims=True))
                elif arm=='permuted':assert np.array_equal(np.sort(phi,1),np.sort(cp,1))
                else:assert cp.shape==(0,101)
                audits.append(dict(task=task.name,seed=seed,arm=arm,contexts=16,kmers=1024,checkpoint_sha256=j['checkpoint_sha256'],elapsed_seconds=j['elapsed_seconds'],conservation_check=True,original_source_row_preserved=True))
        # Fixed simple baselines include native frequency/Phi associations from Candidate A.
        freq={k:0. for k in KMERS};ph={k:0. for k in KMERS};n={k:0 for k in KMERS}
        source=ROOT/'results'/'extracted'/'outputs'/task.name/'seed11'/'aligned'/'attributions.npz';orig=np.load(source)
        for ss,phi in zip(orig['sequence'],orig['conservation']):
            vals=np.lib.stride_tricks.sliding_window_view(phi,5).mean(1)
            for i in range(97):k=ss[i:i+5];freq[k]+=1;ph[k]+=float(vals[i]);n[k]+=1
        ph={k:ph[k]/n[k] for k in KMERS};cons={k:consensus_score(protein,k) for k in KMERS}
        for arm,dictscore in [('frequency',freq),('phi_priority',ph),('known_consensus',cons)]:
            pred=np.array([dictscore[k] for k in KMERS])
            for label,tru in ts[protein].items():
                actual=np.array([tru[k] for k in KMERS]);cc,num,tie,et=concordance(pred,actual,pp)
                rows.append(dict(task=task.name,protein=protein,cell=cell,seed=0,arm=arm,label=label,primary_column=label in primary,spearman_all5mer=float(spearmanr(pred,actual).statistic),within_composition_concordance=cc,n_pairs=num,prediction_ties=tie,experimental_ties=et,n_contexts=16,n_kmers=1024))
    df=pd.DataFrame(rows);df.to_csv(ROOT/'results'/'CONTROLLED_ALL_METRICS.csv',index=False)
    pd.DataFrame(scorerows).to_csv(ROOT/'results'/'CONTROLLED_5MER_SCORES.csv',index=False)
    pd.DataFrame(audits).to_csv(ROOT/'results'/'CONTROLLED_AUDIT.csv',index=False)
    pd.DataFrame([dict(kmer=k,**dict(zip('ACGT',counts)),composition_group=str(counts),group_size=len(members)) for counts,members in groups.items() for k in (KMERS[i] for i in members)]).to_csv(ROOT/'results'/'COMPOSITION_GROUPS.csv',index=False)
    bg=df[df.primary_column].groupby(['protein','cell','arm'],as_index=False)[['spearman_all5mer','within_composition_concordance']].mean();bg.to_csv(ROOT/'results'/'CONTROLLED_BACKGROUND_METRICS.csv',index=False)
    prot=bg.groupby(['protein','arm'],as_index=False)[['spearman_all5mer','within_composition_concordance']].mean();prot.to_csv(ROOT/'results'/'CONTROLLED_PROTEIN_METRICS.csv',index=False)
    comp=[]
    for p in prot.protein.unique():
        sub=prot[prot.protein==p].set_index('arm')
        for ctrl in ['sequence_only','mean','permuted','frequency','phi_priority','known_consensus']:
            comp.append(dict(protein=p,control=ctrl,delta_spearman=float(sub.loc['aligned','spearman_all5mer']-sub.loc[ctrl,'spearman_all5mer']),delta_concordance=float(sub.loc['aligned','within_composition_concordance']-sub.loc[ctrl,'within_composition_concordance'])))
    pd.DataFrame(comp).to_csv(ROOT/'results'/'CONTROLLED_DIFFERENCES.csv',index=False)
    (ROOT/'results'/'CONTROLLED_QA.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(prot.to_string(index=False))

if __name__=='__main__':run()
