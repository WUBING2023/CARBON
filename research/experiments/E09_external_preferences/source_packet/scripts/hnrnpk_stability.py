from literal_tables import read_csv_literal
"""Paired true gene/locus bootstrap; word-level and fair tie diagnostics."""
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits
from evaluate_rules import ROOT,KMERS,ranking,truth,consensus_score

def membership(scores,k=20,eligible=None):
    values=np.asarray(scores,float);ok=np.ones(len(values),bool) if eligible is None else np.asarray(eligible,bool)
    idx=np.flatnonzero(ok);out=np.zeros(len(values));cut=np.sort(values[idx])[::-1][min(k,len(idx))-1]
    above=ok&(values>cut);at=ok&(values==cut);out[above]=1;out[at]=max(0,k-int(above.sum()))/max(1,int(at.sum()))
    return out

def run():
    rng=np.random.default_rng(20261004);B=2000
    coord=read_csv_literal(ROOT/'results'/'locus_mapping'/'HNRNPK_NATIVE_COORDINATES.csv');gkeys=sorted(coord.cluster.unique());G=len(gkeys);assert G==385
    gi={k:i for i,k in enumerate(gkeys)};weights=rng.multinomial(G,np.full(G,1/G),size=B).astype(np.float32)
    gt=read_csv_literal(ROOT/'results'/'CROSS_ASSAY_5MER_TRUTH.csv');gt=gt[gt.protein=='HNRNPK'].set_index('kmer').loc[KMERS];scores=gt.zscore5.to_numpy();gtmem=membership(scores);gt20=np.argsort(-scores,kind='stable')[:20]
    rbns,_=truth();rkeys=next(l for l in rbns['HNRNPK'] if l in set(read_csv_literal(ROOT/'results'/'EXTERNAL_COLUMNS.csv').query('primary==True').label));rbvals=np.array([rbns['HNRNPK'][rkeys][k] for k in KMERS]);rb20=np.argsort(-rbvals,kind='stable')[:20]
    boot={};observed=[];wordrows=[];diagnostics=[]
    for task,part in coord.groupby('task'):
        part=part.sort_values('selected_index');assert list(part.selected_index)==list(range(1024));cluster=np.array([gi[k] for k in part.cluster]);cell=task.split('_')[1]
        for seed in [11,23,47]:
            for arm in ['aligned','sequence_only','mean','permuted']:
                a=np.load(ROOT/'results'/'extracted'/'outputs'/task/f'seed{seed}'/arm/'attributions.npz');v=np.lib.stride_tricks.sliding_window_view(a['position_ig'],5,axis=1).sum(2);win=v.argmax(1);km=[s[i:i+5] for s,i in zip(a['sequence'],win)];ki=np.array([KMERS.index(k) for k in km]);value=v[np.arange(1024),win]
                score=np.zeros((G,1024),np.float32);cnt=np.zeros((G,1024),np.float32);np.add.at(score,(cluster,ki),value);np.add.at(cnt,(cluster,ki),1)
                with threadpool_limits(limits=8):bs=weights@score;bn=weights@cnt
                bs[bn==0]=-np.inf;order=np.argsort(-bs,axis=1,kind='stable')[:,:20];br=np.isin(order,gt20).sum(1)/20;brr=np.isin(order,rb20).sum(1)/20
                boot[(task,seed,arm)]=np.stack([br,brr],1)
                base=score.sum(0);n=cnt.sum(0);top=np.argsort(np.where(n>0,-base,np.inf),kind='stable')[:20]
                ordinary=np.isin(top,gt20).sum()/20;expected=float(np.dot(membership(base,eligible=n>0),gtmem)/20)
                observed.append(dict(task=task,cell=cell,seed=seed,arm=arm,recall20=ordinary,expected_tie_recall20=expected,bootstrap_low=float(np.quantile(br,.025)),bootstrap_high=float(np.quantile(br,.975))))
                for kk in gt20:
                    wordrows.append(dict(task=task,cell=cell,seed=seed,arm=arm,kmer=KMERS[kk],external_Z5=float(scores[kk]),C_count=KMERS[kk].count('C'),in_model_top20=bool(kk in top),primary_model_score=float(base[kk]),n_winning_windows=int(n[kk]),is_polyC4plus=KMERS[kk].count('C')>=4))
        # Frequency/Phi winner simple baselines; all four DNA-count baselines included.
        a=np.load(ROOT/'results'/'extracted'/'outputs'/task/'seed11'/'aligned'/'attributions.npz');freq=np.zeros(1024);phi=np.zeros(1024);pn=np.zeros(1024)
        for ss,ph in zip(a['sequence'],a['conservation']):
            pv=np.lib.stride_tricks.sliding_window_view(ph,5).mean(1);wi=int(np.argmax(pv));kk=KMERS.index(ss[wi:wi+5]);phi[kk]+=float(pv[wi]);pn[kk]+=1
            for offset in range(97):freq[KMERS.index(ss[offset:offset+5])]+=1
        for name,values,eligible in [('frequency',freq,None),('phi_priority',phi,pn>0)]+[(f'{b}_count',np.array([k.count(b) for k in KMERS]),None) for b in 'ACGT']:
            mask=np.ones(1024,bool) if eligible is None else eligible;top=np.argsort(np.where(mask,-values,np.inf),kind='stable')[:20]
            diagnostics.append(dict(task=task,arm=name,lexicographic_recall20=float(np.isin(top,gt20).sum()/20),expected_tie_recall20=float(np.dot(membership(values,eligible=mask),gtmem)/20),tie_rule='uniform among all boundary ties; true ties independent'))
    # Paired bootstrap averages model fits/cell contexts without increasing biological n.
    arrays={arm:np.stack([v for (task,seed,cond),v in boot.items() if cond==arm]).mean(0) for arm in ['aligned','sequence_only','mean','permuted']}
    diff=[]
    for label,i in [('RNAcompete',0),('RBNS',1)]:
        for ctrl in ['sequence_only','mean','permuted']:
            d=arrays['aligned'][:,i]-arrays[ctrl][:,i]
            source=read_csv_literal(ROOT/'results'/'CROSS_ASSAY_PROTEINS.csv').query("protein=='HNRNPK'") if i==0 else read_csv_literal(ROOT/'results'/'PROTEIN_METRICS.csv').query("protein=='HNRNPK'")
            key='IG_recall20' if i==0 else 'recall20';point=float(source.set_index('arm').loc['aligned',key]-source.set_index('arm').loc[ctrl,key])
            diff.append(dict(assay=label,control=ctrl,delta_recall20=point,conditional_gene_locus_bootstrap_low=float(np.quantile(d,.025)),conditional_gene_locus_bootstrap_high=float(np.quantile(d,.975)),bootstrap_positive_fraction=float((d>0).mean()),n_clusters=G,n_bootstrap=B,n_proteins=1,shared_cluster_resampling_across_cells=True))
    pd.DataFrame(observed).to_csv(ROOT/'results'/'HNRNPK_PRIMARY_STABILITY.csv',index=False);pd.DataFrame(wordrows).to_csv(ROOT/'results'/'HNRNPK_ALL_TOP20_WORDS.csv',index=False);pd.DataFrame(diagnostics).to_csv(ROOT/'results'/'HNRNPK_SIMPLE_BASELINE_TIES.csv',index=False);pd.DataFrame(diff).to_csv(ROOT/'results'/'HNRNPK_PAIRED_GENE_BOOTSTRAP.csv',index=False)
    np.savez_compressed(ROOT/'results'/'HNRNPK_BOOTSTRAP_DRAWS.npz',**arrays)
    wd=pd.DataFrame(wordrows);decomp=wd.groupby(['arm','is_polyC4plus'],as_index=False).in_model_top20.mean();decomp.to_csv(ROOT/'results'/'HNRNPK_WORD_FAMILY_DECOMPOSITION.csv',index=False)
    meta=dict(n_bootstrap=B,n_clusters=G,rng_seed=20261004,cluster_definition='Gencode v29 shared gene plus interval overlap; conservative strand-independent union',n_proteins=1,scope='conditional rule-extraction uncertainty within the observed HNRNPK loci; not external biological/population confidence',all_seed_backgrounds_retained=True,do_not_choose_new_endpoint=True)
    (ROOT/'results'/'HNRNPK_BOOTSTRAP_QA.json').write_text(json.dumps(meta,indent=2));print(pd.DataFrame(diff).to_string(index=False));print(pd.DataFrame(diagnostics).groupby('arm')[['lexicographic_recall20','expected_tie_recall20']].mean().to_string())

if __name__=='__main__':run()
