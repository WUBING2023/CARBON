from literal_tables import read_csv_literal
import json,itertools
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from evaluate_rules import ROOT,KMERS,ranking,consensus_score
from evaluate_controlled import pairs,concordance

def run():
    file=ROOT/'sources'/'Zscores_420_origrncmpt.txt'
    with file.open() as f:header=f.readline().split()[1:]
    labels={'RNCMPT00026':'HNRNPK','RNCMPT00818':'RBFOX2'}
    cols=[0]+[header.index(k)+1 for k in labels]
    df=read_csv_literal(file,sep=r'\s+',skiprows=1,header=None,usecols=cols);keys=df.iloc[:,0].astype(str).str.replace('U','T')
    byprotein={};rows=[];sc=[]
    for ii,(label,protein) in enumerate(labels.items(),1):
        total={k:0. for k in KMERS};n={k:0 for k in KMERS}
        for word,val in zip(keys,df.iloc[:,ii]):
            assert np.isfinite(val)
            for offset in range(3):k=word[offset:offset+5];total[k]+=float(val);n[k]+=1
        actual=np.array([total[k]/n[k] for k in KMERS]);byprotein[protein]=(label,actual)
        for k,v in zip(KMERS,actual):sc.append(dict(protein=protein,label=label,kmer=k,zscore5=float(v),available_occurrences=n[k]))
    ig=read_csv_literal(ROOT/'results'/'MODEL_5MER_SCORES.csv');resp=read_csv_literal(ROOT/'results'/'CONTROLLED_5MER_SCORES.csv');pp,_=pairs()
    for key,part in ig.groupby(['task','protein','cell','seed','arm']):
        task,protein,cell,seed,arm=key
        if protein not in byprotein:continue
        label,actual=byprotein[protein];part=part.set_index('kmer').loc[KMERS];rr=resp[(resp.task==task)&(resp.seed==seed)&(resp.arm==arm)].set_index('kmer').loc[KMERS]
        chosen=sorted([k for k in KMERS if part.loc[k,'winner_n']>0],key=lambda k:(-part.loc[k,'winner_sum_ig'],k))[:20]
        actual_top=ranking(dict(zip(KMERS,actual)))[:20]
        pred=rr.mean_delta_log_count.to_numpy();cc,num,tie,et=concordance(pred,actual,pp)
        rows.append(dict(task=task,protein=protein,cell=cell,seed=seed,arm=arm,label=label,IG_recall20=len(set(chosen)&set(actual_top))/20,IG_spearman=float(spearmanr(part.mean_occurrence_ig,actual).statistic),response_spearman=float(spearmanr(pred,actual).statistic),response_composition_concordance=cc,n_pairs=num,prediction_ties=tie,experimental_ties=et))
    for task in sorted(ig.task.unique()):
        protein,cell=task.split('_')[:2]
        if protein not in byprotein:continue
        label,actual=byprotein[protein];actual_top=ranking(dict(zip(KMERS,actual)))[:20]
        ori=np.load(ROOT/'results'/'extracted'/'outputs'/task/'seed11'/'aligned'/'attributions.npz')
        fr={k:0. for k in KMERS};ph={k:0. for k in KMERS};nn={k:0 for k in KMERS};winner={k:0. for k in KMERS};wn={k:0 for k in KMERS}
        for word,pvec in zip(ori['sequence'],ori['conservation']):
            vv=np.lib.stride_tricks.sliding_window_view(pvec,5).mean(1);at=int(np.argmax(vv));winner[word[at:at+5]]+=float(vv[at]);wn[word[at:at+5]]+=1
            for i in range(97):k=word[i:i+5];fr[k]+=1;ph[k]+=float(vv[i]);nn[k]+=1
        ph={k:ph[k]/nn[k] for k in KMERS};con={k:consensus_score(protein,k) for k in KMERS}
        for arm,pri,sec,eligible in [('frequency',fr,fr,KMERS),('known_consensus',con,con,KMERS),('phi_priority',winner,ph,[k for k in KMERS if wn[k]>0])]:
            top=ranking(pri,eligible)[:20];pred=np.array([sec[k] for k in KMERS]);cc,num,tie,et=concordance(pred,actual,pp);rho=float(spearmanr(pred,actual).statistic)
            rows.append(dict(task=task,protein=protein,cell=cell,seed=0,arm=arm,label=label,IG_recall20=len(set(top)&set(actual_top))/20,IG_spearman=rho,response_spearman=rho,response_composition_concordance=cc,n_pairs=num,prediction_ties=tie,experimental_ties=et))
    pd.DataFrame(sc).to_csv(ROOT/'results'/'CROSS_ASSAY_5MER_TRUTH.csv',index=False)
    df=pd.DataFrame(rows);df.to_csv(ROOT/'results'/'CROSS_ASSAY_ALL_METRICS.csv',index=False)
    bg=df.groupby(['protein','cell','arm'],as_index=False)[['IG_recall20','IG_spearman','response_spearman','response_composition_concordance']].mean();bg.to_csv(ROOT/'results'/'CROSS_ASSAY_BACKGROUNDS.csv',index=False)
    prot=bg.groupby(['protein','arm'],as_index=False)[['IG_recall20','IG_spearman','response_spearman','response_composition_concordance']].mean();prot.to_csv(ROOT/'results'/'CROSS_ASSAY_PROTEINS.csv',index=False)
    print(prot.to_string(index=False))

if __name__=='__main__':run()
