from literal_tables import read_csv_literal
"""Fixed occurrence/unique-containing diagnostics; primary analysis remains unchanged."""
import json
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from evaluate_rules import ROOT,KMERS,ranking,consensus_score
from evaluate_controlled import pairs,concordance
from hnrnpk_stability import membership

def run():
    labels={'RNCMPT00026':'HNRNPK','RNCMPT00268':'PTBP1','RNCMPT00269':'PTBP1','RNCMPT00818':'RBFOX2'}
    p=ROOT/'sources'/'Zscores_420_origrncmpt.txt'
    with p.open() as f:header=f.readline().split()[1:]
    cols=[0]+[header.index(k)+1 for k in labels];df=read_csv_literal(p,sep=r'\s+',skiprows=1,header=None,usecols=cols);words=df.iloc[:,0].str.replace('U','T')
    true={};td=[]
    for ii,(label,protein) in enumerate(labels.items(),1):
        for method in ['occurrence','unique_containing']:
            val={k:0. for k in KMERS};n={k:0 for k in KMERS}
            for word,z in zip(words,df.iloc[:,ii]):
                keys=[word[j:j+5] for j in range(3)];keys=set(keys) if method=='unique_containing' else keys
                for k in keys:val[k]+=float(z);n[k]+=1
            x=np.array([val[k]/n[k] for k in KMERS]);true[(label,method)]=x
            for k,v in zip(KMERS,x):td.append(dict(protein=protein,label=label,mapping=method,kmer=k,z5=float(v),n_measured_units=n[k]))
    ig=read_csv_literal(ROOT/'results'/'MODEL_5MER_SCORES.csv');resp=read_csv_literal(ROOT/'results'/'CONTROLLED_5MER_SCORES.csv');pp,_=pairs();rows=[];rankrows=[]
    for label,protein in labels.items():
        orig=true[(label,'occurrence')];alt=true[(label,'unique_containing')];oldtop=np.argsort(-orig,kind='stable')[:20];newtop=np.argsort(-alt,kind='stable')[:20]
        rankrows.append(dict(protein=protein,label=label,mapping_score_spearman=float(spearmanr(orig,alt).statistic),top20_overlap=int(np.isin(oldtop,newtop).sum()),denominator=20))
        for key,part in ig[ig.protein==protein].groupby(['task','cell','seed','arm']):
            task,cell,seed,arm=key;part=part.set_index('kmer').loc[KMERS];pred=part.mean_occurrence_ig.to_numpy();win=part.winner_sum_ig.to_numpy();eligible=part.winner_n.to_numpy()>0;top=np.argsort(np.where(eligible,-win,np.inf),kind='stable')[:20];rr=resp[(resp.task==task)&(resp.seed==seed)&(resp.arm==arm)].set_index('kmer').loc[KMERS].mean_delta_log_count.to_numpy()
            for method in ['occurrence','unique_containing']:
                actual=true[(label,method)];gt=np.argsort(-actual,kind='stable')[:20];mem=membership(actual);cc,n,t,e=concordance(rr,actual,pp)
                rows.append(dict(task=task,protein=protein,cell=cell,seed=seed,arm=arm,label=label,mapping=method,IG_recall20=float(np.isin(top,gt).sum()/20),expected_tie_IG_recall20=float(np.dot(membership(win,eligible=eligible),mem)/20),IG_spearman=float(spearmanr(pred,actual).statistic),response_spearman=float(spearmanr(rr,actual).statistic),within_composition_concordance=cc))
        for base in 'ACGT':
            values=np.array([k.count(base) for k in KMERS]);top=np.argsort(-values,kind='stable')[:20]
            for method in ['occurrence','unique_containing']:
                actual=true[(label,method)];gt=np.argsort(-actual,kind='stable')[:20];cc,n,t,e=concordance(values,actual,pp)
                rows.append(dict(task='universal',protein=protein,cell='none',seed=0,arm=f'{base}_count',label=label,mapping=method,IG_recall20=float(np.isin(top,gt).sum()/20),expected_tie_IG_recall20=float(np.dot(membership(values),membership(actual))/20),IG_spearman=float(spearmanr(values,actual).statistic),response_spearman=float(spearmanr(values,actual).statistic),within_composition_concordance=cc))
    pd.DataFrame(td).to_csv(ROOT/'results'/'MAPPING_SENSITIVITY_TRUTH.csv',index=False);pd.DataFrame(rows).to_csv(ROOT/'results'/'MAPPING_SENSITIVITY_METRICS.csv',index=False);pd.DataFrame(rankrows).to_csv(ROOT/'results'/'MAPPING_SENSITIVITY_TRUTH_RANKS.csv',index=False)
    result=pd.DataFrame(rows).groupby(['protein','mapping','arm'],as_index=False)[['IG_recall20','expected_tie_IG_recall20','IG_spearman','response_spearman','within_composition_concordance']].mean();result.to_csv(ROOT/'results'/'MAPPING_SENSITIVITY_SUMMARY.csv',index=False)
    print(pd.DataFrame(rankrows).to_string(index=False));print(result[result.protein=='HNRNPK'].to_string(index=False))

if __name__=='__main__':run()
