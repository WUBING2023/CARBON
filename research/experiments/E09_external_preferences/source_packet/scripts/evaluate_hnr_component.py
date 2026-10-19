from literal_tables import read_csv_literal
from pathlib import Path
import sys,tarfile,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'scripts'))
from evaluate_rules import KMERS
from evaluate_controlled import pairs,concordance
target=R/'results'/'hnr_component_extracted';target.mkdir(exist_ok=True)
with tarfile.open(R/'results'/'hnrnpk_nogru_external_outputs.tar.gz') as f:
    for item in f:
        assert (target/item.name).resolve().is_relative_to(target.resolve());f.extract(item,target)
truth=read_csv_literal(R/'results'/'EXTERNAL_5MER_SCORES.csv');truth=truth[truth.protein=='HNRNPK']
extra=read_csv_literal(R/'results'/'CROSS_ASSAY_5MER_TRUTH.csv');extra=extra[extra.protein=='HNRNPK'].rename(columns={'zscore5':'score'})
truth=pd.concat([truth[['protein','label','kmer','score']],extra[['protein','label','kmer','score']]],ignore_index=True);pp,_=pairs();rows=[];audit=[];wordrows=[]
for task in sorted((target/'hnrnpk_nogru_external_outputs').iterdir()):
    for seed in [11,23,47]:
        d=task/f'seed{seed}';z=np.load(d/'responses.npz');info=json.loads((d/'summary.json').read_text());old=np.load(R/'results'/'controlled_extracted'/'controlled_outputs'/task.name/f'seed{seed}'/'aligned'/'responses.npz')
        assert list(z['kmers'])==KMERS and z['delta_log_count'].shape==(16,1024) and all(np.array_equal(z[k],old[k]) for k in ['sequence','source_row','chromosome','native_conservation','conditional_conservation'])
        pred=z['delta_log_count'].mean(0)
        for k,v in zip(KMERS,pred):wordrows.append(dict(task=task.name,protein='HNRNPK',cell=task.name.split('_')[1],seed=seed,arm='nogru',kmer=k,score=float(v)))
        for label,tab in truth.groupby('label'):
            actual=tab.set_index('kmer').loc[KMERS,'score'].to_numpy();cc,n,t,e=concordance(pred,actual,pp)
            rows.append(dict(task=task.name,protein='HNRNPK',cell=task.name.split('_')[1],seed=seed,arm='nogru',label=label,spearman_all5mer=float(spearmanr(pred,actual).statistic),within_composition_concordance=cc,n_pairs=n,prediction_ties=t,experimental_ties=e))
        audit.append(dict(task=task.name,seed=seed,same16_backgrounds=True,samePhi=True,n_words=1024,checkpoint_sha256=info['checkpoint_sha256'],elapsed_seconds=info['elapsed_seconds']))
prior=read_csv_literal(R/'results'/'CONTROLLED_ALL_METRICS.csv');prior=prior[(prior.protein=='HNRNPK')&prior.arm.isin(['aligned','sequence_only','mean','permuted','frequency','known_consensus','phi_priority'])]
cross=read_csv_literal(R/'results'/'CROSS_ASSAY_ALL_METRICS.csv');cross=cross[cross.protein=='HNRNPK'].rename(columns={'response_spearman':'spearman_all5mer','response_composition_concordance':'within_composition_concordance'})
df=pd.concat([pd.DataFrame(rows),prior,cross],ignore_index=True);df.to_csv(R/'results'/'HNR_COMPONENT_ALL_METRICS.csv',index=False)
bg=df.groupby(['label','cell','arm'],as_index=False)[['spearman_all5mer','within_composition_concordance']].mean();bg.to_csv(R/'results'/'HNR_COMPONENT_BACKGROUNDS.csv',index=False)
prot=bg.groupby(['label','arm'],as_index=False)[['spearman_all5mer','within_composition_concordance']].mean();prot.to_csv(R/'results'/'HNR_COMPONENT_PROTEIN_BY_ASSAY.csv',index=False)
pd.DataFrame(audit).to_csv(R/'results'/'HNR_COMPONENT_QA.csv',index=False);pd.DataFrame(wordrows).to_csv(R/'results'/'HNR_COMPONENT_WORD_SCORES.csv',index=False)
print(prot.to_string(index=False))
