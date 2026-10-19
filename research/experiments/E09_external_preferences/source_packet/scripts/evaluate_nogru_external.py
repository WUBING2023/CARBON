from literal_tables import read_csv_literal
import pathlib,tarfile,json,itertools,sys
import numpy as np,pandas as pd
from scipy.stats import spearmanr
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent))
from evaluate_controlled import pairs,concordance
out=ROOT/'results'/'nogru_external_extracted';out.mkdir(exist_ok=True)
with tarfile.open(ROOT/'results'/'nogru_external_outputs.tar.gz') as f:
    for item in f:
        assert (out/item.name).resolve().is_relative_to(out.resolve());f.extract(item,out)
truth=read_csv_literal(ROOT/'results'/'EXTERNAL_5MER_SCORES.csv');truth=truth[truth.protein=='PTBP1']
kmers=[''.join(k) for k in itertools.product('ACGT',repeat=5)];pp,_=pairs();rows=[];audit=[];words=[]
for task in sorted((out/'nogru_external_outputs').iterdir()):
    for seed in [11,23,47]:
        folder=task/f'seed{seed}';z=np.load(folder/'responses.npz');j=json.loads((folder/'summary.json').read_text())
        full=np.load(ROOT/'results'/'controlled_extracted'/'controlled_outputs'/task.name/f'seed{seed}'/'aligned'/'responses.npz')
        assert list(z['kmers'])==kmers and z['delta_log_count'].shape==(16,1024)
        assert all(np.array_equal(z[k],full[k]) for k in ['sequence','source_row','chromosome','native_conservation','conditional_conservation'])
        pred=z['delta_log_count'].mean(0)
        for k,v in zip(kmers,pred):words.append(dict(task=task.name,seed=seed,arm='nogru',kmer=k,score=float(v)))
        for label,tab in truth.groupby('label'):
            actual=tab.set_index('kmer').loc[kmers,'score'].to_numpy();cc,n,t,e=concordance(pred,actual,pp)
            rows.append(dict(task=task.name,protein='PTBP1',cell=task.name.split('_')[1],seed=seed,arm='nogru',label=label,spearman_all5mer=float(spearmanr(pred,actual).statistic),within_composition_concordance=cc,n_pairs=n,prediction_ties=t,experimental_ties=e))
        audit.append(dict(task=task.name,seed=seed,checkpoint_sha256=j['checkpoint_sha256'],same_native_backgrounds=True,same_phyloP=True,n_backgrounds=16,n_words=1024,elapsed_seconds=j['elapsed_seconds']))
allprior=read_csv_literal(ROOT/'results'/'CONTROLLED_ALL_METRICS.csv');allprior=allprior[(allprior.protein=='PTBP1') & allprior.arm.isin(['aligned','sequence_only','mean','permuted','frequency','known_consensus','phi_priority'])]
df=pd.concat([pd.DataFrame(rows),allprior],ignore_index=True);df.to_csv(ROOT/'results'/'NOGRU_EXTERNAL_ALL_METRICS.csv',index=False)
bg=df.groupby(['cell','arm'],as_index=False)[['spearman_all5mer','within_composition_concordance']].mean();bg.to_csv(ROOT/'results'/'NOGRU_EXTERNAL_BACKGROUNDS.csv',index=False)
bg.groupby('arm',as_index=False)[['spearman_all5mer','within_composition_concordance']].mean().to_csv(ROOT/'results'/'NOGRU_EXTERNAL_PROTEIN.csv',index=False)
pd.DataFrame(audit).to_csv(ROOT/'results'/'NOGRU_EXTERNAL_QA.csv',index=False);pd.DataFrame(words).to_csv(ROOT/'results'/'NOGRU_EXTERNAL_WORD_SCORES.csv',index=False)
print(bg.to_string(index=False));print(bg.groupby('arm')[['spearman_all5mer','within_composition_concordance']].mean().to_string())
