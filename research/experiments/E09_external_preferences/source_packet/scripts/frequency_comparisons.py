from literal_tables import read_csv_literal
"""Retained additional Full and Sequence versus frequency contrasts."""
from pathlib import Path
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];df=read_csv_literal(R/'results'/'WIDE_PROTEIN_METRICS.csv');new=df[~df.exposed_protein];proteins=sorted(new.protein.unique());assert len(proteins)==25
rng=np.random.default_rng(20261004);rng.integers(0,28,(10000,28));ix=rng.integers(0,25,(10000,25));rows=[]
for metric in ['spearman_all5mer','within_composition_concordance']:
    table=new.pivot(index='protein',columns='arm',values=metric).loc[proteins]
    for arm in ['full','sequence']:
        d=(table[arm]-table.frequency).to_numpy();bs=d[ix].mean(1)
        rows.append(dict(subset='new25',endpoint=metric,model=arm,control='frequency',mean_delta=float(d.mean()),ci_low=float(np.quantile(bs,.025)),ci_high=float(np.quantile(bs,.975)),positive=int((d>0).sum()),ties=int((d==0).sum()),negative=int((d<0).sum()),n_proteins=25,n_bootstrap=10000,rng_seed=20261004))
pd.DataFrame(rows).to_csv(R/'results'/'WIDE_MODELS_VS_FREQUENCY.csv',index=False);print(pd.DataFrame(rows).to_string(index=False))
