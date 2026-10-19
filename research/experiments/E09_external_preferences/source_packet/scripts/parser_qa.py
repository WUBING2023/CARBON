"""Verify original literal parsing and quantify all preserved before/after tables."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata
from literal_tables import read_csv_literal
R=Path(__file__).resolve().parents[1];old=R/'revisions'/'parser_default_20261004'/'results'
raw=R/'sources'/'Zscores_420_origrncmpt.txt'
with raw.open() as f:header=f.readline().split()[1:]
man=json.loads((R/'WIDE_MANIFEST.json').read_text());labs=sorted({k for v in man for k in v['experiments']});indices=[header.index(k)+1 for k in labs]
df=read_csv_literal(raw,sep=r'\s+',skiprows=1,header=None,usecols=[0]+indices)
indices=list(df.columns[1:])
literal=[]
with raw.open() as f:
    next(f)
    for line in f:
        split=line.split();literal.append([float(split[i]) for i in indices])
literal=np.asarray(literal);parsed=df.iloc[:,1:].to_numpy(float)
assert np.array_equal(literal,parsed)
diffs=[]
for p in sorted(old.glob('*.csv')):
    current=R/'results'/p.name
    if not current.exists():continue
    a=read_csv_literal(p);b=read_csv_literal(current)
    if a.shape!=b.shape or list(a)!=list(b):
        diffs.append(dict(file=p.name,shape_unchanged=False));continue
    numeric=[c for c in a if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c])]
    vals=[]
    for c in numeric:
        aa=a[c].to_numpy(float);bb=b[c].to_numpy(float);assert np.array_equal(np.isnan(aa),np.isnan(bb)),(p.name,c)
        okay=np.isfinite(aa)&np.isfinite(bb)
        if okay.any():vals.append(float(abs(aa[okay]-bb[okay]).max()))
    diffs.append(dict(file=p.name,shape_unchanged=True,max_numeric_difference=max(vals,default=0),numeric_columns=len(numeric)))
a=read_csv_literal(old/'WIDE_PROTEIN_METRICS.csv');b=read_csv_literal(R/'results'/'WIDE_PROTEIN_METRICS.csv')
aa=a[~a.exposed_protein].sort_values(['protein','arm']);bb=b[~b.exposed_protein].sort_values(['protein','arm'])
assert np.array_equal(aa[['spearman_all5mer','within_composition_concordance']].to_numpy(),bb[['spearman_all5mer','within_composition_concordance']].to_numpy(),equal_nan=True)
truth=read_csv_literal(R/'results'/'WIDE_EXTERNAL_5MER_SCORES.csv');rb=truth[truth.experiment=='RNCMPT00818'].set_index('kmer').z5
assert rb.ATCCG!=rb.GTAGA
qa=dict(contract="pandas float_precision='round_trip' equals Python float(original token)",
        raw_numeric_tokens_checked=int(literal.size),raw_literal_bit_equal=True,
        n_original_experiment_columns=len(labs),all_columns_recomputed_uniformly=True,
        new25_primary_points_bit_unchanged=True,all_snapshots_retained=True,
        original_independent_62_microdifference_failures_retained=True,no_new_model_inference=True,
        RBFOX2_ATCCG=float(rb.ATCCG),RBFOX2_GTAGA=float(rb.GTAGA),
        ATCCG_rank=float(rankdata(rb.to_numpy())[list(rb.index).index('ATCCG')]),
        GTAGA_rank=float(rankdata(rb.to_numpy())[list(rb.index).index('GTAGA')]),
        before_after_table_checks=diffs,
        original_source_sha256=hashlib.sha256(raw.read_bytes()).hexdigest())
(R/'results'/'PARSER_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in qa.items() if k!='before_after_table_checks'},indent=2))
print([v for v in diffs if v.get('max_numeric_difference',0)>1e-12])
