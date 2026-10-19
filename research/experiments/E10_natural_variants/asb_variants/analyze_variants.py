"""Recompute natural variant diagnostics from saved original arrays. CPU <=4."""
import os
os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['OMP_NUM_THREADS']='4';os.environ['MKL_NUM_THREADS']='4'
import sys
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import json,re,hashlib,datetime,gzip
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr
R=Path(__file__).resolve().parent
P=R/'source_inputs'
if not P.exists():P=Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/expanded_functional_validation_20261003/06_ASB_expansion')
event_file=P/'EVENTS_INCLUDED.json'
E=json.loads(event_file.read_text('utf8') if event_file.exists() else gzip.decompress((P/'EVENTS_INCLUDED.json.gz').read_bytes()).decode('utf8'))
Q=np.load(P/'PREDICTIONS.npz',allow_pickle=False)
rows=[]
for e in E:
    c=e['context']; ref=e['hg38_reference_plus']; maj=e['normalized_major_hg38_plus']; alt=e['normalized_minor_hg38_plus'] if ref==maj else maj
    seq=c['reference_sequence_rna_direction_dna']
    d=dict(event_id=e['event_id'],variant_key=f"{e['chrom']}:{e['pos0']+1}:{ref}:{alt}",chrom=e['chrom'],pos0=e['pos0'],pos1=e['pos0']+1,
      ref=ref,alt=alt,RBP=e['RBP'],cell=e['cell'],task=e['task'],split=e['split'],cluster_id=e['cluster_id'],strand=e['strand'],
      gene_names=';'.join(e['gene_names']),gene_ids=';'.join(e['gene_ids']),ref_preferred=float(ref==maj),
      phi_center=float(c['phi_raw'][50]),phi_window_mean=float(np.mean(c['phi_raw'])),GC101=sum(b in 'GC' for b in seq)/101,
      GC7=sum(b in 'GC' for b in seq[47:54])/7,local5=seq[48:53],reference_sequence_rna_dna=seq,
      major_sequence_rna_dna=c['major_sequence_dna'],minor_sequence_rna_dna=c['minor_sequence_dna'],source_raw_counts=e['source_raw_counts'],
      p_value=e['p_value'],fdr=e['fdr'],source_excel_row=e['source_excel_row'],source_sheet=e['source_sheet'])
    nums=re.findall(r'[ACGT]:(\d+)',e['source_raw_counts']);d['raw_total_depth']=sum(int(z) for z in nums)
    for a in ['A_native','S','A_flat','A_perm0','A_perm1','A_perm2']:
        val=Q['101__'+a];v=float(val[e['major_input_row']]-val[e['minor_input_row']])
        d[a+'_major_minus_minor']=v;d[a+'_ref_minus_alt']=v*(1 if ref==maj else -1)
        d[a+'_credit']=.5 if abs(v)<=1e-7 else float(v>0)
        d[a+'_abs_delta']=abs(v)
    v=np.mean([d[a+'_major_minus_minor'] for a in ['A_perm0','A_perm1','A_perm2']]);d['A_permmean_major_minus_minor']=v;d['A_permmean_ref_minus_alt']=v*(1 if ref==maj else -1);d['A_permmean_abs_delta']=abs(v);d['A_permmean_credit']=.5 if abs(v)<=1e-7 else float(v>0)
    rows.append(d)
d=pd.DataFrame(rows)
assert len(d)==7262 and d.event_id.is_unique
source_file=P/'ASB_ALL_EVENT_RESULTS.csv'
if not source_file.exists():source_file=P.parent/'07_results/ASB_ALL_EVENT_RESULTS.csv'
source=pd.read_csv(source_file)
for a in ['A_native','S','A_flat']:
    old=source.set_index('event_id')[a+'_delta_major_minus_minor'].reindex(d.event_id)
    assert np.allclose(old.to_numpy(),d[a+'_major_minus_minor'],rtol=0,atol=1e-14)
d.to_csv(R/'EVENTS_WITH_CONSERVATION_AND_PREDICTIONS.csv',index=False)
arms=['A_native','S','A_flat','A_permmean','ref']
def credits(f,a):return f.ref_preferred.to_numpy() if a=='ref' else f[a+'_credit'].to_numpy()
def macro(f,a):
    z=pd.DataFrame({'RBP':f.RBP.to_numpy(),'task':f.task.to_numpy(),'credit':credits(f,a)})
    return z.groupby(['RBP','task']).credit.mean().groupby(level='RBP').mean()
def bootstrap_macro(f,a,b=None):
    x=macro(f,a)
    if b is not None:x=x-macro(f,b)
    rng=np.random.default_rng(20261012);v=x.to_numpy();bs=v[rng.integers(len(v),size=(5000,len(v)))].mean(1)
    return float(v.mean()),*np.quantile(bs,[.025,.975]).tolist(),len(v)
summaries=[];quartiles=[];contrasts=[]
for label in ['all','strict_test','validation','training']:
    f=d if label=='all' else d[d.split==label]
    for a in arms:
        mean,lo,hi,n=bootstrap_macro(f,a);summaries.append(dict(cohort=label,arm=a,accuracy=mean,ci_low=lo,ci_high=hi,events=len(f),loci=f[['chrom','pos0']].drop_duplicates().shape[0],RBPs=n))
    for a,b in [('A_native','S'),('A_native','A_flat'),('A_native','A_permmean'),('A_native','ref')]:
        mean,lo,hi,n=bootstrap_macro(f,a,b);contrasts.append(dict(cohort=label,contrast=a+' minus '+b,difference=mean,ci_low=lo,ci_high=hi,events=len(f),RBPs=n))
    if label not in ('all','strict_test'):continue
    loci=f.drop_duplicates(['chrom','pos0']);cuts=np.quantile(loci.phi_center,[.25,.5,.75]);bins=np.searchsorted(cuts,f.phi_center,side='left')+1
    f=f.copy();f['phi_quartile']=bins
    for i in range(1,5):
        s=f[f.phi_quartile==i]
        for a in arms:
            mean,lo,hi,n=bootstrap_macro(s,a);quartiles.append(dict(cohort=label,quartile=i,phi_lower=float(s.phi_center.min()),phi_upper=float(s.phi_center.max()),arm=a,accuracy=mean,ci_low=lo,ci_high=hi,events=len(s),loci=s[['chrom','pos0']].drop_duplicates().shape[0],RBPs=n))
    f[['event_id','phi_quartile']].to_csv(R/f'QUARTILE_ASSIGNMENTS_{label}.csv',index=False)
pd.DataFrame(summaries).to_csv(R/'DIRECTION_SUMMARY.csv',index=False);pd.DataFrame(contrasts).to_csv(R/'DIRECTION_CONTRASTS.csv',index=False);pd.DataFrame(quartiles).to_csv(R/'DIRECTION_BY_CONSERVATION.csv',index=False)

def controls(f,ref_control):
    # Actual local bases are a limited sequence control, NOT an independent RBP motif score.
    cols=[np.ones(len(f)),f.GC101.to_numpy(),f.GC7.to_numpy(),np.log1p(f.raw_total_depth.to_numpy())]
    if ref_control:cols.append(f.ref_preferred.to_numpy())
    for j in range(5):
        for base in 'ACG':cols.append((f.local5.str[j]==base).to_numpy().astype(float))
    cat=pd.get_dummies(f[['task','strand']].astype(str),drop_first=True,dtype=float)
    mutation=pd.get_dummies(f.ref+'>'+f.alt,drop_first=True,dtype=float)
    return np.column_stack(cols+[cat.to_numpy(),mutation.to_numpy()])

associations=[]
for label in ['all','strict_test']:
    f=(d if label=='all' else d[d.split==label]).reset_index(drop=True)
    control=controls(f,True);pinv=np.linalg.pinv(control)
    x=rankdata(f.phi_center.to_numpy());x=(x-x.mean())/x.std();xr=x-control@(pinv@x)
    gene=pd.factorize(f.cluster_id)[0];rbp=pd.factorize(f.RBP)[0];rng=np.random.default_rng(20261013)
    weights=[rng.exponential(size=(2000,gene.max()+1)),rng.exponential(size=(2000,rbp.max()+1))]
    for a in ['A_native','S','A_flat','A_permmean']:
        for outcome in ['abs_delta','credit']:
            y=f[a+'_'+outcome].to_numpy()
            if outcome=='abs_delta':y=rankdata(y);y=(y-y.mean())/y.std()
            yr=y-control@(pinv@y);beta=float(xr@yr/(xr@xr));res=yr-beta*xr
            influence=xr*res/(xr@xr)
            # Positive-exponential TWO-WAY multiplier score approximation.
            values=[]
            for start in range(0,2000,100):
                w=weights[0][start:start+100,gene]*weights[1][start:start+100,rbp]
                values.extend((beta+(w-1)@influence).tolist())
            lo,hi=np.quantile(values,[.025,.975]);rho=float(spearmanr(f.phi_center,f[a+'_'+outcome]).statistic)
            associations.append(dict(cohort=label,arm=a,outcome=outcome,unadjusted_spearman=rho,adjusted_standardized_rank_slope=beta,ci_low=float(lo),ci_high=float(hi),events=len(f),gene_clusters=gene.max()+1,RBPs=rbp.max()+1,interval='two-way exponential multiplier score approximation; exploratory'))
pd.DataFrame(associations).to_csv(R/'CONSERVATION_ASSOCIATIONS.csv',index=False)
manifest=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),events=len(d),strict_events=int(sum(d.split=='strict_test')),new_inference=False,
  old_delta_exactly_recomputed=True,old_input_source=str(P),outcome='ASB direction only; predicted magnitudes are computational outcomes',
  protocol_deviation='For efficient CPU uncertainty, adjusted slope intervals use first-order two-way multiplier score approximation rather than 2000 full weighted nuisance-model refits. No cohort, outcome or comparator changed.',
  files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in R.glob('*.csv')})
(R/'ANALYSIS_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n','utf8')
print(pd.DataFrame(summaries).query("cohort == 'strict_test'").to_string(index=False))
print(pd.DataFrame(associations).query("cohort == 'strict_test'").to_string(index=False))
