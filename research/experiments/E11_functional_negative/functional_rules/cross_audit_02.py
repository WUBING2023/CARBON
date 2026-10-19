"""Independent source/point/selected cluster-bootstrap replay. No worker imports."""
import sys
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import numpy as np,pandas as pd,itertools,json,hashlib,csv
from scipy.stats import spearmanr
R=Path(__file__).resolve().parent;P=R.parent/'02_external_preferences';O=R/'cross_audit_02';O.mkdir(exist_ok=True);km=[''.join(k) for k in itertools.product('ACGT',repeat=5)];ki={k:i for i,k in enumerate(km)};checks=[]
def ck(name,err,tol=1e-12):
 checks.append({'check':name,'max_error':float(err),'pass':bool(err<=tol)})
# Recreate 7->5 occurrence mapping directly from original author-associated matrix.
fp=P/'sources/Zscores_420_origrncmpt.txt';sums=np.zeros(1024);nums=np.zeros(1024,int)
with fp.open() as f:
 head=f.readline().split();col=head.index('RNCMPT00026');n7=0
 for line in f:
  v=line.split();word=v[0].replace('U','T');val=float(v[col]);assert len(word)==7 and np.isfinite(val);n7+=1
  for off in range(3):j=ki[word[off:off+5]];sums[j]+=val;nums[j]+=1
truth=sums/nums;tr=pd.read_csv(P/'results/CROSS_ASSAY_5MER_TRUTH.csv').query("protein=='HNRNPK'").set_index('kmer').loc[km];ck('original420_occurrence_truth',np.max(abs(tr.zscore5-truth)));ck('original420_occurrence_denominator',np.max(abs(tr.available_occurrences-nums)));toptruth=set(np.argsort(-truth,kind='stable')[:20]);coord=pd.read_csv(P/'results/locus_mapping/HNRNPK_NATIVE_COORDINATES.csv');gkeys=sorted(coord.cluster.unique());assert len(gkeys)==385;gidx={k:i for i,k in enumerate(gkeys)};weights=np.random.default_rng(20261004).multinomial(385,np.full(385,1/385),size=2000);saved=np.load(P/'results/HNRNPK_BOOTSTRAP_DRAWS.npz');primary=pd.read_csv(P/'results/CROSS_ASSAY_ALL_METRICS.csv');indep=[];draws={a:[] for a in ['aligned','sequence_only','mean','permuted']}
for task,part in coord.groupby('task'):
 part=part.sort_values('selected_index');assert list(part.selected_index)==list(range(1024));gi=np.array([gidx[k] for k in part.cluster]);
 for seed in [11,23,47]:
  for arm in draws:
   p=P/'results/extracted/outputs'/task/f'seed{seed}'/arm/'attributions.npz';a=np.load(p);ig=a['position_ig'];seq=a['sequence'];assert ig.shape==(1024,101) and len(seq)==1024
   win=np.array([np.convolve(v,np.ones(5),mode='valid').argmax() for v in ig]);value=np.array([ig[i,j:j+5].sum() for i,j in enumerate(win)]);ww=np.array([ki[s[j:j+5]] for s,j in zip(seq,win)]);score=np.bincount(ww,weights=value,minlength=1024);count=np.bincount(ww,minlength=1024);rank=sorted(np.flatnonzero(count),key=lambda j:(-score[j],km[j]))[:20];rec=sum(j in toptruth for j in rank)/20
   row=primary[(primary.task==task)&(primary.seed==seed)&(primary.arm==arm)].iloc[0];ck(f'{task}/{seed}/{arm}/ones_point',abs(rec-row.IG_recall20));indep.append({'task':task,'seed':seed,'arm':arm,'independent_point':rec,'original_point':row.IG_recall20,'n_eligible_5mers':int((count>0).sum())})
   # Direct per-window weights, not worker gene-matrix multiplication; preserves rank per model fit.
   rr=[]
   for wt in weights[:31]:
    sw=np.bincount(ww,weights=value*wt[gi],minlength=1024);nn=np.bincount(ww,weights=wt[gi],minlength=1024);order=sorted(np.flatnonzero(nn),key=lambda j:(-sw[j],km[j]))[:20];rr.append(sum(j in toptruth for j in order)/20)
   draws[arm].append(rr)
for arm,dd in draws.items():ck(arm+'/first31paired_cluster_draws',np.max(abs(np.array(dd).mean(0)-saved[arm][:31,0])))
point=pd.DataFrame(indep).groupby('arm').independent_point.mean();diff=pd.read_csv(P/'results/HNRNPK_PAIRED_GENE_BOOTSTRAP.csv').query("assay=='RNAcompete'");q=[]
for _,r in diff.iterrows():
 d=saved['aligned'][:,0]-saved[r.control][:,0];ck(r.control+'/saved_point',abs(point.aligned-point[r.control]-r.delta_recall20));ck(r.control+'/percentile_low',abs(np.quantile(d,.025)-r.conditional_gene_locus_bootstrap_low));ck(r.control+'/percentile_high',abs(np.quantile(d,.975)-r.conditional_gene_locus_bootstrap_high));q.append({'control':r.control,'point':point.aligned-point[r.control],'bootstrap_mean':d.mean(),'bootstrap_median':np.median(d),'low':np.quantile(d,.025),'high':np.quantile(d,.975),'point_above_high':point.aligned-point[r.control]>np.quantile(d,.975),'strict_positive_fraction':(d>0).mean()})
# Fair expected top20 C-count baseline with all boundary ties, independent labels.
cc=np.array([k.count('C') for k in km]);cut=np.sort(cc)[-20];gt=cc>cut;tie=cc==cut;membership=gt.astype(float);membership[tie]=(20-gt.sum())/tie.sum();expected=sum(membership[list(toptruth)])/20;lex=sum(i in toptruth for i in sorted(range(1024),key=lambda i:(-cc[i],km[i]))[:20])/20
construct={}
with (P/'sources/S1.In_Progress.v6.For_all.DR-2.csv').open(encoding='utf-8-sig') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if any(v=='RNCMPT00026' for v in r.values()):construct=r;break
aa=next((v for k,v in construct.items() if 'Construct_AA_seq' in k),None)
ab=[l.strip() for l in (P/'sources/Fullset_z_scores-setA-vs-setB-comparison.txt').read_text().splitlines() if 'RNCMPT00026' in l]
pd.DataFrame(indep).to_csv(O/'INDEPENDENT_POINTS.csv',index=False);pd.DataFrame(q).to_csv(O/'INDEPENDENT_BOOTSTRAP_SUMMARY.csv',index=False);pd.DataFrame(checks).to_csv(O/'CHECKS.csv',index=False)
qa={'pass':all(c['pass'] for c in checks),'checks':len(checks),'source_7mers':n7,'all_1024_5mers':int((nums>0).sum()),'occurrence_min':int(nums.min()),'occurrence_max':int(nums.max()),'construct_AA_length':len(aa) if aa else None,'construct_metadata':construct,'author_A_B_agreement_lines':ab,'author_AB_raw_matrices_independently_tested':False,'Ccount_lex_recall20':lex,'Ccount_boundary_expected_recall20':expected,'Ccount_all_words_rho':float(spearmanr(cc,truth).statistic),'bootstrap_unit':'385 native gene/locus clusters; two cell backgrounds/3 seeds not independent proteins; one external HNRNPK assay','bootstrap_point_above_upper_explained':'non-smooth top20 rule extraction is biased downward under resampling; direct31-draw independent replay and all24 all-ones points passed if pass=true; percentileCI not calibrated superiority p','source_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [fp,P/'results/HNRNPK_BOOTSTRAP_DRAWS.npz',P/'results/locus_mapping/HNRNPK_NATIVE_COORDINATES.csv',Path(__file__)]}}
(O/'QA.json').write_text(json.dumps(qa,indent=2),encoding='utf8');print(json.dumps({k:v for k,v in qa.items() if k not in ['source_hashes','construct_metadata']}));print(pd.DataFrame(q).to_string(index=False))
