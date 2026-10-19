"""Independent family end-point replay from author matrix and saved raw responses."""
import sys
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import numpy as np,pandas as pd,json,itertools,hashlib,csv,datetime
from scipy.stats import spearmanr
R=Path(__file__).resolve().parent;P=R.parent/'02_external_preferences';O=R/'cross_audit_02';checks=[]
def ck(name,err,tol=1e-12):checks.append({'check':name,'max_error':float(err),'pass':bool(err<=tol)})
km=[''.join(k) for k in itertools.product('ACGT',repeat=5)];ki={k:i for i,k in enumerate(km)}
meta=pd.read_csv(P/'results/WIDE_EXPERIMENT_METADATA.csv');labs=sorted(meta.experiment.unique())
with (P/'sources/Zscores_420_origrncmpt.txt').open() as f:
 header=f.readline().split();ix=[header.index(l) for l in labs];sums=np.zeros((len(labs),1024));counts=np.zeros(1024,int)
 for line in f:
  v=line.split();word=v[0].replace('U','T');z=np.array([float(v[i]) for i in ix])
  for off in range(3):k=ki[word[off:off+5]];sums[:,k]+=z;counts[k]+=1
truth={lab:sums[i]/counts for i,lab in enumerate(labs)};savedtruth=pd.read_csv(P/'results/WIDE_EXTERNAL_5MER_SCORES.csv')
for lab,g in savedtruth.groupby('experiment'):ck(lab+'/author_truth',np.max(abs(g.set_index('kmer').loc[km].z5.to_numpy()-truth[lab])))
groups={}
for i,w in enumerate(km):groups.setdefault(tuple(w.count(c) for c in 'ACGT'),[]).append(i)
pairs=np.array([z for g in groups.values() for z in itertools.combinations(g,2)]);a,b=pairs.T
sourcepred=pd.read_csv(P/'results/WIDE_ALL_5MER_PREDICTIONS.csv');saved=pd.read_csv(P/'results/WIDE_ALL_METRICS.csv');rows=[]
for (task,arm),g in sourcepred.groupby(['task','arm'],sort=True):
 protein=g.protein.iloc[0];cell=g.cell.iloc[0];pred=g.set_index('kmer').loc[km].score.to_numpy()
 if arm in ['full','sequence','frozen_flat','frozen_perm']:
  raw=np.load(P/'results/wide_extracted/wide_outputs'/task/arm/'responses.npz');ck(task+'/'+arm+'/raw_context_mean',np.max(abs(raw['delta_log_count'].mean(0)-pred)))
 for lab in meta.loc[meta.protein==protein,'experiment']:
  y=truth[lab];d=y[a]-y[b];pdiff=pred[a]-pred[b];use=d!=0;con=np.where(pdiff[use]==0,.5,np.sign(pdiff[use])==np.sign(d[use])).mean();rho=spearmanr(pred,y).statistic if np.ptp(pred)>0 else np.nan
  r=saved[(saved.task==task)&(saved.arm==arm)&(saved.experiment==lab)].iloc[0]
  ck(task+'/'+arm+'/'+lab+'/composition',abs(con-r.within_composition_concordance));ck(task+'/'+arm+'/'+lab+'/pair_denominator',abs(int(use.sum())-r.n_pairs))
  if np.isfinite(rho):ck(task+'/'+arm+'/'+lab+'/rho',abs(rho-r.spearman_all5mer))
  rows.append({'task':task,'protein':protein,'cell':cell,'arm':arm,'experiment':lab,'exposed':protein in ['HNRNPK','PTBP1','RBFOX2'],'rho':rho,'concordance':con})
df=pd.DataFrame(rows);bg=df.groupby(['protein','cell','arm','exposed'])[['rho','concordance']].mean().reset_index();prot=bg.groupby(['protein','arm','exposed'])[['rho','concordance']].mean().reset_index();savedprot=pd.read_csv(P/'results/WIDE_PROTEIN_METRICS.csv')
for _,r in prot.iterrows():
 s=savedprot[(savedprot.protein==r.protein)&(savedprot.arm==r.arm)].iloc[0]
 for x,y in [('rho','spearman_all5mer'),('concordance','within_composition_concordance')]:
  if np.isfinite(r[x]):ck(r.protein+'/'+r.arm+'/protein_'+x,abs(r[x]-s[y]))
boots=pd.read_csv(P/'results/WIDE_PAIRED_PROTEIN_BOOTSTRAP.csv');rng=np.random.default_rng(20261004);summary=[]
for name,sub in [('all28',prot),('new25',prot[~prot.exposed]),('exposed3',prot[prot.exposed])]:
 proteins=sorted(sub.protein.unique());indices=rng.integers(0,len(proteins),(10000,len(proteins)))
 for metric,endpoint in [('rho','spearman_all5mer'),('concordance','within_composition_concordance')]:
  t=sub.pivot(index='protein',columns='arm',values=metric).loc[proteins]
  for control in ['sequence','frozen_flat','frozen_perm','frequency','phi_only','A_count','C_count','G_count','T_count']:
   if not np.isfinite(t[control]).all():continue
   d=(t.full-t[control]).to_numpy();values=d[indices].mean(1);r=boots[(boots.subset==name)&(boots.endpoint==endpoint)&(boots.control==control)].iloc[0]
   for label,v,orig in [('point',d.mean(),r.mean_delta),('low',np.quantile(values,.025),r.protein_paired_bootstrap_low),('high',np.quantile(values,.975),r.protein_paired_bootstrap_high)]:ck(name+'/'+metric+'/'+control+'/'+label,abs(v-orig))
   if name=='new25':summary.append({'endpoint':endpoint,'control':control,'n_proteins':len(proteins),'point':d.mean(),'low':np.quantile(values,.025),'high':np.quantile(values,.975),'positive_proteins':int((d>0).sum())})
pd.DataFrame(rows).to_csv(O/'FAMILY_INDEPENDENT_METRICS.csv',index=False);pd.DataFrame(summary).to_csv(O/'FAMILY_INDEPENDENT_CONTRASTS.csv',index=False);pd.DataFrame(checks).to_csv(O/'FAMILY_CHECKS.csv',index=False)
qa={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':all(c['pass'] for c in checks),'checks':len(checks),'max_error':max(c['max_error'] for c in checks),'all28':28,'new25':25,'author_matrix_7mers':16382,'same_composition_candidate_pairs':len(pairs),'bootstrap_n':10000,'bootstrap_unit':'protein; cell backgrounds/assay columns first averaged within protein','scope':'Original RNAcompete truth and raw response means, all protein metrics and paired bootstrap readback; no new model inference','source_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'sources/Zscores_420_origrncmpt.txt',P/'results/WIDE_ALL_5MER_PREDICTIONS.csv',P/'results/WIDE_PAIRED_PROTEIN_BOOTSTRAP.csv',Path(__file__)]}}
(O/'FAMILY_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf8');print(json.dumps({k:v for k,v in qa.items() if k!='source_hashes'}));assert qa['pass']
