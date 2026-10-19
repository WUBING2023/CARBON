"""Diagnose source effect columns without selecting a favorable endpoint."""
from pathlib import Path
import csv,json,math,time
import numpy as np
import openpyxl
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'fr1'
WORKBOOK=ROOT.parents[1]/'expanded_functional_validation_20261003/02_new_public_datasets/sources/hnrnkp_springer_44318_2021_BFEMBJ2020106357_MOESM5_ESM.xlsx'
def num(v):
 try:
  x=float(v);return x if math.isfinite(x) else None
 except (TypeError,ValueError):return None
def rho(a,b):
 def ranks(v):
  _,ii,n=np.unique(v,return_inverse=True,return_counts=True);end=np.cumsum(n);start=end-n
  return ((start+end-1)/2)[ii]
 x=ranks(a);y=ranks(b);x=x-x.mean();y=y-y.mean();d=np.linalg.norm(x)*np.linalg.norm(y)
 return float(x@y/d) if d>1e-14 else None
def write_csv(name,rows):
 with (OUT/name).open('w',encoding='utf8',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));writer.writeheader();writer.writerows(rows)

start=time.perf_counter()
w=openpyxl.load_workbook(WORKBOOK,read_only=True,data_only=True)
it=iter(w['NucLibB'].values);next(it);h=next(it)
raw={r[0]:dict(zip(h,r)) for r in it if r[0]}
labels=list(csv.DictReader((OUT/'LABELS_ELIGIBLE.csv').open(encoding='utf8')))
rows=[];summary=[]
for r in labels:
 z=raw[r['construct_id'].split(':',1)[1]]
 nt=np.array([num(z[f'siNT Nuc/Cyto (log2) rep{j}']) for j in [1,2,3]],float)
 for protein in ['HNRNPK','SLTM','SNRNP70']:
  kd=np.array([num(z[f'si{protein} Nuc/Cyto (log2) rep{j}']) for j in [1,2,3]],float)
  a=num(z['Effect of si'+protein+' on Nuc/Cyto (log2) ratios'])
  nc=num(z['Nuc/Cyto (log2) average'])
  row={'construct_id':r['construct_id'],'gene':r['gene'],'protein':protein,'author_effect':a,
       'condition_mean_KD_minus_NT':float(np.nanmean(kd)-np.nanmean(nt)),
       'condition_median_KD_minus_NT':float(np.nanmedian(kd)-np.nanmedian(nt)),
       'paired_replicate_median_KD_minus_NT':float(np.nanmedian(kd-nt)),
       'all_nine_cross_replicate_median_KD_minus_NT':float(np.nanmedian((kd[:,None]-nt[None,:]).ravel())),
       'KD_mean_minus_baseline2018':float(np.nanmean(kd)-nc) if nc is not None else None,
       'log2_mean_ratio_KD_minus_NT':float(np.log2(np.nanmean(2**kd))-np.log2(np.nanmean(2**nt))),
       'all_three_reps_finite':bool(np.isfinite(nt).all() and np.isfinite(kd).all())}
  rows.append(row)
for gene in ['ALL']+sorted(set(r['gene'] for r in rows)):
 for protein in ['HNRNPK','SLTM','SNRNP70']:
  rs=[r for r in rows if r['protein']==protein and (gene=='ALL' or r['gene']==gene)]
  for formula in ['condition_mean_KD_minus_NT','condition_median_KD_minus_NT','paired_replicate_median_KD_minus_NT','all_nine_cross_replicate_median_KD_minus_NT','KD_mean_minus_baseline2018','log2_mean_ratio_KD_minus_NT']:
   for sign in [1,-1]:
    ok=[r for r in rs if r['author_effect'] is not None and r[formula] is not None]
    aa=np.array([r['author_effect'] for r in ok]);bb=np.array([r[formula]*sign for r in ok]);diff=aa-bb
    summary.append({'gene':gene,'protein':protein,'formula':formula,'sign':sign,'n':len(ok),'rho':rho(aa,bb),
         'mean_diff':float(diff.mean()),'SD_diff':float(diff.std()),'RMSE':float(np.sqrt(np.mean(diff**2))),
         'max_abs_diff':float(np.abs(diff).max()),'n_equal_1e_8':int((np.abs(diff)<1e-8).sum()),
         'direction_agreement':float(np.mean(np.sign(aa)==np.sign(bb)))})
write_csv('EFFECT_AXIS_PER_TILE.csv',rows);write_csv('EFFECT_AXIS_ALL_RNA_FORMULAS.csv',summary)
(OUT/'EFFECT_AXIS_STATUS.json').write_text(json.dumps({'tested_formulas_all_genes':True,'any_exact_formula':any(r['gene']=='ALL' and r['n']==r['n_equal_1e_8'] for r in summary),
    'unknown_author_model_or_normalization_not_invented':True,'wall_seconds':time.perf_counter()-start},indent=2),encoding='utf8')
for r in summary:
 if r['gene']=='ALL' and r['sign']==1:print(json.dumps(r))
