"""Strict source-based hg19->hg38 input census; no model scores loaded."""
from pathlib import Path
import csv,gzip,json,hashlib,datetime,collections
import numpy as np
R=Path(__file__).resolve().parent;R.mkdir(exist_ok=True)
S=R.parent/'ptbp_splicing_feasibility'; CHROMS={'chr1','chr8','chr15','chrX','chrY'}
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
bins=collections.defaultdict(list);h=None;t=q=0
with gzip.open(S/'hg19ToHg38.over.chain.gz','rt') as f:
 for line in f:
  z=line.split()
  if not z: continue
  if z[0]=='chain':
   h=z;t=int(z[5]);q=int(z[10]);continue
  size=int(z[0])
  if h[2] in CHROMS:
   assert h[4]=='+'
   b=(t,t+size,h[7],int(h[8]),h[9],q,h[12])
   for key in range(t//100000,(t+size-1)//100000+1):bins[(h[2],key)].append(b)
  if len(z)==3:t+=size+int(z[1]);q+=size+int(z[2])
def point(ch,pos):
 hits=set()
 for a,b,d,n,st,qp,cid in bins[(ch,pos//100000)]:
  if a<=pos<b:
   x=qp+pos-a
   if st=='-':x=n-x-1
   hits.add((d,x,st))
 return next(iter(hits)) if len(hits)==1 else None
def interval(ch,a,b):
 v=[point(ch,p) for p in range(a,b)]
 if any(x is None for x in v):return None
 if len({(x[0],x[2]) for x in v})!=1:return None
 step=1 if v[0][2]=='+' else -1
 if any(v[i+1][1]-v[i][1]!=step for i in range(len(v)-1)):return None
 return (v[0][0],min(x[1] for x in v),max(x[1] for x in v)+1,v[0][2])
allrows=[];windows=[];seen=set();sources=[];signerr=0
for cell,acc in [('K562','ENCFF061DFK'),('HepG2','ENCFF044JDS')]:
 p=S/(acc+'_SE.MATS.JunctionCountOnly.txt');sources.append({'path':str(p),'sha256':sha(p)})
 for rowno,d in enumerate(csv.DictReader(p.open(encoding='utf8'),delimiter='\t')):
  if d['chr'] not in CHROMS:continue
  gene=d['GeneID'].strip('"').split('.')[0]
  rec={'cell':cell,'source_accession':acc,'source_row':rowno,'source_ID':d['ID'],'gene':gene,'geneSymbol':d['geneSymbol'].strip('"'),'hg19_chrom':d['chr'],'hg19_strand':d['strand'],'scope':'chr1_primary' if d['chr']=='chr1' else 'other_holdouts_secondary'}
  try:
   kd=np.array([float(x) for x in d['IncLevel1'].split(',')]);ct=np.array([float(x) for x in d['IncLevel2'].split(',')]);assert len(kd)==len(ct)==2
   cov=np.array([[float(x) for x in d[k].split(',')] for k in ['IJC_SAMPLE_1','SJC_SAMPLE_1','IJC_SAMPLE_2','SJC_SAMPLE_2']]); assert cov.shape==(4,2)
   delta=float(d['IncLevelDifference']);fdr=float(d['FDR']);finite=np.isfinite(kd).all() and np.isfinite(ct).all() and np.isfinite(delta) and np.isfinite(fdr)
   err=abs(delta-(kd.mean()-ct.mean())) if finite else np.nan
   if finite:signerr=max(signerr,err)
   rec.update(delta_PSI=delta,FDR=fdr,KD_PSI=float(kd.mean()),control_PSI=float(ct.mean()),coverage_min=float(np.min([cov[0]+cov[1],cov[2]+cov[3]])),coverage_mean=float(cov.sum()/4),rounded_delta_error=err)
  except (ValueError,AssertionError):finite=False
  if not finite:rec['exclusion']='nonfinite_or_invalid_source';allrows.append(rec);continue
  # IncLevelDifference is rounded to 3 decimals in the public rMATS table.
  assert err<=.001001,(acc,rowno,err)
  rec['response_class']='responsive' if fdr<=.05 and abs(delta)>=.05 else ('nonresponsive' if fdr>=.5 and abs(delta)<=.02 else 'intermediate')
  if rec['coverage_min']<10:rec['exclusion']='coverage_below10';allrows.append(rec);continue
  es,ee,ues,uee,des,dee=[int(d[k]) for k in ['exonStart_0base','exonEnd','upstreamES','upstreamEE','downstreamES','downstreamEE']]
  rec.update(hg19_exon_start=es,hg19_exon_end=ee,hg19_upstreamES=ues,hg19_upstreamEE=uee,hg19_downstreamES=des,hg19_downstreamEE=dee)
  if es-uee<101 or des-ee<101:rec['exclusion']='intron_below101';allrows.append(rec);continue
  ends=[point(d['chr'],v) for v in [ues,uee-1,es,ee-1,des,dee-1]]
  if any(v is None for v in ends) or len({(v[0],v[2]) for v in ends})!=1 or ends[0][0]!=d['chr']:
   rec['exclusion']='six_boundary_nonunique_or_wrong_chr';allrows.append(rec);continue
  flip=ends[0][2];newstrand=d['strand'] if flip=='+' else ('-' if d['strand']=='+' else '+')
  newes=min(ends[2][1],ends[3][1]);newee=max(ends[2][1],ends[3][1])+1
  rngs=[(es-101,es),(ee,ee+101)] if d['strand']=='+' else [(ee,ee+101),(es-101,es)]
  vals=[interval(d['chr'],a,b) for a,b in rngs]
  if any(v is None or v[0]!=d['chr'] or v[2]-v[1]!=101 or v[3]!=flip for v in vals):
   rec['exclusion']='flank101_nonunique_indel_or_wrong_chr';allrows.append(rec);continue
  if newstrand=='+':expected=[(newes-101,newes),(newee,newee+101)]
  else:expected=[(newee,newee+101),(newes-101,newes)]
  if [(v[1],v[2]) for v in vals]!=expected:rec['exclusion']='mapped_junction_gap';allrows.append(rec);continue
  newbounds=tuple(v[1] for v in ends);dedup=(cell,gene,d['chr'],newstrand,newbounds)
  if dedup in seen:rec['exclusion']='duplicate_gene_three_exons';allrows.append(rec);continue
  seen.add(dedup);ei=len(allrows);rec.update(exclusion='',event_index=ei,hg38_chrom=d['chr'],hg38_strand=newstrand,hg38_exon_start=newes,hg38_exon_end=newee,chain_flip=flip)
  for side,v in zip(['upstream','downstream'],vals):
   wi=len(windows);rec[side+'_window_index']=wi
   windows.append({'window_index':wi,'event_index':ei,'cell':cell,'chrom':v[0],'start':v[1],'end':v[2],'strand':newstrand,'side':side,'source_row':rowno,'source_accession':acc,'gene':gene,'identity':acc+'|'+str(rowno)+'|'+side})
  allrows.append(rec)
def csvout(name,rows):
 fields=list(dict.fromkeys(k for r in rows for k in r));f=(R/name).open('w',newline='',encoding='utf8');w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows);f.close()
csvout('QUALIFICATION_ALL.csv',allrows);csvout('WINDOWS.csv',windows)
den=[]
for cell in ['K562','HepG2']:
 for scope in ['chr1_primary','other_holdouts_secondary']:
  z=[r for r in allrows if r['cell']==cell and r['scope']==scope];ok=[r for r in z if r['exclusion']=='']
  den.append({'cell':cell,'scope':scope,'source_rows':len(z),'qualified_events':len(ok),'genes':len(set(r['gene'] for r in ok)),'responsive':sum(r['response_class']=='responsive' for r in ok),'KD_inclusion':sum(r['response_class']=='responsive' and r['delta_PSI']>0 for r in ok),'KD_exclusion':sum(r['response_class']=='responsive' and r['delta_PSI']<0 for r in ok),'nonresponsive':sum(r['response_class']=='nonresponsive' for r in ok),'intermediate':sum(r['response_class']=='intermediate' for r in ok),'exclusions':json.dumps(dict(collections.Counter(r['exclusion'] for r in z if r['exclusion'])))})
csvout('DENOMINATORS_MAP.csv',den)
freeze={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_files':sources,'chain_sha256':sha(S/'hg19ToHg38.over.chain.gz'),'script_sha256':sha(__file__),'protocol_sha256':sha(R.parent/'PROTOCOL_FR4.md'),'addendum_sha256':sha(R.parent/'PROTOCOL_FR4_HOLDOUT_ADDENDUM.md'),'qualification_sha256':sha(R/'QUALIFICATION_ALL.csv'),'windows_sha256':sha(R/'WINDOWS.csv'),'rounded_source_delta_max_error':signerr,'n_windows':len(windows),'model_scores_loaded':False,'strict_all_101_mapping':True,'six_exon_halfopen_boundaries':True}
(R/'MAPPING_FREEZE.json').write_text(json.dumps(freeze,indent=2),encoding='utf8')
print(json.dumps({'n_windows':len(windows),'denominators':den,'source_delta_error_max':signerr}))
