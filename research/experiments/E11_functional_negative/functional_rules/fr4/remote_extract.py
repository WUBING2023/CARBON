from pathlib import Path
import csv,json,time
import numpy as np,pysam,pyBigWig
R=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/03_functional_rules/fr4')
FA=Path('/data/CARBON_PROJECT/数据处理/gencodeV29/GRCh38.p12.genome.fa');BW=Path('/data/CARBON_PROJECT/数据/phyloP/hg38.phyloP100way.bw')
fa=pysam.FastaFile(str(FA));bw=pyBigWig.open(str(BW));rows=list(csv.DictReader((R/'WINDOWS.csv').open()));xx=[];meta=[]
for r in rows:
 a,b=int(r['start']),int(r['end']);s=fa.fetch(r['chrom'],a,b).upper();v=np.array(bw.values(r['chrom'],a,b),dtype=np.float32)
 status='' if len(s)==101 and set(s)<=set('ACGT') and np.isfinite(v).all() else 'noncanonical_sequence_or_incomplete_native_Phi'
 if r['strand']=='-':s=s.translate(str.maketrans('ACGT','TGCA'))[::-1];v=v[::-1].copy()
 x=np.zeros((101,5),np.float32)
 if not status:
  for k,c in enumerate('ACGT'):x[:,k]=np.fromiter((a==c for a in s),dtype=np.float32,count=101)
  x[:,4]=v
 xx.append(x);meta.append({**r,'native_exclusion':status,'sequence_RNA5to3':s,'CT_fraction':sum(c in 'CT' for c in s)/101,'C_fraction':s.count('C')/101,'CTCT_count':sum(s[i:i+4]=='CTCT' for i in range(98)),'mean_Phi':float(v.mean()) if np.isfinite(v).all() else None})
bad_events={r['event_index'] for r in meta if r['native_exclusion']}
for r in meta:
 if r['event_index'] in bad_events:r['native_exclusion']='pair_has_noncanonical_sequence_or_incomplete_native_Phi'
x=np.stack(xx);valid=np.array([not r['native_exclusion'] for r in meta]);np.savez_compressed(R/'INPUTS.npz',X=x,valid=valid)
(R/'NATIVE_METADATA.json').write_text(json.dumps(meta))
