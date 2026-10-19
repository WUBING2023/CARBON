"""Recover actual source coordinates and conservative gene/overlap clusters."""
import sys,re,json,gzip,hashlib
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/02_external_preferences')
MATCHED=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched')
TASKS=['HNRNPK_HepG2_ENCSR828ZID','HNRNPK_K562_ENCSR268ETU']
rows=[];audit=[]
for task in TASKS:
    cache=np.load(ROOT/'outputs'/task/'seed11'/'aligned'/'attributions.npz')
    keep={(int(c),int(r)):i for i,(c,r) in enumerate(zip(cache['chromosome'],cache['source_row']))}
    for chrom in [8,15]:
        p=Path('/data/CARBON_PROJECT/数据/h5data')/task/f'{task}_chr{chrom}.h5'
        df=pd.read_hdf(p,key='df');df.columns=['chr','chromStart','chromEnd','peak','score','strand','signal','signal_ctl','seq']
        filtered=df[~df['seq'].str.contains('N',case=False)].dropna().reset_index(drop=True)
        # Reproduce the original preprocess_5ch.py row filtering exactly.
        checked=0
        for ri,entry in filtered.iterrows():
            if (chrom,int(ri)) not in keep:continue
            ii=keep[(chrom,int(ri))];sequence=str(entry['seq']).upper().replace('U','T');assert sequence==cache['sequence'][ii]
            assert int(entry.chromEnd)-int(entry.chromStart)==101
            rows.append(dict(task=task,cell=task.split('_')[1],selected_index=ii,chromosome=chrom,source_row=int(ri),start=int(entry.chromStart),end=int(entry.chromEnd),strand=entry.strand,sequence_checked=True));checked+=1
        audit.append(dict(task=task,chromosome=chrom,original_rows=len(df),filtered_rows=len(filtered),mapped_selected=checked,source_path=str(p),source_sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
assert len(rows)==2048
genes={8:[],15:[]};gff=Path('/data/CARBON_PROJECT/数据处理/gencodeV29/gencode.v29.annotation.gff3')
with gff.open() as f:
    for line in f:
        if line.startswith('#'):continue
        parts=line.rstrip().split('\t')
        if len(parts)!=9 or parts[2]!='gene' or parts[0] not in ['chr8','chr15']:continue
        attr=dict(v.split('=',1) for v in parts[8].split(';') if '=' in v)
        ident=attr.get('gene_id',attr.get('ID',''));ident=ident[5:] if ident.startswith('gene:') else ident
        genes[int(parts[0][3:])].append((int(parts[3])-1,int(parts[4]),ident,attr.get('gene_name',''),parts[6]))
parent=list(range(len(rows)))
def find(i):
    while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
    return i
def union(i,j):
    i=find(i);j=find(j)
    if i!=j:parent[max(i,j)]=min(i,j)
gene_to_row={}
for i,r in enumerate(rows):
    hits=[g for g in genes[r['chromosome']] if g[0]<r['end'] and g[1]>r['start']]
    r['overlapping_gene_ids']=';'.join(g[2] for g in hits);r['overlapping_gene_names']=';'.join(g[3] for g in hits)
    # Ignore strand conservatively: overlapping antisense loci are not counted apart.
    for g in hits:
        if g[2] in gene_to_row:union(i,gene_to_row[g[2]])
        else:gene_to_row[g[2]]=i
for chrom in [8,15]:
    ids=sorted([i for i,r in enumerate(rows) if r['chromosome']==chrom],key=lambda i:rows[i]['start']);active=[]
    for i in ids:
        active=[j for j in active if rows[j]['end']>rows[i]['start']]
        for j in active:union(i,j)
        active.append(i)
groups={}
for i,r in enumerate(rows):groups.setdefault(find(i),[]).append(i)
for gid,(root,members) in enumerate(sorted(groups.items()),1):
    for i in members:rows[i]['cluster']=f'gene_locus_{gid:04d}';rows[i]['cluster_windows']=len(members)
out=ROOT/'locus_mapping';out.mkdir(exist_ok=True)
pd.DataFrame(rows).sort_values(['task','selected_index']).to_csv(out/'HNRNPK_NATIVE_COORDINATES.csv',index=False)
pd.DataFrame(audit).to_csv(out/'SOURCE_COORDINATE_AUDIT.csv',index=False)
meta=dict(n_selected=2048,n_clusters=len(groups),n_gene_mapped=sum(bool(r['overlapping_gene_ids']) for r in rows),largest_cluster=max(map(len,groups.values())),cluster_rule='union selected intervals sharing any overlapping Gencode v29 gene body, plus overlapping genomic intervals; ignore strand conservatively; shared across cell contexts',reference_gff3=str(gff),gff3_sha256=hashlib.sha256(gff.read_bytes()).hexdigest(),coordinate_zero_based_half_open=True,all_2048_sequences_checked=True,not_independent_protein_replication=True)
(out/'COORDINATE_QA.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta),flush=True)
