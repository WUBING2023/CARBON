"""server15 read-only public-sequence mapping. No weights or model imports."""
from pathlib import Path
import sys
import json
import collections
import hashlib
import numpy as np
import pysam
import pyBigWig

payload=json.loads(sys.stdin.read())
FA=Path('/data/CARBON_PROJECT/数据处理/gencodeV29/GRCh38.p12.genome.fa')
GFF=Path('/data/CARBON_PROJECT/数据处理/gencodeV29/gencode.v29.annotation.gff3')
PHI=Path('/data/CARBON_PROJECT/数据/phyloP/hg38.phyloP100way.bw')
length=payload['sequence_length'];assert length==109
gene_names=set(payload['genes'])
exons=collections.defaultdict(lambda:collections.defaultdict(set))
genes=collections.defaultdict(list)
for line in GFF.open():
    if line.startswith('#'):continue
    fields=line.rstrip().split('\t')
    if len(fields)!=9 or fields[2] not in ['gene','exon']:continue
    attrs=dict(x.split('=',1) for x in fields[8].split(';') if '=' in x)
    name=attrs.get('gene_name',attrs.get('Name',''))
    if name not in gene_names:continue
    if fields[2]=='gene':genes[name].append({'chrom':fields[0],'start0':int(fields[3])-1,'end':int(fields[4]),'strand':fields[6],'gene_id':attrs.get('gene_id',attrs.get('ID',''))})
    else:
        tx=attrs.get('transcript_id',attrs.get('Parent',''))
        exons[name][tx].add((fields[0],int(fields[3])-1,int(fields[4]),fields[6]))
fa=pysam.FastaFile(str(FA));bw=pyBigWig.open(str(PHI))
comp=str.maketrans('ACGTN','TGCAN')
def hash_file(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def edge_metadata(path):
    size=path.stat().st_size
    with path.open('rb') as f:
        first=f.read(1024*1024);f.seek(max(0,size-1024*1024));last=f.read(1024*1024)
    return {'path':str(path),'size':size,'mtime_ns':path.stat().st_mtime_ns,
            'edge_sha256':hashlib.sha256(first+last).hexdigest(),'checksum_scope':'first_and_last1MiB_not_full_file_SHA256'}
def phi_at(chrom,coords,strand):
    values=[];step=1 if strand=='+' else -1;i=0
    while i<len(coords):
        j=i+1
        while j<len(coords) and coords[j]==coords[j-1]+step:j+=1
        part=np.asarray(bw.values(chrom,min(coords[i:j]),max(coords[i:j])+1),dtype=np.float32)
        if strand=='-':part=part[::-1]
        values.extend(part.tolist());i=j
    assert len(values)==length
    return [float(x) if np.isfinite(x) else None for x in values]
out=[];inventory=[]
for gene in payload['genes']:
    txs=[]
    for tx,es in exons[gene].items():
        chroms={x[0] for x in es};strands={x[3] for x in es}
        assert len(chroms)==len(strands)==1
        chrom=next(iter(chroms));strand=next(iter(strands));seq='';coords=[]
        for _,a,b,_ in sorted(es,key=lambda e:e[1],reverse=strand=='-'):
            part=fa.fetch(chrom,a,b).upper();pos=list(range(a,b))
            if strand=='-':part=part.translate(comp)[::-1];pos.reverse()
            seq+=part;coords+=pos
        txs.append({'transcript_id':tx,'chrom':chrom,'strand':strand,'sequence':seq,'coords':coords})
    inventory.append({'gene_name':gene,'annotation_genes':genes[gene],'transcript_count':len(txs),
                      'transcripts':[{'transcript_id':t['transcript_id'],'chrom':t['chrom'],'strand':t['strand'],'length':len(t['sequence']),'sequence_sha256':hashlib.sha256(t['sequence'].encode()).hexdigest()} for t in txs]})
    for record in [r for r in payload['records'] if r['gene_name']==gene]:
        query=record['sequence_dna'];assert len(query)==length
        hits=collections.defaultdict(list)
        for t in txs:
            offset=0
            while True:
                offset=t['sequence'].find(query,offset)
                if offset<0:break
                coord=tuple(t['coords'][offset:offset+length]);assert len(coord)==length
                hits[(t['chrom'],t['strand'],coord)].append({'transcript_id':t['transcript_id'],'transcript_offset0':offset})
                offset+=1
        item={**record,'coordinate_mapping_count':len(hits),'unique_natural_human_RNA_mapping':len(hits)==1,
              'no_model_inference':True,'maps':[]}
        for (chrom,strand,coords),equivalent in hits.items():
            coords=list(coords)
            rebuilt=''.join(fa.fetch(chrom,i,i+1).upper() for i in coords)
            if strand=='-':rebuilt=rebuilt.translate(comp)
            assert rebuilt==query
            phi=phi_at(chrom,coords,strand)
            item['maps'].append({'chrom':chrom,'strand':strand,'coordinates0':coords,'genomic_base_reconstruction':rebuilt,
                'equivalent_isoforms':equivalent,'phi_raw':phi,'phi_valid':[x is not None for x in phi],
                'phi_missing':sum(x is None for x in phi),'genomic_runs':sum(i==0 or coords[i]!=coords[i-1]+(1 if strand=='+' else -1) for i in range(length))})
        out.append(item)
fa.close();bw.close()
result={'records':out,'reference_inventory':inventory,
        'policy':{'exact109nt_matching_all_declared_human_gene_isoforms':True,'coordinate_dedup_key':'chrom,strand,all109_coordinates',
                  'no_inferred_mouse_human_liftover_or_synthetic_adaptor_coordinates':True,
                  'orientation':'RNA5prime_to3prime','missing_Phi':'explicit_null_plus_valid_mask',
                  'source_species_distinct_from_human_reference_compatibility':True,'model_outputs_queried':False},
        'reference':{'fasta':edge_metadata(FA),'fasta_fai':{'path':str(FA)+'.fai','sha256':hash_file(Path(str(FA)+'.fai'))},
                     'annotation':{'path':str(GFF),'size':GFF.stat().st_size,'sha256':hash_file(GFF)},'phylop':edge_metadata(PHI)}}
print(json.dumps(result,ensure_ascii=False,allow_nan=False))
