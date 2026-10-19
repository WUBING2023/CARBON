from pathlib import Path
import csv,json,ast,hashlib
R=Path(__file__).resolve().parent;BASE=R.parent.parent
tables=json.loads((R/'HISTORICAL_EXPOSURE.json').read_text(encoding='utf-8'))['scanned_coordinate_tables']
records=[];names=set();sources=[]
def parse(value):
 try:return json.loads(value)
 except Exception:
  try:return ast.literal_eval(value)
  except Exception:return None
def number(v):
 try:return int(float(v))
 except Exception:return None
for filename in tables:
 p=Path(filename);sources.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 with open(p,encoding='utf-8-sig',newline='') as f:
  for row in csv.DictReader(f):
   chrom=next((row[k] for k in ['chrom','chromosome','chr'] if row.get(k)),None)
   if chrom is None:continue
   chrom='chr'+str(chrom).removeprefix('chr')
   genes=[]
   for k in ['gene','gene_name','source_gene_name','reference_gene_name','source_gene','gene_names']:
    if row.get(k):
     a=parse(row[k]);genes.extend(a if isinstance(a,list) else [row[k]])
   genes=[str(g) for g in genes if g and str(g) not in ['nan','None','[]']];names.update(genes)
   start=end=None;mode=''
   for a,b in [('start','end'),('genomic_101_start0_inferred','genomic_101_end0_inferred'),('context_start0_inferred','context_end0_inferred')]:
    if row.get(a) and row.get(b):start,end=number(row[a]),number(row[b]);mode=a;break
   if start is None:
    coords=parse(row.get('reference_coordinates109_0based_json',''))
    if isinstance(coords,list) and coords:
     vals=[number(v) for v in coords];vals=[v for v in vals if v is not None]
     if vals:start,end=min(vals),max(vals)+1;mode='reference_coordinate_span'
   if start is None:
    runs=parse(row.get('genomic_runs',''))
    if isinstance(runs,list) and runs:
     vals=[]
     for q in runs:
      if isinstance(q,(list,tuple)):vals.extend(number(v) for v in q[:2])
      elif isinstance(q,dict):vals.extend(number(q.get(v)) for v in ['start','end','start0','end0'] if v in q)
     vals=[v for v in vals if v is not None]
     if vals:start,end=min(vals),max(vals)+1;mode='genomic_run_span'
   if start is None:
    for k,offset in [('pos0',0),('hg38_pos1',1)]:
     pos=number(row.get(k,''))
     if pos is not None:start,end=pos-offset-102,pos-offset+103;mode=k+'_plus102';break
   if start is not None and end is not None and end>start:records.append({'source':str(p),'chromosome':chrom,'start':start,'end':end,'mode':mode,'genes':genes})
out={'scope':'Coordinate/gene exposure from the same15historical CSVs. Genomic span is conservative across splice jumps. Broader manuscript benchmarking and other historical analyses are not fully inventoried; ledger absence is not blindness. Names without coordinates also count as previously exposed RNA names.','sources':sources,'gene_names':sorted(names),'intervals':records}
(R/'HISTORY_REPLAY_ALL_CHROMOSOMES.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
out['intervals']=[v for v in records if v['chromosome'] in ['chrX','chrY']]
(R/'phase2_ambiguity/HISTORICAL_EXPOSURE.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'tables':len(sources),'gene_names':len(names),'XY_records':len(out['intervals']),'XY_distinct_spans':len(set((v['chromosome'],v['start'],v['end']) for v in out['intervals']))},indent=2))
