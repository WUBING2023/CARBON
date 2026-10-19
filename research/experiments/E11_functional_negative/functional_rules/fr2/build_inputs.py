from pathlib import Path
import json,csv,hashlib,datetime
import numpy as np
R=Path(__file__).resolve().parent
source=json.loads((R/'ALL_SOURCE_ROWS.json').read_text('utf8'));maps=json.loads((R/'REFERENCE_MAPPING.json').read_text('utf8'));mi={r['construct_id']:r for r in maps['records']}
qual=[];rr=[];xx=[]
def f(v):
 try:v=float(v);return v if np.isfinite(v) else None
 except (TypeError,ValueError):return None
def avg(v):
 vv=[f(x) for x in v];vv=[x for x in vv if x is not None];return float(np.mean(vv)) if vv else None
for r in source:
 m=mi.get(r['construct_id']);q={k:r.get(k,'') for k in ['construct_id','source_row','author_gene','official_gene','source_accession','source_human_design','new_RNA_primary']}
 ok=bool(m and m['unique_natural_human_RNA_mapping'] and m['maps'][0]['phi_missing']==0)
 q['unique109_human_reference']=bool(m and m['unique_natural_human_RNA_mapping']);q['complete_native_Phi']=ok;q['input_index']=len(rr) if ok else ''
 q['reason']='qualified_exact109_completePhi' if ok else ('not_human_design' if not m else 'not_unique_exact_reference_or_Phi_missing')
 qual.append(q)
 if not ok:continue
 z=m['maps'][0];assert z['genomic_base_reconstruction']==r['sequence109'];v=r['values'];seq=r['sequence109'];x=np.zeros((101,5),np.float32);x[np.arange(101),[{'A':0,'C':1,'G':2,'T':3}[b] for b in seq[4:105]]]=1;x[:,4]=z['phi_raw'][4:105];assert np.isfinite(x).all();xx.append(x)
 row={**q,'reference_gene':m['gene_name'],'sequence109':seq,'sequence101':seq[4:105],'chromosome':z['chrom'],'strand':z['strand'],'coordinates109_json':json.dumps(z['coordinates0']),'unique_key':r['official_gene']+'|'+seq,'source_STable4_row':r['source_row']}
 for placement in ['3','5']:
  for endpoint in ['Nuc/Cyto','WCE/Input']:
   values=[f(v[f"{placement}' UTR Replicate {j} {endpoint} (log2)"]) for j in [1,2,3]];name=placement+('_NC' if endpoint=='Nuc/Cyto' else '_expression')
   row[name]=avg(values);row[name+'_nreps']=sum(x is not None for x in values)
   for j,val in enumerate(values,1):row[name+'_rep'+str(j)]=val
 row['C_content']=seq[4:105].count('C')/101
 for motif in ['CCC','CCCC','CCTCCC','GCCTCCC']:row[motif]=sum(seq[4:105][i:i+len(motif)]==motif for i in range(101-len(motif)+1))
 row['mean_Phi']=float(x[:,4].mean());rr.append(row)
def write(n,rs):
 with (R/n).open('w',encoding='utf8',newline='') as h:w=csv.DictWriter(h,fieldnames=list(dict.fromkeys(k for r in rs for k in r)));w.writeheader();w.writerows(rs)
write('QUALIFICATION_ALL_SOURCE.csv',qual);write('INPUT_ROWS.csv',rr)
np.savez_compressed(R/'INPUTS.npz',X=np.stack(xx),construct_ids=np.asarray([r['construct_id'] for r in rr]),new_RNA_primary=np.asarray([r['new_RNA_primary'] for r in rr],bool))
facts={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_processed_rows':len(source),'source_human_design_rows':len(mi),'qualified_all_human_rows':len(rr),'qualified_new_RNA_rows':sum(r['new_RNA_primary'] for r in rr),'RNA_all':sorted({r['official_gene'] for r in rr}),'RNA_new':sorted({r['official_gene'] for r in rr if r['new_RNA_primary']}),'distinct_window_keys':len({r['unique_key'] for r in rr}),'crop':'4:105','model_queries_yet':0,'expected_queries':len(rr)*24,'source_files_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [R/'ALL_SOURCE_ROWS.json',R/'REFERENCE_MAPPING.json',R/'INPUTS.npz',R/'INPUT_ROWS.csv',R.parent/'PROTOCOL_FR2.md',Path(__file__)]}}
(R/'INPUT_FREEZE.json').write_text(json.dumps(facts,indent=2),encoding='utf8');print(json.dumps(facts))
