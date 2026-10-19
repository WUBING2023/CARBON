from pathlib import Path
import csv,json,hashlib,datetime,time,sys,io,base64
import numpy as np,torch,subprocess
R=Path('/data/CARBON_PROJECT/supervised_expansion_20261004/03_functional_rules/fr4');M=Path('/data/CARBON_PROJECT/conservation_upgrade_20260929_matched')
sys.path.insert(0,str(M));from matched_common import make_model,transformed,sha
torch.set_num_threads(4);torch.set_num_interop_threads(1)
assert not (R/'RESULT.json').exists(),'Frozen result already exists'
start=time.time();FA=Path('/data/CARBON_PROJECT/数据处理/gencodeV29/GRCh38.p12.genome.fa');BW=Path('/data/CARBON_PROJECT/数据/phyloP/hg38.phyloP100way.bw')
subprocess.run(['/usr/bin/python3',str(R/'remote_extract.py')],check=True)
inp=np.load(R/'INPUTS.npz');x=inp['X'];valid=inp['valid'];meta=json.loads((R/'NATIVE_METADATA.json').read_text());bad_events={r['event_index'] for r in meta if r['native_exclusion']};sources=[];out={}
for cell,task in [('K562','PTBP1_K562_ENCSR981WKN'),('HepG2','PTBP1_HepG2_ENCSR384KAN')]:
 idx=np.array([i for i,r in enumerate(meta) if r['cell']==cell and valid[i]],dtype=np.int64);out[task+'__input_indices']=idx
 pm=[]
 for i in idx:
  r=meta[i];key=task+'|FR4|'+r['chrom']+'|'+r['identity'];seed=int.from_bytes(hashlib.sha256(key.encode()).digest()[:8],'little');pm.append(np.random.default_rng(seed).permutation(101))
 pm=np.array(pm);out[task+'__permutation']=pm
 for seed in [11,23,47]:
  for arm in ['aligned','sequence_only','mean','permuted']:
   ck=M/'outputs'/task/f'seed{seed}'/arm/'best.pth';model,init=make_model(seed,arm);run=json.loads((ck.parent/'run.json').read_text());assert init['actual_initial_state_sha256']==run['actual_initial_state_sha256'];model.load_state_dict(torch.load(ck,map_location='cpu',weights_only=True),strict=True);model.eval();val=[];tx=transformed(x[idx],pm,arm)
   with torch.no_grad():
    for st in range(0,len(tx),128):val.append(model(torch.from_numpy(np.ascontiguousarray(tx[st:st+128])))['log_total_count'].numpy().ravel())
   v=np.concatenate(val);assert len(v)==len(idx) and np.isfinite(v).all();out[f'{task}__{arm}__seed{seed}']=v
   sources.append({'task':task,'arm':arm,'seed':seed,'path':str(ck),'sha256':sha(ck),'run_sha256':sha(ck.parent/'run.json')})
   print(json.dumps({'task':task,'arm':arm,'seed':seed,'n_windows':len(idx),'elapsed_seconds':time.time()-start}),file=sys.stderr,flush=True);del model
buf=io.BytesIO();np.savez_compressed(buf,**out);p=buf.getvalue();(R/'PREDICTIONS.npz').write_bytes(p)
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.time()-start,'queries':int(valid.sum())*12,'n_source_windows':len(meta),'native_valid_windows':int(valid.sum()),'native_excluded_events':len(bad_events),'device':'CPU','threads':4,'models':sources,'input_sha256':sha(R/'INPUTS.npz'),'windows_sha256':sha(R/'WINDOWS.csv'),'protocol_sha256':sha(R/'PROTOCOL_FR4.md'),'addendum_sha256':sha(R/'PROTOCOL_FR4_HOLDOUT_ADDENDUM.md'),'prediction_sha256':sha(R/'PREDICTIONS.npz'),'script_sha256':sha(__file__),'native_files':[{ 'path':str(f),'bytes':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns} for f in [FA,BW]],'minus_strand_sequence_and_Phi_reversed_together':True,'permutation_identity':'sha256(task|FR4|chrom|accession|original_source_row|side) first8 little-endian seed, PCG64, fixed acrossseed/arms','quantity':'native log_total_count; ensemble mean log output; no profile normalization used'}
(R/'RESULT.json').write_text(json.dumps(res,indent=2));files=['INPUTS.npz','NATIVE_METADATA.json','PREDICTIONS.npz','RESULT.json'];print(json.dumps({'files':{f:base64.b64encode((R/f).read_bytes()).decode() for f in files},'summary':res}))
