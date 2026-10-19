from pathlib import Path
import csv,json,hashlib,datetime
import numpy as np
R=Path(__file__).resolve().parent;B=R.parents[2];A=B/'supervised_interpretability_20261003/02_execution/A'
ll=sorted(list(csv.DictReader((A/'LABELS_ALL150.csv').open(encoding='utf-8-sig'))),key=lambda r:(int(r['position0']),r['ALT_base']))
inp=np.load(A/'INPUTS.npz');w=inp['WT109'];assert np.all(w==w[0]) and w.shape==(4,109,5)
wt=ll[0]['WT_sequence_dna_109'];assert ''.join('ACGT'[i] for i in w[0,:,:4].argmax(1))==wt
xx=[w[0,4:105].copy()];ids=['NucLibC:Jpx_9'];meta=[{'construct_id':ids[0],'is_WT':True,'position0':-1,'ALT_base':''}]
for q in ll:
 x=w[0,4:105].copy();pos=int(q['position0'])-4;alt=q['ALT_sequence_dna_109'];assert sum(a!=b for a,b in zip(wt,alt))==1 and wt[pos+4]!=q['ALT_base'];x[pos,:4]=0;x[pos,'ACGT'.index(q['ALT_base'])]=1;assert ''.join('ACGT'[i] for i in x[:,:4].argmax(1))==alt[4:105];xx.append(x);ids.append(q['variant_construct_id']);meta.append({'construct_id':ids[-1],'is_WT':False,'position0':pos+4,'ALT_base':q['ALT_base']})
np.savez_compressed(R/'INPUTS_MATCHED.npz',X=np.stack(xx),construct_ids=np.asarray(ids));(R/'INPUT_METADATA.json').write_text(json.dumps(meta,indent=2),encoding='utf8')
(R/'MATCHED_INPUT_FREEZE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [A/'INPUTS.npz',A/'LABELS_ALL150.csv',R/'INPUTS_MATCHED.npz',R.parent/'PROTOCOL_FR3.md',Path(__file__)]},'all151_inputs':True,'WT_Phi_identical_ALT':True,'not_ALT_specific_conservation':True},indent=2),encoding='utf8')
print('151 fixed same-WT-Phi inputs built')
