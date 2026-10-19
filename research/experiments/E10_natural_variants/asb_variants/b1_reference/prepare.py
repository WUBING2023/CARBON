"""Freeze paired sequence edits using original outputs, with no inferential recalculation."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent.parent;SRC=ROOT/'sequence_landscape_20261005'/'direction_b';RAW=ROOT/'02_external_preferences'/'results'/'wide_extracted'
d=np.load(SRC/'display_packet.npz');score=d['protein_score'];external=d['external_mean_Z5'];km=d['kmers'].tolist();pp=d['proteins'].tolist();c=json.loads((SRC/'CASE_SELECTION.json').read_text())
(OUT/'display_packet.npz').write_bytes((SRC/'display_packet.npz').read_bytes())
packet={};meta={};sources=[]
def source(p):sources.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
source(SRC/'display_packet.npz')
controls=['UCAAC','CAUAC','CUACA'];words=controls+['CUAAC'];ix=[c['QKI']['words_RNA'].index(w) for w in words]
packet['QKI_words_RNA']=np.array(words);packet['QKI_assay_Z5']=d['QKI_external_Z5'][ix]
packet['QKI_native_frequency']=d['QKI_frequency_native_count'][ix]
for arm in ['full','sequence']:
    z=d['QKI_'+arm+'_delta_log_count'][:,ix];packet['QKI_'+arm+'_word_delta_log_count']=z
    packet['QKI_'+arm+'_paired_edit_delta']=z[:,-1,None].astype('float64')-z[:,:3].astype('float64')
for key in ['context_ids','native_sequences','native_conservation','cells']:packet['QKI_'+key]=d['QKI_'+key]
meta['QKI']={'words':words,'anchor':'CUAAC','controls':controls,'selected':'all three adjacent different-base exchange neighbors of frozen external anchor; no modelbenefit selection','n_native_contexts':32,'native_example_id':str(d['QKI_context_ids'][0]),'native_center_RNA':str(d['QKI_native_sequences'][0])[48:53].replace('T','U')}

words=['UGAGG','UAGGG'];wi=[km.index(w.replace('U','T')) for w in words];pi=pp.index('HNRNPA1')
packet['HNRNPA1_words_RNA']=np.array(words);packet['HNRNPA1_assay_Z5']=external[pi,wi];packet['HNRNPA1_native_frequency']=score[pi,2,wi]
ctx=[];seqs=[];cells=[];phi=[];byarm={'full':[],'sequence':[]}
tasks=['HNRNPA1_HepG2_ENCSR769UEW','HNRNPA1_K562_ENCSR154HRN']
for task in tasks:
    inp=RAW/'wide_inputs'/(task+'.npz');source(inp)
    a=np.load(RAW/'wide_outputs'/task/'full'/'responses.npz');b=np.load(RAW/'wide_outputs'/task/'sequence'/'responses.npz')
    for arm,z in [('full',a),('sequence',b)]:
        fp=RAW/'wide_outputs'/task/arm/'responses.npz';source(fp);source(fp.parent/'summary.json')
        assert np.array_equal(z['kmers'],d['kmers']);byarm[arm].append(z['delta_log_count'][:,wi])
    for key in ['chromosome','source_row','sequence']:assert np.array_equal(a[key],b[key])
    assert np.array_equal(a['native_conservation'],a['conditional_conservation'])
    for j in range(16):
        ctx.append(task+'|chr'+str(int(a['chromosome'][j]))+'|row'+str(int(a['source_row'][j])));seqs.append(str(a['sequence'][j]));cells.append('HepG2' if 'HepG2' in task else 'K562');phi.append(a['native_conservation'][j])
for arm,z in byarm.items():
    z=np.vstack(z);packet['HNRNPA1_'+arm+'_word_delta_log_count']=z;packet['HNRNPA1_'+arm+'_paired_edit_delta']=z[:,1].astype('float64')-z[:,0].astype('float64')
packet['HNRNPA1_context_ids']=np.array(ctx);packet['HNRNPA1_native_sequences']=np.array(seqs);packet['HNRNPA1_cells']=np.array(cells);packet['HNRNPA1_native_conservation']=np.array(phi)
meta['HNRNPA1']={'words':words,'selected':'perprotein largest independent externalZ5 gap among frequencyopposed exactadjacent swaps; root chose HNRNPA1 (second-largest externalgap) for clear sequence illustration and paired twocell data; exploratory, not preregistered; highest SRSF7 retained in full25candidate table','n_native_contexts':32,'native_example_id':ctx[0],'native_center_RNA':seqs[0][48:53].replace('T','U'),'tasks':tasks}
metrics=pd.read_csv(SRC/'FROZEN_PRIMARY_METRICS.csv',float_precision='round_trip');packet['cohort_proteins']=np.array(pp)
packet['cohort_order_accuracy_percent']=np.array([[metrics[(metrics.protein==p)&(metrics.arm==arm)].iloc[0]['within_composition_concordance']*100 for arm in ['full','sequence','frequency']] for p in pp])
(OUT/'FROZEN_PRIMARY_METRICS.csv').write_bytes((SRC/'FROZEN_PRIMARY_METRICS.csv').read_bytes());(OUT/'FROZEN_PRIMARY_CONTRASTS.csv').write_bytes((SRC/'FROZEN_PRIMARY_CONTRASTS.csv').read_bytes())
np.savez_compressed(OUT/'INTERVENTION_PACKET.npz',**packet)
stats=[]
for p in ['QKI','HNRNPA1']:
    lo=packet[p+'_words_RNA'].tolist();z=packet[p+'_assay_Z5'];freq=packet[p+'_native_frequency'];n=3 if p=='QKI' else 1
    for j in range(n):
        r={'protein':p,'low_word':lo[j],'high_word':lo[-1],'independent_delta_Z5':float(z[-1]-z[j]),'native_frequency_delta':float(freq[-1]-freq[j])}
        for arm in ['full','sequence']:
            dx=packet[p+'_'+arm+'_paired_edit_delta'];dx=dx[:,j] if p=='QKI' else dx
            r[arm+'_mean_delta_log_count']=float(dx.mean());r[arm+'_median_delta_log_count']=float(np.median(dx));r[arm+'_positive_context_count']=int((dx>0).sum());r[arm+'_minimum']=float(dx.min());r[arm+'_maximum']=float(dx.max())
        stats.append(r)
pd.DataFrame(stats).to_csv(OUT/'DISPLAYED_EDITS.csv',index=False,float_format='%.17g')
meta['shared']={'center_zero_based':[48,53],'center_one_based':[49,53],'full_Phi_input_fixed_to_original_native_vector':True,'Sequence_has_no_Phi_input':True,'wordpair_equalcomposition_not_nativeequalcomposition':True,'paired_delta':'subtract existing wordoutputs afterfloat64promotion; native baseline cancels; no refit/normalization/bootstrap','model_output_unit':'original delta_log_count; not affinity, raw count ratio or biological mutant response','32_contexts_not_biological_replication':True,'sources':sources,'displayed_edits':stats}
(OUT/'INTERVENTION_SOURCE.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({'HNRNPA1_native_center_RNA':meta['HNRNPA1']['native_center_RNA'],'edits':stats},ensure_ascii=False))
