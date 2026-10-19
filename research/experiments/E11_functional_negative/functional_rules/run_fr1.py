"""Read-only replay of historical model predictions against matched reporter function.
No training, label-selection, remote model calls, or edits to historical sources.
"""
from pathlib import Path
import csv, json, hashlib, math, itertools, time, datetime
import numpy as np
import openpyxl

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parents[1] / 'expanded_functional_validation_20261003'
WORKBOOK = OLD / '02_new_public_datasets/sources/hnrnkp_springer_44318_2021_BFEMBJ2020106357_MOESM5_ESM.xlsx'
PRED = OLD / '08_HNRNPK_SLTM_validation'
OUT = ROOT / 'fr1'
OUT.mkdir(exist_ok=True)

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def num(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError): return None
def mean(v):
    a=[num(x) for x in v];a=[x for x in a if x is not None]
    return float(np.mean(a)) if a else None
def ranks(v):
    v=np.asarray(v,float)
    _,ii,n=np.unique(v,return_inverse=True,return_counts=True)
    end=np.cumsum(n);start=end-n
    return ((start+end-1)/2)[ii]
def corr(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if len(a)<2:return None
    a=a-a.mean();b=b-b.mean();d=np.linalg.norm(a)*np.linalg.norm(b)
    return float(a@b/d) if d>1e-14 else None
def rho(a,b):return corr(ranks(a),ranks(b))
def partial(a,b,c):
    x=ranks(a);y=ranks(b)
    z=np.column_stack([np.ones(len(x))]+[ranks(v) for v in c])
    x=x-z@np.linalg.lstsq(z,x,rcond=None)[0]
    y=y-z@np.linalg.lstsq(z,y,rcond=None)[0]
    return corr(x,y)
def write_csv(name,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/name).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def motif_count(seq,motif):return sum(seq[i:i+len(motif)]==motif for i in range(len(seq)-len(motif)+1))

start=time.perf_counter()
inputs={str(p):sha(p) for p in [WORKBOOK,PRED/'ROWS.csv',PRED/'INPUTS.npz',PRED/'PREDICTIONS.npz',PRED/'PROTOCOL.json',PRED/'RESULT.json',ROOT/'PROTOCOL_FR1.md',Path(__file__)]}
(OUT/'REPLAY_FREEZE.json').write_text(json.dumps({'utc':now(),'before_new_endpoint_analysis':True,'files_sha256':inputs,'data_exposure':'Historical predictions and binding results exposed. New function endpoint is exploratory; no blind claim.'},indent=2),encoding='utf8')

# Read original reporter columns; keep every qualification and author formula discrepancy.
book=openpyxl.load_workbook(WORKBOOK,read_only=True,data_only=True)
s=book['NucLibB'];it=iter(s.values);next(it);head=next(it)
source={}
for rn,v in enumerate(it,3):
    row=dict(zip(head,v));cid=row.get('ID')
    if not cid:continue
    assert cid not in source,cid
    row['source_row']=rn;source[cid]=row
old=list(csv.DictReader((PRED/'ROWS.csv').open(encoding='utf8')))
proto=json.loads((PRED/'PROTOCOL.json').read_text('utf8'));result=json.loads((PRED/'RESULT.json').read_text('utf8'))
assert result['predictions_sha256']==sha(PRED/'PREDICTIONS.npz')
assert result['checkpoint_hashes_verified'] and result['strict_load']
x=np.load(PRED/'INPUTS.npz',allow_pickle=False);pred=np.load(PRED/'PREDICTIONS.npz',allow_pickle=False)
assert x['construct_ids'].tolist()==[r['construct_id'] for r in old]
n=len(old);assert n==proto['rows_per_task']==4232
methods=['A_native','S','A_flat','A_permutation_mean']
tasks=proto['tasks_in_block_order'];scores={}
for ti,task in enumerate(tasks):
    scores[task]={c:pred['101__'+c][ti*n:(ti+1)*n].astype(float) for c in ['A_native','S','A_flat','A_perm0','A_perm1','A_perm2']}
    scores[task]['A_permutation_mean']=np.mean([scores[task]['A_perm'+str(j)] for j in range(3)],axis=0)
qual=[];eligible=[];formula=[]
for i,r in enumerate(old):
    q={'construct_id':r['construct_id'],'input_index':i,'gene':r['reference_gene_name'],'library':r['library'],'is_WT':r['author_WT']=='True'}
    if r['library']!='NucLibB' or r['author_WT']!='True':
        q.update(eligible=False,reason='not_NucLibB_natural_WT');qual.append(q);continue
    z=source[r['author_id']]
    assert z['Tile sequence'].upper()==r['sequence_dna109']
    nt=[z[f'siNT Nuc/Cyto (log2) rep{j}'] for j in [1,2,3]]
    vals={protein:[z[f'si{protein} Nuc/Cyto (log2) rep{j}'] for j in [1,2,3]] for protein in ['HNRNPK','SLTM','SNRNP70']}
    if mean(nt) is None or any(mean(v) is None for v in vals.values()):
        q.update(eligible=False,reason='missing_primary_function_condition');qual.append(q);continue
    seq=r['sequence_dna109'][4:105];assert len(seq)==101
    rr={**q,'sequence_dna109':r['sequence_dna109'],'sequence101':seq,'chromosome':r['chromosome'],'coordinates109':r['reference_coordinates109_0based_json'],
        'source_row':z['source_row'],'source_species_status':r['source_species_status'],'source_proven_human':r['source_species_proven_human']=='True',
        'baseline_NC':num(z['Nuc/Cyto (log2) average']),'siNT_NC':mean(nt),'siNT_rep_count':sum(num(v) is not None for v in nt),
        'measured_HNRNPK_RIP':num(r['HNRNPK_repmean_log2_IP_Input']),'measured_SLTM_RIP':num(r['SLTM_repmean_log2_IP_Input']),
        'C_content':seq.count('C')/101,'CCC':motif_count(seq,'CCC'),'CCCC':motif_count(seq,'CCCC'),'CCTCCC':motif_count(seq,'CCTCCC'),'GCCTCCC':motif_count(seq,'GCCTCCC'),
        'mean_Phi':float(x['X_101'][i,:,4].mean()),'unique_key':r['reference_gene_name']+'|'+r['sequence_dna109']}
    for protein,v in vals.items():
        rr['D_'+protein]=mean(nt)-mean(v)
        rr['si'+protein+'_NC']=mean(v)
        rr['si'+protein+'_rep_count']=sum(num(z) is not None for z in v)
        author=num(z['Effect of si'+protein+' on Nuc/Cyto (log2) ratios'])
        formula.append({'construct_id':r['construct_id'],'protein':protein,'author_effect':author,'computed_KD_minus_siNT':-rr['D_'+protein],
            'abs_difference':abs(author+rr['D_'+protein]) if author is not None else None})
        for j,vv in enumerate(v,1):rr[f'si{protein}_rep{j}']=num(vv)
    for j,vv in enumerate(nt,1):rr['siNT_rep'+str(j)]=num(vv)
    rr['all_four_conditions_three_reps']=all(rr[z]==3 for z in ['siNT_rep_count','siHNRNPK_rep_count','siSLTM_rep_count','siSNRNP70_rep_count'])
    q.update(eligible=True,reason='source_id_and_sequence_match_native_Phi_complete_function_present')
    qual.append(q);eligible.append(rr)
write_csv('QUALIFICATION_ALL_INPUTS.csv',qual);write_csv('LABELS_ELIGIBLE.csv',eligible);write_csv('AUTHOR_KD_FORMULA_CHECK.csv',formula)

# Duplicate windows are merged before analysis; overlap across different windows retained in RNA clusters.
groups={}
for r in eligible:groups.setdefault(r['unique_key'],[]).append(r)
baselines=['C_content','CCC','CCCC','CCTCCC','GCCTCCC','mean_Phi','measured_HNRNPK_RIP','measured_SLTM_RIP']
endpoints=['D_HNRNPK','D_SLTM','D_SNRNP70','baseline_NC']
table=[]
for task,ss in scores.items():
    for key,rs in groups.items():
        row={k:rs[0][k] for k in ['unique_key','gene','sequence101','sequence_dna109','source_proven_human','chromosome']}
        row['construct_ids']=json.dumps([r['construct_id'] for r in rs]);row['source_row_count']=len(rs)
        row['input_indices']=json.dumps([r['input_index'] for r in rs]);row['task']=task;row['RBP'],row['cell']=task.split('_')[:2]
        row['all_four_conditions_three_reps']=all(r['all_four_conditions_three_reps'] for r in rs)
        for c in methods:row[c]=float(np.mean([ss[c][r['input_index']] for r in rs]))
        for c in baselines+endpoints+['siNT_NC']:row[c]=mean([r[c] for r in rs])
        table.append(row)
write_csv('MODEL_FUNCTION_JOIN.csv',table)

perrna=[];summary=[];contrasts=[];partialrows=[]
all_methods=methods+baselines
for task in tasks:
    taskrows=[r for r in table if r['task']==task]
    for cohort in ['all_reference_compatible','source_proven_human','complete_three_rep']:
        trs=[r for r in taskrows if cohort=='all_reference_compatible' or (r['source_proven_human'] if cohort=='source_proven_human' else r['all_four_conditions_three_reps'])]
        for ep in endpoints:
            records=[]
            for gene in sorted(set(r['gene'] for r in trs)):
                rs=[r for r in trs if r['gene']==gene and r[ep] is not None]
                row={'task':task,'cohort':cohort,'endpoint':ep,'gene':gene,'n_unique_tiles':len(rs),'source_rows':sum(r['source_row_count'] for r in rs),'eligible_rna_min5':len(rs)>=5}
                for c in all_methods:
                    ok=[r for r in rs if r[c] is not None]
                    row[c+'__rho']=rho([r[ep] for r in ok],[r[c] for r in ok]) if len(ok)>=5 else None
                    row[c+'__n']=len(ok)
                perrna.append(row);records.append(row)
                if ep=='D_HNRNPK':
                    for control in ['measured_HNRNPK_RIP','siNT_NC']:
                        ok=[r for r in rs if r[control] is not None]
                        p={'task':task,'cohort':cohort,'gene':gene,'endpoint':ep,'control_rank':control,'n_unique_tiles':len(ok)}
                        for c in methods+['C_content','CCC','CCTCCC']:
                            p[c+'__partial_rho']=partial([r[c] for r in ok],[r[ep] for r in ok],[[r[control] for r in ok]]) if len(ok)>=8 else None
                        partialrows.append(p)
            for c in all_methods:
                vv=[r[c+'__rho'] for r in records if r[c+'__rho'] is not None]
                summary.append({'task':task,'cohort':cohort,'endpoint':ep,'method':c,'RNA_total':len(records),'RNA_defined':len(vv),'undefined_genes':json.dumps([r['gene'] for r in records if r[c+'__rho'] is None]),'equal_RNA_mean_rho':float(np.mean(vv)) if vv else None,'current_library_only':True})
            for c in ['S','A_flat','A_permutation_mean','C_content','CCC','CCCC','CCTCCC','GCCTCCC','mean_Phi']:
                rs=[r for r in records if r['A_native__rho'] is not None and r[c+'__rho'] is not None]
                v=np.array([r['A_native__rho']-r[c+'__rho'] for r in rs]);rng=np.random.default_rng(20261004)
                boot=v[rng.integers(0,len(v),(10000,len(v)))].mean(1) if len(v) else np.array([])
                contrasts.append({'task':task,'cohort':cohort,'endpoint':ep,'comparison':'A_native_minus_'+c,'RNA_n':len(v),'paired_genes':json.dumps([r['gene'] for r in rs]),
                    'delta_rho':float(v.mean()) if len(v) else None,'CI95_low':float(np.quantile(boot,.025)) if len(v) else None,'CI95_high':float(np.quantile(boot,.975)) if len(v) else None,
                    'Bonf3_low':float(np.quantile(boot,.05/6)) if len(v) else None,'Bonf3_high':float(np.quantile(boot,1-.05/6)) if len(v) else None,
                    'A_win_RNA':int((v>0).sum()),'A_loss_RNA':int((v<0).sum()),'A_tie_RNA':int((v==0).sum()),'descriptive_shared_study':True})
write_csv('PER_RNA_METRICS.csv',perrna);write_csv('SUMMARY.csv',summary);write_csv('CONTRASTS.csv',contrasts);write_csv('SECONDARY_PARTIAL_RANK.csv',partialrows)
facts={'completed_utc':now(),'cpu_wall_seconds':time.perf_counter()-start,'new_model_queries':0,'historical_queries_replayed':len(table)*4,
    'source_WT_eligible_rows':len(eligible),'unique_windows':len(groups),'reference_RNAs':sorted(set(r['gene'] for r in eligible)),
    'proven_human_RNAs':sorted(set(r['gene'] for r in eligible if r['source_proven_human'])),
    'KD_author_effect_formula_max_error':max(r['abs_difference'] for r in formula if r['abs_difference'] is not None),
    'main_scope':'NucLibB reporter WT only, not endogenous gene-level KD','primary_endpoint':'D_HNRNPK = condition means siNT minus siHNRNPK nuclear/cytoplasmic log2',
    'source_files_sha256':inputs,'model_family':'old paired A/S; flat/permuted are frozen input controls, not trained controls','comparison_cells_all_reported':tasks,
    'biological_unit':'reference RNA; overlapping windows and experimental replicates not independent RNA','results_are_exploratory':True}
(OUT/'EXECUTION.json').write_text(json.dumps(facts,indent=2),encoding='utf8')
for row in summary:
    if row['cohort']=='all_reference_compatible' and row['endpoint']=='D_HNRNPK' and row['method'] in methods+['C_content','CCC','CCTCCC','measured_HNRNPK_RIP']:
        print(json.dumps(row))
print(json.dumps(facts))
