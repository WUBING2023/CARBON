from pathlib import Path
import csv,json,hashlib,datetime,time
import numpy as np,openpyxl
R=Path(__file__).resolve().parent;B=R.parents[1];A=B/'supervised_interpretability_20261003/02_execution/A';D=B/'two_directions_exploration_20261004/direction2';O=R/'fr3';O.mkdir(exist_ok=True)
W=B/'expanded_functional_validation_20261003/02_new_public_datasets/sources/hnrnkp_springer_44318_2021_BFEMBJ2020106357_MOESM5_ESM.xlsx'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def write(n,rs):
 with (O/n).open('w',encoding='utf8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rs for k in r)));w.writeheader();w.writerows(rs)
def ranks(v):
 a=np.asarray(v,float);_,ii,n=np.unique(a,return_inverse=True,return_counts=True);e=np.cumsum(n);return ((e+e-n-1)/2)[ii]
def corr(a,b):
 a=np.asarray(a,float);b=np.asarray(b,float);a=a-a.mean();b=b-b.mean();d=np.linalg.norm(a)*np.linalg.norm(b);return float(a@b/d) if d>1e-12 else None
def rho(a,b):return corr(ranks(a),ranks(b)) if len(a)>3 else None
def partial(a,b,c):
 a=ranks(a);b=ranks(b);z=np.column_stack([np.ones(len(a)),ranks(c)]);return corr(a-z@np.linalg.lstsq(z,a,rcond=None)[0],b-z@np.linalg.lstsq(z,b,rcond=None)[0])
def top_share(s,e,k=10):
 s=np.asarray(s);e=np.asarray(e);t=np.sort(s)[-k];gt=s>t;tie=s==t;w=gt.astype(float);w[tie]=(k-gt.sum())/tie.sum();assert np.isclose(w.sum(),k);return float(w@e/e.sum()) if e.sum()>0 else None
start=time.perf_counter();files=[W,A/'ATTRIBUTIONS.npz',A/'LABELS_ALL150.csv',A/'LABELS_NucLibB_ALL75.csv',D/'VARIANT_EFFECTS_ALL.csv',R/'PROTOCOL_FR3.md',Path(__file__)]
if not (O/'REPLAY_FREEZE.json').exists():
 (O/'REPLAY_FREEZE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files_sha256':{str(p):sha(p) for p in files},'old_binding_exposure':True,'new_NC_endpoint_exploratory':True},indent=2),encoding='utf8')
else:
 initial=json.loads((O/'REPLAY_FREEZE.json').read_text('utf8'))
 assert all(sha(p)==h for p,h in initial['files_sha256'].items() if Path(p)!=Path(__file__))
 (O/'REPLAY_CODE_REPAIR.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'repair':'sha helper now accepts saved string paths; after all result CSV computed, final provenance check failed before first successful completion. No scientific endpoint or source modified. Original freeze retained.','original_code_sha256':initial['files_sha256'][str(Path(__file__))],'current_code_sha256':sha(__file__)},indent=2),encoding='utf8')
book=openpyxl.load_workbook(W,read_only=True,data_only=True);source={}
for lib in ['NucLibB','NucLibC']:
 it=iter(book[lib].values);next(it);head=next(it)
 source[lib]={str(v[0]):dict(zip(head,v)) for v in it if v[0]}
attrs=np.load(A/'ATTRIBUTIONS.npz',allow_pickle=False);tables={};rows=[];metrics=[];contrasts=[];agreements=[];dose=[];den=[]
for lib,label in [('NucLibC','LABELS_ALL150.csv'),('NucLibB','LABELS_NucLibB_ALL75.csv')]:
 labels=sorted(read(A/label),key=lambda r:(int(r['position0']),r['ALT_base']));pos=np.array(sorted({int(r['position0']) for r in labels}));assert len(labels)==3*len(pos)
 wtseq=labels[0]['WT_sequence_dna_109'];wt=source[lib]['Jpx_9'];assert (wt.get('Tile sequence')==wtseq)
 nk=[('rep'+str(j)) if lib=='NucLibC' else ('Nuc/Cyto (log2) rep'+str(j)) for j in range(1,5 if lib=='NucLibC' else 3)]
 wtvals=np.array([float(wt[k]) for k in nk]);altvals=[];rip=[];avgdiff=[]
 for q in labels:
  z=source[lib][q['variant_construct_id'].split(':',1)[1]];assert z['Tile sequence']==q['ALT_sequence_dna_109'];assert sum(a!=b for a,b in zip(wtseq,z['Tile sequence']))==1
  vv=np.array([float(z[k]) for k in nk]);assert np.isfinite(vv).all();delta=float(vv.mean()-wtvals.mean());altvals.append(vv-wtvals)
  rip.append(float(q['HNRNPK_d_repmean']));avgkey='Nuc/Cyto (log2) average';avgdiff.append(abs(float(z[avgkey])-vv.mean()))
  rows.append({'library':lib,'construct_id':q['variant_construct_id'],'position0':q['position0'],'position1':q['position1'],'WT_base':q['WT_base'],'ALT_base':q['ALT_base'],'NC_delta_repmean':delta,'NC_paired_report_delta_replicates_json':json.dumps((vv-wtvals).tolist()),'HNRNPK_RIP_delta_repmean':q['HNRNPK_d_repmean'],'SLTM_RIP_delta_repmean':q['SLTM_d_repmean'],'WT_NC_repmean':float(wtvals.mean()),'ALT_NC_repmean':float(vv.mean()),'NC_n_reported_reps':len(nk),'source_sheet_row':q['source_sheet_row']})
 delta=np.asarray(altvals).mean(1).reshape(len(pos),3);E=np.mean(abs(delta),1);L=-delta.mean(1);BR=np.array(rip).reshape(len(pos),3);EB=np.mean(abs(BR),1)
 tables[lib]={'pos':pos,'E':E,'L':L,'E_RIP':EB,'delta':delta,'delta_reps':np.asarray(altvals).reshape(len(pos),3,len(nk))}
 den.append({'library':lib,'positions':len(pos),'variants':len(labels),'RNA_n':1,'all_ALT_complete':True,'NC_all_reps_finite':True,'NC_reps':len(nk),'source_author_average_max_difference':max(avgdiff),'binding_vs_function_position_rho':rho(EB,E)})
 for task in sorted({k.split('__')[0] for k in attrs.files}):
  keys=[k for k in attrs.files if k.startswith(task+'__crop4__') and k.endswith('__importance')];ss={}
  for k in keys:
   cond=k.split('__')[2];ss[cond]=attrs[k].ravel()[pos-4]
  phi=attrs[task+'__crop4__A_native__condition_input'].reshape(101,5)[:,4];cov=set()
  for motif in ['CCTCCC','GCCTCCC']:
   for i in range(len(wtseq)-len(motif)+1):
    if wtseq[i:i+len(motif)]==motif:cov.update(range(i,i+len(motif)))
  ss.update({'core_coverage':np.array([p in cov for p in pos],float),'C7':np.array([wtseq[max(0,p-3):p+4].count('C')/len(wtseq[max(0,p-3):p+4]) for p in pos]),'CT7':np.array([sum(wtseq[max(0,p-3):p+4].count(b) for b in 'CT')/len(wtseq[max(0,p-3):p+4]) for p in pos]),'Phi':phi[pos-4],'absPhi':abs(phi[pos-4]),'distance_center':-abs(pos-4-50)})
  rows_task=[]
  for cond,s in ss.items():
   for p,i,e,l,b in zip(pos,s,E,L,EB):rows_task.append({'library':lib,'task':task,'condition':cond,'position0':int(p),'importance':float(i),'function_E_NC':float(e),'function_L_NC':float(l),'binding_E_RIP':float(b)})
   rr=rho(s,E);null=[rho(np.roll(s,k),E) for k in range(1,len(s))]
   mm={'library':lib,'task':task,'condition':cond,'n_positions':len(pos),'RNA_n':1,'rho_function_E_NC':rr,'rho_function_L_NC':rho(s,L),'rho_binding_E_RIP':rho(s,EB),'partial_function_E_controlling_binding_E':partial(s,E,EB),'top10_function_E_share':top_share(s,E),'shift_p':(1+sum(n is not None and n>=rr for n in null))/len(pos) if rr is not None else None}
   metrics.append(mm)
  write(lib+'_'+task+'_POSITION_MAP.csv',rows_task)
  a=next(m for m in metrics if m['library']==lib and m['task']==task and m['condition']=='A_native')
  perm=[m['rho_function_E_NC'] for m in metrics if m['library']==lib and m['task']==task and m['condition'] in ['A_perm0','A_perm1','A_perm2']]
  for cc in ['S','A_flat','no_gru','no_se','core_coverage','C7','CT7']:
   cm=next((m for m in metrics if m['library']==lib and m['task']==task and m['condition']==cc),None)
   if cm:contrasts.append({'library':lib,'task':task,'comparison':'A_native_minus_'+cc,'delta_rho_function_E':a['rho_function_E_NC']-cm['rho_function_E_NC'],'RNA_n':1})
  contrasts.append({'library':lib,'task':task,'comparison':'A_native_minus_mean_perm_rho','delta_rho_function_E':a['rho_function_E_NC']-np.mean(perm),'RNA_n':1})
 for rep in range(len(nk)):
  for other in range(rep+1,len(nk)):
   e1=np.mean(abs(tables[lib]['delta_reps'][:,:,rep]),1);e2=np.mean(abs(tables[lib]['delta_reps'][:,:,other]),1);agreements.append({'scope':lib,'endpoint':'position_E_NC','rep1':rep+1,'rep2':other+1,'n':len(pos),'rho':rho(e1,e2)})
agreements.append({'scope':'B_vs_C_common25','endpoint':'position_E_NC','n':25,'rho':rho(tables['NucLibC']['E'][-25:],tables['NucLibB']['E'])})
write('LABELS_FUNCTION_ALL.csv',rows);write('METRICS_OLD_IG_ALL.csv',metrics);write('CONTRASTS_OLD_IG_ALL.csv',contrasts);write('MEASUREMENT_AGREEMENT.csv',agreements);write('DENOMINATORS.csv',den)
(O/'REPLAY_RESULT.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.perf_counter()-start,'new_model_queries':0,'main_function_endpoint':'E_NC, not RIP','scientific_source_sha256_unchanged':all(sha(p)==h for p,h in json.loads((O/'REPLAY_FREEZE.json').read_text('utf8'))['files_sha256'].items() if Path(p)!=Path(__file__)),'execution_code_sha256':sha(__file__),'single_RNA':True},indent=2),encoding='utf8')
for m in metrics:
 if m['library']=='NucLibC' and ('HNRNPK' in m['task']) and m['condition'] in ['A_native','S','A_flat','no_gru','no_se','core_coverage','C7','CT7']:print(json.dumps(m))
print(json.dumps({'denominators':den,'agreement':agreements}))
