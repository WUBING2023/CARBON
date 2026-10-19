"""Independent plotted-data replay: no worker script imported or changed."""
from pathlib import Path
import numpy as np,csv,json,itertools,hashlib,datetime
R=Path(__file__).resolve().parent;P=R.parent/'02_external_preferences';O=R/'cross_audit_02'
j=json.loads((P/'figures/FIGURE_DATA.json').read_text());checks=[]
def ck(name,err,tol=1e-12):checks.append({'check':name,'max_error':float(err),'pass':bool(err<=tol)})
def rows(path):return list(csv.DictReader(path.open(encoding='utf8')))
tab={(r['protein'],r['arm']):r for r in rows(P/'results/WIDE_PROTEIN_METRICS.csv')}
new=sorted({r['protein'] for r in tab.values() if r['exposed_protein']=='False'})
ck('availability_defined_all25',int(new!=j['main_new_proteins']))
for r in j['panels_ab']:
 t=tab[r['protein'],r['condition']];ck(f"{r['panel']}/{r['protein']}/{r['condition']}/point",abs(float(t[r['endpoint']])-r['value']))
 ck(f"{r['panel']}/{r['protein']}/{r['condition']}/shared_row",abs(r['y']-(24-j['main_new_proteins'].index(r['protein']))))
km=[''.join(k) for k in itertools.product('ACGT',repeat=5)];ki={k:i for i,k in enumerate(km)}
sums=np.zeros(1024);nums=np.zeros(1024,int)
with (P/'sources/Zscores_420_origrncmpt.txt').open() as f:
 col=f.readline().split().index('RNCMPT00026')
 for line in f:
  v=line.split();s=v[0].replace('U','T');z=float(v[col])
  for k in range(3):q=ki[s[k:k+5]];sums[q]+=z;nums[q]+=1
truth=sums/nums;words=[km[x] for x in sorted(range(1024),key=lambda i:(-truth[i],km[i]))[:20]]
ck('c/all_external_top20_order',int(words!=j['panel_c_words']))
for w in words:ck('c/'+w+'/original_Z5',abs(truth[ki[w]]-j['panel_c_Z5'][w]))
coord=rows(P/'results/locus_mapping/HNRNPK_NATIVE_COORDINATES.csv');tasks=sorted(set(r['task'] for r in coord));members={a:[] for a in ['aligned','sequence_only','mean','permuted']}
for task in tasks:
 for seed in [11,23,47]:
  for arm in members:
   a=np.load(P/'results/extracted/outputs'/task/f'seed{seed}'/arm/'attributions.npz');ig=a['position_ig'];seq=a['sequence'];winner=np.array([np.convolve(v,np.ones(5),mode='valid').argmax() for v in ig]);value=np.array([ig[i,k:k+5].sum() for i,k in enumerate(winner)]);wi=np.array([ki[s[k:k+5]] for s,k in zip(seq,winner)]);score=np.bincount(wi,weights=value,minlength=1024);count=np.bincount(wi,minlength=1024);top=set(sorted(np.flatnonzero(count),key=lambda i:(-score[i],km[i]))[:20]);members[arm].append([int(ki[w] in top) for w in words])
for arm,values in members.items():
 result=np.asarray(values).mean(0)
 for w,v in zip(words,result):ck('c/'+arm+'/'+w+'/selection_fraction',abs(v-j['panel_c_selection'][arm][w]))
dt={(r['assay'],r['control']):r for r in rows(P/'results/HNRNPK_PAIRED_GENE_BOOTSTRAP.csv')}
for r in j['panel_d']:
 t=dt[r['assay'],r['control']]
 for a,b in [('delta_pp','delta_recall20'),('low_pp','conditional_gene_locus_bootstrap_low'),('high_pp','conditional_gene_locus_bootstrap_high')]:ck('d/'+r['assay']+'/'+r['control']+'/'+a,abs(r[a]-100*float(t[b])))
for path,h in j['source_sha256'].items():
 p=Path(path)
 if not p.is_absolute():p=P/p
 if p.exists():ck('source_hash/'+path,int(hashlib.sha256(p.read_bytes()).hexdigest()!=h))
qa={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':all(c['pass'] for c in checks),'checks':len(checks),'max_error':max(c['max_error'] for c in checks),'independent_ab_points':250,'independent_c_matrix_cells':80,'independent_c_external_values':20,'d_points_and_bounds':18,'source_hash':hashlib.sha256((P/'figures/FIGURE_DATA.json').read_bytes()).hexdigest(),'PNG_sha256':hashlib.sha256((P/'figures/EXTERNAL_RULES.png').read_bytes()).hexdigest(),'no_new_experiments':True,'scientific_scope':'Plotted data correctness does not establish Phi-specific, biological-function, or cross-assay advantage.'}
with (O/'FIGURE_CHECKS.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=['check','max_error','pass']);w.writeheader();w.writerows(checks)
(O/'FIGURE_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf8');print(json.dumps(qa));assert qa['pass']
