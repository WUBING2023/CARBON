from pathlib import Path
import json,hashlib,datetime,zipfile
R=Path(__file__).resolve().parent;CAM=R.parent
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while True:
   b=f.read(4*1024*1024)
   if not b:break
   h.update(b)
 return h.hexdigest()
audits=[CAM/'02_external_preferences/results'/n for n in ['CROSS_AUDIT_01_PC_QA.json','CROSS_AUDIT_01_PC_FINAL_CI_QA.json','CROSS_AUDIT_01_PC_TRAINED_BRANCHES.json','CROSS_AUDIT_01_POSITION_CAPACITY_MODEL.json']]+[CAM/'05_delivery/CONTEXT_FIGURE_QA.json']
for p in audits:assert p.is_file(),p
scope=json.loads(audits[-1].read_text());fqa=json.loads((R/'figure/FINAL_FIGURE_QA.json').read_text());resource=json.loads((R/'RESOURCE_RELEASE.json').read_text())
assert scope['pass'] and fqa['zero_clipped_gene_points'] and not resource['own_active_python_processes']
assert sha(R/'figure/CONTEXT_FINAL.png')==scope['PNG_sha256']
assert sha(R/'figure/FINAL_FIGURE_SOURCE.json')==scope['source_json_sha256']
assert sha(R/'figure/FINAL_PLOTTED_GENES.csv')==scope['plotted_genes_sha256']
qa={'status':'passed, final candidate experiment/figure delivery; official manuscript not overwritten','frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'final_figure':'figure/CONTEXT_FINAL.png','caption':'FINAL_CAPTION.md','results':'RESULTS.md','figure_machine_qa':fqa,'independent_audits':[{'path':str(p),'sha256':sha(p)} for p in audits],'resource_release':resource,'retained_negative_evidence':['root conservation controls','phase2_ambiguity all four arms and seeds','test_8_15_support old exposed rule support','exact-distance/at-most2 diagnostics','native identical17nt+Phi zero cohort','fixed original counterexamples'],'scope_limits':['native two-support signal allocation, not binding probability or IG accuracy','PTBP1 near-capacity position control only; HNRNPK original NoGRU support','no conservation-specific or downstream functional success claimed','model-heldout chr1/X/Y; prior study exposure retained','RBFOX2 proposed full fits not executed for fixed time budget']}
(R/'DELIVERY_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
(R/'PACKAGE_CONTENTS.md').write_text('''DELIVERY_SOURCE_BUNDLE.zip contains source scripts, protocols, English Results/Caption, Chinese report, complete CSV results, metadata/QA/source manifests, final images and original fixed counterexamples. DELIVERY_MANIFEST.json inventories the entire experiment directory including locally retained raw arrays, bootstrap draws and checkpoints. Those large binary replay artifacts remain at their recorded local paths and are not duplicated inside the compact source bundle. The bundle is therefore a review packet, while exact numerical replay uses the full local directory described in REPLAY.md. No external upload or publication occurred. Historical precontrol figures and reports remain explicitly named in the full inventory.
''',encoding='utf-8')
excluded={'DELIVERY_MANIFEST.json','DELIVERY_SOURCE_BUNDLE.zip','BUNDLE_SHA256.txt'}
files=[]
for p in sorted(R.rglob('*')):
 if not p.is_file() or p.name in excluded or '__pycache__' in p.parts or p.suffix in {'.pyc','.tar','.gz','.zip'}:continue
 files.append({'path':p.relative_to(R).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
manifest={'frozen_utc':qa['frozen_utc'],'root':str(R),'files':files,'excluded':['self manifest','transport tar/gz/zip archives','Python bytecode','bundle digest file'],'no_new_compute_pending':True,'final_image_sha256':sha(R/'figure/CONTEXT_FINAL.png'),'final_caption_sha256':sha(R/'FINAL_CAPTION.md'),'results_sha256':sha(R/'RESULTS.md')}
(R/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
bundle=R/'DELIVERY_SOURCE_BUNDLE.zip';included=[]
with zipfile.ZipFile(bundle,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files:
  p=R/f['path'];suffix=p.suffix.lower()
  keep=suffix in {'.py','.md','.json','.csv'} or (p.parent==R/'figure' and p.stem in {'CONTEXT_FINAL','FIXED_COUNTEREXAMPLES'} and suffix in {'.png','.pdf','.svg'})
  if keep:z.write(p,f['path']);included.append(f['path'])
 z.write(R/'DELIVERY_MANIFEST.json','DELIVERY_MANIFEST.json')
(R/'BUNDLE_SHA256.txt').write_text(sha(bundle)+'  DELIVERY_SOURCE_BUNDLE.zip\n',encoding='ascii')
print(json.dumps({'manifest_files':len(files),'inventory_bytes':sum(v['bytes'] for v in files),'bundle_bytes':bundle.stat().st_size,'bundle_entries':len(included)+1,'bundle_sha256':sha(bundle),'PNG_sha256':manifest['final_image_sha256'],'caption_sha256':manifest['final_caption_sha256'],'results_sha256':manifest['results_sha256']}))
