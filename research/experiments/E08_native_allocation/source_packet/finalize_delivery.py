from pathlib import Path
import json,sys,hashlib,datetime
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
import numpy as np,pandas as pd,torch
from pypdf import PdfReader
R=Path(__file__).resolve().parent;O=R/'figure';C=R/'nogru_position_capacity'
out=[]
for p in sorted((C/'remote/outputs').glob('*/seed*/best.pth')):
 s=torch.load(p,map_location='cpu',weights_only=True)
 keys=[k for k in s if 'position_mlp' in k or 'position_profile_bias' in k]
 d={'task':p.parent.parent.name,'seed':p.parent.name,'checkpoint_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'learned_tensors':{k:{'l1':float(s[k].abs().sum()),'l2':float(torch.linalg.vector_norm(s[k])),'max_abs':float(s[k].abs().max())} for k in keys}}
 assert d['learned_tensors']['mixer.position_mlp.2.weight']['l2']>0
 assert d['learned_tensors']['position_profile_bias']['l2']>0
 out.append(d)
(C/'CONTROL_LEARNING_QA.json').write_text(json.dumps({'six_complete_checkpoints':len(out),'all_zero_initialized_private_tensors_trained_nonzero':len(out)==6,'checkpoints':out},indent=2),encoding='utf-8')
plot=pd.read_csv(O/'FINAL_PLOTTED_GENES.csv');q=json.loads((O/'FINAL_FIGURE_QA.json').read_text());lims=q['c_axis_limits']
assert ((plot.loss_delta>=lims[0])&(plot.loss_delta<=lims[1])).all()
assert len(plot)==865+688
pdf=PdfReader(O/'CONTEXT_FINAL.pdf');dim=[float(pdf.pages[0].mediabox.width)*25.4/72,float(pdf.pages[0].mediabox.height)*25.4/72]
assert np.max(np.abs(np.asarray(dim)-[85,125]))<1e-5
src=json.loads((O/'FINAL_FIGURE_SOURCE.json').read_text());assert len(src['examples'])==2
q.update({'required_visual_reinspection_after_this_export':False,'actual_png_visual_review':{'worker':'passed after CTCT/header, label and full-range corrections','parent':'root passed final second visual review; no further layout changes','review_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},'pdf_dimensions_mm':dim,'gene_points':len(plot),'ptbp1_task_gene_points':int((plot.protein=='PTBP1').sum()),'hnrnpk_task_gene_points':int((plot.protein=='HNRNPK').sum()),'zero_clipped_gene_points':True,'model_native_share_endpoint_not_IG':True,'hnrnpk_has_no_near_capacity_control':True})
(O/'FINAL_FIGURE_QA.json').write_text(json.dumps(q,indent=2),encoding='utf-8')
manifest=json.loads((C/'ANALYSIS_MANIFEST.json').read_text());manifest['supporting']='Original three input-defined scenes x five Full controls: Bonferroni 15; three chromosome ranges: Bonferroni 45. Three center-distance scenes remain post hoc diagnostics, with Full-minus-local-position Bonferroni 3. All six scenes x five controls x three ranges additionally receive Bonferroni 90; no diagnostic replaces the main two-contrast family.';manifest['primary']='PTBP1 combined chr1/X/Y, all_core, strict input-adjusted. Main two contrasts: Full-minus-NoGRU_position_capacity and NoGRU_position_capacity-minus-original NoGRU; paired joint-gene cluster bootstrap, Bonferroni 2.'
(C/'ANALYSIS_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps({'control_learning_qa':str(C/'CONTROL_LEARNING_QA.json'),'pdf_dimensions_mm':dim,'gene_points':len(plot),'learned_last_mlp_l2':[v['learned_tensors']['mixer.position_mlp.2.weight']['l2'] for v in out],'bias_l1':[v['learned_tensors']['position_profile_bias']['l1'] for v in out]}))
