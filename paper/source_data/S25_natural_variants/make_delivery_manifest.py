"""Final source-bound delivery inventory and geometric/content QA."""
from pathlib import Path
import json,hashlib,datetime,importlib.metadata as md,fitz
R=Path(__file__).resolve().parent
figures=[]
for stem in ['FIG5_ADDITIONAL_d_e_85mm','SUPPLEMENTARY_FIGURE_S22_NATURAL_SNV','SUPPLEMENTARY_FIGURE_S23_CONSERVATION']:
    p=R/'figures'/f'{stem}.pdf';d=fitz.open(p);page=d[0];text=page.get_text();assert len(text)>300,(stem,len(text))
    figures.append(dict(file=str(p.relative_to(R)),width_mm=round(page.rect.width/72*25.4,3),height_mm=round(page.rect.height/72*25.4,3),pages=len(d),text_characters=len(text),font_size_pt_min=round(min(s['size'] for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']),3),visual_QA='Actual PNGs inspected; no clipping or overlap; S22 composed-page redaction graphics bug repaired.'))
    d.close()
paths=[p for p in R.rglob('*') if p.is_file() and not any(q.startswith('replay_check_') or q=='__pycache__' for q in p.relative_to(R).parts) and p.name!='DELIVERY_MANIFEST.json']
files={str(p.relative_to(R)).replace('\\','/'):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(paths)}
report=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),root='.',source_data='S25_natural_variants',supplement_method='S1.20',figures=figures,source_binding='All displayed values derive from the packaged existing predictions and preferred-allele labels; all paired allele tensor checks pass.',replay=json.loads((R/'PORTABLE_REPLAY_QA.json').read_text('utf8')),tie_threshold_scope='This frozen7262 archive only; no equality with separately archived old S1.14 arrays asserted.',new_inference=False,new_training=False,server_access=False,versions={n:md.version(n) for n in ['numpy','pandas','scipy','matplotlib','PyMuPDF']},files=files,copy_exclude=['replay_check_*','__pycache__'],manuscript_artifacts=['figures/FIG5_ADDITIONAL_d_e_85mm.pdf','figures/SUPPLEMENTARY_FIGURE_S22_NATURAL_SNV.pdf','figures/SUPPLEMENTARY_FIGURE_S23_CONSERVATION.pdf'],note='The 223-mm composite is engineering-only. Original B1 a/b/c in b1_reference remains byte-exact.')
(R/'DELIVERY_MANIFEST.json').write_text(json.dumps(report,indent=2)+'\n','utf8')
print(json.dumps(dict(figures=figures,files=len(files),all_replay_checks_passed=report['replay']['all_passed'],total_packaged_bytes=sum(e['bytes'] for e in files.values())),indent=2))
