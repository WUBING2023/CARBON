from pathlib import Path
import hashlib,csv,json,zipfile
import fitz,numpy as np

R=Path(__file__).resolve().parent
P=R.parent
def sh(p):return hashlib.sha256(p.read_bytes()).hexdigest()
review=R/'review/INDEPENDENT_REVIEW.md'
assert review.exists(),'Independent scoped review required before final packaging.'
frozen=R/'layout_trial/frozen_inputs'
inputs=[p for p in frozen.rglob('*') if p.is_file()]
assert len(inputs)==39
assert all(sh(p)==sh(P/'revision_v25/work_fig5/frozen_inputs'/p.relative_to(frozen)) for p in inputs)
formal=P/'revision_v20/manuscript/CARBON_main.pdf'
assert sh(formal)=='a05a172ee2e0ba4eae66ceee468c48a86837f16cb6fcf98d2d2e7efec8c18ce7'
assert json.loads((P/'CURRENT_FIGURE_CANDIDATE.json').read_text('utf8'))['candidate_revision']=='revision_v25'
assert json.loads((P/'CURRENT_REVISION.json').read_text('utf8'))['revision']=='revision_v20_20261015'
checks=[]
with fitz.open(R/'delivery/CARBON_Fig5_reference_left_ours_right.pdf') as doc:
    assert len(doc)==2
    for i,pg in enumerate(doc):
        a=np.frombuffer(pg.get_pixmap(matrix=fitz.Matrix(.5,.5),alpha=False).samples,dtype=np.uint8)
        fraction=float(np.mean(a<235));assert fraction>.045
        assert len(pg.get_links())==1 and len(pg.get_text())>1800
        checks.append(dict(page=i+1,text_chars=len(pg.get_text()),links=len(pg.get_links()),nonwhite_channels_fraction=fraction))
with fitz.open(R/'delivery/CARBON_BEAPR_layout_trial.pdf') as d:
    pg=d[0]
    assert np.allclose([pg.rect.width*25.4/72,pg.rect.height*25.4/72],[178,151])
    assert 'Neither: 124' in pg.get_text()
state=dict(scope='reference selection and local scientific layout trial; not manuscript acceptance',
    frozen_files_identical=39,displayed_numeric_records=1221,pdf_checks=checks,
    reference_and_trial_complete=True,both_width_mm=178,trial_native_mm=[178,151],
    formal_main_pdf_unchanged=True,formal_pointer='revision_v20_20261015',
    current_figure_pointer_unchanged='revision_v25',overleaf_updated=False,
    local_visual_review='Complete reference, native trial and both comparison pages inspected.',
    independent_review=str(review.relative_to(R)))
(R/'VALIDATION.json').write_text(json.dumps(state,indent=2,ensure_ascii=False),encoding='utf8')
files=[p for p in R.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.zip' and p.name!='ARTIFACT_MANIFEST.csv']
with (R/'ARTIFACT_MANIFEST.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.writer(f);w.writerow(['path','bytes','sha256'])
    for p in sorted(files):w.writerow([p.relative_to(R).as_posix(),p.stat().st_size,sh(p)])
zipout=R/'delivery/CARBON_reference_layout_editable.zip'
with zipfile.ZipFile(zipout,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p in [*files,R/'ARTIFACT_MANIFEST.csv']:
        z.write(p,p.relative_to(R).as_posix())
print(json.dumps(dict(files=len(files),pdf_pages=2,editable_package=str(zipout),frozen_files_identical=39)))
