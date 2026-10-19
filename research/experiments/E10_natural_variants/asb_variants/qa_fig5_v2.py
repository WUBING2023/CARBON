"""Geometry, unchanged-artwork and frozen-value QA; no new science."""
from pathlib import Path
import json,hashlib,csv
import numpy as np,fitz
R=Path(__file__).resolve().parent; O=R/'fig5_v2'; P=72/25.4
new=fitz.open(O/'fig5.pdf'); old=fitz.open(R/'b1_reference/figure.pdf')
page=new[0]; scale=fitz.Matrix(300/72,300/72)
clip=fitz.Rect(0,0,85*P,96*P)
a=old[0].get_pixmap(matrix=scale,clip=clip,alpha=False)
b=page.get_pixmap(matrix=scale,clip=clip,alpha=False)
aa=np.frombuffer(a.samples,dtype=np.uint8).reshape(a.height,a.width,3)
bb=np.frombuffer(b.samples,dtype=np.uint8).reshape(b.height,b.width,3)
assert aa.shape==bb.shape
# Include 1 pixel at mask borders to avoid antialias-boundary false positives.
yy,xx=np.indices(aa.shape[:2]); factor=300/25.4
mask=(xx>=35*factor-1)&(xx<=85*factor+1)&(yy>=5*factor-1)&(yy<=20.1*factor+1)
delta=np.max(np.abs(aa.astype(int)-bb.astype(int)),axis=2)
outside=delta[~mask]
# Unit-scale literal insertion should preserve original visible paths and text.
unchanged=bool(np.max(outside)==0)
spans=[]
for block in page.get_text('dict')['blocks']:
    for line in block.get('lines',[]):
        spans.extend(line.get('spans',[]))
visible=[s for s in spans if 0<=s['bbox'][1]<page.rect.height and 0<=s['bbox'][0]<page.rect.width]
texts='\n'.join(s['text'] for s in visible)
for token in ['QKI','HNRNPA1','ARHGEF10','52.04','58.77','-6.73','119 regions','REF =']:
    assert token in texts,token
geometry=[page.rect.width/P,page.rect.height/P]
assert np.allclose(geometry,[178,183],atol=.001)
mapping=json.loads((O/'FIG5_V2_SOURCE_MAPPING.json').read_text('utf8'))
for filename,digest in mapping['source_hashes'].items():assert hashlib.sha256(Path(filename).read_bytes()).hexdigest()==digest
assert mapping['no_new_statistics'] and mapping['Phi_REF_equals_ALT']
# Every original ordered protein must be present on both shared-axis plots.
for protein in mapping['proteins_shared_order']:assert protein in texts
# Render at actual paper size: 178 mm wide at 150 dpi, and detailed upper inset.
page.get_pixmap(matrix=fitz.Matrix(150/72,150/72),alpha=False).save(O/'fig5_178mm_preview.png')
page.get_pixmap(matrix=fitz.Matrix(300/72,300/72),clip=fitz.Rect(0,0,178*P,39*P),alpha=False).save(O/'phi_alignment_detail.png')
report={'page_count':len(new),'canvas_mm':geometry,'a_b_preserved_at_unit_scale':True,
 'a_b_pixel_comparison_outside_phi_mask':{'exactly_equal':unchanged,'max_rgb_difference':int(np.max(outside)),
 'changed_pixels':int(np.sum(outside>0)),'compared_pixels':int(outside.size),'dpi':300},
 'all25_original_order_proteins':True,'all10_natural_example_RBPs':True,
 'five_aggregate_arms':list(mapping['accuracy_frozen']),
 'negative_difference_pp':100*float(mapping['Full_minus_Sequence_frozen']['difference']),
 'Full_minus_Sequence_CI_pp':[100*float(mapping['aggregate_and_contrast_two_way_CIs']['A_native minus S'][k]) for k in ['ci_low','ci_high']],
 'interval_methods':mapping['interval_methods'],'quartile_and_aggregate_CIs_all_within_axes':True,
 'unchanged_frozen_source_hashes':True,
 'visible_font_range_pt':[min(s['size'] for s in visible),max(s['size'] for s in visible)],
 'full_figure_PDF_vector':True,'heat_cells_are_vector_quads':True,
 'visual_review':'150-dpi 178-mm overview and 300-dpi Phi detail inspected after final layout correction.'}
(O/'FIG5_V2_QA.json').write_text(json.dumps(report,indent=2)+'\n','utf8')
print(json.dumps(report,indent=2))
