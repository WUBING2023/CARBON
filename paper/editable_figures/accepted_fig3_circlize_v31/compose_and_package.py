from pathlib import Path
import csv,json,hashlib,shutil,math
import fitz,numpy as np
W=Path(__file__).resolve().parent;R=W.parent;P=R.parent;MM=72/25.4
OLD=P/'revision_v30/manuscript/figures/fig3.pdf';CIRCLE=W/'circos_carbon_285mm.pdf'
D=R/'delivery';D.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rect=lambda v:fitz.Rect(*[x*MM for x in v])
old=fitz.open(OLD);circle=fitz.open(CIRCLE)
with (W/'circos_carbon_285mm_geometry.csv').open(encoding='utf-8-sig') as f:geometry=list(csv.DictReader(f))
assert len(geometry)==219
doc=fitz.open();page=doc.new_page(width=85*MM,height=98.5*MM)
page.show_pdf_page(rect([.5,.5,46.5,46.5]),circle,0)
# The existing native panels and legend are placed at exactly the same coordinates.
for box in [[47.5,0,85,60],[0,46.5,47.5,60],[0,60,85,98.5],[0,0,4.5,4.1]]:
    page.show_pdf_page(rect(box),old,0,clip=rect(box))
doc.save(W/'fig3_v31.pdf',garbage=4,deflate=True)
doc.close();doc=fitz.open(W/'fig3_v31.pdf');page=doc[0]
checks={}
for name,box in [('b',[48,0,85,60]),('c',[0,61,85,98.5])]:
    a=old[0].get_pixmap(matrix=fitz.Matrix(300/72,300/72),clip=rect(box),alpha=False)
    b=page.get_pixmap(matrix=fitz.Matrix(300/72,300/72),clip=rect(box),alpha=False)
    assert a.width==b.width and a.height==b.height
    av=np.frombuffer(a.samples,dtype='uint8');bv=np.frombuffer(b.samples,dtype='uint8')
    assert np.array_equal(av,bv),name
    checks[name]=dict(pixel_identical=True,dpi=300,position_mm=box)
page.get_pixmap(matrix=fitz.Matrix(3,3),alpha=False).save(D/'CARBON_Fig3_v31.png')
shutil.copy2(W/'fig3_v31.pdf',D/'CARBON_Fig3_v31.pdf')

atlas=fitz.open();ap=atlas.new_page(width=285*MM,height=306*MM)
ap.show_pdf_page(rect([0,0,285,285]),circle,0)
ap.insert_text((2*MM,5*MM),'a',fontsize=11,fontname='hebo')
keys=[('Splicing','#E6211A'),('Stability/decay','#19499C'),('Translation','#6BB72D'),('Stability + translation','#F4A016')]
def rgb(h):return tuple(int(h[i:i+2],16)/255 for i in [1,3,5])
for x,(label,col) in zip([8,68,136,200],keys):
    ap.draw_line((x*MM,290*MM),((x+5)*MM,290*MM),color=rgb(col),width=1.3)
    ap.insert_text(((x+7)*MM,291*MM),label,fontsize=8.5,color=(.15,.18,.22))
for x,label,col in [(8,'RRM','#52C9EC'),(45,'ZnF','#E5D74A'),(77,'KH','#8BC01C'),(107,'No annotation-rule match','#EEEEEE'),(192,'Unknown','#FFFFFF')]:
    ap.draw_rect(rect([x,296,x+2.6,298.6]),fill=rgb(col),color=(.7,.7,.7),width=.3)
    ap.insert_text(((x+4)*MM,298.5*MM),label,fontsize=8.5)
ap.insert_text((8*MM,304*MM),'Pastel sectors: six cuts of the retained tree for visual orientation; not biological classes.',fontsize=7,color=(.4,.45,.5))
atlas.save(D/'CARBON_Fig3a_enlarged_285mm.pdf',garbage=4,deflate=True)
atlas[0].get_pixmap(matrix=fitz.Matrix(.9,.9)).save(D/'CARBON_Fig3a_enlarged_preview.png')

# S21 keeps its five-page overview + overlapping-quadrant structure.
s21=fitz.open();p=s21.new_page(width=178*MM,height=191.1158*MM)
p.show_pdf_page(p.rect,atlas,0)
centre=142.5;overlap=20
boxes=[[0,0,centre+overlap,centre+overlap],[centre-overlap,0,285,centre+overlap],
       [0,centre-overlap,centre+overlap,285],[centre-overlap,centre-overlap,285,285]]
labels={}
for b in circle[0].get_text('dict')['blocks']:
 for line in b.get('lines',[]):
  for span in line['spans']:labels.setdefault(span['text'],[]).append(fitz.Rect(span['bbox']))
units=[]
for row in geometry:
    # Circle's exact user-coordinate transform: canvas -1.025..1.025.
    t=float(row['theta'])*math.pi/180;da=math.pi/219
    points=[(142.5+radius*285/2.05*math.cos(a),142.5-radius*285/2.05*math.sin(a)) for radius in [.49,1.005] for a in [t-da,t,t+da]]
    box=rect([min(x for x,y in points),min(y for x,y in points),max(x for x,y in points),max(y for x,y in points)])
    assert row['label'] in labels,row['label']
    for lb in labels[row['label']]:box|=lb
    complete=[i+1 for i,b in enumerate(boxes) if rect(b).contains(box)]
    assert complete,row['task_id']
    units.append(dict(task=row['task_id'],complete_in_quadrants=complete))
for i,box in enumerate(boxes,1):
    p=s21.new_page(width=210*MM,height=236*MM)
    p.insert_text((8*MM,9*MM),f'Figure S21 | quadrant {i}',fontsize=9,fontname='hebo')
    p.show_pdf_page(rect([5,14,205,214]),circle,0,clip=rect(box))
    p.insert_text((7*MM,223*MM),'Overlapping views retain every full name, motif and annotation in at least one quadrant.',fontsize=7)
    p.insert_text((7*MM,228*MM),'Tree topology, 219-task population and original annotation values are unchanged.',fontsize=7)
s21.save(W/'figS21_motif_atlas_v31.pdf',garbage=4,deflate=True)

comparison=fitz.open();p=comparison.new_page(width=370*MM,height=210*MM)
ref=fitz.open(R/'reference_bundle/circos_plot.pdf')
p.insert_text((8*MM,10*MM),'Provided reference',fontsize=11,fontname='hebo')
p.insert_text((192*MM,10*MM),'CARBON | adapted code, frozen study data',fontsize=11,fontname='hebo')
p.show_pdf_page(rect([5,16,180,191]),ref,0);p.show_pdf_page(rect([190,16,365,191]),circle,0)
p.insert_text((8*MM,200*MM),'Reference data stay separate. CARBON adds no paired/unpaired structure channels.',fontsize=8)
p=comparison.new_page(width=190*MM,height=220*MM)
p.insert_text((8*MM,10*MM),'Figure 3 | existing b/c preserved',fontsize=11,fontname='hebo')
p.show_pdf_page(rect([8,17,178,214]),doc,0)
comparison.save(D/'CARBON_Fig3_reference_comparison.pdf',garbage=4,deflate=True)
info=dict(main_mm=[85,98.5],a_target_mm=[.5,.5,46.5,46.5],old_source_sha256=sha(OLD),
    main_pdf_sha256=sha(W/'fig3_v31.pdf'),B_C=checks,names=219,raster_images=len(page.get_images()),
    all_219_names_in_R_export=True,S21_pages=len(s21),all_219_complete_in_quadrants=True,quadrant_units=units,
    all_original_data_preserved=True,reference_data_imported=False)
(W/'COMPOSITION_QA.json').write_text(json.dumps(info,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in info.items() if k!='quadrant_units'}))
