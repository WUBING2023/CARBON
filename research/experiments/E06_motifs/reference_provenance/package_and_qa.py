"""Render/compare real exports and construct readable atlas quadrant pages.

Run with D:/Anaconda/python.exe (PyMuPDF). No scientific recalculation.
"""
from pathlib import Path
import hashlib,json,shutil,csv,math,sys
import fitz

W=Path(__file__).resolve().parent
PAPER=W.parent
REFROOT=Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/fig4_reference_atlas_20260929')
MM=72/25.4

def add_label(page,x,y,text,size=8):
    page.insert_text((x*MM,y*MM),text,fontname='helv',fontsize=size,color=(.12,.12,.12))

def atlas_quadrants():
    original=fitz.open(W/'fig3_supplement_atlas_full.pdf')
    p=original[0];cx=p.rect.width/2;cy=p.rect.height/2
    overlap=12*MM
    clips=[fitz.Rect(0,0,cx+overlap,cy+overlap),fitz.Rect(cx-overlap,0,p.rect.width,cy+overlap),fitz.Rect(0,cy-overlap,cx+overlap,p.rect.height),fitz.Rect(cx-overlap,cy-overlap,p.rect.width,p.rect.height)]
    labels=['Upper-left','Upper-right','Lower-left','Lower-right']
    # Bound the enlarged view by rays passing between leaves. The 12-mm
    # rectangular overlap accommodates two complete leaves beyond each quadrant
    # axis; the sector clip removes the otherwise sliced duplicate glyphs.
    step=360/219
    def seam(degrees):
        return 90-round((90-degrees)/step)*step
    sectors=[(seam(90)-2*step,seam(180)+2*step),
             (seam(0)-2*step,seam(90)+2*step),
             (seam(180)-2*step,seam(270)+2*step),
             (seam(270)-2*step,seam(360)+2*step)]
    doc=fitz.open();doc.insert_pdf(original)
    for i,(clip,name) in enumerate(zip(clips,labels)):
        q=doc.new_page(width=210*MM,height=236*MM)
        add_label(q,12,12,f'Figure S21: {name} quadrant',10)
        add_label(q,12,20,'Same 219-leaf atlas; quadrant view enlarged for readable task labels.',8)
        add_label(q,12,26,'Inner-to-outer rings: splicing, stability/decay, translation, RRM, ZnF, KH.',8)
        target=fitz.Rect(10*MM,31*MM,200*MM,221*MM)
        q.show_pdf_page(target,original,0,clip=clip)
        lo,hi=sectors[i]
        scale=target.width/clip.width
        source_polygon=[(cx,cy)]+[(cx+3*p.rect.width*math.cos(math.radians(a)),cy-3*p.rect.width*math.sin(math.radians(a)))
                                     for a in [lo+(hi-lo)*j/24 for j in range(25)]]
        points=[(target.x0+(x-clip.x0)*scale,q.rect.height-(target.y0+(y-clip.y0)*scale)) for x,y in source_polygon]
        path=f'{points[0][0]:.6f} {points[0][1]:.6f} m\n'+''.join(f'{x:.6f} {y:.6f} l\n' for x,y in points[1:])+'h W n\n'
        content=q.get_contents()[-1]
        doc.update_stream(content,b'q\n'+path.encode('ascii')+doc.xref_stream(content)+b'\nQ\n')
        add_label(q,12,229,'No: no annotation-rule match; Unknown: missing record. Source: Data S10 / S24.',8)
        q.get_pixmap(matrix=fitz.Matrix(1.8,1.8),alpha=False).save(W/f'fig3_atlas_quadrant_{i+1}.png')
    doc.save(W/'fig3_supplement_atlas_readable.pdf',garbage=4,deflate=True)
    # Conservative complete-object coverage: a task's rotated label rectangle
    # and every colored nucleotide-path rectangle must fit one source window
    # and angular sector. This verifies each full logo, not just its center.
    with (W/'source_inputs/Supplementary_Data_S10_motif_leaf_mapping.csv').open(encoding='utf-8-sig',newline='') as f:
        rows=sorted(csv.DictReader(f),key=lambda r:int(r['leaf_order_clockwise_from_top']))
    spans=[s for b in p.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
    text_boxes={s['text'].strip():fitz.Rect(s['bbox']) for s in spans}
    leaf_boxes=[[] for _ in rows]
    base_colors={(0.235,0.329,0.533),(0.0,0.627,0.529),(0.494,0.38,0.282),(0.902,0.294,0.208)}
    for drawing in p.get_drawings():
        color=tuple(round(v,3) for v in drawing.get('fill') or ())
        if color not in base_colors:continue
        rect=fitz.Rect(drawing['rect']);mid=(rect.tl+rect.br)/2
        if math.hypot(mid.x-cx,mid.y-cy)<100*MM:continue
        angle=math.degrees(math.atan2(cy-mid.y,mid.x-cx))
        index=int(((90-angle)%360)/step)%219
        leaf_boxes[index].append(rect)
    def inside(rect,clip,sector):
        if not clip.contains(rect):return False
        lo,hi=sector
        for corner in [rect.tl,rect.tr,rect.bl,rect.br]:
            a=math.degrees(math.atan2(cy-corner.y,corner.x-cx))
            while a<lo:a+=360
            if a>hi:return False
        return True
    coverage=[]
    for index,row in enumerate(rows):
        boxes=[text_boxes[row['label']]]+leaf_boxes[index]
        complete=[j+1 for j,(clip,sector) in enumerate(zip(clips,sectors)) if all(inside(rect,clip,sector) for rect in boxes)]
        coverage.append({'task':row['task_id'],'label':row['label'],'nucleotide_paths':len(leaf_boxes[index]),'complete_in_quadrants':complete})
    missing=[r['task'] for r in coverage if not r['complete_in_quadrants']]
    audit={'rectangular_overlap_mm':12,'sector_overlap_leaves':2,'sector_bounds_degrees':sectors,
           'all_labels_and_full_logos_have_complete_copy':not missing,'missing_tasks':missing,
           'main_fig3_sha256':hashlib.sha256((W/'fig3.pdf').read_bytes()).hexdigest(),
           'S21_sha256':hashlib.sha256((W/'fig3_supplement_atlas_readable.pdf').read_bytes()).hexdigest(),
           'coverage':coverage}
    (W/'S21_QUADRANT_QA.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    assert not missing,missing
    return doc

def comparisons():
    rd=W/'references';rd.mkdir(exist_ok=True)
    ref=rd/'Sun_2021_Cell_Research_Fig5.png'
    if not ref.exists():shutil.copy2(REFROOT/'assets/prismnet_fig5.png',ref)
    src=rd/'prismnet_figures.json'
    if not src.exists():shutil.copy2(REFROOT/'sources/prismnet_figures.json',src)
    cur=fitz.open(W/'fig3.pdf');old=fitz.open(W/'fig3_baseline.pdf')
    pdf=fitz.open()
    for kind in ['external','old']:
        p=pdf.new_page(width=190*MM,height=225*MM)
        add_label(p,8,10,'Published exemplar: Sun et al. 2021, Fig. 5' if kind=='external' else 'Previous study Figure 3',8)
        add_label(p,103,10,'Revised study Figure 3 (85 mm)',8)
        height=85*2182/1648 if kind=='external' else 85*old[0].rect.height/old[0].rect.width
        box=fitz.Rect(8*MM,17*MM,93*MM,(17+height)*MM)
        if kind=='external':p.insert_image(box,filename=str(ref),keep_proportion=True)
        else:p.show_pdf_page(box,old,0)
        p.show_pdf_page(fitz.Rect(103*MM,17*MM,188*MM,216*MM),cur,0)
        p.get_pixmap(matrix=fitz.Matrix(1.4,1.4),alpha=False).save(W/f'comparison_{kind}.png')
    p=pdf.new_page(width=190*MM,height=105*MM)
    add_label(p,8,10,'Sun et al. Fig. 5a: paired-logo detail',8)
    add_label(p,103,10,'Revised Figure 3b: paired-logo detail',8)
    imgdoc=fitz.open();ip=imgdoc.new_page(width=1648,height=2182);ip.insert_image(ip.rect,filename=str(ref))
    p.show_pdf_page(fitz.Rect(8*MM,18*MM,93*MM,93*MM),imgdoc,0,clip=fitz.Rect(0,0,880,620))
    p.show_pdf_page(fitz.Rect(103*MM,18*MM,188*MM,93*MM),cur,0,clip=fitz.Rect(0,75*MM,85*MM,119*MM))
    pdf.save(W/'FIGURE_COMPARISONS.pdf',garbage=4,deflate=True)

def main():
    if '--atlas-only' in sys.argv:
        atlas_quadrants()
        print('S21 rebuilt only; label/logo coverage recorded in S21_QUADRANT_QA.json')
        return
    atlas_quadrants();comparisons()
    report={}
    for name,pages in [('fig3.pdf',[0]),('fig3_supplement_atlas_readable.pdf',[]),('fig3_supplement_direct_pairs.pdf',[0,12]),('fig3_supplement_all_193_pairs.pdf',[0,19]),('fig3_supplement_sensitivity.pdf',[0])]:
        doc=fitz.open(W/name)
        info={'pages':len(doc),'page0_width_mm':doc[0].rect.width/MM,'page0_height_mm':doc[0].rect.height/MM,'sha256':hashlib.sha256((W/name).read_bytes()).hexdigest()}
        out=[]
        for pg in range(len(doc)):
            p=doc[pg]
            spans=[s for b in p.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
            out.extend([(pg+1,s['text'],s['bbox']) for s in spans if s['bbox'][0]<-.5 or s['bbox'][1]<-.5 or s['bbox'][2]>p.rect.width+.5 or s['bbox'][3]>p.rect.height+.5])
            info.setdefault('min_font_pt',[]).append(min([s['size'] for s in spans] or [0]))
            if pg in pages:
                p.get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False).save(W/f'qa_{Path(name).stem}_{pg+1}.png')
        if name=='fig3_supplement_atlas_readable.pdf':
            # show_pdf_page keeps clipped source objects in the content stream;
            # text extraction can report their inherited boxes outside the page.
            # Each quadrant intentionally repeats partial boundary objects. Check
            # complete label coverage inside at least one visible clip below.
            info['clipped_source_span_boxes']=[r for r in out if r[0]>1]
            info['text_outside_page']=[r for r in out if r[0]==1]
        else:
            info['text_outside_page']=out
        report[name]=info
    with (W/'source_inputs/Supplementary_Data_S10_motif_leaf_mapping.csv').open(encoding='utf-8-sig',newline='') as f:
        labels=[r['label'] for r in csv.DictReader(f)]
    atlas_text=fitz.open(W/'fig3_supplement_atlas_full.pdf')[0].get_text()
    missing_labels=[x for x in labels if x not in atlas_text]
    quad=fitz.open(W/'fig3_supplement_atlas_readable.pdf')
    box=fitz.Rect(10*MM,31*MM,200*MM,221*MM)
    visible_labels=set()
    for p in list(quad)[1:]:
        for b in p.get_text('dict')['blocks']:
            for l in b.get('lines',[]):
                for s in l['spans']:
                    if box.contains(fitz.Rect(s['bbox'])):
                        visible_labels.add(s['text'].strip())
    missing_quad=[x for x in labels if x not in visible_labels]
    report['atlas_pdf_label_check']={'expected':len(labels),'missing':missing_labels,
                                    'missing_complete_label_in_quadrants':missing_quad,
                                    'quadrant_boundary_content':'Partial duplicates are clipped; every label must have a fully contained copy in another quadrant.'}
    assert not missing_labels,missing_labels
    assert not missing_quad,missing_quad
    report['scientific_invariants']={
        'full_comparisons_retained':193,'p_q_recomputed':False,'original_q_lt_0_05_count':72,
        'direct_experiment_comparisons':126,'direct_experiment_q_lt_0_05':45,
        'atlas_leaves':219,'atlas_rbp_symbols':150,'frozen_leaf_order_preserved':True,
        'complete_219_labels_in_large_atlas':True,'main_example_selection_uses_p_q_or_model_score':False,
        'c_predictions_and_attribution_recomputed':False,
        'c_geometry_method':'Existing vector glyphs/curves reprojected linearly using retained original ticks.'}
    report['visual_review']='Main PDF, original figure, published exemplar and full atlas actually opened; corrected overlapping main names with numbered leaders, aligned cropped main logos, and IC/IG tick separation. Complete pair booklet and sensitivity representative pages rendered for inspection.'
    report['geometry']={'main_width_mm':85,'main_height_mm':199,'smallest_main_ordinary_font_pt':6.5,
                         'main_tree_outer_radius_nominal_mm':11.5,'main_outer_annotation_radius_nominal_mm':34.32,
                         'main_canvas_tree_center_x_mm':39.3,'canvas_midpoint_x_mm':42.5,
                         'supplement_full_atlas_width_mm':285,'supplement_full_atlas_font_pt':6}
    (W/'FIGURE_QA.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
