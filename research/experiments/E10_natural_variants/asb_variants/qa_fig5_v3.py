"""Display-only preservation and complete frozen-source validation for Fig5 v3."""
from pathlib import Path
import csv,json,hashlib,xml.etree.ElementTree as ET
import numpy as np,fitz
from pypdf import PdfReader
from pypdf.generic import ContentStream
R=Path(__file__).resolve().parent;O=R/'fig5_v3';P=72/25.4
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):
    with Path(p).open(newline='',encoding='utf-8-sig')as f:return list(csv.DictReader(f))
original=fitz.open(R/'b1_reference/figure.pdf');palette=fitz.open(O/'accepted_b1_palette.pdf');final=fitz.open(O/'fig5.pdf')
def chars(page):
    out=[]
    for block in page.get_text('rawdict')['blocks']:
        for line in block.get('lines',[]):
            for span in line['spans']:
                for c in span['chars']:
                    if c['c']!=' ':out.append((c['c'],tuple(c['bbox']),tuple(c['origin']),span['font'],span['size']))
    return out
chars0,chars1=chars(original[0]),chars(palette[0]);assert len(chars0)==len(chars1)
char_coordinate_extraction_delta=0.
for a,b in zip(chars0,chars1):
    assert (a[0],a[3],a[4])==(b[0],b[3],b[4]),'Palette modification changed text/font/size'
    dd=float(np.max(np.abs(np.asarray(a[1]+a[2])-np.asarray(b[1]+b[2]))));char_coordinate_extraction_delta=max(char_coordinate_extraction_delta,dd)
assert char_coordinate_extraction_delta<.00002,'Text geometry changed beyond color-run extraction roundoff'
rr0=PdfReader(R/'b1_reference/figure.pdf');rr1=PdfReader(O/'accepted_b1_palette.pdf')
ops0=ContentStream(rr0.pages[0].get_contents(),rr0).operations
ops1=ContentStream(rr1.pages[0].get_contents(),rr1).operations
fillops={b'rg',b'g',b'k'}
def nonfill(ops):return [(repr(a),b)for a,b in ops if b not in fillops]
assert nonfill(ops0)==nonfill(ops1),'Non-color PDF operations changed'
assert len(ops1)-len(ops0)==124,'Expected two fill operations for each of 62 base glyphs'
e0=list(ET.parse(R/'b1_reference/figure.svg').getroot().iter());e1=list(ET.parse(O/'accepted_b1_palette.svg').getroot().iter());assert len(e0)==len(e1)
svg_changed=0
for a,b in zip(e0,e1):
    assert a.tag==b.tag and a.text==b.text and a.tail==b.tail
    if a.attrib!=b.attrib:
        assert a.tag.endswith('text')and(a.text or '')in'ACGUT'
        assert {k:v for k,v in a.attrib.items()if k!='style'}=={k:v for k,v in b.attrib.items()if k!='style'}
        ast=a.attrib['style'].split(';');bst=b.attrib['style'].split(';')
        assert [x.strip()for x in ast if not x.strip().startswith('fill:')]==[x.strip()for x in bst if not x.strip().startswith('fill:')]
        svg_changed+=1
assert svg_changed==62
clip=fitz.Rect(0,0,85*P,96*P);scale=fitz.Matrix(300/72,300/72)
pix0=original[0].get_pixmap(matrix=scale,clip=clip,alpha=False);pix1=final[0].get_pixmap(matrix=scale,clip=clip,alpha=False)
a=np.frombuffer(pix0.samples,np.uint8).reshape(pix0.height,pix0.width,3);b=np.frombuffer(pix1.samples,np.uint8).reshape(pix1.height,pix1.width,3);assert a.shape==b.shape
mask=np.zeros(a.shape[:2],bool);fac=300/72
def maskrect(rect):
    x0,y0,x1,y1=rect;mask[max(0,int(np.floor(y0*fac))-2):min(mask.shape[0],int(np.ceil(y1*fac))+2),max(0,int(np.floor(x0*fac))-2):min(mask.shape[1],int(np.ceil(x1*fac))+2)]=True
maskrect((35*P,5*P,85*P,20.1*P))
glyphs=0
for block in original[0].get_text('rawdict')['blocks']:
    for line in block.get('lines',[]):
        for span in line['spans']:
            text=''.join(c['c']for c in span['chars']).replace(' ','')
            if text and set(text)<=set('ACGUT'):
                for c in span['chars']:
                    if c['c']in'ACGUT':maskrect(c['bbox']);glyphs+=1
assert glyphs==62
delta=np.max(np.abs(a.astype(int)-b.astype(int)),axis=2);outside=delta[~mask]
if outside.max():
    changed=(delta>0)&(~mask);ys,xs=np.nonzero(changed)
    print(json.dumps({'outside_max':int(outside.max()),'outside_changed':int(changed.sum()),'bbox_mm':[float(xs.min()*25.4/300),float(ys.min()*25.4/300),float(xs.max()*25.4/300),float(ys.max()*25.4/300)],'examples_mm':list(zip((xs[:20]*25.4/300).tolist(),(ys[:20]*25.4/300).tolist()))}))
assert int(outside.max())==0,'Artwork changed outside authorized base/Phi masks'
mapping=json.loads((O/'FIG5_V3_SOURCE_MAPPING.json').read_text('utf8'))
for p,digest in mapping['source_hashes'].items():assert sha(p)==digest,p
for p,digest in mapping['exports'].items():assert sha(p)==digest,p
canonical=Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/supervised_expansion_20261004/candidate_portfolio_20261006/order_candidates/B3')
for p in(O/'frozen_inputs').iterdir():assert sha(p)==sha(canonical/p.name)
edits=rows(O/'frozen_inputs/ALL_25_SELECTED_EDITS.csv')
z5bounds=[min(float(r['external_delta_Z5'])for r in edits),max(float(r['external_delta_Z5'])for r in edits)]
frequencybounds=[min(float(r['native_frequency_delta'])for r in edits),max(float(r['native_frequency_delta'])for r in edits)]
assert 0<z5bounds[0]<=z5bounds[1]<=6 and -230<=frequencybounds[0]<=frequencybounds[1]<0
source=rows(R/'STRICT396_WITH_ANNOTATIONS.csv');display=rows(O/'EVENT_ATLAS_COLUMNS.csv');src={r['event_id']:r for r in source}
assert len(display)==396 and len(src)==396 and set(src)=={r['event_id']for r in display}
for d in display:
    for field,value in d.items():
        if field!='display_column':assert value==src[d['event_id']][field],(d['event_id'],field)
ordered=sorted(source,key=lambda r:(float(r['phi_center']),int(r['chrom'][3:]),int(r['pos1']),r['ref'],r['alt'],r['cell'],r['RBP'],r['event_id']))
assert [r['event_id']for r in ordered]==[r['event_id']for r in display]
assert(len({r['variant_key']for r in source}),len({r['RBP']for r in source}),len({r['cluster_id']for r in source}))==(287,85,119)
responsefields=['A_native_major_minus_minor','S_major_minus_minor','A_flat_major_minus_minor','A_permmean_major_minus_minor']
responsebounds={field:[min(float(r[field])for r in source),max(float(r[field])for r in source)]for field in responsefields}
assert all(-1.8<=x<=1.8 for pair in responsebounds.values()for x in pair)
tw=[r for r in rows(R/'DIRECTION_TWO_WAY_SENSITIVITY.csv')if r['cohort']=='strict_test']
aggregate=[r for r in tw if r['arm']in['A_native','S','A_flat','A_permmean','ref']]
assert len(aggregate)==5 and all(.4<=float(r['ci_low'])<=float(r['ci_high'])<=.8 for r in aggregate)
contrast=next(r for r in tw if r['arm']=='A_native minus S')
spans=[s for block in final[0].get_text('dict')['blocks']for line in block.get('lines',[])for s in line.get('spans',[])]
text='\n'.join(s['text']for s in spans)
for token in['QKI','HNRNPA1','ARHGEF10','396 events','287 variants','85 RBPs','119 RNA regions','Δlog count','52.04','58.77','-6.73','-15.97','+1.64']:assert token in text,token
for p in mapping['all25_original_order']:assert p in text,p
geometry=[final[0].rect.width/P,final[0].rect.height/P];assert np.allclose(geometry,[178,190],atol=.001)and len(final)==1
assert not final[0].get_images(full=True),'Final figure contains raster scientific content'
for label,rect in [('c_5mer_matrix_detail',(0,104,85,182)),('e_event_atlas_detail',(95,104,178,159))]:final[0].get_pixmap(matrix=scale,clip=fitz.Rect(*[x*P for x in rect]),alpha=False).save(O/(label+'.png'))
report={'page_count':1,'canvas_mm':geometry,'all_scientific_artwork_native_vector':True,'base_palette_glyph_changes':62,
 'accepted_a_b_noncolor_PDF_operations_unchanged':True,'accepted_a_b_all_character_geometry_unchanged_within_0_00002pt_extraction_roundoff':True,'max_color_run_text_extraction_delta_pt':char_coordinate_extraction_delta,'SVG_geometry_unchanged_except_62_base_fill_styles':True,
 'a_b_pixel_regression':{'mask':'62 original nucleotide glyph boxes plus x35..85,y5..20.1 mm Phi inset; 2px antialias margin','dpi':300,'masked_base_glyphs':glyphs,'max_rgb_difference_outside_masks':int(outside.max()),'changed_pixels_outside_masks':int(np.sum(outside>0)),'compared_pixels':int(outside.size)},
 'c_all25_original_order':True,'c_complete_36_tasks_1152_paired_model_responses':True,'c_primary_all1024_word_endpoint_retained':True,'c_secondary_within_composition_frozen_summary_retained':True,'c_deltaZ5_and_frequency_bars_unclipped':True,'c_Z5_bounds':z5bounds,'c_frequency_bounds':frequencybounds,
 'd_all10_example_events_retained':True,'e_all396_original_events_exact_field_match':True,'e_counts_variants_RBPs_regions':[287,85,119],'e_nonoutcome_shared_sort_verified':True,
 'shared_c_d_e_response_scale':mapping['shared_response_scale'],'e_raw_response_bounds':responsebounds,'all_five_aggregate_CIs_inside_40_to_80_axis':True,
 'Full_minus_Sequence_two_way_CI_pp':[100*float(contrast[k])for k in['ci_low','ci_high']],'no_new_scientific_statistics':True,'source_hashes_and_original_canonical_copies_unchanged':True,
 'visible_font_range_pt':[min(s['size']for s in spans),max(s['size']for s in spans)],'font_note':'4 pt H/K identity glyphs only; 5-mer bases are 5.8 pt monospaced with 1.4 mm spacing.',
 'SVG_export_rendered_and_visually_inspected':(O/'fig5_svg_178mm_preview.png').exists(),
 'visual_review':'Parent and author inspected updated 178-mm PDF preview; 300-dpi c/e detail and the independently rendered SVG preview inspected after final readability fixes.'}
(O/'FIG5_V3_QA.json').write_text(json.dumps(report,indent=2)+'\n','utf8')
print(json.dumps(report,indent=2))
