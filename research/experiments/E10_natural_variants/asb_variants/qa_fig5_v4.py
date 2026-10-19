"""Source, rendered-vector and uniform-transformation QA for compact Fig.5."""
from pathlib import Path
import csv,json,hashlib,re,xml.etree.ElementTree as ET
import numpy as np,fitz
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm,AsinhNorm,ListedColormap,Normalize,to_hex
R=Path(__file__).resolve().parent;O=R/'fig5_v4';OLD=R/'fig5_v3';P=72/25.4
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):
    with Path(p).open(newline='',encoding='utf-8-sig')as f:return list(csv.DictReader(f))
m=json.loads((O/'FIG5_V4_SOURCE_MAPPING.json').read_text('utf8'))
for p,digest in m['source_hashes'].items():assert sha(p)==digest,p
for p,digest in m['exports'].items():assert sha(p)==digest,p
final=fitz.open(O/'fig5.pdf');orig=fitz.open(OLD/'fig5.pdf')
assert len(final)==1 and np.allclose([final[0].rect.width/P,final[0].rect.height/P],[178,142],atol=.001)
assert not final[0].get_images(full=True)

# Independently compose only the accepted source under the declared transforms.
expected=fitz.open();epg=expected.new_page(width=178*P,height=142*P)
for pl in m['a_b_placements']:
    dest,clip=[fitz.Rect(*[v*P for v in pl[k]])for k in['destination_mm','source_clip_mm']]
    assert abs(dest.width/clip.width-m['a_b_uniform_scale'])<1e-7
    assert abs(dest.height/clip.height-m['a_b_uniform_scale'])<1e-7
    epg.show_pdf_page(dest,orig,0,clip=clip)
dpi=300;mat=fitz.Matrix(dpi/72,dpi/72)
p0=epg.get_pixmap(matrix=mat,alpha=False);p1=final[0].get_pixmap(matrix=mat,alpha=False)
a=np.frombuffer(p0.samples,np.uint8).reshape(p0.height,p0.width,3)
b=np.frombuffer(p1.samples,np.uint8).reshape(p1.height,p1.width,3)
delta=np.max(np.abs(a.astype(int)-b.astype(int)),axis=2);ab_reports={}
for pl in m['a_b_placements']:
    x0,y0,x1,y1=pl['destination_mm']
    valid=np.zeros(delta.shape,bool)
    valid[int(np.ceil(y0*dpi/25.4))+2:int(np.floor(y1*dpi/25.4))-2,
          int(np.ceil(x0*dpi/25.4))+2:int(np.floor(x1*dpi/25.4))-2]=True
    for xx,yy,ww,hh in m['a_b_only_overlay_masks_mm']:
        valid[max(0,int(np.floor(yy*dpi/25.4))-2):int(np.ceil((yy+hh)*dpi/25.4))+2,
              max(0,int(np.floor(xx*dpi/25.4))-2):int(np.ceil((xx+ww)*dpi/25.4))+2]=False
    dd=delta[valid];assert not np.any(dd),pl['panel']
    ab_reports[pl['panel']]={'dpi':dpi,'compared_pixels':int(dd.size),
        'max_rgb_difference':int(dd.max()),'changed_pixels':int(np.sum(dd>0)),
        'scale':m['a_b_uniform_scale'],'source_crop_mm':pl['source_clip_mm'],
        'destination_mm':pl['destination_mm'],'exception':'panel identifier glyph only; 2px antialias margin'}

# Exactly the same frozen display columns as v3, not only the same counts.
display=rows(O/'EVENT_ATLAS_COLUMNS.csv');source=rows(R/'STRICT396_WITH_ANNOTATIONS.csv')
src={r['event_id']:r for r in source};assert len(display)==len(src)==396
assert sha(O/'EVENT_ATLAS_COLUMNS.csv')==sha(OLD/'EVENT_ATLAS_COLUMNS.csv')
for d in display:
    for k,v in d.items():
        if k!='display_column':assert v==src[d['event_id']][k],(d['event_id'],k)
assert(len({r['variant_key']for r in source}),len({r['RBP']for r in source}),len({r['cluster_id']for r in source}))==(287,85,119)
ordered=sorted(source,key=lambda r:(float(r['phi_center']),int(r['chrom'][3:]),int(r['pos1']),r['ref'],r['alt'],r['cell'],r['RBP'],r['event_id']))
assert [r['event_id']for r in display]==[r['event_id']for r in ordered]
arms=['A_native','S','A_flat','A_permmean']
for r in source:
    for arm in arms:
        expected_sign=1 if float(r['ref_preferred'])==1 else -1
        assert np.isclose(float(r[arm+'_major_minus_minor']),expected_sign*float(r[arm+'_ref_minus_alt']),atol=1e-12)
example=rows(R/'FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv')
assert len(example)==10 and [r['RBP']for r in example]==sorted(r['RBP']for r in example)
assert all(r['variant_key']=='chr8:1844440:C:T'and float(r['ref_preferred'])==1 for r in example)
exdisplay=[r for r in display if r['variant_key']=='chr8:1844440:C:T']
assert {r['event_id']for r in example}=={r['event_id']for r in exdisplay}
assert [int(r['display_column'])for r in exdisplay]==list(range(*[m['c_example_e_columns_1based'][0],m['c_example_e_columns_1based'][1]+1]))

# Verify actual exported overlay QuadMesh cell colors against every source value.
BASE=json.loads((R.parent/'coordination/PALETTE_V3.json').read_text('utf8'))['nucleotides']
cm=LinearSegmentedColormap.from_list('signed',['#3C5488','#FFFFFF','#E64B35']);cm.set_bad('#D8DFE5')
phi=TwoSlopeNorm(vmin=-7,vcenter=0,vmax=10);response=AsinhNorm(linear_width=.02,vmin=-1.8,vmax=1.8);rho=Normalize(-1,1)
import gzip
events=json.loads(gzip.decompress((R/'source_inputs/EVENTS_INCLUDED.json.gz').read_bytes()).decode('utf8'))
ev=next(e for e in events if e['event_id']==example[0]['event_id'])
assert ev['context']['reference_sequence_rna_direction_dna'][50]=='C'
assert ev['context']['minor_sequence_dna'][50]=='T'
edits=rows(OLD/'frozen_inputs/ALL_25_SELECTED_EDITS.csv')
task=rows(OLD/'frozen_inputs/ALL_36_TASKS.csv');responses=rows(OLD/'frozen_inputs/ALL_1152_PAIRED_MODEL_RESPONSES.csv')
assert len(edits)==25 and len(task)==36 and len(responses)==1152
npz=np.load(OLD/'frozen_inputs/SOURCE_PACKET.npz')
coeff={(r['protein'],r['arm']):float(r['spearman_all5mer'])for r in rows(R.parent/'manuscript/source_data/S23_external_order_preferences/FROZEN_PRIMARY_METRICS.csv')}
fm=np.full((50,16),np.nan);sm=fm.copy();rh=[]
for i,r in enumerate(edits):
    p=r['protein'];assert p==m['all25_original_order'][i]
    assert len(r['lower_assay_word'])==len(r['higher_assay_word'])==5
    assert sorted(r['lower_assay_word'])==sorted(r['higher_assay_word'])
    assert 0<float(r['external_delta_Z5'])<=6 and -230<=float(r['native_frequency_delta'])<0
    rh.append([coeff[(p,arm)]for arm in['full','sequence','frequency']])
    for j,cell in enumerate(['HepG2','K562']):
        inds=np.flatnonzero(npz[p+'_cells']==cell)
        if len(inds):assert len(inds)==16;fm[2*i+j]=npz[p+'_full_paired_delta'][inds];sm[2*i+j]=npz[p+'_sequence_paired_delta'][inds]
assert np.isfinite(fm).sum()+np.isfinite(sm).sum()==1152
regions=['intronic','UTR','coding','splice','other_annotated','unavailable']
region_cm=ListedColormap(['#8491B4','#91D1C2','#7E6148','#F39B7F','#B09C85','#DFE3E7'])
nm={b:i for i,b in enumerate('ACGT')};ncm=ListedColormap([BASE[b]for b in'ACGT'])
checks=[
 ('c_phi',np.asarray(ev['context']['phi_raw'])[45:56],phi,cm),
 ('c_10x4',[[float(r[arm+'_ref_minus_alt'])for arm in arms]for r in example],response,cm),
 ('d_all25x3_primary_rho',rh,rho,cm),('d_Full_50x16',fm,response,cm),('d_Sequence_50x16',sm,response,cm),
 ('d_rho_colorbar',np.linspace(-1,1,80),rho,cm),('d_response_colorbar',response.inverse(np.linspace(0,1,160)),response,cm),
 ('e_phi',[float(r['phi_center'])for r in ordered],phi,cm),
 ('e_REF',[nm[r['ref']]for r in ordered],Normalize(-.5,3.5),ncm),
 ('e_ALT',[nm[r['alt']]for r in ordered],Normalize(-.5,3.5),ncm),
 ('e_BEAPR',[float(r['ref_preferred'])for r in ordered],Normalize(0,1),ListedColormap(['#F4A016','#19499C'])),
 ('e_region',[regions.index(r['NCBI_annotation_region'])for r in ordered],Normalize(-.5,5.5),region_cm)]
checks += [('e_'+arm,[float(r[arm+'_major_minus_minor'])for r in ordered],response,cm)for arm in arms]
checks += [('e_response_colorbar',response.inverse(np.linspace(0,1,160)),response,cm)]
root=ET.parse(O/'fig5_overlay.svg').getroot()
quad={int(el.attrib['id'].split('_')[-1]):el for el in root.iter()if el.attrib.get('id','').startswith('QuadMesh_')}
assert len(quad)==len(checks)==17
meshqa=[]
for idx,(name,data,norm,colormap)in enumerate(checks,1):
    paths=[el for el in quad[idx].iter()if el.tag.endswith('path')]
    painted=[]
    for path in paths:
        match=re.search(r'fill:\s*([^;]+)',path.attrib.get('style',''))
        painted.append(match.group(1).lower()if match else'#000000')
    values=np.ma.masked_invalid(np.asarray(data,dtype=float)).ravel()
    expected_colors=[to_hex(color).lower()for color in colormap(norm(values))]
    assert len(painted)==len(expected_colors)and painted==expected_colors,name
    meshqa.append({'mesh':name,'cells':len(painted),'source_to_export_exact_color_match':True})
assert meshqa[6]['cells']==meshqa[-1]['cells']==160
assert checks[6][2]is checks[-1][2]is checks[1][2]

tw={r['arm']:r for r in rows(R/'DIRECTION_TWO_WAY_SENSITIVITY.csv')if r['cohort']=='strict_test'}
summ={r['arm']:r for r in rows(R/'DIRECTION_SUMMARY.csv')if r['cohort']=='strict_test'}
aggregate={arm:{'accuracy_pct':100*float(summ[arm]['accuracy']),
 'ci_pct':[100*float(tw[arm][k])for k in['ci_low','ci_high']]}for arm in arms+['ref']}
assert all(40<=v['ci_pct'][0]<=v['ci_pct'][1]<=80 for v in aggregate.values())
contrast=[100*float(tw['A_native minus S'][k])for k in['ci_low','ci_high']]
assert np.allclose(contrast,[-15.970237369747828,1.6411083571633116])
for name,clip in[('FIG5_V4_TOP_178mm_300dpi',(0,0,178,50)),('FIG5_V4_D_DETAIL_300dpi',(0,50,87,142)),('FIG5_V4_E_DETAIL_300dpi',(95,50,178,142))]:
    final[0].get_pixmap(matrix=mat,clip=fitz.Rect(*[v*P for v in clip]),alpha=False).save(O/(name+'.png'))
report={'canvas_mm':[178,142],'page_count':1,'all_scientific_artwork_native_vector':True,
 'old_to_new_panels':m['old_to_new_panels'],'a_b_uniform_transform_pixel_regression':ab_reports,
 'a_b_original_data_geometry_ticks_and_response_points_preserved_under_uniform_scale':True,
 'a_b_original_axis_font_pt_6_5_becomes_pt':6.5*m['a_b_uniform_scale'],
 'a_original_fixed_Phi_inset_font_pt_5_2_becomes_pt':5.2*m['a_b_uniform_scale'],
 'd_H_K_identity_font_pt':4.,'d_5mer_letters_font_pt':5.8,'d_5mer_letter_spacing_mm':1.4,
 'c_all10_actual_events_retained':True,'c_response_direction':'C reference minus U alternate; positive favors REF',
 'e_response_direction':'BEAPR-observed preferred minus other; positive agrees with BEAPR; REF can be either allele',
 'all396_event_fields_exactly_match_frozen_source':True,'all396_display_csv_byte_identical_to_v3':True,
 'e_non_output_sort_unchanged':True,'counts_25_RBPs_36_tasks_1152_responses_396_events_287_variants_85_RBPs_119_regions_verified':True,
 'all17_native_heatmaps_5690_cells_source_to_export_exact':meshqa,
 'shared_c_d_e_response_mapping':m['shared_c_d_e_response_scale'],
 'common_response_colorbars_exactly_identical':True,
 'd_all1024_primary_and_within_composition_secondary_unchanged':True,
 'd_all25_edits_and_negative_response_cells_retained':True,'all5_accuracy_arms_unclipped':aggregate,
 'e_quartile_counts':m['e_quartile_counts'],'Full_minus_Sequence_two_way_CI_pp':contrast,
 'source_and_export_sha256_all_verified':True,'no_new_scientific_statistics_or_inference':True,
 'final_pdf_sha256':sha(O/'fig5.pdf'),'final_svg_sha256':sha(O/'fig5.svg'),
 'SVG_independent_Cairo_render_exists':(O/'fig5_svg_178mm_preview.png').exists(),
 'visual_review':'Author inspected full 178-mm PDF and independently rendered SVG preview; parent inspected early 178-mm PNG and accepted this compact arrangement. Detail inspection is recorded in the caption/mapping note.'}
(O/'FIG5_V4_QA.json').write_text(json.dumps(report,indent=2)+'\n','utf8')
print(json.dumps({'passed':True,'heatmap_cells':sum(v['cells']for v in meshqa),'a_b_regression':ab_reports,'final_pdf_sha256':report['final_pdf_sha256']},indent=2))
