"""One 178-mm Fig.5, pure display of frozen inputs and results.

Accepted a/b are inserted literally as a vector PDF/SVG crop. Only a's Phi
inset is covered and replaced; all original values remain, now as heat cells.
No inference, statistics, scientific selection or data acquisition occurs.
"""
from pathlib import Path
import os, csv, json, gzip, hashlib, shutil, xml.etree.ElementTree as ET
os.environ['OPENBLAS_NUM_THREADS']='4'
os.environ['OMP_NUM_THREADS']='4'
import numpy as np
import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

R=Path(__file__).resolve().parent
V2=R.parent
O=R/'fig5_v2'; O.mkdir(exist_ok=True)
B=R/'b1_reference'
M=V2/'manuscript/source_data/S23_external_order_preferences/FROZEN_PRIMARY_METRICS.csv'
W,H,PT=178.,183.,72/25.4
FULL,SEQ,FREQ,FLAT,PERM='#E64B35','#4DBBD5','#3C5488','#00A087','#8491B4'
INK,GRAY,PALE='#263238','#7B8790','#C7CED3'
BASE={'A':FREQ,'C':FLAT,'G':'#7E6148','U':FULL,'T':FULL}
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.labelsize':6.5,
 'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.55,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'none'})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='none')

def rows(p):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def t(x,y,z,fs=6.5,**kw):return fig.text(x/W,1-y/H,z,fontsize=fs,va='center',**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def line(x0,y0,x1,y1,**kw):
    obj=Line2D([x0/W,x1/W],[1-y0/H,1-y1/H],transform=fig.transFigure,**kw)
    fig.add_artist(obj);return obj
def rect(x,y,w,h,**kw):
    obj=Rectangle((x/W,1-(y+h)/H),w/W,h/H,transform=fig.transFigure,**kw)
    fig.add_artist(obj);return obj
def clean(a,keep_left=False):
    a.spines[['top','right']].set_visible(False)
    if not keep_left:a.spines['left'].set_visible(False)
    a.tick_params(length=2,width=.5,pad=1.4)
    if not keep_left:a.set_yticks([])
def letter_row(x,y,seq,step=2.7,center=5):
    for j,b in enumerate(seq):t(x+j*step,y,b,fs=7,color=BASE[b],ha='center',fontweight='bold' if j==center else 'normal')

cmap=LinearSegmentedColormap.from_list('raw_phi',[FREQ,'#FFFFFF',FULL])
norm=TwoSlopeNorm(vmin=-7,vcenter=0,vmax=3)
def heat(x,y,w,h,values,central=None):
    vals=np.asarray(values,dtype=float)
    a=ax(x,y,w,h);a.set_axis_off()
    # All cells are native vector quads in PDF and SVG.
    a.pcolormesh(np.arange(len(vals)+1),[0,1],vals[None,:],cmap=cmap,norm=norm,shading='flat',rasterized=False)
    a.set_xlim(0,len(vals));a.set_ylim(0,1)
    if central is not None:
        lo,hi=central;a.add_patch(Rectangle((lo,0),hi-lo,1,fill=False,ec=INK,lw=.65,clip_on=False))
    return a
def scale(x,y,w=22,h=1.15,fs=5.5):
    heat(x,y,w,h,np.linspace(-7,3,120))
    for value in [-7,0,3]:
        # TwoSlopeNorm allocates each sign half of the color scale.
        xx=x+w*norm(value);t(xx,y+h+1.8,f'{value:+d}' if value>0 else str(value),fs=fs,ha='center',color=INK)

# Native artwork retained in output compositing: left top 85 x 96 mm.
# Only the old Phi inset is masked, avoiding every native sequence/response.
rect(35,5,50,15.1,fc='white',ec='none',zorder=0)
packet=np.load(B/'INTERVENTION_PACKET.npz')
qphi=packet['QKI_native_conservation'][0]
t(61,6.7,'Φ fixed reference input · PhyloP',fs=6.1,ha='center',color=INK)
heat(40,9.1,42,2.65,qphi,central=(48,53))
t(40,13.4,'1',fs=5.6,ha='center',color=GRAY)
t(61.0,13.4,'49-53',fs=5.6,ha='center',color=INK)
t(82,13.4,'101 nt',fs=5.6,ha='right',color=GRAY)
t(44,17,'Raw Φ',fs=5.7,color=INK)
scale(59.5,15.9,22,1.1)

# c: shared protein ordering for both primary and secondary endpoints.
primary=rows(M)
proteins=packet['cohort_proteins'].tolist()
metric={(r['protein'],r['arm']):r for r in primary}
assert len(primary)==75 and len(proteins)==25
accuracy=np.array([[100*float(metric[(p,a)]['within_composition_concordance']) for a in ['full','sequence','frequency']] for p in proteins])
assert np.allclose(accuracy,packet['cohort_order_accuracy_percent'],atol=1e-10)
t(2.5,99,'c',fs=9,fontweight='bold');t(8,99,'External preferences · all 25 RBPs',fs=7.1,fontweight='bold')
for x,label,color,marker in [(25,'Full',FULL,'o'),(45,'Sequence',SEQ,'|'),(70,'Frequency',FREQ,'o')]:
    line(x-2.5,104,x-2.5,104,color=color,lw=0,marker=marker,markersize=3 if marker=='o' else 5,markeredgewidth=1)
    t(x,104,label,fs=6.2,color=color)
t(35.5,109,'All-word Spearman ρ',fs=6.3,ha='center')
t(68.5,109,'Within-composition\naccuracy (%)',fs=6.1,ha='center',linespacing=1.1)
cy,ch=113.,63.
sp=ax(22,cy,27,ch);ac=ax(55,cy,28,ch)
for aa in [sp,ac]:
    aa.set_ylim(24.5,-.5);clean(aa)
    for j in range(25):
        if j%2==0:aa.axhspan(j-.5,j+.5,fc='#F3F5F6',ec='none',zorder=0)
sp.set_xlim(-.7,.85);sp.set_xticks([-.5,0,.5]);sp.axvline(0,color=GRAY,lw=.5,ls=(0,(2,2)),zorder=1)
ac.set_xlim(39,73.5);ac.set_xticks([40,50,60,70]);ac.axvline(50,color=GRAY,lw=.5,ls=(0,(2,2)),zorder=1)
for j,p in enumerate(proteins):
    t(20.5,cy+(j+.5)/25*ch,p,fs=6.0,ha='right')
    vv=[float(metric[(p,a)]['spearman_all5mer']) for a in ['full','sequence','frequency']]
    for aa,val in [(sp,vv),(ac,accuracy[j])]:
        aa.plot([min(val),max(val)],[j,j],color='#BEC6CB',lw=.6,zorder=1)
        aa.scatter([val[2]],[j],s=8.2,c=FREQ,zorder=3)
        aa.scatter([val[0]],[j],s=10,c=FULL,zorder=4)
        aa.scatter([val[1]],[j],s=22,c=SEQ,marker='|',lw=1.05,zorder=5)

# d: previously frozen coverage-first natural SNV, all ten RBP events.
ex=rows(R/'FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv')
assert len(ex)==10 and all(r['variant_key']=='chr8:1844440:C:T' for r in ex)
assert [r['RBP'] for r in ex]==sorted(r['RBP'] for r in ex)
events=json.loads(gzip.decompress((R/'source_inputs/EVENTS_INCLUDED.json.gz').read_bytes()).decode('utf8'))
event=next(e for e in events if e['event_id']==ex[0]['event_id'])
ref=event['context']['reference_sequence_rna_direction_dna'].replace('T','U')
alt=event['context']['minor_sequence_dna'].replace('T','U')
phi=np.asarray(event['context']['phi_raw'])
assert ref[50]=='C' and alt[50]=='U' and ref[:50]+ref[51:]==alt[:50]+alt[51:]
t(95.5,4,'d',fs=9,fontweight='bold');t(101,4,'Natural SNV · ARHGEF10 / K562',fs=7.1,fontweight='bold')
t(101,8,'GRCh38 chr8:1,844,440 C>T',fs=6.5)
t(101,14.4,'REF',fs=6.2);t(101,18.8,'ALT',fs=6.2)
sx,step=115.,2.7
rect(sx+5*step-step/2,12.55,step,8.1,fc='#F7E4DF',ec='none',zorder=0)
letter_row(sx,14.4,ref[45:56],step=step);letter_row(sx,18.8,alt[45:56],step=step)
heat(sx-step/2,22.0,step*11,2.65,phi[45:56],central=(3,8))
t(sx-step/2,26.2,'46',fs=5.6,ha='center',color=GRAY)
t(sx+5*step,26.2,'51',fs=5.6,ha='center',color=INK)
t(sx+10.5*step,26.2,'56 nt',fs=5.6,ha='right',color=GRAY)
t(162.5,13,'101-nt reference Φ',fs=5.8,ha='center')
heat(149,15.1,27,2.65,phi,central=(48,53))
t(149,19.45,'1',fs=5.5,ha='center',color=GRAY)
t(162.5,19.45,'49-53',fs=5.5,ha='center')
t(176,19.45,'101',fs=5.5,ha='right',color=GRAY)
t(162.5,23.3,'Center Φ = -3.50',fs=6,ha='center')
scale(152,26.0,23,1.1)
t(101,32.3,'Φ REF = Φ ALT · same fixed input',fs=6.3,color=INK)
t(101,36.7,'RBP',fs=6.5)
t(126.4,35.8,'BEAPR',fs=6.2,ha='center');t(126.4,38,'preferred',fs=5.8,ha='center')
t(143.5,36.7,'Full',fs=6.5,color=FULL);t(161,36.7,'Sequence',fs=6.5,color=SEQ)
dx,dy,dw,dh=136.,40.,39.,45.
da=ax(dx,dy,dw,dh);da.set_xlim(-.04,.04);da.set_ylim(9.5,-.5);clean(da)
da.set_xticks([-.04,0,.04]);da.axvline(0,color=GRAY,lw=.65,ls=(0,(2,2)),zorder=0)
for j,r in enumerate(ex):
    y=dy+(j+.5)/10*dh
    t(101,y,r['RBP'],fs=6.5)
    assert float(r['ref_preferred'])==1
    t(126.4,y,'C',fs=6.5,ha='center',fontweight='bold')
    aa,bb=float(r['A_native_ref_minus_alt']),float(r['S_ref_minus_alt'])
    da.plot([aa,bb],[j,j],c='#BFC7CC',lw=.65,zorder=1)
    da.scatter([aa],[j],s=11.5,c=FULL,zorder=3)
    da.scatter([bb],[j],s=22,c=SEQ,marker='|',lw=1.1,zorder=4)
t(155.5,91.8,'Pred. Δlog count (REF - ALT)',fs=6.4,ha='center')
t(101,93.6,'All 10 RBP events · BEAPR prefers REF (C)',fs=6.1,color='#53616A')

# e: direct quartile relation and complete frozen five-arm cohort summary.
quart=rows(R/'DIRECTION_BY_CONSERVATION.csv')
quart=[r for r in quart if r['cohort']=='strict_test']
summary={r['arm']:r for r in rows(R/'DIRECTION_SUMMARY.csv') if r['cohort']=='strict_test'}
contrast=next(r for r in rows(R/'DIRECTION_CONTRASTS.csv') if r['cohort']=='strict_test' and r['contrast']=='A_native minus S')
two_way={r['arm']:r for r in rows(R/'DIRECTION_TWO_WAY_SENSITIVITY.csv') if r['cohort']=='strict_test'}
t(95.5,99,'e',fs=9,fontweight='bold');t(101,99,'Natural allele preferences',fs=7.1,fontweight='bold')
t(101,103.5,'396 events · 287 variants · 85 RBPs · 119 regions',fs=6.2)
t(101,108.5,'Accuracy by Φ quartile',fs=6.5)
qe=ax(109,114,64,26);qe.set_xlim(-.38,3.38);qe.set_ylim(30,80);clean(qe,keep_left=True)
qe.set_yticks([40,60,80]);qe.set_xticks([0,1,2,3],['Q1','Q2','Q3','Q4'])
qe.axhline(50,color=GRAY,lw=.5,ls=(0,(2,2)),zorder=0)
for arm,color,off,marker in [('A_native',FULL,-.08,'o'),('S',SEQ,.08,'s')]:
    qr=sorted([r for r in quart if r['arm']==arm],key=lambda r:int(r['quartile']))
    means=np.array([float(r['accuracy'])*100 for r in qr])
    low=np.array([float(r['ci_low'])*100 for r in qr]);high=np.array([float(r['ci_high'])*100 for r in qr])
    assert np.min(low)>=30 and np.max(high)<=80, 'quartile CI clipped'
    xx=np.arange(4)+off
    qe.plot(xx,means,color=color,lw=.8,alpha=.9)
    qe.errorbar(xx,means,yerr=[means-low,high-means],fmt=marker,color=color,ms=3.0,lw=.75,capsize=1.6,capthick=.6,zorder=3)
qe.set_ylabel('Accuracy (%)',labelpad=3)
for x,label,color,marker in [(145,'Full',FULL,'o'),(160,'Sequence',SEQ,'s')]:
    line(x-2.3,109,x-2.3,109,color=color,lw=0,marker=marker,markersize=3)
    t(x,109,label,fs=6.1,color=color)
# Fixed, source-defined ranges; no quartile fitting or uncertainty recomputation.
for j in range(4):
    q=next(r for r in quart if int(r['quartile'])==j+1)
    xx=109+(j+.38)/3.76*64
    t(xx,145.3,f"n={q['events']}",fs=5.7,ha='center',color=GRAY)
t(141,149.2,'Increasing center Φ  →',fs=6.1,ha='center')
t(101,154,'Overall direction accuracy (%) · 95% CIs',fs=6.3)
sa=ax(121,158,43,17.5);sa.set_xlim(40,80);sa.set_ylim(4.5,-.5);clean(sa)
sa.set_xticks([40,60,80]);sa.axvline(50,color=GRAY,lw=.5,ls=(0,(2,2)),zorder=0)
arms=['A_native','S','A_flat','A_permmean','ref']
for j,(arm,label,color) in enumerate(zip(arms,['Full','Sequence','Flat Φ','Permuted Φ','REF rule'],[FULL,SEQ,FLAT,PERM,FREQ])):
    rr=summary[arm];ci=two_way[arm];yy=158+(j+.5)/5*17.5
    assert 40<=100*float(ci['ci_low'])<=100*float(ci['ci_high'])<=80, 'aggregate CI clipped'
    t(101,yy,label,fs=6.1,color=INK)
    sa.plot([100*float(ci['ci_low']),100*float(ci['ci_high'])],[j,j],c=color,lw=.85)
    sa.scatter([100*float(rr['accuracy'])],[j],s=12,c=color,zorder=3)
    t(176,yy,f"{100*float(rr['accuracy']):.2f}",fs=6.1,ha='right',color=color)
delta_ci=two_way['A_native minus S']
t(101,181.1,f"Full - Sequence: {100*float(contrast['difference']):.2f} pp [95% CI {100*float(delta_ci['ci_low']):.2f}, {100*float(delta_ci['ci_high']):+.2f}]",fs=6.0,color=INK)
line(89,2,89,181,color='#D9DFE3',lw=.6)

# Overlay PDF/SVG are transparent apart from scoped masks and chart backgrounds.
overlay=O/'fig5_new_panels_overlay'
fig.savefig(overlay.with_suffix('.pdf'),transparent=True)
fig.savefig(overlay.with_suffix('.svg'),transparent=True)
plt.close(fig)

# Literal PDF a/b crop, then scoped overlay. No scaling of accepted artwork.
doc=fitz.open();page=doc.new_page(width=W*PT,height=H*PT)
original=fitz.open(B/'figure.pdf')
page.show_pdf_page(fitz.Rect(0,0,85*PT,96*PT),original,0,clip=fitz.Rect(0,0,85*PT,96*PT))
new=fitz.open(overlay.with_suffix('.pdf'));page.show_pdf_page(page.rect,new,0)
pdf=O/'fig5.pdf';doc.save(pdf,garbage=4,deflate=True)
page.get_pixmap(matrix=fitz.Matrix(600/72,600/72),alpha=False).save(O/'fig5.png')
page.get_pixmap(matrix=fitz.Matrix(150/72,150/72),alpha=False).save(O/'fig5_178mm_preview.png')
doc.close();original.close();new.close()

# Native SVG: insert original artwork with an exact 85 x 96 mm clipping region.
NS='http://www.w3.org/2000/svg';ET.register_namespace('',NS)
root=ET.Element('{'+NS+'}svg',{'width':'178mm','height':'183mm','viewBox':f'0 0 {W*PT} {H*PT}'})
ET.SubElement(root,'{'+NS+'}rect',{'x':'0','y':'0','width':str(W*PT),'height':str(H*PT),'fill':'white'})
defs=ET.SubElement(root,'{'+NS+'}defs');clip=ET.SubElement(defs,'{'+NS+'}clipPath',{'id':'accepted_ab_clip'})
ET.SubElement(clip,'{'+NS+'}rect',{'x':'0','y':'0','width':str(85*PT),'height':str(96*PT)})
def import_svg(path,parent,prefix):
    src=ET.parse(path).getroot()
    ids={el.attrib['id']:prefix+el.attrib['id'] for el in src.iter() if 'id' in el.attrib}
    for el in src.iter():
        if 'id' in el.attrib:el.attrib['id']=ids[el.attrib['id']]
        for attr,val in list(el.attrib.items()):
            if attr.endswith('href') and val.startswith('#'):el.attrib[attr]='#'+ids.get(val[1:],val[1:])
            elif 'url(#' in val:
                for old,newid in ids.items():val=val.replace('url(#'+old+')','url(#'+newid+')')
                el.attrib[attr]=val
    for child in list(src):parent.append(child)
oldgroup=ET.SubElement(root,'{'+NS+'}g',{'clip-path':'url(#accepted_ab_clip)'})
import_svg(B/'figure.svg',oldgroup,'b1_')
newgroup=ET.SubElement(root,'{'+NS+'}g')
import_svg(overlay.with_suffix('.svg'),newgroup,'v2_')
ET.ElementTree(root).write(O/'fig5.svg',encoding='utf-8',xml_declaration=True)

sources=[B/'figure.pdf',B/'figure.svg',B/'INTERVENTION_PACKET.npz',M,R/'FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv',
 R/'DIRECTION_BY_CONSERVATION.csv',R/'DIRECTION_SUMMARY.csv',R/'DIRECTION_CONTRASTS.csv',R/'DIRECTION_TWO_WAY_SENSITIVITY.csv',R/'source_inputs/EVENTS_INCLUDED.json.gz']
mapping={
 'canvas_mm':[W,H],'one_combined_figure':True,'left_right_width_mm':85,'gutter_mm':8,
 'a_b_preserved_crop_mm':[0,0,85,96],'a_b_scale':1,
 'a_only_changed_inset_mask_mm':[35,5,50,15.1],
 'phi_scale':{'raw_PhyloP_range':[-7,3],'center':0,'colors':[FREQ,'#FFFFFF',FULL],'shared_between_a_and_d':True,'no_clipping':True},
 'proteins_shared_order':proteins,'primary_metric':'spearman_all5mer','secondary_metric':'within_composition_concordance',
 'natural_variant':'chr8:1844440:C:T','natural_RBP_events':10,'Phi_REF_equals_ALT':True,
 'accuracy_frozen':{arm:summary[arm] for arm in arms},'Full_minus_Sequence_frozen':contrast,
 'aggregate_and_contrast_two_way_CIs':two_way,
 'interval_methods':{'quartile':'5000 paired RBP bootstrap, seed20261012',
 'aggregate_and_Full_minus_Sequence':'2000 two-way RNA-region/RBP positive-exponential multiplier, seed20261015'},
 'no_new_statistics':True,'source_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
 'exports':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [pdf,O/'fig5.svg',O/'fig5.png']}}
(O/'FIG5_V2_SOURCE_MAPPING.json').write_text(json.dumps(mapping,indent=2)+'\n','utf8')
# Authorized single target copied for parent's one-figure manuscript integration.
shutil.copyfile(pdf,V2/'manuscript/figures/fig5.pdf')
print('Generated native Fig.5, 178 x 183 mm, retained original a/b at unit scale.')
print(O/'fig5_178mm_preview.png')
