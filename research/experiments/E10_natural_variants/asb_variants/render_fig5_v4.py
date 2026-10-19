"""Compact Fig.5 v4: frozen evidence, three illustrative panels over two atlases."""
from pathlib import Path
import os,csv,json,gzip,hashlib,xml.etree.ElementTree as ET
os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['OMP_NUM_THREADS']='4'
import numpy as np,fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm,AsinhNorm,ListedColormap,Normalize
R=Path(__file__).resolve().parent;V=R.parent;O=R/'fig5_v4';O.mkdir(exist_ok=True)
OLD=R/'fig5_v3';INPUT=OLD/'frozen_inputs'
S=V/'manuscript/source_data/S23_external_order_preferences';B=R/'b1_reference'
W,H,PT=178.,142.,72/25.4
FULL,SEQ,FREQ,FLAT,PERM='#E64B35','#4DBBD5','#3C5488','#00A087','#8491B4'
INK,GRAY,PALE='#263238','#74808B','#D8DFE5'
BASE=json.loads((V/'coordination/PALETTE_V3.json').read_text('utf8'))['nucleotides']
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.labelsize':6.5,
 'xtick.labelsize':6.0,'ytick.labelsize':6.0,'axes.linewidth':.55,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'none'})
def rows(p):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='none')
def t(x,y,z,fs=6.5,**kw):return fig.text(x/W,1-y/H,z,fontsize=fs,va='center',**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def line(x0,y0,x1,y1,**kw):
    ob=Line2D([x0/W,x1/W],[1-y0/H,1-y1/H],transform=fig.transFigure,**kw);fig.add_artist(ob);return ob
def rect(x,y,w,h,**kw):
    ob=Rectangle((x/W,1-(y+h)/H),w/W,h/H,transform=fig.transFigure,**kw);fig.add_artist(ob);return ob
def clean(a):
    a.spines[['top','right','left']].set_visible(False);a.set_yticks([]);a.tick_params(length=2,width=.5,pad=1.2)
def letter_row(x,y,seq,step=1.4,center=-1,fs=5.8,mono=True):
    for j,b in enumerate(seq):t(x+j*step,y,b,fs=fs,color=BASE[b],ha='center',
        fontweight='bold'if j==center else'normal',fontfamily='DejaVu Sans Mono'if mono else'Arial')
cmap=LinearSegmentedColormap.from_list('signed',[FREQ,'#FFFFFF',FULL]);cmap.set_bad(PALE)
phinorm=TwoSlopeNorm(vmin=-7,vcenter=0,vmax=10)
resp_norm=AsinhNorm(linear_width=.02,vmin=-1.8,vmax=1.8);rhonorm=Normalize(-1,1)
def matrix(x,y,w,h,values,norm=resp_norm,colors=cmap,outline=False):
    data=np.ma.masked_invalid(np.asarray(values,dtype=float))
    if data.ndim==1:data=data[None,:]
    a=ax(x,y,w,h);a.set_axis_off();nr,nc=data.shape
    a.pcolormesh(np.arange(nc+1),np.arange(nr+1),data,cmap=colors,norm=norm,shading='flat',rasterized=False)
    a.set_xlim(0,nc);a.set_ylim(nr,0)
    if outline:a.add_patch(Rectangle((0,0),nc,nr,fill=False,ec=PALE,lw=.4,clip_on=False))
    return a
def phi_heat(x,y,w,h,values,central=None):
    a=matrix(x,y,w,h,values,norm=phinorm)
    if central:a.add_patch(Rectangle((central[0],0),central[1]-central[0],1,fill=False,ec=INK,lw=.65,clip_on=False))
    return a
def response_scale(x,y,w=24,h=1.1):
    matrix(x,y,w,h,resp_norm.inverse(np.linspace(0,1,160)))
    for v in [-1.8,0,1.8]:t(x+w*resp_norm(v),y+h+1.8,f'{v:+g}'if v>0 else f'{v:g}',fs=5.3,ha='center')

# a/b are literal v3 accepted artwork under uniform translation/scaling.
abscale=57.5/85
placements=[
 {'panel':'a','source_clip_mm':[0,0,85,60],'destination_mm':[0,1,57.5,1+60*abscale]},
 {'panel':'b','source_clip_mm':[0,60,85,95.4],'destination_mm':[59.25,1,116.75,1+35.4*abscale]}]
# Panel identifiers are independently readable; these masks contain no data.
for x in [0,59.25]:
    rect(x+1.2,1.9,3.1,3.0,fc='white',ec='none',zorder=0)
t(2.5,3.5,'a',fs=9,fontweight='bold');t(61.75,3.5,'b',fs=9,fontweight='bold')
# Essential decoder below the shorter b graphic.
t(64.0,29.0,'Full: sequence + Φ',fs=6.1,color=FULL)
t(64.0,32.5,'Sequence: sequence only',fs=6.1,color=SEQ)
t(64.0,36.0,'Φ is fixed during the paired edit',fs=5.9)
t(64.0,39.5,'● HepG2     ▲ K562',fs=5.9,color=GRAY)

# c: literal natural alleles, unchanged Φ and all ten frozen RBP events.
example=rows(R/'FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv')
events=json.loads(gzip.decompress((R/'source_inputs/EVENTS_INCLUDED.json.gz').read_bytes()).decode('utf8'))
event=next(e for e in events if e['event_id']==example[0]['event_id'])
ref=event['context']['reference_sequence_rna_direction_dna'].replace('T','U')
alt=event['context']['minor_sequence_dna'].replace('T','U');phi=np.asarray(event['context']['phi_raw'])
assert len(example)==10 and ref[50]=='C'and alt[50]=='U'
t(119,3.5,'c',fs=9,fontweight='bold');t(124.3,3.5,'ARHGEF10 · K562',fs=7.0,fontweight='bold')
t(124.3,7.5,'chr8:1,844,440 C>T · GRCh38',fs=6.0)
t(121.1,11.2,'REF',fs=5.6);t(121.1,14.4,'ALT',fs=5.6)
sx,step=134.5,2.1
rect(sx+5*step-step/2,9.6,step,6.4,fc='#F7E4DF',ec='none',zorder=0)
letter_row(sx,11.2,ref[45:56],step=step,fs=6.1,center=5)
letter_row(sx,14.4,alt[45:56],step=step,fs=6.1,center=5)
phi_heat(sx-step/2,16.6,step*11,1.9,phi[45:56],(3,8))
t(121.1,17.55,'Φ',fs=6.0)
t(159.2,11.2,'nt 46–56',fs=5.3,color=GRAY)
t(159.2,14.4,'Φ51 = −3.50',fs=5.3)
t(121.1,20.7,'Φ REF = Φ ALT; observed preference: C',fs=5.7)
darms=['A_native','S','A_flat','A_permmean']
mx,mw,my,mh=135.0,41.,25.0,19.5
for xx,label,color in zip([140.125,150.375,160.625,170.875],
                          ['Full','Sequence','Flat Φ','Perm. Φ'],[FULL,SEQ,FLAT,PERM]):
    t(xx,23.6,label,fs=5.1,ha='center',color=color)
t(121.1,23.6,'RBP',fs=5.5)
dvalues=np.array([[float(r[a+'_ref_minus_alt'])for a in darms]for r in example])
matrix(mx,my,mw,mh,dvalues,outline=True)
for i,r in enumerate(example):
    yy=my+(i+.5)*mh/10;t(121.1,yy,r['RBP'],fs=5.3)
    assert float(r['ref_preferred'])==1
    for j in range(4):t(mx+(j+.5)*mw/4,yy,f'{dvalues[i,j]:+.3f}',fs=5.0,ha='center')
    line(121.1,my+(i+1)*mh/10,176,my+(i+1)*mh/10,color='#E2E7EB',lw=.25)
t(121.1,47.0,'Δlog count (C − U); all 10 events marked in e',fs=5.6)

# d: all 25 frozen, external-first/frequency-opposed sequence-order edits.
edits=rows(INPUT/'ALL_25_SELECTED_EDITS.csv');tasks=rows(INPUT/'ALL_36_TASKS.csv')
response=rows(INPUT/'ALL_1152_PAIRED_MODEL_RESPONSES.csv');ep=np.load(INPUT/'SOURCE_PACKET.npz')
packet=np.load(B/'INTERVENTION_PACKET.npz');proteins=packet['cohort_proteins'].tolist()
primary=rows(S/'FROZEN_PRIMARY_METRICS.csv')
coeff={(r['protein'],r['arm']):float(r['spearman_all5mer'])for r in primary}
assert len(edits)==25 and [r['protein']for r in edits]==proteins and len(tasks)==36 and len(response)==1152
cy,ch=65.,62.5
t(2.5,53,'d',fs=9,fontweight='bold');t(8,53,'Sequence-order preferences',fs=7.1,fontweight='bold')
t(8,57.2,'25 RBPs · 1024 words · 36 tasks × 16 contexts',fs=6.2)
headers=[(2.7,'RBP'),(23.8,'Same-count swap'),(36.25,'ΔZ5'),(42.25,'Δfreq'),(50.5,'All-word ρ'),(63.9,'Full'),(77.5,'Sequence')]
for x,label in headers:t(x,60.0,label,fs=5.8,ha='left'if x==2.7 else'center',color=FULL if label=='Full'else SEQ if label=='Sequence'else INK)
t(23.8,62.7,'low → high external Z5',fs=5.2,ha='center')
for x,label in [(47.5,'F'),(50.5,'S'),(53.5,'Fr'),(63.9,'16 contexts'),(77.5,'16 contexts')]:t(x,62.7,label,fs=5.1,ha='center')
rh=ch/25;fmat=np.full((50,16),np.nan);smat=fmat.copy();rhomat=np.zeros((25,3))
zax=ax(34,cy,4.5,ch);fax=ax(40,cy,4.5,ch)
for aa in[zax,fax]:aa.set_ylim(25,0);aa.set_axis_off()
zax.set_xlim(0,6);fax.set_xlim(-230,0)
for ri,r in enumerate(edits):
    p=r['protein'];yy=cy+(ri+.5)*rh;t(2.7,yy,p,fs=5.8)
    lo,hi=r['lower_assay_word'],r['higher_assay_word']
    assert sorted(lo)==sorted(hi)and float(r['external_delta_Z5'])>0 and float(r['native_frequency_delta'])<0
    letter_row(16.5,yy,lo);t(23.8,yy,'→',fs=5.5,ha='center');letter_row(25.5,yy,hi)
    zax.barh(ri+.5,float(r['external_delta_Z5']),height=.55,color=PERM)
    fax.barh(ri+.5,float(r['native_frequency_delta']),height=.55,color=FREQ)
    rhomat[ri]=[coeff[(p,a)]for a in['full','sequence','frequency']]
    for ci,cell in enumerate(['HepG2','K562']):
        for lx in[56.6,70.9]:t(lx,cy+(2*ri+ci+.5)*rh/2,'H'if ci==0 else'K',fs=4.0,ha='center',color=GRAY)
        inds=np.flatnonzero(ep[p+'_cells']==cell)
        if len(inds):
            assert len(inds)==16
            fmat[2*ri+ci]=ep[p+'_full_paired_delta'][inds];smat[2*ri+ci]=ep[p+'_sequence_paired_delta'][inds]
    line(2.7,cy+(ri+1)*rh,83,cy+(ri+1)*rh,color='#E4E8EB',lw=.25)
matrix(46,cy,9,ch,rhomat,norm=rhonorm,outline=True)
for vals,x,width in[(fmat,58,11.8),(smat,72,11)]:
    ma=matrix(x,cy,width,ch,vals,outline=True)
    for yy in range(2,50,2):ma.axhline(yy,color='white',lw=.55)
    for yy in range(1,50,2):ma.axhline(yy,color='white',lw=.18)
assert np.sum(np.isfinite(fmat))+np.sum(np.isfinite(smat))==1152
for xx,label in[(36.25,'0–6'),(42.25,'−230–0')]:t(xx,129.3,label,fs=4.7,ha='center',color=GRAY)
matrix(46,129.5,9,1.1,np.linspace(-1,1,80),norm=rhonorm)
t(46,132.5,'-1',fs=5,ha='center');t(55,132.5,'+1',fs=5,ha='center')
response_scale(61,129.5,21)
t(2.7,132.5,'H: HepG2; K: K562; gray: no task',fs=5.5,color=GRAY)
t(2.7,136.3,'Full − Frequency: ρ +0.146 [95% CI +0.038, +0.248]',fs=6.0)
t(2.7,139.5,'Within-composition: +3.78 pp [95% CI +0.90, +6.73]',fs=6.0)

# e: complete natural-allele atlas with one fixed non-outcome order.
cohort=rows(R/'STRICT396_WITH_ANNOTATIONS.csv')
cohort.sort(key=lambda r:(float(r['phi_center']),int(r['chrom'][3:]),int(r['pos1']),r['ref'],r['alt'],r['cell'],r['RBP'],r['event_id']))
assert len(cohort)==396 and len({r['variant_key']for r in cohort})==287 and len({r['RBP']for r in cohort})==85
qr={r['event_id']:int(r['phi_quartile'])for r in rows(R/'QUARTILE_ASSIGNMENTS_strict_test.csv')}
quart=[r for r in rows(R/'DIRECTION_BY_CONSERVATION.csv')if r['cohort']=='strict_test']
summary={r['arm']:r for r in rows(R/'DIRECTION_SUMMARY.csv')if r['cohort']=='strict_test'}
two={r['arm']:r for r in rows(R/'DIRECTION_TWO_WAY_SENSITIVITY.csv')if r['cohort']=='strict_test'}
t(95.5,53,'e',fs=9,fontweight='bold');t(101,53,'Natural allele preferences',fs=7.1,fontweight='bold')
t(101,57.2,'396 events · 287 variants · 85 RBPs · 119 RNA regions',fs=6.2)
excols=np.array([i for i,r in enumerate(cohort)if r['variant_key']=='chr8:1844440:C:T'])
assert len(excols)==10 and np.array_equal(excols,np.arange(excols[0],excols[-1]+1))
exstart,exend=int(excols[0]),int(excols[-1]+1)
cx,cw=118.,58.
counts=[sum(qr[r['event_id']]==i for r in cohort)for i in[1,2,3,4]];assert counts==[116,101,86,93]
starts=np.r_[0,np.cumsum(counts)]
for i in range(4):t(cx+cw*(starts[i]+counts[i]/2)/396,62.0,f'Q{i+1} · n={counts[i]}',fs=5.3,ha='center')
phi_vals=[float(r['phi_center'])for r in cohort];assert min(phi_vals)>=-7 and max(phi_vals)<=10
phi_heat(cx,64,cw,3.3,phi_vals)
nucmap={b:i for i,b in enumerate(['A','C','G','T'])};nuc_cmap=ListedColormap([BASE[b]for b in['A','C','G','T']]);nnorm=Normalize(-.5,3.5)
matrix(cx,68.3,cw,2.1,[nucmap[r['ref']]for r in cohort],norm=nnorm,colors=nuc_cmap)
matrix(cx,71.,cw,2.1,[nucmap[r['alt']]for r in cohort],norm=nnorm,colors=nuc_cmap)
matrix(cx,73.7,cw,2.2,[float(r['ref_preferred'])for r in cohort],norm=Normalize(0,1),colors=ListedColormap(['#F4A016','#19499C']))
regions=['intronic','UTR','coding','splice','other_annotated','unavailable']
regioncolors=['#8491B4','#91D1C2','#7E6148','#F39B7F','#B09C85','#DFE3E7']
matrix(cx,76.5,cw,2.2,[regions.index(r['NCBI_annotation_region'])for r in cohort],norm=Normalize(-.5,5.5),colors=ListedColormap(regioncolors))
for label,y,h in[('Φ',64,3.3),('REF',68.3,2.1),('ALT',71,2.1),('BEAPR',73.7,2.2),('Region',76.5,2.2)]:t(101,y+h/2,label,fs=6.0)
for arm,label,col,y in zip(darms,['Full','Sequence','Flat Φ','Perm. Φ'],[FULL,SEQ,FLAT,PERM],[80.2,85.2,90.2,95.2]):
    vv=np.array([float(r[arm+'_major_minus_minor'])for r in cohort]);assert min(vv)>=-1.8 and max(vv)<=1.8
    matrix(cx,y,cw,4,vv,outline=True);t(101,y+2,label,fs=6.1,color=col)
for boundary in starts[1:-1]:line(cx+cw*boundary/396,63.5,cx+cw*boundary/396,99.5,color='white',lw=.8)
rect(cx+cw*exstart/396,63.5,cw*(exend-exstart)/396,36,fill=False,ec='#B98932',lw=.7)
xx=cx+cw*(exstart+exend)/2/396
line(xx,99.5,xx,100.7,color='#B98932',lw=.6);t(xx,102.3,'c',fs=6.4,fontweight='bold',ha='center',color='#9B7023')
t(128,102.3,'Δlog count',fs=5.7);response_scale(149,100.6,26)
for i,(label,col)in enumerate(zip(['intron','UTR','coding','splice','other','unknown'],regioncolors)):
    px=101+i*12.7;rect(px,105.7,1.5,1.25,fc=col,ec='none');t(px+2.,106.4,label,fs=5.2)
t(101,109.8,'BEAPR:',fs=5.5);t(111,109.8,'REF',fs=5.5,color='#19499C');t(118,109.8,'ALT',fs=5.5,color='#F4A016')
t(134,109.8,'Φ increases →',fs=5.7)
t(136.3,114.0,'RBP-equal accuracy (%)',fs=5.6,ha='center');t(165,114.0,'Overall · 95% CI',fs=5.6,ha='center')
for j in range(4):t(124.2+6.5*j,117.0,f'Q{j+1}',fs=5.4,ha='center')
arms=['A_native','S','A_flat','A_permmean','ref'];labels=['Full','Sequence','Flat Φ','Perm. Φ','REF rule'];colors=[FULL,SEQ,FLAT,PERM,FREQ]
sy,sh=120.,13.5
sa=ax(152.5,sy,18.5,sh);sa.set_xlim(40,80);sa.set_ylim(4.5,-.5);clean(sa);sa.set_xticks([40,60,80]);sa.axvline(50,color=GRAY,lw=.4,ls=(0,(2,2)))
for i,(arm,label,col)in enumerate(zip(arms,labels,colors)):
    yy=sy+(i+.5)*sh/5;t(101,yy,label,fs=6.)
    if i%2==0:rect(120.5,sy+i*sh/5,26,sh/5,fc='#F2F5F7',ec='none',zorder=0)
    for j in range(4):
        q=next(r for r in quart if r['arm']==arm and int(r['quartile'])==j+1)
        t(124.2+6.5*j,yy,f"{100*float(q['accuracy']):.1f}",fs=5.8,ha='center')
    rr=summary[arm];ci=two[arm];assert 40<=100*float(ci['ci_low'])<=100*float(ci['ci_high'])<=80
    sa.plot([100*float(ci['ci_low']),100*float(ci['ci_high'])],[i,i],color=col,lw=.8)
    sa.scatter([100*float(rr['accuracy'])],[i],c=col,s=8,zorder=3)
    t(176.5,yy,f"{100*float(rr['accuracy']):.2f}",fs=5.7,ha='right',color=col)
t(101,139.5,'Full − Sequence: −6.73 pp [95% CI −15.97, +1.64]',fs=6.)
line(89,51,89,140.5,color='#D9DFE3',lw=.6)
line(1.5,49.5,176.5,49.5,color='#D9DFE3',lw=.6)
line(118,2,118,47,color='#D9DFE3',lw=.5)

overlay=O/'fig5_overlay';fig.savefig(overlay.with_suffix('.pdf'),transparent=True);fig.savefig(overlay.with_suffix('.svg'),transparent=True);plt.close(fig)
doc=fitz.open();page=doc.new_page(width=W*PT,height=H*PT);original=fitz.open(OLD/'fig5.pdf')
for pl in placements:page.show_pdf_page(fitz.Rect(*[v*PT for v in pl['destination_mm']]),original,0,
    clip=fitz.Rect(*[v*PT for v in pl['source_clip_mm']]))
new=fitz.open(overlay.with_suffix('.pdf'));page.show_pdf_page(page.rect,new,0)
pdf=O/'fig5.pdf';doc.save(pdf,garbage=4,deflate=True)
page.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False).save(O/'fig5.png')
page.get_pixmap(matrix=fitz.Matrix(180/72,180/72),alpha=False).save(O/'fig5_178mm_preview.png')
doc.close();original.close();new.close()
ns='http://www.w3.org/2000/svg';ET.register_namespace('',ns)
root=ET.Element('{'+ns+'}svg',{'width':'178mm','height':'142mm','viewBox':f'0 0 {W*PT} {H*PT}'})
ET.SubElement(root,'{'+ns+'}rect',{'x':'0','y':'0','width':str(W*PT),'height':str(H*PT),'fill':'white'})
defs=ET.SubElement(root,'{'+ns+'}defs')
def import_svg(path,parent,prefix):
    src=ET.parse(path).getroot();ids={el.attrib['id']:prefix+el.attrib['id']for el in src.iter()if'id'in el.attrib}
    for el in src.iter():
        if'id'in el.attrib:el.attrib['id']=ids[el.attrib['id']]
        for attr,val in list(el.attrib.items()):
            if attr.endswith('href')and val.startswith('#'):el.attrib[attr]='#'+ids.get(val[1:],val[1:])
            elif'url(#'in val:
                for old,nn in ids.items():val=val.replace('url(#'+old+')','url(#'+nn+')')
                el.attrib[attr]=val
    for child in list(src):parent.append(child)
for pl in placements:
    x0,y0,x1,y1=pl['source_clip_mm'];dx0,dy0,dx1,dy1=pl['destination_mm']
    cid=pl['panel']+'_clip';cp=ET.SubElement(defs,'{'+ns+'}clipPath',{'id':cid})
    ET.SubElement(cp,'{'+ns+'}rect',{'x':str(x0*PT),'y':str(y0*PT),'width':str((x1-x0)*PT),'height':str((y1-y0)*PT)})
    trans=f'translate({dx0*PT} {dy0*PT-y0*PT*abscale}) scale({abscale})'
    gg=ET.SubElement(root,'{'+ns+'}g',{'transform':trans})
    cg=ET.SubElement(gg,'{'+ns+'}g',{'clip-path':'url(#'+cid+')'})
    import_svg(OLD/'fig5.svg',cg,pl['panel']+'_v3_')
import_svg(overlay.with_suffix('.svg'),ET.SubElement(root,'{'+ns+'}g'),'v4_')
ET.ElementTree(root).write(O/'fig5.svg',encoding='utf-8',xml_declaration=True)
fields=[{k:r[k]for k in['event_id','variant_key','RBP','cell','phi_center','ref','alt','ref_preferred','NCBI_annotation_region']+[a+'_major_minus_minor'for a in darms]}for r in cohort]
with(O/'EVENT_ATLAS_COLUMNS.csv').open('w',newline='',encoding='utf8')as f:
    w=csv.DictWriter(f,fieldnames=['display_column']+list(fields[0]));w.writeheader()
    for i,r in enumerate(fields):w.writerow({'display_column':i+1,**r})
sources=[OLD/'fig5.pdf',OLD/'fig5.svg',B/'figure.pdf',B/'figure.svg',B/'INTERVENTION_PACKET.npz',
 S/'FROZEN_PRIMARY_METRICS.csv',S/'FROZEN_PRIMARY_CONTRASTS.csv',R/'STRICT396_WITH_ANNOTATIONS.csv',
 R/'FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv',R/'QUARTILE_ASSIGNMENTS_strict_test.csv',
 R/'DIRECTION_BY_CONSERVATION.csv',R/'DIRECTION_SUMMARY.csv',R/'DIRECTION_CONTRASTS.csv',
 R/'DIRECTION_TWO_WAY_SENSITIVITY.csv',R/'source_inputs/EVENTS_INCLUDED.json.gz']+list(INPUT.iterdir())
mapping={'canvas_mm':[W,H],'old_to_new_panels':{'a':'a','b':'b','c':'d','d':'c','e':'e'},
 'a_b_uniform_scale':abscale,'a_b_vector_source':'copied accepted v3 Fig5; same base palette and fixed-Phi inset',
 'a_b_placements':placements,'a_b_only_overlay_masks_mm':[[1.2,1.9,3.1,3.0],[60.45,1.9,3.1,3.0]],
 'a_b_axes_points_and_labels_rescaled_together':True,'all25_original_order':proteins,
 'd_complete_context_response_cells':1152,'d_existing_tasks':36,
 'd_primary':'exact frozen spearman_all5mer, all1024-word protein estimates',
 'd_secondary':'exact frozen within_composition_concordance summary; complete original table retained',
 'd_edits':'same external-first frequency-opposed selection, all25; no selection changed',
 'c_example_all10':True,'c_example_e_columns_1based':[exstart+1,exend],
 'e_complete_event_columns':396,'e_variant_count':287,'e_RBPs':85,'e_region_clusters':119,
 'e_sort':['phi_center','chromosome','pos1','ref','alt','cell','RBP','event_id'],
 'e_response_fields':[a+'_major_minus_minor'for a in darms],
 'shared_c_d_e_response_scale':{'range':[-1.8,1.8],'display_transform':'asinh','linear_width':.02,'model_specific_normalization':False},
 'shared_phi_raw_scale':[-7,10],'e_quartile_counts':counts,
 'e_quartile_summary':'frozen RBP-equal means; complete RBP-bootstrap CIs retained in supplement/source',
 'aggregate_CI_method':'2000 two-way RNA-region/RBP multiplier, seed20261015',
 'no_new_scientific_statistics':True,'source_hashes':{str(p):sha(p)for p in sources},
 'exports':{str(p):sha(p)for p in[pdf,O/'fig5.svg',O/'fig5.png',O/'fig5_178mm_preview.png']}}
(O/'FIG5_V4_SOURCE_MAPPING.json').write_text(json.dumps(mapping,indent=2)+'\n','utf8')
print('Fig5 v4 early render: 178 x142mm; original a/b uniformly scaled, all frozen cohorts retained.')
