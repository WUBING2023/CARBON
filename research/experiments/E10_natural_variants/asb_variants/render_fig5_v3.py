"""Linked full-cohort Fig5 v3; only frozen-data display and palette changes."""
from pathlib import Path
import os,csv,json,gzip,hashlib,shutil,re,xml.etree.ElementTree as ET
os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['OMP_NUM_THREADS']='4'
import numpy as np,fitz
from pypdf import PdfReader,PdfWriter
from pypdf.generic import ContentStream,FloatObject,DecodedStreamObject
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm,AsinhNorm,ListedColormap,Normalize

R=Path(__file__).resolve().parent;V3=R.parent;O=R/'fig5_v3';O.mkdir(exist_ok=True)
S=V3/'manuscript/source_data/S23_external_order_preferences'
B=R/'b1_reference'
EXT=Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/supervised_expansion_20261004/candidate_portfolio_20261006/order_candidates/B3')
INPUT=O/'frozen_inputs';INPUT.mkdir(exist_ok=True)
for name in ['SOURCE_PACKET.npz','ALL_25_SELECTED_EDITS.csv','ALL_36_TASKS.csv','ALL_1152_PAIRED_MODEL_RESPONSES.csv']:
    shutil.copyfile(EXT/name,INPUT/name)
W,H,PT=178.,190.,72/25.4
FULL,SEQ,FREQ,FLAT,PERM='#E64B35','#4DBBD5','#3C5488','#00A087','#8491B4'
INK,GRAY,PALE='#263238','#74808B','#D8DFE5'
BASE=json.loads((V3/'coordination/PALETTE_V3.json').read_text('utf8'))['nucleotides']
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.labelsize':6.5,
 'xtick.labelsize':6.0,'ytick.labelsize':6.0,'axes.linewidth':.55,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'none'})
def rows(p):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

# Scoped base-letter PDF fill changes preserve every drawing/text coordinate.
writer=PdfWriter(clone_from=B/'figure.pdf');page0=writer.pages[0]
stream=ContentStream(page0.get_contents(),writer)
fill=([FloatObject(0)],b'g');stack=[];base_count=0;base_fills=[]
for args,op in stream.operations:
    if op==b'q':stack.append(fill)
    elif op==b'Q':fill=stack.pop()
    elif op in [b'rg',b'g',b'k']:fill=(args,op)
    chars=[]
    if op==b'TJ':chars=[str(z) for z in args[0] if isinstance(z,str)]
    elif op==b'Tj':chars=[str(args[0])]
    if len(chars)==1 and chars[0] in BASE:
        base_fills.append((chars[0],fill));base_count+=1
# Preserve every original content-stream byte: generic PDF reserialization can
# round existing model-color operands and alter 8-bit raster colors by one unit.
raw=page0.get_contents().get_data();matches=list(re.finditer(rb'\[\s*\(\x00([ACGUT])\)\s*\]\s*TJ',raw))
assert len(matches)==base_count
chunks=[];cursor=0
for match,(base,prior)in zip(matches,base_fills):
    assert match.group(1).decode()==base
    hx=BASE[base].lstrip('#');color=' '.join(format(int(hx[i:i+2],16)/255,'.15g')for i in[0,2,4])+' rg\n'
    restore='\n'+' '.join(format(float(v),'.15g')for v in prior[0])+' '+prior[1].decode()
    chunks.extend([raw[cursor:match.start()],color.encode(),match.group(),restore.encode()]);cursor=match.end()
chunks.append(raw[cursor:]);modified=DecodedStreamObject();modified.set_data(b''.join(chunks));page0.replace_contents(modified)
with (O/'accepted_b1_palette.pdf').open('wb') as f:writer.write(f)
ns='http://www.w3.org/2000/svg';ET.register_namespace('',ns)
svg=ET.parse(B/'figure.svg');svg_bases=0
for el in svg.getroot().iter():
    if el.tag.endswith('text') and (el.text or '') in BASE:
        style=el.attrib.get('style','')
        if re.search(r'fill:\s*[^;]+',style):style=re.sub(r'fill:\s*[^;]+','fill: '+BASE[el.text],style)
        else:style='fill: '+BASE[el.text]+'; '+style
        el.set('style',style);svg_bases+=1
svg.write(O/'accepted_b1_palette.svg',encoding='utf-8',xml_declaration=True)
assert base_count==svg_bases

fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='none')
def t(x,y,z,fs=6.5,**kw):return fig.text(x/W,1-y/H,z,fontsize=fs,va='center',**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def line(x0,y0,x1,y1,**kw):
    obj=Line2D([x0/W,x1/W],[1-y0/H,1-y1/H],transform=fig.transFigure,**kw);fig.add_artist(obj);return obj
def rect(x,y,w,h,**kw):
    obj=Rectangle((x/W,1-(y+h)/H),w/W,h/H,transform=fig.transFigure,**kw);fig.add_artist(obj);return obj
def clean(a,keep_left=False):
    a.spines[['top','right']].set_visible(False)
    if not keep_left:a.spines['left'].set_visible(False);a.set_yticks([])
    a.tick_params(length=2,width=.5,pad=1.2)
def letter_row(x,y,seq,step=2.7,center=5,fs=7,mono=False):
    for j,b in enumerate(seq):t(x+j*step,y,b,fs=fs,color=BASE[b],ha='center',fontweight='bold' if j==center else 'normal',fontfamily='DejaVu Sans Mono' if mono else 'Arial')
cmap=LinearSegmentedColormap.from_list('signed',[FREQ,'#FFFFFF',FULL]);cmap.set_bad('#D8DFE5')
phinorm=TwoSlopeNorm(vmin=-7,vcenter=0,vmax=10)
resp_norm=AsinhNorm(linear_width=.02,vmin=-1.8,vmax=1.8)
rhonorm=Normalize(-1,1)
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
def phi_scale(x,y,w=22,h=1.1):
    # A uniform color-ramp coordinate, labelled with exact inverse raw values.
    vals=phinorm.inverse(np.linspace(0,1,120));phi_heat(x,y,w,h,vals)
    for v in [-7,0,10]:t(x+w*phinorm(v),y+h+1.8,str(v)if v<=0 else '+'+str(v),fs=5.2,ha='center')
def response_scale(x,y,w=24,h=1.1):
    matrix(x,y,w,h,resp_norm.inverse(np.linspace(0,1,160)))
    for v in [-1.8,0,1.8]:t(x+w*resp_norm(v),y+h+1.8,f'{v:+g}'if v>0 else f'{v:g}',fs=5.3,ha='center')

# a/b: literal accepted vector crop with declared nucleotide-only recoloring.
packet=np.load(B/'INTERVENTION_PACKET.npz');qphi=packet['QKI_native_conservation'][0]
rect(35,5,50,15.1,fc='white',ec='none',zorder=0)
t(61,6.7,'Φ fixed reference input · PhyloP',fs=6.1,ha='center')
phi_heat(40,9.1,42,2.65,qphi,(48,53))
for xx,label in [(40,'1'),(61,'49-53'),(82,'101 nt')]:t(xx,13.4,label,fs=5.6,ha='center'if xx<82 else'right',color=GRAY)
t(43,17,'Raw Φ',fs=5.6);phi_scale(59.5,15.9)

# c: complete linked edit/context atlas, all-word primary coefficients.
edits=rows(INPUT/'ALL_25_SELECTED_EDITS.csv');tasks=rows(INPUT/'ALL_36_TASKS.csv')
response=rows(INPUT/'ALL_1152_PAIRED_MODEL_RESPONSES.csv');ep=np.load(INPUT/'SOURCE_PACKET.npz')
primary=rows(S/'FROZEN_PRIMARY_METRICS.csv');coeff={(r['protein'],r['arm']):float(r['spearman_all5mer'])for r in primary}
proteins=packet['cohort_proteins'].tolist()
assert len(edits)==25 and [r['protein']for r in edits]==proteins
assert len(tasks)==36 and len(response)==1152
cy,ch=113.5,62.5
t(2.5,99,'c',fs=9,fontweight='bold');t(8,99,'External sequence-order preferences',fs=7.1,fontweight='bold')
t(8,103.2,'25 RBPs · 1024 words · 36 tasks × 16 contexts',fs=6.3)
headers=[(2.7,'RBP'),(23.8,'Same-count swap'),(36.25,'ΔZ5'),(42.25,'Δfreq'),(50.5,'All-word ρ'),(63.9,'Full'),(77.5,'Sequence')]
for x,label in headers:t(x,108,label,fs=5.8,ha='left'if x==2.7 else'center',color=FULL if label=='Full'else SEQ if label=='Sequence'else INK)
t(23.8,110.7,'low → high external Z5',fs=5.2,ha='center')
for x,label in [(47.5,'F'),(50.5,'S'),(53.5,'Fr'),(63.9,'16 contexts'),(77.5,'16 contexts')]:t(x,110.7,label,fs=5.1,ha='center')
rh=ch/25;fmat=np.full((50,16),np.nan);smat=fmat.copy();rhomat=np.zeros((25,3))
zax=ax(34.0,cy,4.5,ch);fax=ax(40.0,cy,4.5,ch)
for aa in [zax,fax]:aa.set_ylim(25,0);aa.set_axis_off()
zax.set_xlim(0,6);fax.set_xlim(-230,0)
for ri,r in enumerate(edits):
    p=r['protein'];yy=cy+(ri+.5)*rh
    t(2.7,yy,p,fs=5.8)
    lo,hi=r['lower_assay_word'],r['higher_assay_word']
    assert sorted(lo)==sorted(hi) and float(r['external_delta_Z5'])>0 and float(r['native_frequency_delta'])<0
    letter_row(16.5,yy,lo,step=1.4,center=-1,fs=5.8,mono=True)
    t(23.8,yy,'→',fs=5.5,ha='center');letter_row(25.5,yy,hi,step=1.4,center=-1,fs=5.8,mono=True)
    zax.barh(ri+.5,float(r['external_delta_Z5']),height=.55,color=PERM)
    fax.barh(ri+.5,float(r['native_frequency_delta']),height=.55,color=FREQ)
    rhomat[ri]=[coeff[(p,a)]for a in ['full','sequence','frequency']]
    for ci,cell in enumerate(['HepG2','K562']):
        for lx in [56.6,70.9]:t(lx,cy+(2*ri+ci+.5)*rh/2,'H' if ci==0 else 'K',fs=4.0,ha='center',color=GRAY)
        inds=np.flatnonzero(ep[p+'_cells']==cell)
        if len(inds):
            assert len(inds)==16
            fmat[2*ri+ci]=ep[p+'_full_paired_delta'][inds]
            smat[2*ri+ci]=ep[p+'_sequence_paired_delta'][inds]
    line(2.7,cy+(ri+1)*rh,83,cy+(ri+1)*rh,color='#E4E8EB',lw=.25)
matrix(46,cy,9,ch,rhomat,norm=rhonorm,outline=True)
for vals,x,width in [(fmat,58,11.8),(smat,72,11)]:
    ma=matrix(x,cy,width,ch,vals,outline=True)
    for yy in range(2,50,2):ma.axhline(yy,color='white',lw=.55)
    for yy in range(1,50,2):ma.axhline(yy,color='white',lw=.18)
assert np.sum(np.isfinite(fmat))+np.sum(np.isfinite(smat))==1152
for xx,label in [(36.25,'0–6'),(42.25,'−230–0')]:t(xx,177.8,label,fs=4.7,ha='center',color=GRAY)
matrix(46,178.0,9,1.1,np.linspace(-1,1,80),norm=rhonorm)
t(46,181.0,'-1',fs=5.0,ha='center');t(55,181.0,'+1',fs=5.0,ha='center')
response_scale(61,178.0,21)
t(2.7,181,'H: HepG2; K: K562; gray: no task',fs=5.6,color=GRAY)
t(2.7,184.8,'Full - Frequency: ρ +0.146 [95% CI +0.038, +0.248]',fs=6.0)
t(2.7,188.0,'Within-composition: +3.78 pp [95% CI +0.90, +6.73]',fs=6.0)

# d: natural allele identity, coordinate input and complete four-arm responses.
example=rows(R/'FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv')
assert len(example)==10 and [r['RBP']for r in example]==sorted(r['RBP']for r in example)
events=json.loads(gzip.decompress((R/'source_inputs/EVENTS_INCLUDED.json.gz').read_bytes()).decode('utf8'))
event=next(e for e in events if e['event_id']==example[0]['event_id'])
ref=event['context']['reference_sequence_rna_direction_dna'].replace('T','U')
alt=event['context']['minor_sequence_dna'].replace('T','U');phi=np.asarray(event['context']['phi_raw'])
assert ref[50]=='C'and alt[50]=='U'
t(95.5,4,'d',fs=9,fontweight='bold');t(101,4,'Natural SNV · ARHGEF10 / K562',fs=7.1,fontweight='bold')
t(101,8,'GRCh38 chr8:1,844,440 C>T',fs=6.5)
t(101,14.4,'REF',fs=6.2);t(101,18.8,'ALT',fs=6.2)
sx,step=115.,2.7
rect(sx+5*step-step/2,12.55,step,8.1,fc='#F7E4DF',ec='none',zorder=0)
letter_row(sx,14.4,ref[45:56],step=step);letter_row(sx,18.8,alt[45:56],step=step)
phi_heat(sx-step/2,22.0,step*11,2.65,phi[45:56],(3,8))
for xx,label in [(sx-step/2,'46'),(sx+5*step,'51'),(sx+10.5*step,'56 nt')]:t(xx,26.2,label,fs=5.5,ha='center'if label!='56 nt'else'right',color=GRAY)
t(162.5,13,'101-nt reference Φ',fs=5.8,ha='center');phi_heat(149,15.1,27,2.65,phi,(48,53))
for xx,label in [(149,'1'),(162.5,'49-53'),(176,'101')]:t(xx,19.5,label,fs=5.4,ha='center')
t(162.5,23.3,'Center Φ = -3.50',fs=6,ha='center');phi_scale(152,26,23)
t(101,32.3,'Φ REF = Φ ALT · same fixed input',fs=6.3)
t(101,37,'RBP',fs=6.4);t(125,36.1,'BEAPR',fs=5.9,ha='center');t(125,38.2,'preferred',fs=5.4,ha='center')
darms=['A_native','S','A_flat','A_permmean']
for xx,label,color in zip([137.8,148.7,159.6,170.5],['Full','Sequence','Flat Φ','Perm. Φ'],[FULL,SEQ,FLAT,PERM]):t(xx,37,label,fs=5.9,ha='center',color=color)
dvalues=np.array([[float(r[a+'_ref_minus_alt'])for a in darms]for r in example])
matrix(132.5,40,43.5,44,dvalues,outline=True)
for i,r in enumerate(example):
    yy=40+(i+.5)*4.4;t(101,yy,r['RBP'],fs=6.5)
    assert float(r['ref_preferred'])==1;t(125,yy,'C',fs=6.5,ha='center',fontweight='bold',color=BASE['C'])
    for j,xx in enumerate([137.8,148.7,159.6,170.5]):t(xx,yy,f'{dvalues[i,j]:+.3f}',fs=6.0,ha='center',color=INK)
    line(101,40+(i+1)*4.4,176,40+(i+1)*4.4,color='#E2E7EB',lw=.35)
t(101,88.4,'Predicted Δlog count (C - U); positive agrees with REF',fs=6.2)
t(101,93.1,'All 10 RBP events are retained and located in e',fs=6.2,color=GRAY)

# e: one fixed, non-output ordering for every event/input/model row.
cohort=rows(R/'STRICT396_WITH_ANNOTATIONS.csv')
cohort.sort(key=lambda r:(float(r['phi_center']),int(r['chrom'][3:]),int(r['pos1']),r['ref'],r['alt'],r['cell'],r['RBP'],r['event_id']))
assert len(cohort)==396 and len({r['variant_key']for r in cohort})==287 and len({r['RBP']for r in cohort})==85
qr={r['event_id']:int(r['phi_quartile'])for r in rows(R/'QUARTILE_ASSIGNMENTS_strict_test.csv')}
quart=rows(R/'DIRECTION_BY_CONSERVATION.csv');quart=[r for r in quart if r['cohort']=='strict_test']
summary={r['arm']:r for r in rows(R/'DIRECTION_SUMMARY.csv')if r['cohort']=='strict_test'}
two={r['arm']:r for r in rows(R/'DIRECTION_TWO_WAY_SENSITIVITY.csv')if r['cohort']=='strict_test'}
contr=next(r for r in rows(R/'DIRECTION_CONTRASTS.csv')if r['cohort']=='strict_test'and r['contrast']=='A_native minus S')
t(95.5,99,'e',fs=9,fontweight='bold');t(101,99,'Complete allele-response atlas',fs=7.1,fontweight='bold')
t(101,103.2,'396 events · 287 variants · 85 RBPs · 119 RNA regions',fs=6.2)
excols=np.array([i for i,r in enumerate(cohort)if r['variant_key']=='chr8:1844440:C:T'])
assert len(excols)==10 and np.array_equal(excols,np.arange(excols[0],excols[-1]+1))
exstart,exend=int(excols[0]),int(excols[-1]+1)
cx,cw=118.,58.
counts=[sum(qr[r['event_id']]==i for r in cohort)for i in [1,2,3,4]]
assert counts==[116,101,86,93]
starts=np.r_[0,np.cumsum(counts)]
for i in range(4):
    mid=cx+cw*(starts[i]+counts[i]/2)/396
    t(mid,108,f'Q{i+1} · n={counts[i]}',fs=5.3,ha='center')
features=[('Φ',110,3.3),('REF',114.3,2.1),('ALT',117.0,2.1),('BEAPR',119.7,2.2),('Region',122.5,2.2)]
phi_vals=[float(r['phi_center'])for r in cohort];assert min(phi_vals)>=-7 and max(phi_vals)<=10
phi_heat(cx,110,cw,3.3,phi_vals)
nucmap={b:i for i,b in enumerate(['A','C','G','T'])};nuc_cmap=ListedColormap([BASE[b]for b in ['A','C','G','T']]);nnorm=Normalize(-.5,3.5)
matrix(cx,114.3,cw,2.1,[nucmap[r['ref']]for r in cohort],norm=nnorm,colors=nuc_cmap)
matrix(cx,117,cw,2.1,[nucmap[r['alt']]for r in cohort],norm=nnorm,colors=nuc_cmap)
matrix(cx,119.7,cw,2.2,[float(r['ref_preferred'])for r in cohort],norm=Normalize(0,1),colors=ListedColormap(['#F4A016','#19499C']))
regions=['intronic','UTR','coding','splice','other_annotated','unavailable']
regioncolors=['#8491B4','#91D1C2','#7E6148','#F39B7F','#B09C85','#DFE3E7']
matrix(cx,122.5,cw,2.2,[regions.index(r['NCBI_annotation_region'])for r in cohort],norm=Normalize(-.5,5.5),colors=ListedColormap(regioncolors))
for label,y,h in features:t(101,y+h/2,label,fs=6.0)
for arm,label,color,y in zip(darms,['Full','Sequence','Flat Φ','Perm. Φ'],[FULL,SEQ,FLAT,PERM],[126.2,131.2,136.2,141.2]):
    vv=np.array([float(r[arm+'_major_minus_minor'])for r in cohort]);assert np.min(vv)>=-1.8 and np.max(vv)<=1.8
    matrix(cx,y,cw,4.0,vv,outline=True);t(101,y+2,label,fs=6.1,color=color)
for boundary in starts[1:-1]:line(cx+cw*boundary/396,109.5,cx+cw*boundary/396,145.5,color='white',lw=.8)
rect(cx+cw*exstart/396,109.5,cw*(exend-exstart)/396,36,fill=False,ec='#B98932',lw=.7)
xx=cx+cw*(exstart+exend)/2/396
line(xx,145.5,xx,146.7,color='#B98932',lw=.6)
t(xx,148.3,'d',fs=6.4,fontweight='bold',ha='center',color='#9B7023')
t(128,148.3,'Δlog count',fs=5.7)
response_scale(149,146.6,26)
# Annotation legend is factual context; no clinical/pathogenicity categories.
for i,(label,col)in enumerate(zip(['intron','UTR','coding','splice','other','unknown'],regioncolors)):
    px=101+i*12.7;rect(px,152.2,1.5,1.25,fc=col,ec='none');t(px+2.0,152.9,label,fs=5.2)
t(101,156.3,'BEAPR:',fs=5.5);t(111,156.3,'REF',fs=5.5,color='#19499C');t(118,156.3,'ALT',fs=5.5,color='#F4A016')
t(134,156.3,'Φ increases →',fs=5.7)
t(137.3,161.0,'RBP-equal accuracy by Φ quartile',fs=5.6,ha='center')
t(165,161.0,'Overall · 95% CI',fs=5.6,ha='center')
for j in range(4):t(124.2+6.5*j,164.0,f'Q{j+1}',fs=5.4,ha='center')
arms=['A_native','S','A_flat','A_permmean','ref'];labels=['Full','Sequence','Flat Φ','Perm. Φ','REF rule'];colors=[FULL,SEQ,FLAT,PERM,FREQ]
sa=ax(152.5,167,18.5,16.5);sa.set_xlim(40,80);sa.set_ylim(4.5,-.5);clean(sa);sa.set_xticks([40,60,80]);sa.axvline(50,color=GRAY,lw=.4,ls=(0,(2,2)))
for i,(arm,label,col)in enumerate(zip(arms,labels,colors)):
    yy=167+(i+.5)*3.3;t(101,yy,label,fs=6.0)
    if i%2==0:rect(120.5,167+i*3.3,26,3.3,fc='#F2F5F7',ec='none',zorder=0)
    for j in range(4):
        q=next(r for r in quart if r['arm']==arm and int(r['quartile'])==j+1)
        t(124.2+6.5*j,yy,f"{100*float(q['accuracy']):.1f}",fs=5.8,ha='center')
    rr=summary[arm];ci=two[arm];assert 40<=100*float(ci['ci_low'])<=100*float(ci['ci_high'])<=80
    sa.plot([100*float(ci['ci_low']),100*float(ci['ci_high'])],[i,i],color=col,lw=.8)
    sa.scatter([100*float(rr['accuracy'])],[i],c=col,s=8,zorder=3)
    t(176.5,yy,f"{100*float(rr['accuracy']):.2f}",fs=5.7,ha='right',color=col)
t(101,188.0,'Full - Sequence: -6.73 pp [95% CI -15.97, +1.64]',fs=6.0)
line(89,2,89,189,color='#D9DFE3',lw=.6)

overlay=O/'fig5_overlay';fig.savefig(overlay.with_suffix('.pdf'),transparent=True);fig.savefig(overlay.with_suffix('.svg'),transparent=True);plt.close(fig)
doc=fitz.open();page=doc.new_page(width=W*PT,height=H*PT)
original=fitz.open(O/'accepted_b1_palette.pdf');page.show_pdf_page(fitz.Rect(0,0,85*PT,96*PT),original,0,clip=fitz.Rect(0,0,85*PT,96*PT))
new=fitz.open(overlay.with_suffix('.pdf'));page.show_pdf_page(page.rect,new,0)
pdf=O/'fig5.pdf';doc.save(pdf,garbage=4,deflate=True)
page.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False).save(O/'fig5.png')
page.get_pixmap(matrix=fitz.Matrix(180/72,180/72),alpha=False).save(O/'fig5_178mm_preview.png')
doc.close();original.close();new.close()
root=ET.Element('{'+ns+'}svg',{'width':'178mm','height':'190mm','viewBox':f'0 0 {W*PT} {H*PT}'})
ET.SubElement(root,'{'+ns+'}rect',{'x':'0','y':'0','width':str(W*PT),'height':str(H*PT),'fill':'white'})
defs=ET.SubElement(root,'{'+ns+'}defs');clip=ET.SubElement(defs,'{'+ns+'}clipPath',{'id':'ab_clip'});ET.SubElement(clip,'{'+ns+'}rect',{'width':str(85*PT),'height':str(96*PT)})
def import_svg(path,parent,prefix):
    src=ET.parse(path).getroot();ids={el.attrib['id']:prefix+el.attrib['id']for el in src.iter()if'id'in el.attrib}
    for el in src.iter():
        if'id'in el.attrib:el.attrib['id']=ids[el.attrib['id']]
        for attr,val in list(el.attrib.items()):
            if attr.endswith('href')and val.startswith('#'):el.attrib[attr]='#'+ids.get(val[1:],val[1:])
            elif'url(#'in val:
                for old,newid in ids.items():val=val.replace('url(#'+old+')','url(#'+newid+')')
                el.attrib[attr]=val
    for child in list(src):parent.append(child)
oldgroup=ET.SubElement(root,'{'+ns+'}g',{'clip-path':'url(#ab_clip)'});import_svg(O/'accepted_b1_palette.svg',oldgroup,'b1_')
import_svg(overlay.with_suffix('.svg'),ET.SubElement(root,'{'+ns+'}g'),'v3_');ET.ElementTree(root).write(O/'fig5.svg',encoding='utf-8',xml_declaration=True)
fields=[{k:r[k]for k in ['event_id','variant_key','RBP','cell','phi_center','ref','alt','ref_preferred','NCBI_annotation_region']+ [a+'_major_minus_minor'for a in darms]}for r in cohort]
with(O/'EVENT_ATLAS_COLUMNS.csv').open('w',newline='',encoding='utf8')as f:
    w=csv.DictWriter(f,fieldnames=['display_column']+list(fields[0]));w.writeheader()
    for i,r in enumerate(fields):w.writerow({'display_column':i+1,**r})
sources=[B/'figure.pdf',B/'figure.svg',B/'INTERVENTION_PACKET.npz',S/'FROZEN_PRIMARY_METRICS.csv',S/'FROZEN_PRIMARY_CONTRASTS.csv',R/'STRICT396_WITH_ANNOTATIONS.csv',R/'FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv',R/'QUARTILE_ASSIGNMENTS_strict_test.csv',R/'DIRECTION_BY_CONSERVATION.csv',R/'DIRECTION_SUMMARY.csv',R/'DIRECTION_CONTRASTS.csv',R/'DIRECTION_TWO_WAY_SENSITIVITY.csv',R/'source_inputs/EVENTS_INCLUDED.json.gz']+list(INPUT.iterdir())
mapping={'canvas_mm':[W,H],'base_palette':BASE,'original_geometry_retained':True,'base_glyph_fill_changes':base_count,
 'only_other_a_b_change_phi_mask_mm':[35,5,50,15.1],'all25_original_order':proteins,
 'c_complete_context_response_cells':1152,'c_existing_tasks':36,'c_primary':'spearman_all5mer, exact all1024-word frozen protein estimates',
 'c_secondary':'within_composition_concordance, original frozen summary and full table retained',
 'c_selected_edits':'external-first frequency-opposed; all25; no selection changed',
 'e_complete_event_columns':396,'e_variant_count':287,'e_RBPs':85,'e_region_clusters':119,
 'e_sort':['phi_center','chromosome','pos1','ref','alt','cell','RBP','event_id'],
 'e_response_fields':[a+'_major_minus_minor'for a in darms],'shared_response_scale':{'range':[-1.8,1.8],'display_transform':'asinh','linear_width':.02,'model_specific_normalization':False},
 'shared_phi_raw_scale':[-7,10],'d_example_e_columns_1based':[exstart+1,exend],
 'e_quartile_counts':counts,'e_quartile_summary':'exact frozen RBP-equal means; full RBP-bootstrap intervals retained in source/supplement',
 'aggregate_CI_method':'2000 two-way RNA-region/RBP multiplier, seed20261015',
 'no_new_scientific_statistics':True,'source_hashes':{str(p):sha(p)for p in sources},'exports':{str(p):sha(p)for p in [pdf,O/'fig5.svg',O/'fig5.png']}}
(O/'FIG5_V3_SOURCE_MAPPING.json').write_text(json.dumps(mapping,indent=2)+'\n','utf8')
shutil.copyfile(pdf,V3/'manuscript/figures/fig5.pdf')
print('First native full-atlas Fig5 v3 rendered, 178 x190 mm; '+str(base_count)+' nucleotide glyph fills changed.')
