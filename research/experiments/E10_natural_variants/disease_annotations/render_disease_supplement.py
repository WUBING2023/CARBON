"""Native supplementary atlas from complete disease/source joins."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

R=Path(__file__).resolve().parent
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
es=rows(R/'DISEASE_ALL13_EVENT_RESULTS.csv');vs=rows(R/'DISEASE_ALL8_VARIANT_RESULTS.csv')
allv=rows(R/'ALL4669_VARIANT_CONTEXT.csv');cv=[r for r in allv if r['ClinVar_exact_match']=='exact_GRCh38_REF_ALT']
af=[r for r in allv if r['ClinVar_AF_TGP']!='']
assert len(es)==13 and len(vs)==8 and len(cv)==520 and len(af)==428
FULL,SEQ,INK,FREQ,PERM,GOLD='#E64B35','#4DBBD5','#263238','#3C5488','#8491B4','#B58A31'
W,H=178.,192.
plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':6.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.5})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');values=[]
def txt(x,y,s,fs=6.5,col=INK,ha='left',weight='normal'):
    return fig.text(x/W,1-y/H,s,fontsize=fs,color=col,ha=ha,va='center',weight=weight)
def axis(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def tag(x,y,l,s):txt(x,y,l,10,weight='bold');txt(x+5,y,s,7.7,weight='bold')
def clean(a,y=True):
    a.spines[['top','right']].set_visible(False)
    if not y:a.spines['left'].set_visible(False)
    a.tick_params(length=2,pad=1.3,width=.45,labelsize=6.)
def record(panel,key,field,value):values.append(dict(panel=panel,id=key,field=field,value=value))

tag(1.5,3,'a','Whole existing variant population')
txt(2,10,'4,669 exact GRCh38 REF/ALT variants · 7,262 previously exposed ASB events',7.)
txt(2,17,'ClinVar',6.8,FREQ,weight='bold');txt(24,17,'520 exact matches',7.,weight='bold')
txt(63,17,'P/LP3 · B/LB453 · VUS4 · conflicts7 · other1 · no germline class52',6.3)
txt(2,24,'TCGA',6.8,FREQ,weight='bold');txt(24,24,'6 exact matches',7.,weight='bold')
txt(63,24,'6 singleton case observations · 8-site union · 13 ASB events',6.3)
txt(2,31,'ClinVar unmatched4,149; TCGA unmatched4,663. Successful absence in source is not biological absence.',6.2,'#74808B')

tag(1.5,39,'b','All disease-associated events')
for x,s in [(2,'Gene / exact allele'),(58,'RBP / cell'),(83,'Database'),(102,'ASB'),(170,'Φ')]:txt(x,45,s,6.4,weight='bold')
labels=[r['event_id'] for r in es]
bs=axis(122,49,38,52);bs.set_xscale('asinh',linear_width=.02)
bs.set(xlim=(-1,1),ylim=(12.7,-.7),xticks=[-1,-.1,0,.1,1],yticks=[])
bs.set_xticklabels(['−1','−.1','0','.1','1']);bs.axvline(0,color='#BDC8D2',lw=.65);clean(bs,False)
vmap={v['variant_key']:v for v in vs};last=None
for i,e in enumerate(es):
    k=e['variant_key'];v=vmap[k];y=49+(i+.7)/13.4*52
    if k!=last:
        txt(2,y,v['gene'],6.6,weight='bold')
        txt(17,y,k.rsplit(':',2)[0]+' '+k.split(':')[2]+'>'+k.split(':')[3],5.8)
        db=('CV' if int(v['ClinVar_PLP']) else '')+('/' if int(v['ClinVar_PLP']) and int(v['TCGA']) else '')+v['TCGA_projects'].replace('TCGA-','')
        txt(83,y,db,6.,FREQ)
        txt(176,y,f"{float(v['reference_PhyloP']):.2f}",6.1,PERM,'right')
    txt(58,y,e['RBP']+' / '+e['cell'],5.9)
    txt(102,y,'REF',6.,FREQ)
    fv=float(e['full_ALT_minus_REF']);sv=float(e['sequence_ALT_minus_REF'])
    bs.plot([fv,sv],[i,i],color='#BBC6D1',lw=.7)
    bs.scatter(fv,i,s=15,c=FULL,edgecolor='white',lw=.2,zorder=3)
    bs.scatter(sv,i,s=14,c=SEQ,marker='^',edgecolor='white',lw=.2,zorder=3)
    record('b',e['event_id'],'Full_ALT_minus_REF',fv);record('b',e['event_id'],'Sequence_ALT_minus_REF',sv)
    last=k
bs.set_xlabel('ALT − REF Δlog count · asinh',labelpad=1.5)
txt(2,106,'CV: aggregate P/LP, 2 stars. Cancer abbreviations: TCGA project; one case/site.',6.1)
txt(2,111,'Source prefers REF all13 · Full7/13, Sequence10/13 · 11 training /2 validation',6.1)
txt(123,111,'● Full',6.3,FULL);txt(145,111,'▲ Sequence',6.3,SEQ)

tag(1.5,120,'c','Clinical class and reference conservation')
groups=['Benign/Likely benign','Uncertain significance','Conflicting classifications','No germline classification','Other classification','Pathogenic/Likely pathogenic']
short=['B/LB','VUS','Conflict','No germline','Other','P/LP']
ca=axis(13,131,70,43)
rng=np.random.default_rng(20261014)
for i,(g,label) in enumerate(zip(groups,short)):
    rr=[r for r in cv if r['ClinVar_group']==g];ys=np.array([float(r['reference_PhyloP']) for r in rr]);jitter=rng.uniform(-.22,.22,len(rr))
    col=GOLD if g=='Pathogenic/Likely pathogenic' else '#8494A7'
    ca.scatter(i+jitter,ys,s=13 if col==GOLD else 4,color=col,alpha=1 if col==GOLD else .5,edgecolor='none')
    ca.plot([i-.25,i+.25],[np.median(ys)]*2,color=INK,lw=1.)
    txt(13+(i+.55)/6.1*70,128.5,'n='+str(len(rr)),5.9,ha='center')
    for r,val in zip(rr,ys):record('c',r['variant_key'],'reference_PhyloP',float(val))
    if col==GOLD:
        for r,j,val in zip(rr,jitter,ys):
            name=vmap[r['variant_key']]['gene'];ca.annotate(name,(i+j,val),xytext=(4,0),textcoords='offset points',fontsize=5.8,color=GOLD,va='center')
ca.set(xlim=(-.55,5.7),xticks=range(6),xticklabels=short)
ca.tick_params(axis='x',labelrotation=35);ca.set_ylabel('Reference PhyloP100way',labelpad=2.);clean(ca)
ca.axhline(0,color='#C7D0D8',lw=.55,zorder=0)

tag(92,120,'d','Population ALT frequency')
da=axis(103,131,72,43)
xs=np.array([float(r['ClinVar_AF_TGP']) for r in af]);ys=np.array([float(r['reference_PhyloP']) for r in af])
assert np.all((xs>=0)&(xs<=1))
da.set_xscale('symlog',linthresh=.0001)
da.scatter(xs,ys,s=6,color='#8494A7',alpha=.58,edgecolor='none')
da.axhline(0,color='#C7D0D8',lw=.55,zorder=0)
for r,x,y in zip(af,xs,ys):
    record('d',r['variant_key'],'ClinVar_TGP_ALT_AF',float(x));record('d',r['variant_key'],'reference_PhyloP',float(y))
    if r['ClinVar_group']=='Pathogenic/Likely pathogenic':
        da.scatter([x],[y],s=25,color=GOLD,edgecolor='white',lw=.4,zorder=5)
        da.annotate(vmap[r['variant_key']]['gene'],(x,y),xytext=(5,2),textcoords='offset points',fontsize=6.3,color=GOLD)
da.set(xlim=(-.00002,1.2),xticks=[0,.0001,.01,1],xticklabels=['0','10⁻⁴','10⁻²','1'])
da.set_xlabel('ClinVar TGP ALT AF · symlog',labelpad=1.5);da.set_ylabel('Reference PhyloP100way',labelpad=1.8);clean(da)
txt(92,128.5,'428 available · missing AF remains unavailable',5.9)
txt(2,188,'Clinical class and tumor occurrence annotate context; they do not validate a disease-specific RBP mechanism.',6.2,'#74808B')

fig.canvas.draw();outside=[];rend=fig.canvas.get_renderer()
for t in fig.findobj(matplotlib.text.Text):
    if not t.get_visible() or not t.get_text():continue
    b=t.get_window_extent(rend)
    if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:outside.append(t.get_text())
assert not outside,outside
for ext in ('pdf','svg','png'):fig.savefig(R/f'figS27.{ext}',dpi=220,facecolor='white')
with (R/'S27_PLOTTED_VALUES.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(values[0]));w.writeheader();w.writerows(values)
(R/'S27_RENDER_QA.json').write_text(json.dumps({'native_mm':[W,H],'unique_ClinVar_points':520,'population_AF_points':428,
 'all_disease_events':13,'out_of_page_text':outside,'record_count':len(values),'no_p_values_or_CIs_generated':True},indent=2),encoding='utf-8')
print(json.dumps({'native_mm':[W,H],'records':len(values),'bounds':outside}))
