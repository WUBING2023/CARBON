"""Reference-faithful QKI phenotype landscape; drawing frozen data only."""
from pathlib import Path
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[k]='2'
import json, time, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
from matplotlib.text import Text

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent.parent/'order_candidates'/'B2'
START=time.perf_counter()
W,H,MM=85,165,1/25.4
FULL,SEQ,K,HG,ASSAY,INK='#E64B35','#4DBBD5','#3C5488','#F39B7F','#8491B4','#252A30'
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.labelsize':6.5,
    'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.5,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
def axmm(fig,x,y,w,h):return fig.add_axes([x/W,y/H,w/W,h/H])
def txt(fig,x,y,s,**kw):return fig.text(x/W,y/H,s,**kw)
def clean(a):
    a.spines['top'].set_visible(False);a.spines['right'].set_visible(False)
    a.tick_params(length=1.7,width=.4,pad=1.2)

p=np.load(SOURCE/'SOURCE_PACKET.npz')
words=[str(s).replace('T','U') for s in p['QKI_words']]
z=np.asarray(p['QKI_external_Z5'],dtype=float)
f=np.asarray(p['QKI_full_delta_log_count'],dtype=float)
s=np.asarray(p['QKI_sequence_delta_log_count'],dtype=float)
edges=pd.read_csv(SOURCE/'QKI_ALL_48_ADJACENT_EDGES.csv',float_precision='round_trip')
assert f.shape==s.shape==(32,30) and len(edges)==48 and len(set(words))==30
assert list(p['QKI_cells'])==['HepG2']*16+['K562']*16
def rank(w):return sum('ACU'.index(w[i])>'ACU'.index(w[j]) for i in range(5) for j in range(i+1,5))
order=sorted(words,key=lambda w:(rank(w),w));ids={w:i+1 for i,w in enumerate(order)}
idx={w:words.index(w) for w in words}
literal_edges=set()
for w in words:
    assert sorted(w)==list('AACCU')
    for j in range(4):
        if w[j]!=w[j+1]:literal_edges.add(tuple(sorted((w,w[:j]+w[j+1]+w[j]+w[j+2:]))))
assert len(literal_edges)==48
seen=set()
for e in edges.itertuples(index=False):
    lo,hi=idx[e.lower_assay_word],idx[e.higher_assay_word]
    seen.add(tuple(sorted((e.lower_assay_word,e.higher_assay_word))))
    assert np.isclose(z[hi]-z[lo],e.external_delta_Z5,rtol=0,atol=1e-12)
    assert np.isclose((f[:,hi]-f[:,lo]).mean(),e.Full_mean_delta_log_count,rtol=0,atol=1e-12)
    assert np.isclose((s[:,hi]-s[:,lo]).mean(),e.Sequence_mean_delta_log_count,rtol=0,atol=1e-12)
assert seen==literal_edges
coords={};project_x={}
for r in range(9):
    group=[w for w in order if rank(w)==r]
    for j,w in enumerate(group):
        off=j-(len(group)-1)/2
        coords[w]=(44+11.3*off+.8*(r-4),112.5+2.2*r)
        project_x[w]=r+.10*off
lift=8.8
raised={w:(coords[w][0],coords[w][1]+lift*z[idx[w]]) for w in words}
zmap=LinearSegmentedColormap.from_list('external',['#F9FAFC',ASSAY])
znorm=Normalize(0,2.7)
response=LinearSegmentedColormap.from_list('absolute_response',[K,'#FFFFFF',FULL])
rnorm=Normalize(-1.2,1.2)
fig=plt.figure(figsize=(W*MM,H*MM))
txt(fig,1.6,163.2,'a',fontsize=9,fontweight='bold',va='top')
txt(fig,8,160.8,'QKI  A2C2U',fontsize=7.5)
txt(fig,82,160.8,'30 words · 48 exchanges',ha='right')
land=axmm(fig,0,109.5,85,48.5);land.set_xlim(0,85);land.set_ylim(109.5,158);land.axis('off')
for e in edges.itertuples(index=False):
    lo,hi=e.lower_assay_word,e.higher_assay_word
    land.add_patch(FancyArrowPatch(coords[lo],coords[hi],arrowstyle='-',
        connectionstyle='arc3,rad=.08',color='#E0E5EA',lw=.35,zorder=0))
for w in words:
    x,y=coords[w];xr,yr=raised[w]
    land.plot([x,xr],[y,yr],color='#D5DBE2',lw=.35,zorder=1)
    land.plot(x,y,'o',ms=1.8,color='#DFE4EB',mec='none',zorder=2)
for e in edges.itertuples(index=False):
    lo,hi=e.lower_assay_word,e.higher_assay_word
    anchor='CUAAC' in (lo,hi)
    land.add_patch(FancyArrowPatch(raised[lo],raised[hi],arrowstyle='-',
        connectionstyle='arc3,rad=.10',color=INK if anchor else '#BAC3D0',
        lw=.95 if anchor else .55,alpha=1 if anchor else .68,zorder=3))
for w in words:
    x,y=raised[w];zz=z[idx[w]]
    land.add_patch(Circle((x,y),.95+.90*zz/2.7,facecolor=zmap(znorm(zz)),
        edgecolor=INK if w=='CUAAC' else '#5E6A7D',lw=.9 if w=='CUAAC' else .45,zorder=4))
# Source Fig3A uses short genotype identifiers. Every current node has its exact
# 01–30 identifier; panel c is the complete literal decoder. Only labels move.
# Candidate positions are bounded inside panel a and avoid other literal labels.
boxes=[];labels=[]
candidates=[(0,3.2),(0,-3.2),(4.9,0),(-4.9,0),(4.2,3.1),(-4.2,3.1),
    (4.2,-3.1),(-4.2,-3.1),(0,5.8),(0,-5.8),(6.5,3.7),(-6.5,3.7),
    (6.5,-3.7),(-6.5,-3.7),(0,8.2),(0,-8.2),(7,6),(-7,6),(7,-6),(-7,-6)]
candidates += [(dx,dy) for radius in (10,13,16) for dx,dy in [(radius,0),(-radius,0),(0,radius),(0,-radius),(radius*.7,radius*.7),(-radius*.7,radius*.7),(radius*.7,-radius*.7),(-radius*.7,-radius*.7)]]
for w in sorted(words,key=lambda w:(-z[idx[w]],w)):
    x,y=raised[w];best=None
    for dx,dy in candidates:
        bx,by=x+dx,y+dy
        if not (10<bx<79 and 110<by<157):continue
        box=(bx-1.9,by-1.4,bx+1.9,by+1.4)
        overlap=sum(max(0,min(box[2],b[2])-max(box[0],b[0]))*
            max(0,min(box[3],b[3])-max(box[1],b[1])) for b in boxes)
        nodehit=sum(max(0,2.5-((bx-nx)**2+(by-ny)**2)**.5) for ww,(nx,ny) in raised.items() if ww!=w)
        score=100*overlap+15*nodehit+np.hypot(dx,dy)
        if best is None or score<best[0]:best=(score,bx,by,box)
    _,bx,by,box=best;boxes.append(box)
    land.plot([x,bx],[y,by],color='#A8B1BE',lw=.28,zorder=5)
    t=land.text(bx,by,f'{ids[w]:02}',ha='center',va='center',fontfamily='Arial',
        fontsize=6.5,color=INK,zorder=7,bbox={'facecolor':'white','edgecolor':'none','pad':.35})
    labels.append(t)
xanchor,yanchor=raised['CUAAC']
land.annotate('CUAAC (22)',xy=(xanchor,yanchor),xytext=(xanchor,min(156,yanchor+8)),
    fontsize=6.5,ha='center',arrowprops={'arrowstyle':'-','color':INK,'lw':.5},zorder=8)
# Geometric z-axis has exactly the same 8.8 mm per Z5 lift as every node.
land.plot([5.4,5.4],[113,113+2.7*lift],color=INK,lw=.6)
for zz in (0,1,2):
    yy=113+zz*lift;land.plot([4.4,6.4],[yy,yy],color=INK,lw=.5)
    land.text(3.7,yy,str(zz),ha='right',va='center',fontsize=6.5)
land.text(5.4,139.2,'Z5',ha='center',fontsize=6.5)
txt(fig,9,109.1,'External RNAcompete',fontsize=6.5)
txt(fig,81,109.1,'Height and saturation: Z5',ha='right')

txt(fig,1.6,107.3,'b',fontsize=9,fontweight='bold',va='top')
means=[f.mean(0),s.mean(0)]
for j,(data,color,name,count) in enumerate(zip(means,[FULL,SEQ],['Full','Sequence'],[31,30])):
    a=axmm(fig,10+39*j,87.8,32,15.6)
    for e in edges.itertuples(index=False):
        lo,hi=e.lower_assay_word,e.higher_assay_word
        a.plot([project_x[lo],project_x[hi]],[data[idx[lo]],data[idx[hi]]],
            color=color,alpha=.30,lw=.55,zorder=1)
    a.scatter([project_x[w] for w in order],[data[idx[w]] for w in order],
        s=9,facecolor=color,edgecolor='white',linewidth=.3,zorder=3)
    # Same literal anchor as a; annotation is orientation only, not favorable selection.
    w='CUAAC';a.annotate(str(ids[w]),(project_x[w],data[idx[w]]),xytext=(4,-11),
        textcoords='offset points',fontsize=6.5,
        arrowprops={'arrowstyle':'-','color':'#818A95','lw':.4})
    a.axhline(0,color='#858D95',lw=.5,ls=(0,(2,2)),zorder=0)
    a.set_xlim(-.45,8.45);a.set_ylim(-.07,.31)
    a.set_xticks([0,4,8]);a.set_yticks([0,.15,.3]);clean(a)
    if j==1:a.set_yticklabels([])
    a.set_title(f'{name}   {count}/48 positive',fontsize=6.5,color=color,pad=2)
txt(fig,2.0,95.6,'Mean Δ log count',rotation=90,ha='center',va='center')
txt(fig,46,82.9,'Literal order rank (A<C<U inversions)',ha='center')

txt(fig,1.6,80.3,'c',fontsize=9,fontweight='bold',va='top')
txt(fig,8,76.9,'Word',ha='left')
mx0,mw,my,mh=22,60,10.5,65.2
matrix=np.concatenate([f[:,[idx[w] for w in order]].T,s[:,[idx[w] for w in order]].T],axis=1)
a=axmm(fig,mx0,my,mw,mh)
a.imshow(matrix,cmap=response,norm=rnorm,aspect='auto',interpolation='nearest',extent=(0,64,30,0))
a.set_xticks([]);a.set_yticks([])
for v in (16,32,48):a.axvline(v-.01,color='white',lw=1.3 if v==32 else .8)
for v in range(1,30):a.axhline(v,color='white',lw=.35)
for sp in a.spines.values():sp.set_linewidth(.4);sp.set_edgecolor('#BDC4CB')
for i,w in enumerate(order):
    yy=my+mh-(i+.5)*mh/30
    txt(fig,3.1,yy,f'{ids[w]:02}',va='center',color='#6D7783')
    txt(fig,8.4,yy,w,va='center',fontfamily='DejaVu Sans Mono')
for j,(start,color,cell) in enumerate([(0,HG,'H'),(16,K,'K'),(32,HG,'H'),(48,K,'K')]):
    x=mx0+mw*start/64
    fig.add_artist(Rectangle((x/W,76.0/H),mw*16/64/W,.65/H,transform=fig.transFigure,
        facecolor=color,edgecolor='none'))
    txt(fig,x+mw*8/64,77.0,cell,ha='center')
txt(fig,mx0+mw*.25,80.0,'Full',ha='center',color=FULL)
txt(fig,mx0+mw*.75,80.0,'Sequence',ha='center',color=SEQ)
for start in (0,16,32,48):
    txt(fig,mx0+mw*(start+8)/64,8.3,'1–16',ha='center')
ca=axmm(fig,8.5,4.3,25,1.45)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=rnorm,cmap=response),cax=ca,orientation='horizontal',ticks=[-1.2,0,1.2])
cb.ax.tick_params(labelsize=6.5,width=.4,length=1.5,pad=1);cb.outline.set_linewidth(.4)
txt(fig,40,4.9,'Δ log count vs native',va='center')


nodes=pd.DataFrame([{'node_id':ids[w],'word':w,'literal_inversion_rank':rank(w),
    'external_Z5':z[idx[w]],'base_x_mm':coords[w][0],'base_y_mm':coords[w][1],
    'raised_y_mm':raised[w][1],'Full_mean_delta_log_count':means[0][idx[w]],
    'Sequence_mean_delta_log_count':means[1][idx[w]]} for w in order])
nodes.to_csv(HERE/'ALL_30_NODE_PLOT_DATA.csv',index=False,float_format='%.17g')
display=edges.copy();display.insert(0,'exchange_id',np.arange(1,49))
display.insert(1,'lower_node_id',[ids[w] for w in display.lower_assay_word])
display.insert(2,'higher_node_id',[ids[w] for w in display.higher_assay_word])
display.to_csv(HERE/'ALL_48_EDGE_PLOT_DATA.csv',index=False,float_format='%.17g')
rows=[]
for i,w in enumerate(order):
    for arm,values in [('Full',f),('Sequence',s)]:
        for bi,v in enumerate(values[:,idx[w]]):rows.append({'node_id':ids[w],'word':w,
            'arm':arm,'context_id':str(p['QKI_context_ids'][bi]),'cell':str(p['QKI_cells'][bi]),
            'context_within_cell':bi%16+1,'delta_log_count_vs_native':v})
pd.DataFrame(rows).to_csv(HERE/'ALL_1920_RESPONSE_PLOT_DATA.csv',index=False,float_format='%.17g')
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
label_boxes=[t.get_window_extent(renderer) for t in labels]
label_overlaps=[]
for i in range(len(labels)):
    for j in range(i):
        if label_boxes[i].overlaps(label_boxes[j]):label_overlaps.append([labels[i].get_text(),labels[j].get_text()])
texts=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
outside=[]
for t in texts:
    b=t.get_window_extent(renderer)
    if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.x1+.5 or b.y1>fig.bbox.y1+.5:outside.append(t.get_text())
assert min(t.get_fontsize() for t in texts)>=6.5
assert not outside,outside
for ext in ('pdf','svg','png'):fig.savefig(HERE/f'figure.{ext}',dpi=600)
fig.savefig(HERE/'PRINT_85MM.png',dpi=180)
qa={'physical_width_mm':W,'physical_height_mm':H,'minimum_ordinary_font_pt':6.5,
    'panel_font_pt':9,'nodes':30,'literal_edges':48,'response_cells':len(rows),
    'projected_graphs':2,'edges_per_projected_graph':48,'graph_base_layout_depends_only_on_literal_words':True,
    'raised_node_height_mm_per_external_Z5':lift,'all_literal_edge_means_match_raw_arrays':True,
    'Full_positive_edge_means':int((edges.Full_mean_delta_log_count>0).sum()),
    'Sequence_positive_edge_means':int((edges.Sequence_mean_delta_log_count>0).sum()),
    'node_id_label_bbox_overlaps':label_overlaps,'all_literal_word_identities_retained_in_panel_c':True,'text_outside_canvas':outside,
    'raw_response_range':[float(matrix.min()),float(matrix.max())],
    'common_response_scale':[-1.2,1.2],'source_packet_sha256':hashlib.sha256((SOURCE/'SOURCE_PACKET.npz').read_bytes()).hexdigest(),
    'thread_limit':2,'new_inference':False,'new_training':False,'new_bootstrap':False,
    'elapsed_seconds':time.perf_counter()-START}
(HERE/'RENDER_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf8')
print(json.dumps(qa,indent=2))
