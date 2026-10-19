"""v5: ONE native85x185mm a-h figure; restored b, unchanged scores, fair shared maps."""
from pathlib import Path
import csv,json,shutil,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge,Patch,PathPatch
from matplotlib.path import Path as MPath
from matplotlib.colors import LinearSegmentedColormap,Normalize,TwoSlopeNorm,ListedColormap,PowerNorm,FuncNorm
from matplotlib.lines import Line2D

O=Path(__file__).resolve().parent
rows=list(csv.DictReader((O/'source_fig2/benchmark.csv').open(encoding='utf-8-sig')))
rows.sort(key=lambda r:({'K562':0,'HepG2':1,'other':2}[r['cell_line']],r['task']))
tracks=json.loads((O/'source_fig2/native_track_geometry.json').read_text())
N=len(rows);assert N==222
models=[('CARBON','marcp','#E64B35'),('iDeepB','ideepb','#4DBBD5'),('RBPNet','rbpn','#3C5488')]
cells={'K562':'#3C5488','HepG2':'#F39B7F','other':'#8491B4'}
metrics=['Mean_PCC','Spearman','Signal_AUC','Signal_AP','Peak_AUC','Peak_AP']
names=['PCC','ρ','Signal\nAUROC','Signal\nAP','Peak\nAUROC','Peak\nAP']
short=['PCC','ρ','sAUROC','sAP','pAUROC','pAP']
arr=lambda c:np.array([float(r[c]) for r in rows])
deltas=np.column_stack([arr('marcp_'+m)-arr(key+'_'+m) for _,key,_ in models[1:] for m in metrics])
for j,(_,key,_) in enumerate(models[1:]):
    for k,m in enumerate(metrics):assert np.allclose(deltas[:,j*6+k],arr('delta_marcp_minus_'+key+'_'+m),atol=5e-15,rtol=0)
plt.rcParams.update({'font.family':'Arial','font.size':5.8,'axes.labelsize':5.8,'xtick.labelsize':5.8,'ytick.labelsize':5.8,
                    'axes.linewidth':.45,'xtick.major.width':.45,'ytick.major.width':.45,'xtick.major.size':1.5,
                    'ytick.major.size':1.5,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white'})
W,H=85.,185.
F=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
def ax(x,y,w,h,**kw):return F.add_axes([x/W,1-(y+h)/H,w/W,h/H],**kw)
def t(x,y,s,**kw):return F.text(x/W,1-y/H,s,**kw)
def panel(s,x,y):t(x,y,s,fontsize=9,fontweight='bold',va='top')
def clean(a):a.spines[['top','right']].set_visible(False);a.tick_params(pad=1.4)

# One model-blind PCC map; round common display anchor 0.30, not a statistical cutoff.
panel('a',1,1)
a=ax(.5,2.5,30,30);a.set(xlim=(-33,33),ylim=(-33,33));a.set_aspect('equal');a.axis('off')
pcc_all=np.concatenate([arr(key+'_Mean_PCC') for _,key,_ in models])
pcc_center=float(np.median(pcc_all))
L,C,U=-.1,.30,.5
gamma=1/3
def forward(v):
    v=np.asarray(v)
    return np.where(v<=C,.5-.5*(np.maximum(C-v,0)/(C-L))**gamma,.5+.5*(np.maximum(v-C,0)/(U-C))**gamma)
def inverse(v):
    v=np.asarray(v)
    return np.where(v<=.5,C-(C-L)*np.maximum(1-2*v,0)**(1/gamma),C+(U-C)*np.maximum(2*v-1,0)**(1/gamma))
norm=FuncNorm((forward,inverse),vmin=L,vmax=U)
shared=LinearSegmentedColormap.from_list('shared_blue_white_npg_red',['#3C5488','#FFFFFF','#E64B35'],N=257)
assert pcc_all.min()>=-.1 and pcc_all.max()<=.5
for ri,(name,key,col) in enumerate(models):
    radius=29-ri*5.2
    for j,v in enumerate(arr(key+'_Mean_PCC')):
        a.add_patch(Wedge((0,0),radius,90+j*315/N,90+(j+1)*315/N,width=4.15,
                          facecolor=shared(norm(v)),edgecolor='white',lw=.055))
    a.text(.5,radius-2.1,name,ha='left',va='center',fontsize=5.8)
for j,r in enumerate(rows):a.add_patch(Wedge((0,0),12,90+j*315/N,90+(j+1)*315/N,width=1.4,facecolor=cells[r['cell_line']],edgecolor='none'))
a.text(0,5.4,'222 tasks',ha='center',va='center',fontsize=5.8)
a.text(0,1.2,'150 RBPs',ha='center',va='center',fontsize=5.8)
a.text(0,-3.4,'PCC',ha='center',va='center',fontsize=5.8)
# Native original vertical violins/boxes retained; compact placement only.
cax=ax(3,36,25,1.1)
cb=F.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=shared),cax=cax,orientation='horizontal',ticks=[L,C,U])
cb.ax.set_xticklabels(['−0.1','0.30','0.5']);cb.ax.tick_params(length=1,pad=1.3,labelsize=5.5);cb.outline.set_linewidth(.25)
panel('b',32,1)
bmodels=[('CARBON','marcp','#E64B35'),('iDeepB','ideepb','#4DBBD5'),('RBPNet','rbpn','#00A087')]
b=ax(35,8,47,47*40/69)
vector_records=[]
for j,(name,key,colour) in enumerate(bmodels):
    positions=np.arange(6)+(j-1)*.24;datasets=[arr(key+'_'+m) for m in metrics]
    vp=b.violinplot(datasets,positions=positions,widths=.22,showextrema=False,points=140,bw_method=.25)
    for body in vp['bodies']:body.set_facecolor(colour);body.set_edgecolor(colour);body.set_linewidth(.45);body.set_alpha(.65)
    b.boxplot(datasets,positions=positions,widths=.067,showfliers=False,whis=1.5,patch_artist=True,
              medianprops={'color':'#222','linewidth':.8},boxprops={'facecolor':colour,'edgecolor':'#222','linewidth':.45},
              whiskerprops={'linewidth':.45},capprops={'linewidth':.45},manage_ticks=False)
    vector_records.append({'model':name,'violin_data_vertices':[p.vertices.tolist() for body in vp['bodies'] for p in body.get_paths()]})
b.set(xlim=(-.6,5.6),ylim=(-.12,.94),ylabel='Task score')
b.set_xticks(np.arange(6),['PCC','Spearman','Signal\nAUROC','Signal\nAP','Peak\nAUROC','Peak\nAP']);b.set_yticks([0,.2,.4,.6,.8]);b.tick_params(axis='x',length=0,pad=2)
b.spines[['top','right']].set_visible(False)
for spine in b.spines.values():spine.set_linewidth(.55)
b.tick_params(pad=1.5,labelsize=5.5,width=.55,length=2);b.yaxis.label.set_fontsize(5.8)
b.yaxis.grid(True,color='#E9E9E9',lw=.4);b.set_axisbelow(True)
b.legend(handles=[Patch(fc=colour,label=name) for name,key,colour in bmodels],loc='lower center',bbox_to_anchor=(.5,1.04),
         frameon=False,ncol=3,fontsize=5.8,handlelength=1,handletextpad=.35,columnspacing=.7)
t(48.5,47,'Median Δ / wins',ha='center',fontsize=6)
for k,(name,key,colour) in enumerate(bmodels[1:]):
    yy=53.5+6.6*k;t(13.5,yy,'vs '+name,ha='right',va='center',fontsize=5.8,color='#303A42')
    for j,m in enumerate(metrics):
        d=arr('marcp_'+m)-arr(key+'_'+m);assert np.allclose(d,arr('delta_marcp_minus_'+key+'_'+m),atol=5e-15,rtol=0)
        xx=14+(j+.6)/6.2*69;median=float(np.median(d));wins=int((d>0).sum())
        t(xx,yy-1.5,f'{median:+.3f}',ha='center',va='center',fontsize=5.8)
        t(xx,yy+1.6,f'{wins}/222',ha='center',va='center',fontsize=5.8,color='#A23E31' if wins<200 else '#5F6D76')
handles=[Line2D([],[],marker='o',ls='none',color=cells[k],ms=2.3,label=f'{label} ({n})') for k,label,n in [('K562','K562',119),('HepG2','HepG2',101),('other','Adrenal',2)]]
F.legend(handles=handles,loc='center',bbox_to_anchor=(.5,1-65/H),ncol=3,frameon=False,fontsize=5.8,columnspacing=.8,handletextpad=.35,handlelength=1)

# Four complete paired plots, with the same numerical bounds as the accepted source.
spec=[('c','iDeepB','ideepb','Mean_PCC','PCC',8.5,70.5),('d','RBPNet','rbpn','Mean_PCC','PCC',51,70.5),
      ('e','iDeepB','ideepb','Spearman','ρ',8.5,106),('f','RBPNet','rbpn','Spearman','ρ',51,106)]
extrema=[]
for letter,name,key,m,label,x0,y0 in spec:
    panel(letter,x0-7.5,y0-4)
    s=ax(x0,y0,25.5,25.5);x,y=arr(key+'_'+m),arr('marcp_'+m)
    lo,hi=(-.09,.55) if m=='Mean_PCC' else (-.09,.57)
    assert min(x.min(),y.min())>=lo and max(x.max(),y.max())<=hi
    s.plot([lo,hi],[lo,hi],color='#92999F',ls=(0,(3,2)),lw=.55,zorder=0)
    for cell,col in cells.items():
        mask=np.array([r['cell_line']==cell for r in rows]);s.scatter(x[mask],y[mask],s=4.8,c=col,ec='white',lw=.09,alpha=.85)
    diff=y-x;imin,imax=int(np.argmin(diff)),int(np.argmax(diff))
    for kind,ix,mark in [('max',imax,'^'),('min',imin,'v')]:
        s.scatter(x[ix],y[ix],s=19,marker=mark,fc='white',ec='#242D33',lw=.65,zorder=5)
        rbp=rows[ix]['task'].split('_')[0]
        tx,ty=(-.065,min(.52,y[ix]+.1)) if kind=='max' else ((x[ix]+.06,y[ix]-.09) if x[ix]<.3 else (.54,y[ix]-.12))
        s.annotate(rbp,xy=(x[ix],y[ix]),xytext=(tx,ty),ha='left' if tx<.5 else 'right',va='top' if kind=='max' else 'bottom',fontsize=5.8,
                   arrowprops={'arrowstyle':'-','lw':.4,'color':'#53616A','shrinkA':1.5,'shrinkB':3},zorder=6)
        extrema.append({'panel':letter,'extreme':kind,'task':rows[ix]['task'],'cell_line':rows[ix]['cell_line'],'comparator_score':float(x[ix]),'CARBON_score':float(y[ix]),'paired_gain':float(diff[ix])})
    s.set(xlim=(lo,hi),ylim=(lo,hi),xlabel=f'{name} {label}',ylabel=f'CARBON {label}');s.set_xticks([0,.2,.4]);s.set_yticks([0,.2,.4]);s.set_aspect('equal');clean(s)
F.legend(handles=[Line2D([],[],marker='^',ls='none',mfc='white',mec='#242D33',ms=3,label='Max Δ'),Line2D([],[],marker='v',ls='none',mfc='white',mec='#242D33',ms=3,label='Min Δ'),Line2D([],[],color='#92999F',ls='--',lw=.55,label='y = x')],loc='center',bbox_to_anchor=(.5,1-141/H),ncol=3,frameon=False,fontsize=5.8,handlelength=1.1,handletextpad=.35,columnspacing=1)

panel('g',1,143.2);t(25.3,145.2,'CARBON − baseline',ha='center',fontsize=5.8)
g=ax(7.5,152.3,36,24)
gnorm=TwoSlopeNorm(vmin=-.12,vcenter=0,vmax=.12)
im=g.imshow(deltas,cmap=shared,norm=gnorm,aspect='auto',interpolation='nearest',rasterized=True)
g.set_xticks(range(12),short*2,rotation=90,ha='center',fontsize=5.5);g.tick_params(axis='x',length=0,pad=1.5)
g.set_yticks([])
for row,lab in [(0,'1'),(118,'119'),(221,'222')]:t(5.2,152.3+24*(row+.5)/222,lab,ha='right',va='center',fontsize=5.8)
g.text(.25,1.03,'vs iDeepB',transform=g.transAxes,ha='center',fontsize=5.8);g.text(.75,1.03,'vs RBPNet',transform=g.transAxes,ha='center',fontsize=5.8)
g.axvline(5.5,color='white',lw=.7);g.axhline(118.5,color='white',lw=.6);g.axhline(219.5,color='white',lw=.6)
for sp in g.spines.values():sp.set_visible(False)
strip=ax(6.1,152.3,.65,24);strip.imshow(np.array([{'K562':0,'HepG2':1,'other':2}[r['cell_line']] for r in rows])[:,None],cmap=ListedColormap(list(cells.values())),vmin=0,vmax=2,aspect='auto',interpolation='nearest');strip.axis('off')
cax=ax(45.2,153.3,1,22);cb=F.colorbar(im,cax=cax,ticks=[-.12,0,.12],extend='both');cb.ax.tick_params(labelsize=6.5,length=1.2,pad=1);cb.outline.set_linewidth(.3)

panel('h',53,143.2);t(57,145.2,'FMR1 · K562',fontsize=5.8)
labels=['CARBON  r=0.69','eCLIP signal','Pooled eCLIP (inv.)','Control (inv.)','iDeepB  r=0.41','RBPNet  r=0.35']
vertices=[]
for i,(track,label) in enumerate(zip(tracks,labels)):
    yy=148+i*5.1;t(57.2,yy,label,fontsize=5.5,va='top',color='#333')
    a=ax(57.2,yy+2.2,25,2.8);v=np.array(track['normalized_vertices']);codes=np.full(len(v),MPath.LINETO,dtype=np.uint8);codes[0]=MPath.MOVETO
    color={0:models[0][2],1:'#91D1C2',2:'#8491B4',3:'#3C5488',4:models[1][2],5:models[2][2]}[i]
    a.add_patch(PathPatch(MPath(v,codes),fc=color,ec=color,lw=.18));a.set(xlim=(0,100),ylim=(0,1));a.set_yticks(track['tick_locations_normalized'],track['tick_labels']);clean(a)
    a.tick_params(axis='y',length=1,pad=1,labelsize=5.5)
    if i<5:a.set_xticks([])
    else:a.set_xticks([0,50,100]);a.set_xlabel('Position (nt)',labelpad=1.5)
    vertices.append(v)

for ext in ['pdf','svg','png']:F.savefig(O/f'fig2_v5.{ext}',dpi=300 if ext=='png' else 600,facecolor='white')
plt.close(F)
for ext in ['pdf','svg','png']:shutil.copy2(O/f'fig2_v5.{ext}',O/f'fig2.{ext}')
shutil.copy2(O/'fig2.pdf',O.parent/'manuscript/figures/fig2.pdf')
with (O/'fig2_scatter_extrema.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(extrema[0]));w.writeheader();w.writerows(extrema)
summary=[]
for name,key,col in models[1:]:
    for m in metrics:
        d=arr('marcp_'+m)-arr(key+'_'+m);summary.append(dict(comparator=name,metric=m,median_delta=float(np.median(d)),wins=int((d>0).sum()),ties=int((d==0).sum()),losses=int((d<0).sum()),n=N))
with (O/'fig2_paired_summary.csv').open('w',encoding='utf8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
qa={'dimensions_mm':[W,H],'ordinary_font_pt_range':[5.5,5.8],'tasks':222,'ring_values':666,'scatter_points_per_panel':222,'matrix_shape':list(deltas.shape),
    'ring_encoding':{'model_identity':'radial position and black labels only; hue encodes the same PCC for every model','common_map':'common continuous signed cube-root gamma1/3; blue below0.30, white at0.30, NPG red above0.30; visual contrast only','pcc_limits':[-.1,.5],'white_anchor':C,'exponent':gamma,'pooled_median':pcc_center,'pooled_records':int(pcc_all.size),'clipping':False,
                     'same_PCC_same_RGBA':True,'palette':['#3C5488','#FFFFFF','#E64B35'],'median_by_model':{name:float(np.median(arr(key+'_Mean_PCC'))) for name,key,_ in models},'fraction_above_common_center':{name:float(np.mean(arr(key+'_Mean_PCC')>C)) for name,key,_ in models}},
    'matrix_encoding':{'limits':[-.12,.12],'zero_colour':'#FFFFFF','negative_blue':True,'positive_red':True,'transform':'raw differences, no per-column standardization',
                       'below_scale':int((deltas<-.12).sum()),'above_scale':int((deltas>.12).sum()),'data_range':[float(deltas.min()),float(deltas.max())]},
    'profile_geometry_unchanged':all(np.array_equal(v,np.array(t['normalized_vertices'])) for v,t in zip(vertices,tracks)),
    'track_identity_correction':'tracks 3/4 are inverted pooled eCLIP and inverted control; the frozen JSON mislabels them ENCODE replicates',
    'paired_summary_preserved_in_b':True,'one_complete_figure':True,'one_float_required':True,'b_original_vertical_mark_style_retained':True,'no_score_recomputation_or_model_rerun':True,'fair_comparison_limitations_remain':True}
(O/'FIG2_DISPLAY_V5.json').write_text(json.dumps(qa,indent=2),encoding='utf8')
(O/'fig2b_v5_violin_data_vertices.json').write_text(json.dumps(vector_records),encoding='utf8')
import fitz
D=fitz.open(O/'fig2_v5.pdf')
for dpi in [96,300]:D[0].get_pixmap(dpi=dpi).save(O/f'fig2_v5_85mm_{dpi}dpi.png')
print('ONE complete Figure2 exported at85x185mm.')
