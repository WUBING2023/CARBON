"""85-mm native context figures. Uses accepted coordinate grammar and frozen arrays."""
from pathlib import Path
import json, shutil, hashlib
import fitz
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, FuncNorm, PowerNorm, TwoSlopeNorm, Normalize, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, ConnectionPatch
from matplotlib.backends.backend_pdf import PdfPages

O=Path(__file__).resolve().parent
D=json.loads((O/'frozen_inputs/fig4/SOURCE_MAP.json').read_text())
E=json.loads((O/'frozen_inputs/fig4/native_examples_all.json').read_text())
Z={k:v for k,v in np.load(O/'frozen_inputs/fig4/display_packet.npz').items()}
delta=np.load(O/'frozen_inputs/fig4/gene_adjusted_error_difference.npy')
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.labelsize':6.5,'xtick.labelsize':6.5,
                    'ytick.labelsize':6.5,'axes.linewidth':.45,'pdf.fonttype':42,'svg.fonttype':'none'})
col=dict(observed='#333333',full='#E64B35',nogru='#3C5488',local='#4DBBD5',input='#A9B0B6',phi='#8491B4')
basecol={'A':'#3C5488','C':'#00A087','G':'#7E6148','T':'#E64B35','U':'#E64B35','N':'#9A9A9A'}
W=85

def setup(H):
    F=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
    def axmm(x,y,w,h):return F.add_axes([x/W,y/H,w/W,h/H])
    def t(x,y,s,**kw):return F.text(x/W,y/H,s,**kw)
    return F,axmm,t

def native(F,axmm,t,H,gene,head,letter):
    e=E[gene];r=e['pair_row'];off=head-118.2;x,w=14.,68.
    pos=np.r_[np.arange(11),np.arange(15,26)]
    t(1,head+.7,letter,fontweight='bold',fontsize=9,va='top')
    t(x,head,gene,fontstyle='italic',fontweight='bold',fontsize=7,va='center')
    t(x+max(12,len(gene)*1.9),head,e['cell'],fontsize=6.5,va='center',color='#555')
    t(82,head,'0–100 nt',fontsize=6.5,ha='right',va='center',color='#666')
    def simple(a):
        a.set_yticks([]);a.spines[['left','top','right']].set_visible(False);a.tick_params(length=1.5,pad=1.1)
    overview=axmm(x,112.5+off,w,3.5);yy=np.array(e['observed_signal']);inp=np.array(e['input_signal'])
    for c in [e['c1'],e['c2']]:overview.axvspan(c-5.5,c+5.5,color='#E8EDEE',lw=0)
    overview.fill_between(np.arange(101),0,yy,color=col['observed'],alpha=.55,lw=0)
    overview.plot(np.arange(101),yy,color=col['observed'],lw=.35)
    overview.plot(np.arange(101),inp,color=col['input'],lw=.7)
    top=max(float(yy.max()),float(inp.max()));overview.set(xlim=(-.5,100.5),ylim=(0,top*1.12),xticks=[]);simple(overview)
    overview.set_yticks([top],[f'{top:g}']);t(9.5,114.1+off,'Reads',ha='right',va='center',color='#666')
    bar=axmm(x,111.2+off,w,.65)
    bar.imshow(np.array([to_rgb(basecol[b]) for b in e['native_sequence']])[None,:,:],aspect='auto',interpolation='nearest',extent=(-.5,100.5,0,1));bar.axis('off')
    for j,c in enumerate([e['c1'],e['c2']]):
        target=x+w*((5.5 if j==0 else 20.5)+.5)/26
        F.add_artist(ConnectionPatch((c+.5,0),(target/W,(110.2+off)/H),coordsA=bar.transData,coordsB=F.transFigure,color='#A7AFB5',lw=.5,ls=(0,(2,2))))
    for left,right in [(x,x+w*11/26),(x+w*15/26,x+w)]:
        F.add_artist(Line2D([left/W,right/W],[(110.2+off)/H]*2,transform=F.transFigure,color='#B2BBC1',lw=.4))
    for xx,b in zip(pos,Z['pair_sequence'][r]):
        core=2<=xx<=8 or 17<=xx<=23
        t(x+w*(xx+.5)/26,108.6+off,str(b),fontfamily='DejaVu Sans Mono',fontsize=6.5,fontweight='bold' if core else 'normal',color=basecol[str(b)],ha='center',va='center')
    for key,y,h in [('raw_observed',104.4,2.1),('raw_input',102.,1.),('prediction',96.,4.6),('native_phi',91.,3.5)]:
        a=axmm(x,y+off,w,h)
        for left in [3.5,18.5]:a.axvspan(left,left+4,color='#EDF1F3',lw=0)
        if key in ['raw_observed','raw_input']:
            val=Z['pair_'+key][r]
            for sl in [slice(0,11),slice(11,22)]:a.fill_between(pos[sl],0,val[sl],color=col['observed'] if key=='raw_observed' else col['input'],alpha=.55,lw=.25,edgecolor=col['observed'] if key=='raw_observed' else col['input'])
            a.set_ylim(0,max(1,float(val.max()))*1.12)
            if key=='raw_observed':a.set_yticks([float(val.max())],[f'{val.max():g}'])
        elif key=='prediction':
            top=max(Z[f'pair_native_{m}'][r].max() for m in ['full','nogru','local'])
            for m,ls in [('full','-'),('nogru','--'),('local',':')]:
                for sl in [slice(0,11),slice(11,22)]:a.plot(pos[sl],Z[f'pair_native_{m}'][r][sl],color=col[m],ls=ls,lw=.8)
            a.set_ylim(0,top*1.16);a.set_yticks([top],[f'{top:.2g}'])
        else:
            phi_values=Z['pair_native_phi'][r]
            assert np.isfinite(phi_values).all() and phi_values.min() >= -3 and phi_values.max() <= 6
            a.set_ylim(-3,6);a.axhline(0,color='#D7DDE0',lw=.35)
            for sl in [slice(0,11),slice(11,22)]:a.plot(pos[sl],Z['pair_native_phi'][r][sl],color=col['phi'],lw=.65)
            a.set_yticks([-2,5],['−2','5'])
        a.set_xlim(-.5,25.5);ticks,labels=a.get_yticks(),[v.get_text() for v in a.get_yticklabels()];simple(a)
        if key!='raw_input':a.set_yticks(ticks,labels)
        a.spines['bottom'].set_visible(False)
        for left,right in [(-.5,10.5),(14.5,25.5)]:a.plot([left,right],[0,0],transform=a.get_xaxis_transform(),lw=.45,color='#333',clip_on=False)
        if key=='native_phi':
            c1,c2=e['c1'],e['c2'];a.set_xticks([0,5,10,15,20,25],[str(c1-5),str(c1),str(c1+5),str(c2-5),str(c2),str(c2+5)])
        else:a.set_xticks([])
        if key=='prediction':t(2.5,y+off+h/2,'p(nt)',rotation=90,ha='center',va='center',color='#666')
        else:
            local_input=Z['pair_raw_input'][r].sum()
            keylabel={'raw_observed':'eCLIP','raw_input':'Input (0)' if local_input==0 else 'Input','native_phi':'PhyloP'}[key]
            t(9.5 if key!='raw_input' else 12,y+off+h/2,keylabel,ha='right',va='center',color='#666')
    for arm,fx in [('observed',.12),('aligned',.36),('nogru',.60),('local_pos',.84)]:
        share=e['plotted_array_recomputed_shares'][arm]
        shown=f'{share:.4f}' if share!=.5 and abs(share-.5)<.0005 else f'{share:.3f}'
        t(x+w*fx,86.5+off,shown,ha='center',va='center',color=col[{'aligned':'full','local_pos':'local'}.get(arm,arm)])
    t(12,86.5+off,'q1 / p1',ha='right',va='center',color='#666')
    for xx in [47.2,48.6]:F.add_artist(Line2D([(xx-.30)/W,(xx+.30)/W],[(90.7+off)/H,(91.3+off)/H],transform=F.transFigure,color='#8E999F',lw=.55))

def legend(F,H,y):
    handles=[Patch(facecolor=col['observed'],label='eCLIP'),Line2D([],[],color=col['input'],lw=.7,label='Input'),
             Line2D([],[],color=col['full'],lw=.8,label='Full'),Line2D([],[],color=col['nogru'],ls='--',lw=.8,label='NoGRU'),
             Line2D([],[],color=col['local'],ls=':',lw=.8,label='Local+pos')]
    F.legend(handles=handles,loc='center',bbox_to_anchor=(48/W,y/H),ncol=5,frameon=False,fontsize=6.5,handlelength=1.2,handletextpad=.3,columnspacing=.7)

def save(F,stem):
    for ext in ['pdf','svg','png']:F.savefig(O/f'{stem}.{ext}',dpi=600,facecolor='white')

palette=json.loads((O/'frozen_inputs/PALETTE_V3.json').read_text(encoding='utf8'))
col['observed']=palette['observed'];col['input']=palette['observed']
basecol={**palette['nucleotides'],'N':'#9A9A9A'}
H=139.;F,axmm,t=setup(H)
native(F,axmm,t,H,'ZMYM4',137.2,'a');native(F,axmm,t,H,'OCRL',99.7,'b');legend(F,H,102.8)
plt.rcParams.update({'font.size':5.5,'axes.labelsize':5.5,'xtick.labelsize':5.5,'ytick.labelsize':5.5})
t(1,65.2,'c',fontweight='bold',fontsize=9,va='top')
t(14,64.0,'PTBP1',fontweight='bold',fontsize=6.5,va='center')
t(82,64.0,'865 cell–gene rows',fontsize=6.5,ha='right',va='center')
t(14,61.5,'O: Observed   F: Full   G: NoGRU   L: Local+pos',fontsize=5.5,va='center')
keys=['observed','full','nogru','local'];short=['O','F','G','L']
pool=np.concatenate([Z['gene_'+key].ravel() for key in keys])
C=float(np.median(pool));U=float(np.quantile(pool,.99));gamma=.5
def forward(v):
    v=np.asarray(v)
    return np.where(v<=C,.5-.5*(np.maximum(C-v,0)/C)**gamma,.5+.5*(np.maximum(v-C,0)/(U-C))**gamma)
def inverse(v):
    v=np.asarray(v)
    return np.where(v<=.5,C-C*np.maximum(1-2*v,0)**(1/gamma),C+(U-C)*np.maximum(2*v-1,0)**(1/gamma))
norm=FuncNorm((forward,inverse),vmin=0,vmax=U)
# Apply one common endpoint lightening to every observed/model matrix and raw-share strip.
heat_lightening=.24
def lighten(c):
    return tuple((1-heat_lightening)*v+heat_lightening for v in to_rgb(c))
cm=LinearSegmentedColormap.from_list('shared_absolute_mass_light',[lighten('#3C5488'),'#FFFFFF',lighten('#E64B35')],N=257)
share_norm=TwoSlopeNorm(vmin=0,vcenter=.5,vmax=1)
shares={key:Z['gene_'+key].reshape(865,2,11).sum(2)[:,0] for key in keys}
assert all(Z['gene_'+k].shape==(865,22) for k in keys)
assert all(np.allclose(Z['gene_'+k].sum(1),1,rtol=0,atol=1e-14) for k in keys)
axes=[];mass_images=[];y0,hh=9.,45.
for support,x0 in [(0,14.),(1,34.)]:
    t(x0+8.8,58.6,f'Support {support+1}',fontsize=5.8,ha='center',va='center')
    for j,(key,label) in enumerate(zip(keys,short)):
        x=x0+j*4.4;w=4.1
        a=axmm(x,y0,w,hh)
        val=Z['gene_'+key][:,support*11:(support+1)*11]
        im=a.imshow(val,cmap=cm,norm=norm,aspect='auto',origin='upper',interpolation='nearest',rasterized=True)
        a.axis('off');a.axhline(411.5,color='white',lw=.65)
        t(x+w/2,55.5,label,color=col[key],fontsize=5.5,ha='center',va='center')
        axes.append(a);mass_images.append(im)

# Raw input fraction stays a separate grey quantity with its own 0-1 key.
inp=Z['gene_input'].reshape(865,2,11).sum(2)
icm=LinearSegmentedColormap.from_list('input',['#FFFFFF',lighten('#A9B0B6'),lighten('#333333')])
a=axmm(52.1,y0,1.5,hh)
input_im=a.imshow(inp,cmap=icm,norm=PowerNorm(gamma=.5,vmin=0,vmax=1),aspect='auto',origin='upper',interpolation='nearest',rasterized=True)
a.axis('off');a.axhline(411.5,color='white',lw=.65);a.axvline(.5,color='white',lw=.3);axes.append(a)
t(52.85,58.1,'Input',fontsize=5.5,rotation=90,ha='center',va='center',color='#777')

# These two stripes are direct sums of the leftmost source arrays, not errors.
share_im=None
for key,x,label in [('observed',55.,'O'),('full',57.2,'F')]:
    a=axmm(x,y0,1.7,hh)
    share_im=a.imshow(shares[key][:,None],cmap=cm,norm=share_norm,aspect='auto',origin='upper',interpolation='nearest',rasterized=True)
    a.axis('off');a.axhline(411.5,color='white',lw=.65);axes.append(a)
    t(x+.85,55.5,label,fontsize=5.5,color=col[key],ha='center',va='center')
t(56.95,58.4,'Raw\nq / p',fontsize=5.5,ha='center',va='center',linespacing=.85)

for cell,lo,hi in [('K562',0,412),('HepG2',412,865)]:
    t(4,y0+hh*(1-(lo+hi)/2/865),cell,fontsize=5.8,rotation=90,ha='center',va='center',color='#777')
for gene,letter in [('ZMYM4','a'),('OCRL','b')]:
    row=E[gene]['gene_row'];yy=y0+hh*(1-(row+.5)/865)
    t(12.5,yy,letter,fontsize=6.5,ha='right',va='center',fontweight='bold')
    for a in axes:a.axhline(row,color='#27364B',lw=.25,alpha=.7)

# Independent cell-specific x/y plots, not heatmap row coordinates.
scatter_records=[]
for cell,lo,hi,yy in [('K562',0,412,35.8),('HepG2',412,865,9.8)]:
    a=axmm(65.0,yy,17.0,17.0)
    a.plot([0,1],[0,1],color='#92999F',lw=.5,ls=(0,(2,2)),zorder=0)
    a.scatter(shares['observed'][lo:hi],shares['full'][lo:hi],s=1.8,color=col['full'],lw=0,alpha=.6,rasterized=True,zorder=2)
    a.set(xlim=(-.035,1.035),ylim=(-.035,1.035),xticks=[0,.5,1],yticks=[0,.5,1],xlabel='Observed q',ylabel='Full p')
    a.set_xticklabels(['0','0.5','1']);a.set_yticklabels(['0','0.5','1'])
    a.xaxis.label.set_color(col['observed']);a.yaxis.label.set_color(col['full'])
    a.xaxis.labelpad=.7;a.yaxis.labelpad=.8;a.tick_params(length=1,pad=.7,labelsize=5.5)
    a.spines[['top','right']].set_visible(False)
    t(73.5,yy+18.2,cell,fontsize=5.5,ha='center',va='center')
    for gene,letter in [('ZMYM4','a'),('OCRL','b')]:
        row=E[gene]['gene_row']
        if lo<=row<hi:
            q,p=shares['observed'][row],shares['full'][row]
            a.scatter([q],[p],s=10,fc='white',ec=col['observed'],lw=.6,zorder=6)
            a.annotate(letter,(q,p),xytext=(-6,-7) if q>.8 else (6,-7),textcoords='offset points',fontsize=6,fontweight='bold',ha='center',va='center')
    scatter_records.append({'cell':cell,'rows':hi-lo,'first_row':lo,'last_row':hi-1,'q_sha256':hashlib.sha256(shares['observed'][lo:hi].tobytes()).hexdigest(),'p_sha256':hashlib.sha256(shares['full'][lo:hi].tobytes()).hexdigest()})

cb=F.colorbar(mass_images[0],cax=axmm(14,3.2,19,1.25),orientation='horizontal',ticks=[0,C,U],extend='max')
cb.ax.set_xticklabels(['0','0.040','0.2']);cb.ax.tick_params(length=1,pad=.4,labelsize=5.5);cb.outline.set_linewidth(.3)
t(23.5,6.5,'Mass / nt',fontsize=5.5,ha='center',va='center')
cbp=F.colorbar(share_im,cax=axmm(36,3.2,13,1.25),orientation='horizontal',ticks=[0,.5,1])
cbp.ax.tick_params(length=1,pad=.4,labelsize=5.5);cbp.outline.set_linewidth(.3)
t(42.5,6.5,'Support 1 share',fontsize=5.5,ha='center',va='center')
cbi=F.colorbar(input_im,cax=axmm(52,3.2,7,1.25),orientation='horizontal',ticks=[0,1])
cbi.ax.tick_params(length=1,pad=.4,labelsize=5.5);cbi.outline.set_linewidth(.3)
t(55.5,6.5,'Input',fontsize=5.5,ha='center',va='center')
save(F,'fig4');plt.close(F)


Dout=fitz.open(O/'fig4.pdf')
for dpi in [96,300]:Dout[0].get_pixmap(dpi=dpi).save(O/f'fig4_v32_85mm_{dpi}dpi.png')
report={'dimensions_mm':[85,139],'one_complete_figure':True,'native_geometry_values_line_colours_preserved':True,'native_fill_alpha':.55,'common_heatmap_white_blend':heat_lightening,
        'matrix_shape_each':[865,22],'all_four_arrays_preserved':True,'row_order_frozen':True,'865_rows_preserved':True,'1851_pairs_preserved':True,
        'support_grouping':'Support1/2 each11nt; O/F/G/L adjacent within each support; no column interpolation',
        'mass_mapping':{'gamma':gamma,'anchors_L_C_U':[0,C,U],'pooled_records':int(pool.size),'pooled_median_white':True,'pooled_q99_high':True,'same_mass_same_RGBA_for_all_models':True,'saturation_count':int((pool>U).sum()),'saturation_fraction':float(np.mean(pool>U)),'upper_extension':True,'scale_unchanged_from_v3':True},
        'raw_share_mapping':{'operation':'sum11nt of first support from frozen22nt matrix','anchors':[0,.5,1],'linear_common_scale':True,'observed_full_columns_all865':True,'no_error_or_adjusted_q_in_raw_atlas':True},
        'scatter':{'x':'raw observed first-support mass q','y':'Full first-support mass p','diagonal':'y=x','all865_points':True,'cell_panels':scatter_records,'not_heatmap_row_coordinates':True},
        'selected_cell_gene_rows':{gene:{'row':int(E[gene]['gene_row']),'q_raw':float(shares['observed'][E[gene]['gene_row']]),'p_full':float(shares['full'][E[gene]['gene_row']]),'a_b_are_pair_specific':True,'c_is_gene_average':True} for gene in ['ZMYM4','OCRL']},
        'no_training_or_inference_or_bootstrap_rerun':True,'no_per_model_normalization':True,'control_arrays_retained':True,'soft_prediction_not_forced_to_match_sparse_observed_saturation':True}
(O/'FIG4_DISPLAY_V32.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
print('Figure4 v32: lighter fills, same norm and data; 865 rows, 1851 pairs, failure retained.')
