"""Reference-guided full native context -> literal support zoom -> full cohort.

No scientific endpoint or native packet change. Uniform nucleotide identity only.
"""
from pathlib import Path
import json,hashlib,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,PowerNorm,to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Patch

O=Path(__file__).resolve().parent
D=json.loads((O/'SOURCE_MAP.json').read_text())
Z={k:v for k,v in np.load(O/'display_packet.npz').items()}
W,H=85.,124.
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.labelsize':6.5,
    'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.45,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
F=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
col=dict(observed='#19499C',full='#E64B35',nogru='#3C5488',local='#4DBBD5',input='#A9B0B6',phi='#8491B4')
basecol={'A':'#6BB72D','C':'#19499C','G':'#F4A016','T':'#E6211A','U':'#E6211A','N':'#9A9A9A'}
def axmm(x,y,w,h):return F.add_axes([x/W,y/H,w/W,h/H])
def t(x,y,s,**kw):return F.text(x/W,y/H,s,**kw)
def simple(ax):
    ax.set_yticks([]);ax.spines[['left','top','right']].set_visible(False);ax.tick_params(length=1.5,pad=1.1)
pos=np.r_[np.arange(11),np.arange(15,26)]
slice_checks={}
for gene,head,letter in [('LEPR',122.2,'a'),('OCRL',84.7,'b')]:
    e=D['examples'][gene];r=e['pair_row'];off=head-118.2;x,w=14.,68.
    fullkeys={'observed_signal':'raw_observed','input_signal':'raw_input','full_profile':'native_full',
              'nogru_profile':'native_nogru','local_position_profile':'native_local','phi':'native_phi'}
    idx=np.array(e['support_positions'],dtype=int)
    errors={k:float(np.max(np.abs(np.array(e[k])[idx]-Z['pair_'+v][r]))) for k,v in fullkeys.items()}
    assert all(len(e[k])==101 for k in fullkeys) and max(errors.values())==0
    assert ''.join(np.array(list(e['native_sequence']))[idx])==''.join(Z['pair_sequence'][r])
    slice_checks[gene]=dict(lengths={k:len(e[k]) for k in fullkeys},errors=errors,
                            support_positions=idx.tolist(),gene_row=e['gene_row'],gene_row_pairs=e['gene_row_pairs'])
    t(1,head+.7,letter,fontweight='bold',fontsize=9,va='top')
    t(x,head,gene,fontstyle='italic',fontweight='bold',fontsize=7,va='center')
    t(x+10,head,e['cell'],fontsize=6.5,va='center',color='#555555')
    t(82,head,'0–100 nt',fontsize=6.5,ha='right',va='center',color='#666666')
    # Whole-window measured context. Both stored arrays are literal counts.
    overview=axmm(x,112.5+off,w,3.5)
    yy=np.array(e['observed_signal']);inp=np.array(e['input_signal'])
    for c in [e['c1'],e['c2']]:overview.axvspan(c-5.5,c+5.5,color='#E8EDEE',lw=0)
    overview.fill_between(np.arange(101),0,yy,color=col['observed'],lw=0)
    overview.plot(np.arange(101),inp,color='#A9B0B6',lw=.7)
    overview.set(xlim=(-.5,100.5),ylim=(0,max(yy.max(),inp.max())*1.12),xticks=[]);simple(overview)
    t(12,114.1+off,'Reads',ha='right',va='center',fontsize=6.5,color='#666666')
    overview.text(.012,.94,f'{max(yy.max(),inp.max()):g}',transform=overview.transAxes,color='#666666',fontsize=6.5,va='top')
    bar=axmm(x,111.2+off,w,.65)
    base=np.array([to_rgb(basecol[b]) for b in e['native_sequence']])[None,:,:]
    bar.imshow(base,aspect='auto',interpolation='nearest',extent=(-.5,100.5,0,1));bar.axis('off')
    # Leader centers use actual four-base CTCT midpoint: floor center +0.5.
    for j,c in enumerate([e['c1'],e['c2']]):
        target=x+w*((5.5 if j==0 else 20.5)+.5)/26
        F.add_artist(ConnectionPatch((c+.5,0),(target/W,(110.2+off)/H),
            coordsA=bar.transData,coordsB=F.transFigure,color='#A7AFB5',lw=.5,ls=(0,(2,2))))
        for left,right in [(x+w*(0/26),x+w*(11/26)),(x+w*(15/26),x+w*(26/26))]:
            F.add_artist(Line2D([left/W,right/W],[(110.2+off)/H,(110.2+off)/H],transform=F.transFigure,color='#B2BBC1',lw=.4))
    for xx,b in zip(pos,Z['pair_sequence'][r]):
        core=4<=xx<=7 or 19<=xx<=22
        t(x+w*(xx+.5)/26,108.6+off,str(b),fontfamily='DejaVu Sans Mono',fontsize=6.5,
          fontweight='bold' if core else 'normal',color=basecol[str(b)],ha='center',va='center')
    tracks=[('raw_observed',104.4,2.1),('raw_input',102.,1.),('prediction',96.,4.6),('native_phi',91.,3.5)]
    for key,y,h in tracks:
        ax=axmm(x,y+off,w,h)
        for left in [3.5,18.5]:ax.axvspan(left,left+4,color='#EDF1F3',lw=0)
        if key in ['raw_observed','raw_input']:
            val=Z['pair_'+key][r]
            for sl in [slice(0,11),slice(11,22)]:ax.fill_between(pos[sl],0,val[sl],color=col['observed'] if key=='raw_observed' else col['input'],lw=0)
            ax.set_ylim(0,max(1,float(val.max()))*1.12)
            if key=='raw_observed':ax.set_yticks([float(val.max())],[f'{val.max():g}'])
        elif key=='prediction':
            top=max(Z[f'pair_native_{m}'][r].max() for m in ['full','nogru','local'])
            for m,ls in [('full','-'),('nogru','--'),('local',':')]:
                for sl in [slice(0,11),slice(11,22)]:ax.plot(pos[sl],Z[f'pair_native_{m}'][r][sl],color=col[m],ls=ls,lw=.8)
            ax.set_ylim(0,top*1.16)
            ax.set_yticks([top],[f'{top:.2g}'])
        else:
            ax.set_ylim(-3,2);ax.axhline(0,color='#D7DDE0',lw=.35)
            for sl in [slice(0,11),slice(11,22)]:ax.plot(pos[sl],Z['pair_native_phi'][r][sl],color=col['phi'],lw=.65)
            ax.set_yticks([-2,2],['−2','2'])
        ax.set_xlim(-.5,25.5)
        ytickvals,yticklabels=ax.get_yticks(),[v.get_text() for v in ax.get_yticklabels()]
        simple(ax)
        if key!='raw_input':ax.set_yticks(ytickvals,yticklabels)
        ax.spines['bottom'].set_visible(False)
        for left,right in [(-.5,10.5),(14.5,25.5)]:ax.plot([left,right],[0,0],transform=ax.get_xaxis_transform(),lw=.45,color='#333333',clip_on=False)
        if key=='native_phi':
            c1,c2=e['c1'],e['c2'];ax.set_xticks([0,5,10,15,20,25],
                [str(c1-5),str(c1),str(c1+5),str(c2-5),str(c2),str(c2+5)])
        else:ax.set_xticks([])
        if key=='prediction':t(2.5,y+off+h/2,'p(nt)',rotation=90,ha='center',va='center',fontsize=6.5,color='#666666')
        else:t(9.5 if key!='raw_input' else 12,y+off+h/2,{'raw_observed':'eCLIP','raw_input':'Input (0)','native_phi':'PhyloP'}[key],ha='right',va='center',fontsize=6.5,color='#666666')
    for arm,fx in [('observed',.12),('aligned',.36),('nogru',.60),('local_pos',.84)]:
        t(x+w*fx,86.5+off,f"{e['plotted_array_recomputed_shares'][arm]:.2f}",ha='center',va='center',fontsize=6.5,
          color=col[{'aligned':'full','local_pos':'local'}.get(arm,arm)])
    t(12,86.5+off,'q1',ha='right',va='center',fontsize=6.5,color='#666666')
    # The central blank is intentional broken genomic support coordinates.
    for xx in [47.2,48.6]:
        F.add_artist(Line2D([(xx-.30)/W,(xx+.30)/W],[(90.7+off)/H,(91.3+off)/H],transform=F.transFigure,color='#8E999F',lw=.55))

# Population retains every frozen gene row with identical independent ordering.
t(1,51.2,'c',fontweight='bold',fontsize=9,va='top')
t(14,50.3,'PTBP1',fontsize=6.5,fontweight='bold',va='center')
t(82,50.3,'865 cell–gene rows',fontsize=6.5,ha='right',va='center')
cm=LinearSegmentedColormap.from_list('npg_density',[(0,'#FFFFFF'),(.12,'#D8EBE9'),(.25,'#91D1C2'),(.40,'#4DBBD5'),(.65,'#3C5488'),(1,'#24395B')])
icm=LinearSegmentedColormap.from_list('input',['#FFFFFF','#A9B0B6','#333333'])
layout=[('observed','Observed',14,15.),('input','Input',30,2.5),('full','Full',34,14.8),('nogru','NoGRU',49.8,14.8),('local','Local+pos',65.6,14.8)]
images={}
for key,label,x,w in layout:
    v=Z[f'gene_{key}'];v=v.reshape(865,2,11).sum(axis=2) if key=='input' else v
    ax=axmm(x,6.5,w,37.)
    im=ax.imshow(v,aspect='auto',origin='upper',interpolation='nearest',rasterized=True,
                extent=(-.5,v.shape[1]-.5,864.5,-.5),cmap=icm if key=='input' else cm,
                norm=PowerNorm(gamma=.5 if key=='input' else .25,vmin=0,vmax=1,clip=False))
    images[key]=im;ax.axhline(411.5,color='white',lw=.6);ax.axvline(.5 if key=='input' else 10.5,color='white',lw=.6)
    ax.set(xticks=[],yticks=[],ylim=(864.5,-.5));[sp.set_visible(False) for sp in ax.spines.values()]
    t(x+w/2,47.5,label,ha='center',va='center',fontsize=6.5,color=col[key])
    if key!='input':
        for j in range(2):t(x+w*(.25+j*.5),44.9,str(j+1),ha='center',va='center',fontsize=6.5,color='#777777')
    for gene,letter in [('LEPR','a'),('OCRL','b')]:ax.plot(-.8,D['examples'][gene]['gene_row'],marker='>',ms=1.7,color='#333333',clip_on=False)
for cell,lo,hi in [('K562',0,412),('HepG2',412,865)]:t(4,6.5+37*(1-(lo+hi)/2/865),cell,rotation=90,ha='center',va='center',fontsize=6.5,color='#777777')
for gene,letter in [('LEPR','a'),('OCRL','b')]:t(12.5,6.5+37*(1-(D['examples'][gene]['gene_row']+.5)/865),letter,ha='right',va='center',fontsize=6.5,fontweight='bold')
for key,x,w in [('observed',14,36),('input',57,24)]:
    cb=F.colorbar(images[key],cax=axmm(x,2.8,w,1.),orientation='horizontal',ticks=[0,.25,1]);cb.ax.tick_params(length=1.2,pad=.2,labelsize=6.5);cb.set_ticklabels(['0','0.25','1']);cb.outline.set_linewidth(.25)
    t(x+w/2,5.,'Support fraction' if key=='observed' else 'Window input fraction',ha='center',va='center',fontsize=6.5)
handles=[Patch(facecolor=col['observed'],label='eCLIP'),Line2D([],[],color=col['input'],lw=.7,label='Input'),
         Line2D([],[],color=col['full'],lw=.8,label='Full'),Line2D([],[],color=col['nogru'],ls='--',lw=.8,label='NoGRU'),
         Line2D([],[],color=col['local'],ls=':',lw=.8,label='Local+pos')]
F.legend(handles=handles,loc='center',bbox_to_anchor=(48/W,87.8/H),ncol=5,frameon=False,
         fontsize=6.5,handlelength=1.3,handletextpad=.3,columnspacing=.8)
for ext in ['png','pdf','svg']:
    F.savefig(O/f'figure.{ext}',dpi=600,facecolor='white');shutil.copy2(O/f'figure.{ext}',O/f'FIGURE_A.{ext}')
D['reference_polish_display']={'dimensions_mm':[W,H],'source_101_to_original22_checks':slice_checks,
    'literal_base_palette':basecol,'identity_glyphs_uniform_height':True,'full_window_counts_only_no_fake_attribution':True,
    'coordinate_domain':'overview original0–100; detail two separated11nt supports atoriginal0-based coords; no intervening signal invented',
    'local_core_artist_indices':'4..7 and19..22 atdisplay0..10,15..25, exactly floor-relative−1..+2',
    'atlas_height_mm':37,'atlas_min_png_height_px':37*600/25.4,'all_865_rows_retained':True,
    'source_packet_sha256':hashlib.sha256((O/'display_packet.npz').read_bytes()).hexdigest(),
    'references':'PrismNet native identity→input/response coordinate linkage; RBPNet whole context→literal local enlargement; no transferred assays or IG'}
if 'display' in D and 'historical_display' not in D:D['historical_display']=D['display']
if 'revision_provenance' in D and 'historical_revision_provenance' not in D:D['historical_revision_provenance']=D.pop('revision_provenance')
D['display']=D['reference_polish_display']
D['revision_provenance']={'current_revision':str(O),'display_reorganized_not_pixel_identical':True,
                         'native_packet_unchanged':True,'source_packet_sha256':D['display']['source_packet_sha256'],
                         'current_png_sha256':hashlib.sha256((O/'figure.png').read_bytes()).hexdigest()}
(O/'SOURCE_MAP.json').write_text(json.dumps(D,indent=2),encoding='utf-8')
(O/'DISPLAY_MAPPING.json').write_text(json.dumps(D['reference_polish_display'],indent=2),encoding='utf-8')
print(json.dumps({'png':str(O/'figure.png'),'size_mm':[W,H],'max_101_to22_error':0,'atlas_rows':865}))
