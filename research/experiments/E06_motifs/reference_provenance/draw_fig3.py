"""Recompose the retained evidence at 85 mm, without rerunning scientific analyses.

Execute with D:/Anaconda/envs/rbp/python.exe -X utf8.
Editable outputs: SVG plus native vector PDF; source matrices/geometry kept locally.
"""
from pathlib import Path
import json, csv, shutil, math
import sys
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Wedge
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

W = Path(__file__).resolve().parent
PAPER = W.parent
MOTIF = Path('M:/RBP/26年数据RBP/可解释性分析/Motif')
FIGA = Path('M:/RBP/26年数据RBP/图/figure5/figureA_data')
SDIR = PAPER/'manuscript/source_data/S1_S14_original'
plt.rcParams.update({'font.family':'Arial','font.size':6,'pdf.fonttype':42,'svg.fonttype':'none',
                     'axes.linewidth':.4,'axes.labelsize':6,'xtick.labelsize':6,'ytick.labelsize':6,
                     'xtick.major.size':1.7,'ytick.major.size':1.7,'xtick.major.width':.4,'ytick.major.width':.4})
BASES={'A':'#3C5488','C':'#00A087','G':'#7E6148','U':'#E64B35'}
SIGNAL={'eCLIP':'#91D1C2','Pred. WT':'#E64B35','Pred. mut.':'#3C5488'}
STATUS={'Yes':'#EFA493','Not':'#A5D5E6','Unknown':'#FFFFFF'}
GLYPHS={}
for b in 'ACGU':
    p=TextPath((0,0),b,size=1,prop=FontProperties(family='Arial',weight='bold'))
    bb=p.get_extents()
    GLYPHS[b]=p.transformed(Affine2D().translate(-bb.x0,-bb.y0).scale(1/bb.width,1/bb.height))

def mm_text(fig,x,y,s,fontsize=6,**kw):
    fw,fh=fig.get_size_inches()*25.4
    return fig.text(x/fw,1-y/fh,s,fontsize=fontsize,va='top',**kw)

def mm_axes(fig,x,y,width,height):
    fw,fh=fig.get_size_inches()*25.4
    return fig.add_axes([x/fw,1-(y+height)/fh,width/fw,height/fh])

def ic(ppm):
    p=np.asarray(ppm,float);p=p/p.sum(axis=1,keepdims=True)
    h=2+np.sum(p*np.log2(np.maximum(p,1e-300)),axis=1)
    return p*np.maximum(h,0)[:,None]

def logo(ax,ppm,start=0,extent=None):
    vals=ic(ppm)
    for i,v in enumerate(vals):
        y=0
        for j in np.argsort(v):
            h=float(v[j]);b='ACGU'[j]
            if h<.002:continue
            tr=Affine2D().scale(.90,h).translate(start+i-.45,y)
            ax.add_patch(PathPatch(GLYPHS[b],transform=tr+ax.transData,color=BASES[b],lw=0))
            y+=h
    ax.set_ylim(0,2.08)
    ax.set_xlim(*(extent if extent else (-.5,len(ppm)-.5)))
    ax.set_axis_off()

def aligned_mats(row,mats,compact=False):
    q=mats[row.comparison_id+'_query'];r=mats[row.comparison_id+'_reference']
    if row.Orientation=='-':r=r[::-1, :][:,[3,2,1,0]]
    rs=-int(row.Optimal_offset)
    extent=(min(0,rs)-.5,max(len(q),rs+len(r))-.5)
    if compact:
        # Main-panel display crop is centered on the stored aligned reference.
        # Full query/reference inputs, offset, orientation and p/q stay intact.
        width=max(9,len(r)+2)
        mid=rs+(len(r)-1)/2
        extent=(mid-width/2,mid+width/2)
    return q,r,rs,extent

def fmt(value):
    v=float(value)
    return f'{v:.2e}' if v<.001 else f'{v:.4f}' if v<.1 else f'{v:.3f}'

def load_atlas():
    inp=W/'source_inputs';inp.mkdir(exist_ok=True)
    for p in [SDIR/'Supplementary_Data_S10_motif_leaf_mapping.csv',FIGA/'figureA_dist_matrix.csv']:
        dest=inp/p.name
        if not dest.exists():shutil.copy2(p,dest)
    meta=pd.read_csv(inp/'Supplementary_Data_S10_motif_leaf_mapping.csv').sort_values('leaf_order_clockwise_from_top')
    dm=pd.read_csv(inp/'figureA_dist_matrix.csv',index_col=0)
    cache=W/'atlas_source_matrices.npz'
    if not cache.exists():
        arrays={}
        for row in meta.itertuples():
            with h5py.File(MOTIF/'outputs'/row.task_id/'tfmodisco/modisco_results.h5','r') as hf:
                arrays[row.task_id]=hf[row.h5_group+'/sequence'][()][int(row.crop_start_0based):int(row.crop_end_exclusive)]
        np.savez_compressed(cache,**arrays)
    arrays=np.load(cache)
    z=linkage(squareform(dm.values,checks=False),method='average')
    positions={t:i for i,t in enumerate(meta.task_id)}
    members={i:[positions[t]] for i,t in enumerate(dm.index)}
    n=len(dm)
    for j,r in enumerate(z):
        members[n+j]=members[int(r[0])]+members[int(r[1])]
        assert max(members[n+j])-min(members[n+j])+1==len(members[n+j]),'Frozen leaf order inconsistent with retained average linkage'
    return meta,dm,z,members,arrays

def atlas(ax,meta,dm,z,members,arrays,main,large=False):
    n=len(meta)
    # The center tree outer radius changes explicitly (not only track thickness).
    tree_r,root_r,label_r,logo_r,logo_length,ring_r,ring_h=(11.5,1.6,12.2,23.0,7.5,30.9,.57) if not large else (73,8,74.5,108.0,21.0,130.0,1.35)
    rmax=ring_r+6*ring_h
    theta=np.pi/2-2*np.pi*(np.arange(n)+.5)/n
    height=z[-1,2]
    node_angle={i:theta[members[i][0]] for i in range(n)}
    node_radius={i:tree_r for i in range(n)}
    for j,row in enumerate(z):
        k=n+j;a,b=int(row[0]),int(row[1]);node_angle[k]=(node_angle[a]+node_angle[b])/2
        node_radius[k]=tree_r-(tree_r-root_r)*row[2]/height
        aa=np.linspace(node_angle[a],node_angle[b],max(3,int(abs(node_angle[a]-node_angle[b])*30)))
        rr=node_radius[k]
        ax.plot(rr*np.cos(aa),rr*np.sin(aa),color='#222222',lw=.35 if not large else .45)
        for c in [a,b]:
            tc=node_angle[c];rc=node_radius[c]
            ax.plot([rr*np.cos(tc),rc*np.cos(tc)],[rr*np.sin(tc),rc*np.sin(tc)],color='#222222',lw=.35 if not large else .45)
    names=set(main.task)
    display_ids={row.task:int(row.main_order) for row in main.itertuples()}
    # Frozen display positions avoid squeezing nearby representatives together.
    main_number_xy={1:(5.0,20.5),2:(12.7,14.8),3:(-2.2,20.0),4:(17.3,9.7),5:(11.0,18.4),6:(0,-19.5)}
    for i,row in enumerate(meta.itertuples()):
        ang=theta[i]
        # Full quantitative logo ring, using the same stored 9-nt crops.
        vals=ic(arrays[row.task_id]);step=logo_length/9
        tangential=.78*(logo_r+logo_length/2)*2*np.pi/n
        for pos,v in enumerate(vals):
            base=0
            for j in np.argsort(v):
                hv=float(v[j])*tangential/2
                if hv<.0007:continue
                # Sequence progresses outward; information stack occupies the wedge's tangent.
                trans=(Affine2D().scale(step*.94,hv).translate(logo_r+pos*step,-tangential/2+base).rotate(ang))
                ax.add_patch(PathPatch(GLYPHS['ACGU'[j]],transform=trans+ax.transData,facecolor=BASES['ACGU'[j]],lw=0))
                base+=hv
        if large:
            label=row.label
            deg=np.degrees(ang)
            if np.cos(ang)>=0:rot=deg;ha='left'
            else:rot=deg+180;ha='right'
            ax.text(label_r*np.cos(ang),label_r*np.sin(ang),label,rotation=rot,rotation_mode='anchor',ha=ha,va='center',
                    fontsize=6 if large else 6,color='#222222')
        elif row.task_id in names:
            number=display_ids[row.task_id];tx,ty=main_number_xy[number]
            lx,ly=tree_r*np.cos(ang),tree_r*np.sin(ang)
            # Numbers correspond to fixed row identities in panel b, not significance.
            ax.plot([lx,lx*1.18,tx],[ly,ly*1.18,ty],color='#555555',lw=.45)
            ax.text(tx,ty,str(number),ha='center',va='center',fontsize=6.5,color='#222222',
                    bbox=dict(boxstyle='circle,pad=.13',facecolor='white',edgecolor='#AAAAAA',linewidth=.35))
        for k,col in enumerate(['splicing','stability_decay','translation','RRM','ZnF','KH']):
            value=getattr(row,col)
            good={'splicing':'Yes','stability_decay':'Yes','translation':'Yes','RRM':'RRM','ZnF':'ZnF','KH':'KH'}[col]
            status='Yes' if value==good else 'Not' if value=='Not' else 'Unknown'
            d0=90-360*(i+1)/n;d1=90-360*i/n
            ax.add_patch(Wedge((0,0),ring_r+(k+1)*ring_h,d0,d1,width=ring_h,facecolor=STATUS[status],edgecolor='#C9C9C9',lw=.12))
    ax.set_xlim(-rmax-2,rmax+2);ax.set_ylim(-rmax-2,rmax+2);ax.set_aspect('equal');ax.set_axis_off()

def glyph_path(record):
    verts=[];codes=[];prev=None
    for item in record['items']:
        if prev!=item[1]:verts.append(item[1]);codes.append(MPath.MOVETO)
        if item[0]=='l':verts.append(item[2]);codes.append(MPath.LINETO)
        elif item[0]=='c':verts.extend(item[2:]);codes.extend([MPath.CURVE4]*3)
        prev=item[-1]
    if record.get('closePath') and verts:verts.append(verts[0]);codes.append(MPath.CLOSEPOLY)
    return MPath(verts,codes)

def c_tracks(fig,geometry,start=125,step=23,large=False):
    for i,panel in enumerate(geometry['panels']):
        sy=start+i*step
        parts=panel['task'].split('_')
        label=f'{parts[0]} · {parts[1]}'
        if large:label+='   '+panel['locus']
        mm_text(fig,12.5,sy,label,fontsize=6.5 if not large else 9)
        left,width=12.5,71
        if large:left,width=22,160
        for tr,dy,hh,lim,ticks in [('ic',2.3,4.6,(0,2.15),[2] if i==2 and not large else [0,2]),('ig',8.8,5.2,[(-.7,3.8),(-.75,2.8),(-.85,1.15)][i],[[0,3],[0,2],[-.5,1]][i])]:
            if large:dy*=2;hh*=2.5
            ax=mm_axes(fig,left,sy+dy,width,hh)
            for g in panel['glyphs'][tr]:ax.add_patch(PathPatch(glyph_path(g),facecolor=BASES[g['base']],lw=0))
            ax.set_xlim(.5,50.5 if tr=='ic' else 101.5);ax.set_ylim(*lim);ax.set_yticks(ticks)
            ax.axhline(0,color='#555555',lw=.3)
            ax.tick_params(axis='x',bottom=False,labelbottom=False,pad=1.2)
            ax.tick_params(axis='y',pad=1.2,labelsize=6.5 if not large else 8)
            for s in ['top','right','bottom']:ax.spines[s].set_visible(False)
            title='IC\n(bits)' if tr=='ic' else 'IG ×\ninput'
            ax.set_ylabel(title,fontsize=6.5 if not large else 8,rotation=0,ha='right',va='center',labelpad=3)
        sdy=15.0 if not large else 33
        sh=6.6 if not large else 22
        ax=mm_axes(fig,left,sy+sdy,width,sh)
        ax.axvspan(*panel['shade_positions'],color='#4DBBD5',alpha=.18,lw=0)
        ax.axvline(panel['mutation_position'],color='#555555',ls='--',lw=.55)
        for cr in panel['curves']:
            xy=np.asarray(cr['points'])
            ax.plot(xy[:,0],xy[:,1],color=SIGNAL[cr['name']],lw=.8 if not large else 1.3,ls='--' if cr['name']=='Pred. mut.' else '-')
        ax.set_xlim(.5,101.5);ax.set_ylim(-.04,1.07);ax.set_yticks([0,1]);ax.set_xticks([1,26,51,76,101])
        ax.tick_params(pad=1.2,labelsize=6.5 if not large else 8)
        for s in ['top','right']:ax.spines[s].set_visible(False)
        ax.set_ylabel('Signal',fontsize=6.5 if not large else 8,rotation=0,ha='right',va='center',labelpad=3)
        if i<2 and not large:ax.tick_params(labelbottom=False)
        if large:ax.set_xlabel('Position in 101-nt window',fontsize=8)
        if i==2 and not large:ax.set_xlabel('Position',fontsize=6.5,labelpad=1.2)

def build_main(meta,dm,z,members,arrays,main,mats,geometry):
    fig=plt.figure(figsize=(85/25.4,199/25.4),facecolor='white')
    ax=mm_axes(fig,3.3,1.0,72,72)
    atlas(ax,meta,dm,z,members,arrays,main)
    mm_text(fig,1,1,'a',fontsize=9,fontweight='bold')
    handles=[plt.Rectangle((0,0),1,1,facecolor=STATUS[v],edgecolor='#AAAAAA',lw=.4) for v in ['Yes','Not','Unknown']]
    fig.legend(handles,['Yes','No','Unknown'],loc='upper center',bbox_to_anchor=(.5,1-70.7/199),ncol=3,frameon=False,
               fontsize=6.5,handlelength=.8,handletextpad=.5,columnspacing=1.6,borderpad=0)
    mm_text(fig,1,76.0,'b',fontsize=9,fontweight='bold')
    # Six examples, with frozen score-independent eligibility/choice; paired alignments.
    for j,row in enumerate(main.itertuples()):
        col=j%2;rr=j//2;x=3+col*42;y=79.0+rr*13.9
        cell=row.task.split('_')[-2]
        mm_text(fig,x,y,f'{j+1}  {row.cisbp_RBP_Name} · {cell}',fontsize=6.5)
        mm_text(fig,x+25.1,y,row.Target_ID.split('_')[0],fontsize=6.5)
        q,r,rs,extent=aligned_mats(row,mats,compact=True)
        logo(mm_axes(fig,x,y+3.4,18.5,4.3),q,extent=extent)
        logo(mm_axes(fig,x+22,y+3.4,18.5,4.3),r,start=rs,extent=extent)
        stats=main.iloc[j]
        mm_text(fig,x,y+8.5,f'p = {fmt(stats["p_value"])}  q = {fmt(stats["q-value"])}',fontsize=6.5)
    for col in range(2):
        x=3+col*42
        mm_text(fig,x+9.3,76.2,'CARBON',fontsize=6.5,ha='center')
        mm_text(fig,x+31.3,76.2,'CISBP',fontsize=6.5,ha='center')
    mm_text(fig,1,120.0,'c',fontsize=9,fontweight='bold')
    handles=[Line2D([0],[0],color=SIGNAL[k],lw=.9,ls='--' if k=='Pred. mut.' else '-') for k in SIGNAL]
    fig.legend(handles,list(SIGNAL),loc='upper center',bbox_to_anchor=(.61,1-119.5/199),ncol=3,frameon=False,
               fontsize=6.5,handlelength=1.15,handletextpad=.35,columnspacing=.8,borderpad=0)
    c_tracks(fig,geometry)
    fig.savefig(W/'fig3.pdf');fig.savefig(W/'fig3.svg');fig.savefig(W/'fig3.png',dpi=600)
    sizes=[float(t.get_fontsize()) for t in fig.texts]
    for aa in fig.axes:
        sizes.extend(t.get_fontsize() for t in aa.texts)
        sizes.extend(t.get_fontsize() for t in aa.get_yticklabels()+aa.get_xticklabels() if t.get_visible())
        sizes.extend([aa.yaxis.label.get_fontsize(),aa.xaxis.label.get_fontsize()])
    print('main dimensions mm',[85,199],'smallest ordinary text pt',min(sizes))
    plt.close(fig)

def build_large_atlas(meta,dm,z,members,arrays,main):
    fig=plt.figure(figsize=(285/25.4,285/25.4),facecolor='white')
    ax=mm_axes(fig,1,1,283,283);atlas(ax,meta,dm,z,members,arrays,main,large=True)
    fig.savefig(W/'fig3_supplement_atlas_full.pdf');fig.savefig(W/'fig3_supplement_atlas_full.svg');fig.savefig(W/'fig3_supplement_atlas_full.png',dpi=180)
    plt.close(fig)

def build_pair_supplement(frame,mats,filename):
    perpage=10
    with PdfPages(W/filename) as pdf:
        for pg in range(math.ceil(len(frame)/perpage)):
            f=plt.figure(figsize=(210/25.4,297/25.4),facecolor='white')
            mm_text(f,12,10,'CARBON motif / same-RBP CISBP reference',fontsize=11,fontweight='bold')
            mm_text(f,12,17,'Stored Tomtom alignments and p/q values; no search or significance recalculation.',fontsize=8)
            for j,(_,row) in enumerate(frame.iloc[pg*perpage:(pg+1)*perpage].iterrows()):
                col=j%2;rr=j//2;x=12+col*97;y=30+rr*50
                q,r,rs,extent=aligned_mats(row,mats)
                mm_text(f,x,y,f'{row["comparison_id"]}  {row["task"]}',fontsize=8)
                mm_text(f,x,y+5,f'{row["Query_ID"]} / {row["Target_ID"]}',fontsize=8)
                logo(mm_axes(f,x,y+12,40,9),q,extent=extent)
                logo(mm_axes(f,x+48,y+12,40,9),r,start=rs,extent=extent)
                mm_text(f,x+20,y+22,'CARBON',fontsize=7,ha='center');mm_text(f,x+68,y+22,'CISBP',fontsize=7,ha='center')
                mm_text(f,x,y+27,f'p = {fmt(row["p_value"])}   q = {fmt(row["q-value"])}',fontsize=8)
                label=row['reference_class'].replace('_',' ')
                mm_text(f,x,y+33,f'{label}; status {row["metadata_all_statuses"]}; {row["Motif_Type"]}',fontsize=7)
                mm_text(f,x,y+38,f'{row["MSource_Author"]} {row["MSource_Year"]}; PMID {row["PMID"]}',fontsize=7)
            mm_text(f,198,287,f'{pg+1} / {math.ceil(len(frame)/perpage)}',fontsize=8,ha='right')
            pdf.savefig(f);plt.close(f)

def build_c_supplement(geometry):
    f=plt.figure(figsize=(210/25.4,270/25.4),facecolor='white')
    mm_text(f,15,8,'Retained attribution and prediction-sensitivity examples',fontsize=11,fontweight='bold')
    c_tracks(f,geometry,start=21,step=75,large=True)
    f.savefig(W/'fig3_supplement_sensitivity.pdf');f.savefig(W/'fig3_supplement_sensitivity.svg');f.savefig(W/'fig3_supplement_sensitivity.png',dpi=180);plt.close(f)

def main():
    meta,dm,z,members,arrays=load_atlas()
    frame=pd.read_csv(W/'motif_comparisons_all_193.csv')
    main=frame[frame.main_selected].sort_values('main_order')
    mats=np.load(W/'paired_logo_matrices.npz')
    geo=json.loads((W/'retained_c_geometry.json').read_text(encoding='utf-8'))
    build_main(meta,dm,z,members,arrays,main,mats,geo)
    if '--main-only' in sys.argv:
        build_c_supplement(geo)
        return
    build_large_atlas(meta,dm,z,members,arrays,main)
    build_pair_supplement(frame[frame.reference_class=='direct_experiment'].sort_values(['task','Query_ID','Target_ID']),mats,'fig3_supplement_direct_pairs.pdf')
    build_pair_supplement(frame.sort_values(['task','Query_ID','Target_ID']),mats,'fig3_supplement_all_193_pairs.pdf')
    build_c_supplement(geo)

if __name__=='__main__':main()
