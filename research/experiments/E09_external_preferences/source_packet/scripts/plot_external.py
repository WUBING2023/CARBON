from literal_tables import read_csv_literal
"""Pure plotting of retained, frozen external-rule results."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'results';OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.labelsize':6.5,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.5,'xtick.major.size':2,'ytick.major.size':0,'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})
colors={'full':'#E64B35','sequence':'#4DBBD5','frozen_flat':'#8491B4','frozen_perm':'#00A087','frequency':'#3C5488'}
names={'full':'Full','sequence':'Sequence','frozen_flat':'Flat*','frozen_perm':'Perm*','frequency':'Frequency'}
markers={'full':'o','sequence':'s','frozen_flat':'^','frozen_perm':'v','frequency':'D'}
table=read_csv_literal(DATA/'WIDE_PROTEIN_METRICS.csv');table=table[~table.exposed_protein];proteins=sorted(table.protein.unique());assert len(proteins)==25
metrics=['spearman_all5mer','within_composition_concordance'];arms=list(colors)
fig=plt.figure(figsize=(85/25.4,132/25.4))
axA=fig.add_axes([.19,.52,.352,.41]);axB=fig.add_axes([.59,.52,.36,.41],sharey=axA)
points=[]
for ax,metric,letter in [(axA,metrics[0],'a'),(axB,metrics[1],'b')]:
    tab=table.pivot(index='protein',columns='arm',values=metric).loc[proteins]
    for i,p in enumerate(proteins):
        y=len(proteins)-1-i
        ax.axhline(y,color='#ECEEEF',lw=.35,zorder=0)
        ax.plot([tab.loc[p,'sequence'],tab.loc[p,'full']],[y,y],color='#BDC3C5',lw=.7,zorder=1)
        for arm in arms:
            val=float(tab.loc[p,arm]);off={'full':.12,'sequence':-.12,'frozen_flat':.03,'frozen_perm':-.03,'frequency':0}[arm]
            ax.scatter(val,y+off,s=8 if arm in ['full','sequence'] else 5,marker=markers[arm],facecolors=colors[arm] if arm in ['full','sequence','frequency'] else 'white',edgecolors=colors[arm],linewidths=.45,zorder=3 if arm=='full' else 2)
            points.append(dict(panel=letter,protein=p,condition=arm,endpoint=metric,value=val,y=y))
    ax.set_ylim(-.65,24.65);ax.spines[['top','right','left']].set_visible(False);ax.set_yticks(range(25),proteins[::-1]);ax.text(-.12,1.035,letter,transform=ax.transAxes,fontweight='bold',fontsize=8)
    if letter=='a':ax.set_xlim(-.85,.9);ax.set_xticks([-.5,0,.5]);ax.axvline(0,color='#AAB0B3',lw=.55,ls=(0,(2,2)));ax.set_xlabel('Spearman ρ')
    else:ax.set_xlim(.36,.83);ax.set_xticks([.4,.5,.6,.7,.8]);ax.axvline(.5,color='#AAB0B3',lw=.55,ls=(0,(2,2)));ax.set_xlabel('Same-composition order');ax.tick_params(labelleft=False)
legend=[Line2D([],[],marker=markers[a],color='none',markerfacecolor=colors[a] if a in ['full','sequence','frequency'] else 'white',markeredgecolor=colors[a],markeredgewidth=.6,markersize=3,label=names[a]) for a in arms]
fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.988),ncol=5,frameon=False,columnspacing=.65,handletextpad=.25,handlelength=.65,labelspacing=.35,fontsize=6.25)

# One complete experimentally selected twenty-word set, not selected model winners.
word=read_csv_literal(DATA/'HNRNPK_ALL_TOP20_WORDS.csv');tru=word.groupby('kmer').external_Z5.first().sort_values(ascending=False,kind='stable');kmers=list(tru.index);assert len(kmers)==20
rank=word.groupby(['kmer','arm']).in_model_top20.mean().unstack().loc[kmers,['aligned','sequence_only','mean','permuted']]
axC=fig.add_axes([.19,.105,.34,.30]);axC.text(-.12,1.06,'c',transform=axC.transAxes,fontweight='bold',fontsize=8)
axC.set_xlim(-.5,4.5);axC.set_ylim(19.5,-.5);axC.set_yticks(range(20),[k.replace('T','U') for k in kmers]);axC.set_xticks(range(5),['Z5','Full','Seq','Mean','Perm']);axC.xaxis.tick_top();axC.tick_params(axis='x',length=0,labelsize=6)
cmgray=LinearSegmentedColormap.from_list('truth',['#F4F5F6','#475A6A']);arr=tru.to_numpy()[:,None];axC.imshow(arr,extent=(-.5,.5,19.5,-.5),aspect='auto',cmap=cmgray,vmin=0,vmax=float(tru.max()))
for j,arm in enumerate(['aligned','sequence_only','mean','permuted']):
    cm=LinearSegmentedColormap.from_list('selection',['#F5F6F7','#3C5488']);axC.imshow(rank[[arm]].to_numpy(),extent=(j+.5,j+1.5,19.5,-.5),aspect='auto',cmap=cm,vmin=0,vmax=1)
for yy in np.arange(-.5,20,.999999):axC.axhline(yy,color='white',lw=.35)
for xx in np.arange(-.5,5,1):axC.axvline(xx,color='white',lw=.6)
axC.spines[:].set_visible(False)
fig.text(.19,.445,'HNRNPK · RNAcompete',fontsize=6.5)
# Compact keys explicitly distinguish independent measured Z5 from model selection rates.
keyT=fig.add_axes([.19,.065,.07,.012]);keyT.imshow(np.linspace(0,1,100)[None],aspect='auto',cmap=cmgray);keyT.axis('off');fig.text(.19,.045,'Z5 0–3.27',fontsize=6)
keyR=fig.add_axes([.32,.065,.10,.012]);keyR.imshow(np.linspace(0,1,100)[None],aspect='auto',cmap=LinearSegmentedColormap.from_list('selection',['#F5F6F7','#3C5488']));keyR.axis('off');fig.text(.32,.045,'Selected 0–100%',fontsize=6)

axD=fig.add_axes([.715,.18,.235,.217]);ci=read_csv_literal(DATA/'HNRNPK_PAIRED_GENE_BOOTSTRAP.csv');controls=['sequence_only','mean','permuted'];dpoints=[]
for i,ctrl in enumerate(controls):
    yy=2-i
    for assay,color,marker,off in [('RNAcompete','#E64B35','o',.12),('RBNS','#3C5488','s',-.12)]:
        r=ci[(ci.assay==assay)&(ci.control==ctrl)].iloc[0];v=r.delta_recall20*100;lo=r.conditional_gene_locus_bootstrap_low*100;hi=r.conditional_gene_locus_bootstrap_high*100
        # Point may lie outside a percentile CI for discrete top-k; draw CI independently.
        axD.plot([lo,hi],[yy+off]*2,color=color,lw=.9);axD.scatter(v,yy+off,s=10,marker=marker,color=color,zorder=3)
        dpoints.append(dict(assay=assay,control=ctrl,delta_pp=float(v),low_pp=float(lo),high_pp=float(hi)))
axD.axvline(0,color='#AAB0B3',lw=.55,ls=(0,(2,2)));axD.set_xlim(-7,11);axD.set_xticks([-5,0,5,10]);axD.set_ylim(-.5,2.5);axD.set_yticks([2,1,0],['Full−Seq','Full−Mean','Full−Perm'],fontsize=6.3);axD.spines[['top','right','left']].set_visible(False);axD.set_xlabel('Δ recall (pp)');axD.text(-.17,1.1,'d',transform=axD.transAxes,fontweight='bold',fontsize=8)
axD.legend(handles=[Line2D([],[],marker='o',color='#E64B35',markersize=3,linewidth=.8,label='RNAcompete'),Line2D([],[],marker='s',color='#3C5488',markersize=3,linewidth=.8,label='RBNS')],loc='upper left',bbox_to_anchor=(-.16,-.37),frameon=False,fontsize=6.3,handlelength=.9,handletextpad=.4)

for ext in ['png','pdf','svg']:fig.savefig(OUT/f'EXTERNAL_RULES.{ext}',dpi=600)
figuredata=dict(width_mm=85,height_mm=132,main_new_proteins=proteins,panels_ab=points,panel_c_words=kmers,panel_c_Z5=tru.to_dict(),panel_c_selection=rank.to_dict(),panel_d=dpoints,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [DATA/'WIDE_PROTEIN_METRICS.csv',DATA/'HNRNPK_ALL_TOP20_WORDS.csv',DATA/'HNRNPK_PAIRED_GENE_BOOTSTRAP.csv']})
(OUT/'FIGURE_DATA.json').write_text(json.dumps(figuredata,indent=2),encoding='utf-8')
fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
for text in fig.findobj(matplotlib.text.Text):
    if not text.get_visible() or not text.get_text():continue
    b=text.get_window_extent(renderer)
    if b.x0<-.5 or b.y0<-.5 or b.x1>fig.bbox.width+.5 or b.y1>fig.bbox.height+.5:outside.append(dict(text=text.get_text(),bounds=list(b.bounds)))
(OUT/'LAYOUT_QA.json').write_text(json.dumps(dict(outside_canvas=outside,width_mm=85,height_mm=132,min_font_pt=6.0),indent=2))
print('saved',OUT,'outside text',outside)
