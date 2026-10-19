"""Complete conservation strata, without selecting favorable arms or bins."""
import os,sys
os.environ['OPENBLAS_NUM_THREADS']='4';sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;F=R/'figures'
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.labelsize':6.5,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.55,'pdf.fonttype':42,'svg.fonttype':'none'})
q=pd.read_csv(R/'DIRECTION_BY_CONSERVATION.csv');a=pd.read_csv(R/'CONSERVATION_ASSOCIATIONS.csv')
colors={'A_native':'#E64B35','S':'#4DBBD5','A_flat':'#00A087','A_permmean':'#8491B4','ref':'#3C5488'}
names={'A_native':'Full','S':'Sequence','A_flat':'Flat Φ','A_permmean':'Perm. Φ','ref':'REF rule'}
W,H=85,137;fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
def t(x,y,s,**kw):return fig.text(x/W,1-y/H,s,va='center',fontsize=kw.pop('fontsize',6.5),**kw)
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
for j,(cohort,note,y) in enumerate([('strict_test','396 events · chr8/15',8),('all','7,262 events · all splits',55)]):
    t(2.5,y-4,chr(97+j),fontweight='bold',fontsize=9);t(8,y-4,note)
    p=ax(15,y,67,30);p.set_ylim(20,100);p.set_xlim(.6,4.4);p.set_xticks([1,2,3,4],['Q1','Q2','Q3','Q4']);p.set_yticks([25,50,75,100]);p.spines[['top','right']].set_visible(False)
    p.set_ylabel('Direction accuracy (%)');p.set_xlabel('Center PhyloP quartile');p.tick_params(length=2,width=.5,pad=1.5)
    for arm,offset in [('A_native',-.18),('S',0),('ref',.18)]:
        d=q[(q.cohort==cohort)&(q.arm==arm)].sort_values('quartile');x=d.quartile.to_numpy()+offset
        p.errorbar(x,100*d.accuracy,yerr=np.vstack((100*(d.accuracy-d.ci_low),100*(d.ci_high-d.accuracy))),fmt='o',ms=3,color=colors[arm],elinewidth=.6,capsize=1.2,label=names[arm]);p.plot(x,100*d.accuracy,c=colors[arm],lw=.65)
    if j==0:p.legend(loc='upper left',frameon=False,fontsize=6,ncol=3,columnspacing=.9,handletextpad=.3,borderpad=.1)
    # Different RBP coverage by stratum is kept explicit in the source table/caption.
t(2.5,99,'c',fontweight='bold',fontsize=9);t(8,99,'396 events · predicted response magnitude')
p=ax(24,105,58,23);p.set_xlim(-.3,.3);p.set_ylim(-.5,3.5);p.set_yticks([]);p.set_xticks([-.3,0,.3]);p.spines[['top','right','left']].set_visible(False);p.axvline(0,c='#9CA6AD',lw=.6,ls=(0,(2,2)));p.tick_params(length=2,width=.5,pad=1.5)
for i,arm in enumerate(['A_native','S','A_flat','A_permmean']):
    r=a[(a.cohort=='strict_test')&(a.arm==arm)&(a.outcome=='abs_delta')].iloc[0];y=3-i
    t(8,105+(3.5-y)/4*23,names[arm]);p.plot([r.ci_low,r.ci_high],[y,y],c=colors[arm],lw=.8);p.scatter([r.adjusted_standardized_rank_slope],[y],s=12,c=colors[arm])
t(53,134,'Adjusted standardized-rank Φ slope',ha='center')
for ext in ['pdf','png','svg']:fig.savefig(F/('SUPPLEMENTARY_FIGURE_S23_CONSERVATION.'+ext),dpi=600,facecolor='white')
plt.close(fig)
# S22 is the separately usable literal-SNV/whole-cohort addition. Relabel e→b in native SVG/PDF.
import fitz,xml.etree.ElementTree as ET
src=fitz.open(F/'FIG5_ADDITIONAL_a_e_85mm.pdf');page=src[0]
span=next(s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text']=='e')
page.add_redact_annot(fitz.Rect(span['bbox']),fill=(1,1,1));page.apply_redactions(images=0,graphics=0)
from matplotlib import font_manager
fontfile=font_manager.findfont(font_manager.FontProperties(family='Arial',weight='bold'))
page.insert_font(fontname='ReplayBold',fontfile=fontfile)
page.insert_text(fitz.Point(*span['origin']),'b',fontname='ReplayBold',fontsize=span['size'])
src.save(F/'SUPPLEMENTARY_FIGURE_S22_NATURAL_SNV.pdf');page.get_pixmap(matrix=fitz.Matrix(600/72,600/72),alpha=False).save(F/'SUPPLEMENTARY_FIGURE_S22_NATURAL_SNV.png');src.close()
svg=(F/'FIG5_ADDITIONAL_a_e_85mm.svg').read_text('utf8');svg=svg.replace('>e</text>','>b</text>');(F/'SUPPLEMENTARY_FIGURE_S22_NATURAL_SNV.svg').write_text(svg,'utf8')
# Same bootstrap draws give paired adjusted association contrasts, without a new fitted endpoint.
z=np.load(R/'WEIGHTED_REFIT_SAMPLES_200.npz');samples=z['samples'];point=z['point'];out=[]
for ai,arm in [(2,'S'),(4,'A_flat'),(6,'A_permmean')]:
    for offset,outcome in [(0,'abs_delta'),(1,'credit')]:
        v=samples[:,offset]-samples[:,ai+offset];lo,hi=np.quantile(v,[.025,.975]);out.append(dict(cohort='strict_test',contrast='A_native minus '+arm,outcome=outcome,rank_slope_difference=point[offset]-point[ai+offset],ci_low=lo,ci_high=hi,replicates=200,interval='paired full weighted nuisance refits; secondary sensitivity'))
pd.DataFrame(out).to_csv(R/'ADJUSTED_ASSOCIATION_PAIRED_CONTRASTS.csv',index=False)
print(pd.DataFrame(out).to_string(index=False))
