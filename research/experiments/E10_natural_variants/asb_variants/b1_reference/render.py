"""Reference-led literal sequence intervention/readout/frozen-cohort figure.

All plotted numerical values come from INTERVENTION_PACKET.npz. Pairwise
contrasts are descriptive computational interventions, not new inferential tests.
Canvas dimensions are fixed; no tight bounding-box resizing is used.
"""
from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

OUT = Path(__file__).resolve().parent
D = np.load(OUT / 'INTERVENTION_PACKET.npz')
W, H = 85., 132.
FULL, SEQ, FREQ = '#E64B35', '#4DBBD5', '#3C5488'
BASE = {'A': '#3C5488', 'C': '#00A087', 'G': '#7E6148', 'U': '#E64B35', 'T': '#E64B35', 'N': '#9A9A9A'}
plt.rcParams.update({'font.family':'Arial', 'font.size':6.5, 'axes.labelsize':6.5,
                    'xtick.labelsize':6.5, 'ytick.labelsize':6.5,
                    'axes.linewidth':.55, 'pdf.fonttype':42, 'ps.fonttype':42,
                    'svg.fonttype':'none', 'savefig.facecolor':'white'})
fig = plt.figure(figsize=(W/25.4, H/25.4), facecolor='white')

def text(x,y,s,fs=6.5,ha='left',va='center',**kw):
    return fig.text(x/W,1-y/H,s,fontsize=fs,ha=ha,va=va,**kw)

def line(x0,y0,x1,y1,**kw):
    l=Line2D([x0/W,x1/W],[1-y0/H,1-y1/H],transform=fig.transFigure,**kw)
    fig.add_artist(l);return l

def rect(x,y,w,h,**kw):
    p=Rectangle((x/W,1-(y+h)/H),w/W,h/H,transform=fig.transFigure,**kw)
    fig.add_artist(p);return p

def axmm(x,y,w,h):
    return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])

def letters(x,y,word,step=2.0,fs=7.0,bold=False):
    for j,b in enumerate(word):
        text(x+j*step,y,b,fs=fs,ha='center',color=BASE[b],fontweight='bold' if bold else 'normal')

def edit(x,y,lo,hi,step=2.0):
    letters(x,y,lo,step=step)
    text(x+5*step+1.0,y,'→',fs=8,color='#42494C',ha='center')
    letters(x+5*step+4.0,y,hi,step=step,bold=True)

def native(protein,x,y):
    seq=str(D[protein+'_native_sequences'][0]).replace('T','U')
    # The literal native central 5-mer is shown between its actual 3-nt flanks.
    letters(x,y,seq[45:56],step=2.1,fs=7.0)
    line(x+3*2.1-.7,y+2.1,x+7*2.1+.7,y+2.1,color='#323D43',lw=.7)
    text(x+5*2.1,y+4.2,'49–53',ha='center')

def clean(ax):
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='x',length=2,width=.5,pad=1.5)
    ax.set_yticks([])

def paired_strip(ax,protein,col,center,spread=.22,arm_gap=.18):
    a=np.asarray(D[protein+'_full_paired_edit_delta'])
    b=np.asarray(D[protein+'_sequence_paired_edit_delta'])
    if a.ndim==2:a,b=a[:,col],b[:,col]
    # Fixed context-ID offsets, shared across both models. No values sort points.
    offsets=(np.arange(32)%16-7.5)/15*spread
    offsets += np.where(np.arange(32)<16,-.055,.055)
    yf=center+offsets+arm_gap/2;ys=center+offsets-arm_gap/2
    for j in range(32):
        ax.plot([a[j],b[j]],[yf[j],ys[j]],color='#BAC0C4',lw=.32,alpha=.7,zorder=1)
    for marker,sl in [('o',slice(0,16)),('^',slice(16,32))]:
        ax.scatter(a[sl],yf[sl],s=5.4,c=FULL,marker=marker,lw=.15,edgecolors='white',zorder=3)
        ax.scatter(b[sl],ys[sl],s=5.4,c=SEQ,marker=marker,lw=.15,edgecolors='white',zorder=3)
    # Large diamonds are descriptive means, not CIs or biological replicates.
    ax.scatter([a.mean()],[center+arm_gap/2],s=17,marker='|',c='#A82F1E',lw=1.1,zorder=4)
    ax.scatter([b.mean()],[center-arm_gap/2],s=17,marker='|',c='#238EAB',lw=1.1,zorder=4)

# a: dominant concrete operation, independent assay, paired native responses.
text(2.5,4,'a',fs=9,fontweight='bold');text(8,4,'QKI',fs=8,fontweight='bold')
native('QKI',8.5,10.7)
text(9,7.0,'Native · 101 nt',fs=6.5)
phi=np.asarray(D['QKI_native_conservation'][0])
pa=axmm(40,8.6,42,9.1)
pa.plot(np.arange(1,102),phi,color='#8491B4',lw=.65)
pa.axvspan(48.5,53.5,color='#DAE1E6',alpha=.55,lw=0)
pa.set_xlim(1,101);pa.set_ylim(min(0,float(phi.min())),float(phi.max()))
pa.set_xticks([]);pa.set_yticks([min(0,float(phi.min())),float(phi.max())])
pa.set_yticklabels([f'{min(0,float(phi.min())):.1f}',f'{float(phi.max()):.1f}'])
pa.tick_params(axis='y',length=1.5,pad=1,width=.4)
pa.spines[['top','right','bottom']].set_visible(False)
text(55,6.9,'Φ fixed',color='#8491B4',ha='center')
line(19,17.8,19,19.2,color='#B5BEC4',lw=.5)
line(19,19.2,3.0,19.2,color='#B5BEC4',lw=.5)
line(3.0,19.2,3.0,47.5,color='#B5BEC4',lw=.5)
text(8,22,'Swap',fs=6.5)
text(36.5,22,'Assay ΔZ5',ha='center')
text(63.5,22,'Model response',ha='center')
qa=axmm(30.5,25,12.5,27)
qa.set_xlim(0,2.4);qa.set_ylim(-.5,2.5);clean(qa)
qa.set_xticks([0,1,2]);qa.axvline(0,color='#C9CDCF',lw=.5)
qr=axmm(47.5,25,34.5,27)
qr.set_xlim(-.1,1.17);qr.set_ylim(-.5,2.5);clean(qr)
qr.set_xticks([0,.5,1]);qr.axvline(0,color='#939DA4',lw=.55,ls=(0,(2,2)))
for j in range(3):
    yc=2-j;ym=25+(2.5-yc)/3*27
    line(3.0,ym,4.3,ym,color='#B5BEC4',lw=.5)
    edit(5.3,ym,str(D['QKI_words_RNA'][j]),'CUAAC',step=2.0)
    dz=float(D['QKI_assay_Z5'][-1]-D['QKI_assay_Z5'][j])
    qa.plot([0,dz],[yc,yc],color='#A7B0B5',lw=.8)
    qa.scatter([dz],[yc],s=12,color='#3B464C',zorder=3)
    paired_strip(qr,'QKI',j,yc,spread=.27,arm_gap=.25)
text(64.5,56,'Δlog count',ha='center')
text(10,57.5,'● HepG2   ▲ K562',fs=6.5,color='#4D565C')
text(50,58,'Full',color=FULL);text(64,58,'Sequence',color=SEQ)

# b: a frequency-opposed word-order change. Frequency is a separate scalar.
text(2.5,64,'b',fs=9,fontweight='bold');text(8,64,'HNRNPA1',fs=8,fontweight='bold')
native('HNRNPA1',8.5,70.4)
text(37,69.2,'Assay ΔZ5 +5.05',fs=6.5)
text(37,73,'Native frequency −101',fs=6.5,color=FREQ)
edit(5.3,82,'UGAGG','UAGGG',step=2.0)
hr=axmm(34.5,78,47.5,11)
hr.set_xlim(-.06,.225);hr.set_ylim(-.5,.5);clean(hr)
hr.set_xticks([0,.1,.2]);hr.axvline(0,color='#939DA4',lw=.55,ls=(0,(2,2)))
paired_strip(hr,'HNRNPA1',0,0,spread=.34,arm_gap=.27)
text(60.5,93.5,'Δlog count',ha='center')

# c: every frozen new25 protein and all three endpoints on matched axes.
text(2.5,98,'c',fs=9,fontweight='bold')
text(8,98,'25 RBPs',fs=7.3,fontweight='bold')
text(76,98,'Freq.',color=FREQ,ha='center')
line(70.5,98,70.5,98,color=FREQ,lw=0,marker='o',markersize=2.2)
protein=D['cohort_proteins'].tolist();values=D['cohort_order_accuracy_percent']
assert values.shape==(25,3)
for start,n,x in [(0,13,18.5),(13,12,59.8)]:
    ax=axmm(x,101.5,22.4,25.2)
    ax.set_xlim(39,73.5);ax.set_ylim(n-.5,-.5);clean(ax)
    ax.set_xticks([40,50,60,70]);ax.axvline(50,color='#939DA4',lw=.5,ls=(0,(2,2)),zorder=0)
    for local,i in enumerate(range(start,start+n)):
        row=values[i]
        text(x-1.5,101.5+(local+.5)/n*25.2,protein[i],ha='right',fs=6.5)
        ax.plot([row[2],row[0]],[local,local],color='#BFC6CA',lw=.55,zorder=1)
        ax.scatter([row[2]],[local],s=8,c=FREQ,zorder=3)
        ax.scatter([row[0]],[local],s=9.5,c=FULL,zorder=3)
        ax.scatter([row[1]],[local],s=20,marker='|',c=SEQ,lw=1.1,zorder=4)
text(50,130.8,'Within-composition order accuracy (%)',ha='center')

# Pixel outputs plus native vector exports have the exact physical canvas.
for ext in ['png','pdf','svg']:
    fig.savefig(OUT/('figure.'+ext),dpi=600,facecolor='white')
    shutil.copyfile(OUT/('figure.'+ext),OUT/('FIGURE_B.'+ext))
    v=OUT/'candidate_v1';v.mkdir(exist_ok=True)
    if not (v/('figure.'+ext)).exists():shutil.copyfile(OUT/('figure.'+ext),v/('figure.'+ext))
plt.close(fig)
print(str(OUT/'figure.png'))
