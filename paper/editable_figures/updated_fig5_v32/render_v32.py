"""Complete Figure5 v25: explicit population, rule, natural case and clinical audit."""
from data import *
from matplotlib.colors import Normalize
from matplotlib.ticker import NullLocator

W,H=178.,178.
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white');records=[]
def ys(y):return y-6 if y>=74 else y

def ax(x,y,w,h):return fig.add_axes([x/W,1-(ys(y)+h)/H,w/W,h/H])
def txt(x,y,s,fs=6.,col=None,**kw):return fig.text(x/W,1-ys(y)/H,s,fontsize=max(fs,5.5),color=INK if col is None else col,va='center',**kw)
def line(x,y,x2,y2,**kw):fig.add_artist(Line2D([x/W,x2/W],[1-ys(y)/H,1-ys(y2)/H],transform=fig.transFigure,**kw))
def tag(x,y,s,label):txt(x,y,s,9,'#000000',weight='bold');txt(x+5,y,label,6.5,'#000000')
def letters(x,y,s,step=1.6,fs=6.6):
    for i,b in enumerate(s):txt(x+i*step,y,b,fs,BASE[b],ha='center',weight='bold',fontfamily='DejaVu Sans Mono')
def clean(a,left=False):
    for sp in ['top','right']+([] if left else ['left']):a.spines[sp].set_visible(False)
    a.tick_params(length=1.6,pad=1.3,width=.55,colors='#000000')
    for sp in a.spines.values():sp.set_linewidth(.55);sp.set_color('#000000')
def rec(p,k,id,arm,q,v):records.append(dict(panel=p,kind=k,id=str(id),arm=arm,quantity=q,value=float(v)))


# Figure 3 supplies the typography/axis grammar; semantic series colours
# remain the accepted NPG bindings. Annotation grey is a neutral text tier.
INK,GRAY='#000000','#62686E'
plt.rcParams.update({'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':INK,
    'xtick.color':INK,'ytick.color':INK})

tag(1,3.5,'a','RNAcompete · 25 proteins')
txt(9,10.6,'All 5-mers',6.3,INK)
txt(62,10.6,'Same composition',6.3,INK)
# This jitter is identity-stable geometry, not biological data or a new sample.
jitter=np.array([(int(hashlib.sha256(p.encode()).hexdigest()[:8],16)/0xffffffff-.5)*.17 for p in proteins])
for metric,x,mult,limits,ticks,dticks,rawlabel in [
    ('spearman_all5mer',9.,1,(-.7,.86),[-.5,0,.5],[-.5,0,.5],'Spearman ρ'),
    ('within_composition_concordance',62.,100,(35.,80.),[40,60,80],[-10,0,10,20],'Order agreement (%)')]:
    freq=np.array([float(lookup[(p,'frequency')][metric])*mult for p in proteins])
    full=np.array([float(lookup[(p,'full')][metric])*mult for p in proteins])
    delta=np.array([(float(lookup[(p,'full')][metric])-float(lookup[(p,'frequency')][metric]))*mult for p in proteins])
    mf,mm=float(np.mean(freq)),float(np.mean(full))
    c=next(r for r in contrasts if r['endpoint']==metric and r['control']=='frequency')
    mean,lo,hi=[float(c[k])*mult for k in ['mean_delta','protein_paired_bootstrap_low','protein_paired_bootstrap_high']]
    assert np.isclose(mm-mf,mean,atol=1e-12)
    assert min(freq.min(),full.min())>=limits[0] and max(freq.max(),full.max())<=limits[1]
    aa=ax(x,18.,22.,37.)
    for i,(p,ff,vv,dd) in enumerate(zip(proteins,freq,full,delta)):
        aa.plot([0+jitter[i],1+jitter[i]],[ff,vv],color='#BAC5CE',alpha=.57,lw=.48,zorder=1)
        aa.plot(jitter[i],ff,marker='s',ms=2.6,mfc=FREQ,mec='white',mew=.2,linestyle='none',zorder=3)
        aa.plot(1+jitter[i],vv,marker='o',ms=2.8,mfc=FULL,mec='white',mew=.2,linestyle='none',zorder=3)
        rec('a','protein',p,'full',metric+'_minus_frequency',dd)
        for arm,value in [('full',vv),('frequency',ff)]:
            rec('a','protein',p,arm,metric,lookup[(p,arm)][metric])
            if mult!=1:rec('a','protein',p,arm,metric+'_percent',value)
    aa.plot([-.21,.21],[mf,mf],color=FREQ,lw=1.5,zorder=5)
    aa.plot([.79,1.21],[mm,mm],color=FULL,lw=1.5,zorder=5)
    aa.set(xlim=(-.4,1.4),ylim=limits,yticks=ticks,xticks=[0,1],xticklabels=['Frequency','Full'])
    clean(aa,True);aa.tick_params(labelsize=6.0)
    aa.tick_params(axis='x',length=0,pad=3)
    for t in aa.get_xticklabels():t.set_color(INK)
    aa.set_ylabel(rawlabel,fontsize=6.0,labelpad=1.1)
    da=ax(x+27.,18.,12.,37.)
    da.set(ylim=(limits[0]-mf,limits[1]-mf),xlim=(-.15,1.15),yticks=dticks,xticks=[.5],xticklabels=['Δρ' if mult==1 else 'Δorder (pp)'])
    assert delta.min()>=limits[0]-mf and delta.max()<=limits[1]-mf
    da.axhline(0,color='#A9B6C0',lw=.55,zorder=0)
    # Gray marks are the actual 25 paired protein differences. No density or
    # unobserved bootstrap samples are fabricated from the saved CI endpoints.
    da.scatter(.26+jitter*1.6,delta,s=5.2,color='#A8B8C5',alpha=.85,edgecolors='none',zorder=2)
    da.errorbar(.79,mean,yerr=[[mean-lo],[hi-mean]],fmt='o',ms=4.0,color=FULL,capsize=2.1,lw=1.25,zorder=4)
    for sp in ['left','top']:da.spines[sp].set_visible(False)
    da.yaxis.tick_right();da.yaxis.set_label_position('right')
    da.tick_params(axis='y',labelsize=5.7,length=1.8,pad=1.4);da.tick_params(axis='x',labelsize=5.7,length=0,pad=3)
    txt(x+33.,14.0,'Full − Freq',4.9,GRAY,ha='center')
    # Cross-axis guide: zero is aligned with the Frequency mean, while the
    # mean difference is aligned with the Full mean, as in an estimation plot.
    y0=18+37*(limits[1]-mf)/(limits[1]-limits[0])
    y1=18+37*(limits[1]-mm)/(limits[1]-limits[0])
    freq_tick_end=x+(.21+.4)/1.8*22.
    full_tick_end=x+(1.21+.4)/1.8*22.
    delta_mean_x=x+27.+(.79+.15)/1.3*12.
    line(freq_tick_end,y0,x+39.,y0,color='#C8D1D8',lw=.45)
    line(full_tick_end,y1,delta_mean_x,y1,color='#E8ADA3',lw=.45)
    for q,value in [('mean',mean),('CI_low',lo),('CI_high',hi)]:rec('a','macro',metric,'Full_minus_Frequency',q,value)
    for arm,value in [('full',mm),('frequency',mf)]:rec('a','macro',metric,arm,'absolute_mean' if mult==1 else 'absolute_mean_percent',value)

# B: true lifted adjacency remains; duplicate floor edges are not data layers.
tag(107,3.5,'b','QKI')
nodes=rows(HERE/'frozen_inputs/qki_complete_graph/ALL_30_NODE_PLOT_DATA.csv')
edges=rows(HERE/'frozen_inputs/qki_complete_graph/ALL_48_EDGE_PLOT_DATA.csv')
nw={r['word']:r for r in nodes};highlights={'CUAAC','UCAAC','CAUAC','CUACA'}
bx=np.array([float(r['base_x_mm']) for r in nodes]);by=np.array([float(r['base_y_mm']) for r in nodes])
base={w:(110+64*(float(r['base_x_mm'])-bx.min())/np.ptp(bx),47+7*(float(r['base_y_mm'])-by.min())/np.ptp(by)) for w,r in nw.items()}
lift={w:(p[0],p[1]-12*float(nw[w]['external_Z5'])) for w,p in base.items()}
def hull(points):
    pts=sorted(set(points))
    def cross(o,a,b):return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    low=[]
    for p in pts:
        while len(low)>=2 and cross(low[-2],low[-1],p)<=0:low.pop()
        low.append(p)
    high=[]
    for p in reversed(pts):
        while len(high)>=2 and cross(high[-2],high[-1],p)<=0:high.pop()
        high.append(p)
    return low[:-1]+high[:-1]
floor=hull(list(base.values()))
for p,q in zip(floor,floor[1:]+floor[:1]):line(*p,*q,color='#E2E9EF',lw=.4)
gold_edges=0
for r in edges:
    lo,hi=r['lower_assay_word'],r['higher_assay_word'];gold=hi=='CUAAC' and lo in highlights
    line(*lift[lo],*lift[hi],color='#B58A31' if gold else '#A6B5C4',lw=1.05 if gold else .47,alpha=1 if gold else .62)
    gold_edges+=gold
assert len(nodes)==30 and len(edges)==48 and gold_edges==3
for w in highlights:line(*base[w],*lift[w],color='#D9E1E8',lw=.4)
signed=LinearSegmentedColormap.from_list('response',[FREQ,'white','#EFC000'])
effectnorm=AsinhNorm(linear_width=.02,vmin=-1.8,vmax=1.8)
for w,r in nw.items():
    x,y=lift[w];f=float(r['Full_mean_delta_log_count']);s=float(r['Sequence_mean_delta_log_count'])
    fig.add_artist(Line2D([x/W],[1-y/H],transform=fig.transFigure,marker='o',ms=4.7,fillstyle='left',markerfacecolor=signed(effectnorm(f)),markerfacecoloralt=signed(effectnorm(s)),markeredgecolor='#B58A31' if w in highlights else '#8F9FAD',markeredgewidth=.6,linestyle='none'))
    rec('b','word',w,'assay','Z5',r['external_Z5'])
    for arm,value in [('full',f),('sequence',s)]:rec('b','word',w,arm,'mean_delta_log_count_vs_native',value)
for w,(dx,dy) in {'CUAAC':(-4.5,-3.7),'UCAAC':(-6.1,-3.2),'CAUAC':(-9.6,2.4),'CUACA':(2,4)}.items():
    letters(lift[w][0]+dx,lift[w][1]+dy,w,1.6,7.2 if w=='CUAAC' else 6.6)
txt(108,11,r'RNAcompete $Z_5$',6,INK);line(108.4,49,108.4,22.5,color=INK,lw=.55)
for z in [0,1,2]:line(107.9,49-12*z,109.1,49-12*z,color=INK,lw=.55);txt(107.2,49-12*z,str(z),5.7,INK,ha='right')
txt(148,9.8,'Full',6,INK,ha='center');txt(155.5,9.8,'|',6,INK,ha='center');txt(165,9.8,'Sequence',6,INK,ha='center')
txt(156.5,14,'Node halves: left | right',5.5,GRAY,ha='center')

txt(110,60.5,'Δlog count vs native',5.7,INK)
cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=effectnorm,cmap=signed),cax=ax(145,60,30,1.4),orientation='horizontal')
cb.set_ticks([-1,0,1]);cb.set_ticklabels(['−1','0','1']);cb.ax.tick_params(labelsize=5.5,length=1,pad=.3);cb.ax.minorticks_off();cb.outline.set_visible(False)

# C: one allele object, primary models, grouped fixed-input controls, raw assay.
RP=HERE/'frozen_inputs/reporter_validation'
case=json.loads((RP/'SELECTED_CASES.json').read_text('utf8'))['K562:UPF1:chr12_52388891_+']
cp=case['prediction'];ce=case['event']
assert cp['variant_key']=='chr12:51995107:A:G' and cp['split']=='training'
tag(1,78,'c','ACVR1B · chr12:51995107 A>G')
txt(1,82.7,'UPF1 · K562 ASB: A preferred',5.7,INK)
txt(68,82.7,'Training locus',5.5,GRAY)
cr=ce['context']['reference_sequence_rna_direction_dna'].replace('T','U')
ca=ce['context']['minor_sequence_dna'].replace('T','U')
assert cr[50]=='A' and ca[50]=='G'
for yy,allele,seq in [(88.1,'A',cr),(92.,'G',ca)]:
    txt(1,yy,allele,6.0,INK,weight='bold');letters(7,yy,seq[42:59],1.5,6.1)
fig.add_artist(Rectangle(((7+8*1.5-.78)/W,1-ys(93.25)/H),1.56/W,6.4/H,transform=fig.transFigure,fc='none',ec='#B58A31',lw=.6))
txt(41,90.0,r'Reference $\Phi_0$ −2.29',5.7,INK)

txt(1,98,'K562 binding prediction',6.0,INK)
txt(67,98,'HeLa RNA reporter',6.0,INK)
cm=ax(18,103,36,22)
mainitems=[('Full','full','A_native_major_minus_minor',FULL,0),('Sequence','sequence','S_major_minus_minor',SEQ,1),('Flat Φ','flat','A_flat_major_minus_minor',FLAT,2.7),('Perm Φ mean','perm','A_permmean_major_minus_minor',PERM,3.7)]
for lab,arm,key,col,yy in mainitems:
    v=float(cp[key]);cm.barh(yy,v,height=.64 if yy<2 else .56,color=col,alpha=1 if yy<2 else .73,zorder=2)
    rec('c','event',cp['event_id'],arm,'A_minus_G_log_count',v)
for j,marker in enumerate(['o','^','D']):
    arm=f'perm{j}';v=float(cp[f'A_perm{j}_major_minus_minor'])
    cm.plot(v,3.7+(j-1)*.16,marker=marker,ms=3.8,mfc='white',mec=PERM,mew=.8,linestyle='none',zorder=4)
    rec('c','event',cp['event_id'],arm,'A_minus_G_log_count',v)
cm.axvline(0,color=INK,lw=.6)
cm.set(xlim=(-.02,.02),ylim=(4.45,-.65),yticks=[0,1,2.7,3.7],yticklabels=[x[0] for x in mainitems],xticks=[-.02,0,.02])
clean(cm);cm.tick_params(axis='y',length=0,pad=3,labelsize=5.6);cm.tick_params(axis='x',labelsize=5.8)
for t in cm.get_yticklabels():t.set_color(INK)
cm.set_xlabel('A − G Δlog count',fontsize=6.0,labelpad=.7)
for x,marker,lab in [(3,'o','0'),(12,'^','1'),(21,'D','2')]:
    fig.add_artist(Line2D([x/W],[1-ys(133.5)/H],transform=fig.transFigure,marker=marker,ms=3.0,mfc='white',mec=PERM,mew=.7,linestyle='none'));txt(x+2.3,133.5,lab,4.9,GRAY)
txt(29,133.5,'Fixed Φ permutations',4.9,GRAY)
ra=ax(76,103,25,22)
for j,allele in enumerate(['A','G']):
    vals=[r for r in case['reporter'] if r['RNA_allele']==allele];assert len(vals)==3
    vv=[float(r['value']) for r in vals]
    ra.scatter(j+np.array([-.14,0,.14]),vv,s=17,facecolor='white',edgecolor=FREQ,lw=.9,zorder=3)
    ra.plot([j-.25,j+.25],[np.mean(vv)]*2,color=FREQ,lw=1.5)
    for r in vals:rec('c','reporter',cp['variant_key']+'|'+allele+'|'+str(r['biological_replicate']),'published_HeLa',r['endpoint'],r['value'])
ra.set(xlim=(-.5,1.5),ylim=(0,1.65),xticks=[0,1],xticklabels=['A','G'],yticks=[0,.5,1,1.5])
clean(ra,True);ra.tick_params(labelsize=5.5);ra.set_ylabel('Normalized RNA',fontsize=6.0,labelpad=1)
txt(67,133.5,'qRT-PCR · n = 3 per allele',5.5,GRAY)

# D: event-level geometry dominates, all five original summaries remain.
tag(107,78,'d','Native ASB · 396 events')
for x,col,label in [(108,FULL,'Full only correct'),(141,SEQ,'Sequence only correct')]:
    fig.add_artist(Line2D([x/W],[1-ys(84.2)/H],transform=fig.transFigure,marker='o',ms=2.7,mfc=col,mec='none',linestyle='none'));txt(x+2,84.2,label,5.5,INK)
sc=ax(114,90.5,34.5,34.5)
sc.set_xscale('asinh',linear_width=.02);sc.set_yscale('asinh',linear_width=.02)
fvals=np.array([float(r['A_native_major_minus_minor']) for r in cohort]);svals=np.array([float(r['S_major_minus_minor']) for r in cohort])
assert np.all(fvals!=0) and np.all(svals!=0)
for ss,ff,col in [(False,False,'#AABBC9'),(True,True,'#AABBC9'),(False,True,FULL),(True,False,SEQ)]:
    ids=np.where(((svals>0)==ss)&((fvals>0)==ff))[0]
    sc.scatter(svals[ids],fvals[ids],s=6.4 if ss!=ff else 4.5,c=col,alpha=.9 if ss!=ff else .58,edgecolors='none',zorder=3 if ss!=ff else 2)
census=np.zeros((2,2),int)
for r,f,s in zip(cohort,fvals,svals):
    census[int(s<0),int(f<0)]+=1
    for arm,value in [('full',f),('sequence',s)]:rec('d','event',r['event_id'],arm,'source_preferred_minus_other_log_count',value)
    rec('d','event',r['event_id'],'full_sequence','Full_direction_agreement',int(f>0));rec('d','event',r['event_id'],'full_sequence','Sequence_direction_agreement',int(s>0))
assert census.tolist()==[[168,61],[43,124]]
sc.axhline(0,color='#ABB8C2',lw=.55);sc.axvline(0,color='#ABB8C2',lw=.55)
sc.set(xlim=(-.2,1.8),ylim=(-.2,1.8),xticks=[-.1,0,.1,1],yticks=[-.1,0,.1,1])
sc.set_xticklabels(['−.1','0','.1','1']);sc.set_yticklabels(['−.1','0','.1','1'])
sc.xaxis.set_minor_locator(NullLocator());sc.yaxis.set_minor_locator(NullLocator());clean(sc,True)
sc.tick_params(labelsize=5.7);sc.set_xlabel('Sequence Δlog count',fontsize=5.8,labelpad=1);sc.set_ylabel('Full Δlog count',fontsize=5.8,labelpad=1)
for x,y,lab,col,ha in [(.025,.97,'43',FULL,'left'),(.97,.97,'168',INK,'right'),(.025,.025,'124',GRAY,'left'),(.97,.025,'61',SEQ,'right')]:
    sc.text(x,y,lab,transform=sc.transAxes,ha=ha,va='top' if y>.5 else 'bottom',fontsize=6.2,color=col,weight='bold',bbox={'fc':'white','ec':'none','pad':.25,'alpha':.9})
for i in range(2):
    for j in range(2):rec('d','census',str(i)+str(j),'full_sequence','event_count',census[i,j])
txt(153,88.9,'RBP-equal',5.0,GRAY);txt(153,92.1,'accuracy (%)',5.0,GRAY)
acc=ax(162,96,14.8,27.)
for j,(arm,col) in enumerate(zip(ARMS,COLORS)):
    v,lo,hi=[float((summary if q=='accuracy' else two)[arm][q])*100 for q in ['accuracy','ci_low','ci_high']]
    acc.errorbar(v,j,xerr=[[v-lo],[hi-v]],fmt='o',ms=2.6,color=col,capsize=1.,lw=.7)
    for q,val in [('accuracy_percent',v),('two_way_ci_low_percent',lo),('two_way_ci_high_percent',hi)]:rec('d','macro',arm,arm,q,val)
acc.axvline(50,color='#CBD4DB',lw=.5)
acc.set(xlim=(40,78),ylim=(4.5,-.5),yticks=range(5),yticklabels=['Full','Sequence','Flat Φ','Perm Φ','REF rule'],xticks=[40,60])
clean(acc);acc.tick_params(axis='y',length=0,pad=2,labelsize=5.6);acc.tick_params(axis='x',labelsize=5.5)
for t in acc.get_yticklabels():t.set_color(INK)
acc.set_xlabel('95% CI',fontsize=5.5,labelpad=1)
txt(108,134.6,'85 RBPs · 287 sites · development-exposed',5.5,GRAY)

# E: primary site identity and wider, jointly aligned quantitative readouts.
tag(1,139,'e','Disease-annotated variants · 8 sites / 13 events')
DS=HERE/'frozen_inputs/disease_overlap';dr=rows(DS/'DISEASE_ALL13_EVENT_RESULTS.csv');vr=rows(DS/'DISEASE_ALL8_VARIANT_RESULTS.csv')
txt(5,145.1,'Site / gene',5.2,GRAY);txt(30,145.1,'REF>ALT',5.1,GRAY);txt(46,145.1,r'$\Phi_0$',5.6,GRAY,ha='center')
txt(54,145.1,'RBP / cell',5.2,GRAY);txt(76,145.1,'ASB',5.5,INK,ha='center')
earms=['full','sequence','flat','perm'];ekeys=['A_native_ref_minus_alt','S_ref_minus_alt','A_flat_ref_minus_alt','A_permmean_ref_minus_alt']
eclrs=[FULL,SEQ,FLAT,PERM];elabs=['Full','Sequence','Flat Φ','Perm Φ'];ex=[87.,110.,133.,156.]
body_top,body_h,gap=148.8,30.5,.55
pitch=(body_h-gap*7)/13
eaxes=[]
for x,col,lab,key in zip(ex,eclrs,elabs,ekeys):
    ea=ax(x,body_top,20.,body_h);ea.set_xscale('asinh',linear_width=.02)
    ea.axvline(0,color='#9BAAB6',lw=.55,zorder=1)
    ea.set(xlim=(-.12,.88),ylim=(body_h,0),yticks=[],xticks=[-.1,0,.1,.8])
    ea.set_xticklabels(['−.1','0','.1','.8']);ea.xaxis.set_minor_locator(NullLocator());clean(ea)
    ea.tick_params(labelsize=5.5,pad=.8,length=1.3);eaxes.append(ea)
    txt(x+10,144.3,lab,6.0,INK,ha='center');txt(x+10,147.3,f'{sum(float(r[key])>0 for r in dr)}/13',5.0,GRAY,ha='center')
cursor=0.;index=[];ordered=[]
for sid,v in enumerate(vr,1):
    er=[r for r in dr if r['variant_key']==v['variant_key']]
    cy=body_top+cursor+len(er)*pitch/2
    for flag,x,marker,col in [('ClinVar_PLP',1.2,'o','#B58A31'),('TCGA',2.6,'^',FREQ)]:
        if int(v[flag]):fig.add_artist(Line2D([x/W],[1-ys(cy)/H],transform=fig.transFigure,marker=marker,ms=2.6,mfc=col,mec=col,linestyle='none'))
    txt(5,cy,str(sid),5.0,GRAY);txt(9,cy,v['gene'],5.9,weight='bold')
    txt(30,cy,v['population_ref']+'>'+v['population_alt'],5.6,weight='bold')
    txt(46,cy,f'{float(v["reference_PhyloP"]):.2f}',5.3,ha='center')
    rec('e','variant',v['variant_key'],'reference','Phi_center',v['reference_PhyloP'])
    index.append(dict(display_site=sid,**v))
    for k,r in enumerate(er):
        yy=cursor+(k+.5)*pitch
        txt(54,body_top+yy,r['RBP']+(' H' if r['cell']=='HepG2' else ' K'),5.3)
        assert r['observed_ASB_preference']=='REF'
        txt(76,body_top+yy,r['ref'],5.8,INK,ha='center',weight='bold')
        for ea,key,arm,col in zip(eaxes,ekeys,earms,eclrs):
            v0=float(r[key]);agree=v0>0
            ea.barh(yy,v0,height=1.26,color=col if agree else 'white',edgecolor=col,linewidth=.55,hatch=None if agree else '///',zorder=3)
            rec('e','event',r['event_id'],arm,'REF_minus_ALT_log_count',v0)
            rec('e','event',r['event_id'],arm,'source_direction_agreement',int(agree))
        ordered.append(r['event_id'])
    cursor+=len(er)*pitch+(gap if sid<8 else 0)
assert len(ordered)==13 and len(set(ordered))==13 and np.isclose(cursor,body_h)
writecsv(HERE/'CLINICAL_SITE_INDEX.csv',index)
for x,marker,lab,col in [(126,'o','ClinVar P/LP','#B58A31'),(153,'^','TCGA',FREQ)]:
    fig.add_artist(Line2D([x/W],[1-ys(139)/H],transform=fig.transFigure,marker=marker,ms=2.7,mfc=col,mec=col,linestyle='none'));txt(x+2,139,lab,4.9,GRAY)
txt(5,182.7,'11 training / 2 validation',4.9,GRAY);txt(48,182.7,'REF rule: 13/13',5.5,INK)
txt(99,182.7,'REF − ALT Δlog count · asinh',5.5,GRAY)


fig.canvas.draw();ren=fig.canvas.get_renderer();bad=[]
for t in fig.findobj(matplotlib.text.Text):
    if not t.get_visible() or not t.get_text():continue
    bb=t.get_window_extent(ren).transformed(fig.transFigure.inverted())
    if bb.x0<-.005 or bb.y0<-.005 or bb.x1>1.005 or bb.y1>1.005:bad.append(t.get_text())
assert not bad,bad
assert len({(r['panel'],r['kind'],r['id'],r['arm'],r['quantity']) for r in records})==len(records)
baseline=rows(HERE/'PLOTTED_VALUES_V29.csv')
key=lambda r:(r['panel'],r['kind'],r['id'],r['arm'],r['quantity'])
old={key(r):float(r['value']) for r in baseline};new={key(r):float(r['value']) for r in records}
assert old.keys()==new.keys()
assert all(old[k]==new[k] for k in old)
assert len(new)==2028
writecsv(HERE/'PLOTTED_VALUES_V32.csv',records)
for ext in ['pdf','svg','png']:fig.savefig(HERE/f'fig5_v32.{ext}',dpi=400,facecolor='white')
mapping=dict(version='v32',accepted_baseline='v29',native_mm=[W,H],
    panels=['a','b','c','d','e'],all_plotted_values_identical_to_accepted_v29=True,
    plotted_records=len(records),out_of_page_text=bad,
    a_proteins=25,b_nodes=30,b_edges=48,c_model_control_responses=7,c_reporter_measurements=6,
    d_events=396,d_RBPs=85,d_sites=287,d_intervals=5,e_sites=8,e_events=13,e_effects=52,
    e_agreement=[7,10,12,9],constant_REF=13,
    new_analysis=False,new_experiments=False,manuscript_changed=False,Overleaf_changed=False)
(HERE/'SOURCE_MAPPING_V32.json').write_text(json.dumps(mapping,indent=2),encoding='utf8')
print(json.dumps(mapping))
