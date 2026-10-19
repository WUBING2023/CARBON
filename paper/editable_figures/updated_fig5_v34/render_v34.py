"""Figure 5 v34: repeated distribution modules and aligned evidence tracks.

Only display geometry/encoding changes. Every saved numeric record is checked
against the immutable local v32 plotted-value table before export.
"""
from data import *
from matplotlib.ticker import NullLocator
from matplotlib.transforms import ScaledTranslation

W, H = 178., 168.
INK, GRAY, LIGHT = '#111111', '#686E73', '#D9DFE3'
OBS = '#19499C'
BASE = {'A':'#6BB72D','C':'#19499C','G':'#F4A016','U':'#E6211A','T':'#E6211A'}
plt.rcParams.update({'font.family':'Arial','font.size':6.3,
    'axes.labelsize':6.3,'xtick.labelsize':6.3,'ytick.labelsize':6.3,
    'axes.linewidth':.45,'xtick.major.width':.45,'ytick.major.width':.45,
    'xtick.major.size':1.5,'ytick.major.size':1.5,
    'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':INK,
    'xtick.color':INK,'ytick.color':INK,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig = plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
records=[]
coordinate_audit=[]

def ax(x,y,w,h):
    return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],facecolor='none')

def txt(x,y,s,fs=6.3,col=INK,**kw):
    return fig.text(x/W,1-y/H,s,fontsize=fs,color=col,va='center',**kw)

def line(x,y,x2,y2,**kw):
    fig.add_artist(Line2D([x/W,x2/W],[1-y/H,1-y2/H],transform=fig.transFigure,**kw))

def tag(x,y,s,label=None):
    txt(x,y,s,9,'#000000',weight='bold')
    if label is not None:txt(x+5,y,label,6.5)

def letters(x,y,s,step=1.6,fs=6.5):
    for i,b in enumerate(s):
        txt(x+i*step,y,b,fs,BASE[b],ha='center',weight='bold',fontfamily='DejaVu Sans Mono')

def clean(a,left=False):
    for sp in ['top','right']+([] if left else ['left']):a.spines[sp].set_visible(False)
    a.tick_params(length=1.5,pad=1.1,width=.45,colors=INK)
    for sp in a.spines.values():sp.set_linewidth(.45);sp.set_color(INK)

def rec(p,k,id,arm,q,v):
    records.append(dict(panel=p,kind=k,id=str(id),arm=arm,quantity=q,value=float(v)))

def check_values(name,values,limits):
    values=np.asarray(values,dtype=float)
    outside=int(((values<limits[0])|(values>limits[1])|~np.isfinite(values)).sum())
    coordinate_audit.append(dict(name=name,n=int(values.size),minimum=float(values.min()),
        maximum=float(values.max()),axis_limits=list(limits),outside=outside))
    assert outside==0,coordinate_audit[-1]

def frozen_manifest():
    return {str(p.relative_to(HERE)).replace('\\','/'):sha(p)
            for p in sorted((HERE/'frozen_inputs').rglob('*')) if p.is_file()}

input_before=frozen_manifest()
baseline_path=HERE/'PLOTTED_VALUES_V32.csv'
baseline_sha=sha(baseline_path)

# a: two repeated metric modules, actual absolute distributions above paired
# difference strips. No inferred bootstrap distribution is drawn.
tag(1,3.5,'a','RNAcompete · 25 proteins')
jitter=np.array([(int(hashlib.sha256(p.encode()).hexdigest()[:8],16)/0xffffffff-.5)*.17 for p in proteins])
a_audit=[]
for metric,x,mult,limits,ticks,dlimits,dticks,title,rawlabel,dlabel in [
    ('spearman_all5mer',9.,1,(-.7,.86),[-.5,0,.5],(-.7,.75),[-.5,0,.5],
     'All 5-mers','Spearman ρ','Δρ'),
    ('within_composition_concordance',61.,100,(35.,80.),[40,60,80],(-12.,25.),[-10,0,10,20],
     'Same composition','Order agreement (%)','Δorder (pp)')]:
    freq=np.array([float(lookup[(p,'frequency')][metric])*mult for p in proteins])
    full=np.array([float(lookup[(p,'full')][metric])*mult for p in proteins])
    delta=np.array([(float(lookup[(p,'full')][metric])-float(lookup[(p,'frequency')][metric]))*mult for p in proteins])
    mf,mm=float(np.mean(freq)),float(np.mean(full))
    c=next(r for r in contrasts if r['endpoint']==metric and r['control']=='frequency')
    mean,lo,hi=[float(c[k])*mult for k in ['mean_delta','protein_paired_bootstrap_low','protein_paired_bootstrap_high']]
    assert np.isclose(mm-mf,mean,atol=1e-12)
    assert min(freq.min(),full.min())>=limits[0] and max(freq.max(),full.max())<=limits[1]
    assert delta.min()>=dlimits[0] and delta.max()<=dlimits[1]
    check_values('a '+metric+' all raw scores',np.r_[freq,full],limits)
    check_values('a '+metric+' paired differences and saved mean/CI',np.r_[delta,mean,lo,hi],dlimits)
    txt(x+17.5,10.2,title,6.5,ha='center')
    aa=ax(x,15.,35.,22.5)
    for pos,vals,col in [(0,freq,FREQ),(1,full,FULL)]:
        vp=aa.violinplot([vals],positions=[pos],widths=.52,showextrema=False,
                         points=140,bw_method=.3)
        for body in vp['bodies']:
            body.set_facecolor(col);body.set_edgecolor(col);body.set_alpha(.3);body.set_linewidth(.65)
    for i,(p,ff,vv,dd) in enumerate(zip(proteins,freq,full,delta)):
        aa.plot(jitter[i],ff,marker='o',ms=2.6,mfc=FREQ,mec='white',mew=.2,linestyle='none',zorder=3)
        aa.plot(1+jitter[i],vv,marker='o',ms=2.6,mfc=FULL,mec='white',mew=.2,linestyle='none',zorder=3)
        rec('a','protein',p,'full',metric+'_minus_frequency',dd)
        for arm,value in [('full',vv),('frequency',ff)]:
            rec('a','protein',p,arm,metric,lookup[(p,arm)][metric])
            if mult!=1:rec('a','protein',p,arm,metric+'_percent',value)
    aa.plot([-.20,.20],[mf,mf],color=FREQ,lw=.8,zorder=5)
    aa.plot([.80,1.20],[mm,mm],color=FULL,lw=.8,zorder=5)
    aa.set(xlim=(-.52,1.52),ylim=limits,yticks=ticks,xticks=[0,1],xticklabels=['Frequency','Full'])
    clean(aa,True);aa.tick_params(axis='x',length=0,pad=2)
    aa.set_ylabel(rawlabel,labelpad=1.3)
    aa.yaxis.grid(True,color='#ECEFF1',lw=.35);aa.set_axisbelow(True)
    txt(x+17.5,43.1,'Full − Frequency',6.3,GRAY,ha='center')
    da=ax(x,46.0,35.,6.8)
    da.axvline(0,color='#A6ADB3',lw=.4,zorder=0)
    da.scatter(delta,.25+jitter*1.7,s=2.6**2,color='#9BA9B5',edgecolors='white',linewidths=.15,zorder=2)
    da.errorbar(mean,.94,xerr=[[mean-lo],[hi-mean]],fmt='o',ms=3.,color=FULL,capsize=1.5,lw=.75,zorder=4)
    da.set(xlim=dlimits,ylim=(-.10,1.32),xticks=dticks,yticks=[])
    clean(da);da.set_xlabel(dlabel,labelpad=.5)
    for q,value in [('mean',mean),('CI_low',lo),('CI_high',hi)]:rec('a','macro',metric,'Full_minus_Frequency',q,value)
    for arm,value in [('full',mm),('frequency',mf)]:rec('a','macro',metric,arm,'absolute_mean' if mult==1 else 'absolute_mean_percent',value)
    a_audit.append(dict(metric=metric,n=25,absolute_axis=list(limits),delta_axis=list(dlimits),
                        negative_differences=int((delta<0).sum()),violin_source='actual 25 per arm',
                        mean=mean,CI=[lo,hi]))

# b: complete frozen lifted geometry and common response scale retained.
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
for p,q in zip(floor,floor[1:]+floor[:1]):line(*p,*q,color='#E0E5E9',lw=.4)
highlight_edges=0
for r in edges:
    lo,hi=r['lower_assay_word'],r['higher_assay_word'];emphasized=hi=='CUAAC' and lo in highlights
    line(*lift[lo],*lift[hi],color=OBS if emphasized else '#A8B4BF',
         lw=.75 if emphasized else .4,alpha=1 if emphasized else .65)
    highlight_edges+=emphasized
assert len(nodes)==30 and len(edges)==48 and highlight_edges==3
for w in highlights:line(*base[w],*lift[w],color='#D9E1E8',lw=.4)
signed=LinearSegmentedColormap.from_list('navy_white_red',[FREQ,'white',FULL],N=257)
effectnorm=AsinhNorm(linear_width=.02,vmin=-1.8,vmax=1.8)
for w,r in nw.items():
    x,y=lift[w];f=float(r['Full_mean_delta_log_count']);s=float(r['Sequence_mean_delta_log_count'])
    fig.add_artist(Line2D([x/W],[1-y/H],transform=fig.transFigure,marker='o',ms=3.2,fillstyle='left',
        markerfacecolor=signed(effectnorm(f)),markerfacecoloralt=signed(effectnorm(s)),
        markeredgecolor=OBS if w in highlights else '#8F9AA5',markeredgewidth=.65 if w in highlights else .45,linestyle='none'))
    rec('b','word',w,'assay','Z5',r['external_Z5'])
    for arm,value in [('full',f),('sequence',s)]:rec('b','word',w,arm,'mean_delta_log_count_vs_native',value)
check_values('b node-half responses',
    [float(r[k]) for r in nodes for k in ['Full_mean_delta_log_count','Sequence_mean_delta_log_count']],[-1.8,1.8])
for w,(dx,dy) in {'CUAAC':(-4.5,-3.7),'UCAAC':(-6.1,-3.2),'CAUAC':(-9.6,2.4),'CUACA':(2,4)}.items():
    letters(lift[w][0]+dx,lift[w][1]+dy,w,1.6,6.5)
txt(108,11,r'RNAcompete $Z_5$',6.3)
line(108.4,49,108.4,22.5,color=INK,lw=.45)
for z in [0,1,2]:
    line(107.9,49-12*z,109.1,49-12*z,color=INK,lw=.45)
    txt(107.2,49-12*z,str(z),6.3,ha='right')
txt(148,9.8,'Full',6.3,ha='center');txt(155.5,9.8,'|',6.3,ha='center');txt(165,9.8,'Sequence',6.3,ha='center')
txt(156.5,14,'Node halves: left | right',6.3,GRAY,ha='center')
txt(110,56.3,'Δlog count vs native',6.3)
cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=effectnorm,cmap=signed),cax=ax(146,55.6,29,1.4),orientation='horizontal')
cb.set_ticks([-1,0,1]);cb.set_ticklabels(['−1','0','1']);cb.ax.tick_params(labelsize=6.3,length=1,pad=.5)
cb.ax.minorticks_off();cb.outline.set_visible(False)

# c: compact locus header and two aligned readout plots, distinct endpoints.
RP=HERE/'frozen_inputs/reporter_validation'
case=json.loads((RP/'SELECTED_CASES.json').read_text('utf8'))['K562:UPF1:chr12_52388891_+']
cp=case['prediction'];ce=case['event']
assert cp['variant_key']=='chr12:51995107:A:G' and cp['split']=='training'
tag(1,62.5,'c')
txt(6,62.5,'ACVR1B',7,weight='bold',style='italic')
txt(26,62.5,'chr12:51995107 A>G',6.5)
txt(6,66.5,'UPF1 · K562 ASB: A preferred',6.3)
cr=ce['context']['reference_sequence_rna_direction_dna'].replace('T','U')
ca=ce['context']['minor_sequence_dna'].replace('T','U')
assert cr[50]=='A' and ca[50]=='G'
for yy,allele,seq in [(71.,'A',cr),(74.6,'G',ca)]:
    txt(6,yy,allele,6.3,weight='bold');letters(11,yy,seq[42:59],1.5,6.5)
fig.add_artist(Rectangle(((11+8*1.5-.78)/W,1-75.95/H),1.56/W,6.3/H,
                         transform=fig.transFigure,fc='none',ec=OBS,lw=.65))
txt(43,71.0,'Training locus',6.3,GRAY)
txt(43,74.6,r'Reference $\Phi_0$ −2.29',6.3,GRAY)
txt(1,79.7,'K562 binding prediction',6.5)
txt(67,79.7,'HeLa RNA reporter',6.5)
cm=ax(20,84.,33.,21.5)
mainitems=[('Full','full','A_native_major_minus_minor',FULL,0),
           ('Sequence','sequence','S_major_minus_minor',SEQ,1),
           ('Flat Φ','flat','A_flat_major_minus_minor',FLAT,2),
           ('Perm Φ mean','perm','A_permmean_major_minus_minor',PERM,3)]
for lab,arm,key,col,yy in mainitems:
    v=float(cp[key]);cm.barh(yy,v,height=.60,color=col,edgecolor=col,linewidth=.65,zorder=2)
    rec('c','event',cp['event_id'],arm,'A_minus_G_log_count',v)
for j,marker in enumerate(['o','^','D']):
    arm=f'perm{j}';v=float(cp[f'A_perm{j}_major_minus_minor'])
    cm.plot(v,3+(j-1)*.18,marker=marker,ms=3.,mfc='white',mec=PERM,mew=.7,linestyle='none',zorder=4)
    rec('c','event',cp['event_id'],arm,'A_minus_G_log_count',v)
cm.axvline(0,color='#969FA6',lw=.45,zorder=1)
cm.set(xlim=(-.02,.02),ylim=(3.6,-.6),yticks=[0,1,2,3],yticklabels=[x[0] for x in mainitems],xticks=[-.02,0,.02])
check_values('c all seven model/control responses',
    [float(cp[k]) for _,_,k,_,_ in mainitems]+[float(cp[f'A_perm{j}_major_minus_minor']) for j in range(3)],[-.02,.02])
clean(cm);cm.tick_params(axis='y',length=0,pad=2)
cm.set_xlabel('A − G Δlog count',labelpad=.6)
for t in cm.get_yticklabels():t.set_color(INK)
for x,marker,lab in [(2,'o','0'),(12,'^','1'),(22,'D','2')]:
    fig.add_artist(Line2D([x/W],[1-112.5/H],transform=fig.transFigure,marker=marker,ms=2.8,
                          mfc='white',mec=PERM,mew=.7,linestyle='none'))
    txt(x+2,112.5,lab,6.3,GRAY)
txt(29,112.5,'Fixed Φ permutations',6.3,GRAY)
ra=ax(77,84.,23.,21.5)
for j,allele in enumerate(['A','G']):
    vals=[r for r in case['reporter'] if r['RNA_allele']==allele];assert len(vals)==3
    vv=[float(r['value']) for r in vals]
    ra.plot(j+np.array([-.14,0,.14]),vv,'o',ms=2.8,mfc=OBS,mec='white',mew=.2,zorder=3)
    ra.plot([j-.25,j+.25],[np.mean(vv)]*2,color=OBS,lw=.8,zorder=4)
    for r in vals:rec('c','reporter',cp['variant_key']+'|'+allele+'|'+str(r['biological_replicate']),'published_HeLa',r['endpoint'],r['value'])
ra.set(xlim=(-.5,1.5),ylim=(0,1.65),xticks=[0,1],xticklabels=['A','G'],yticks=[0,.5,1,1.5])
check_values('c all six reporter measurements',[float(r['value']) for r in case['reporter']],[0,1.65])
clean(ra,True);ra.set_ylabel('Normalized RNA',labelpad=1.1)
txt(67,112.5,'qRT-PCR · n = 3 per allele',6.3,GRAY)

# d: one aligned scatter/summary frame with the full exposed population.
tag(107,62.5,'d','Native ASB · 396 events')
for x,col,label in [(108,FULL,'Full only correct'),(141,SEQ,'Sequence only correct')]:
    fig.add_artist(Line2D([x/W],[1-67.3/H],transform=fig.transFigure,marker='o',ms=2.6,mfc=col,mec='none',linestyle='none'))
    txt(x+2,67.3,label,6.3)
sc=ax(114,77.,28.5,28.5)
sc.set_xscale('asinh',linear_width=.02);sc.set_yscale('asinh',linear_width=.02)
fvals=np.array([float(r['A_native_major_minus_minor']) for r in cohort])
svals=np.array([float(r['S_major_minus_minor']) for r in cohort])
assert np.all(fvals!=0) and np.all(svals!=0)
check_values('d all 396 Full scatter coordinates',fvals,[-.2,1.8])
check_values('d all 396 Sequence scatter coordinates',svals,[-.2,1.8])
sc.plot([-.2,1.8],[-.2,1.8],color='#BAC0C5',lw=.4,ls=(0,(2,2)),zorder=0)
for ss,ff,col in [(False,False,'#AAB8C4'),(True,True,'#AAB8C4'),(False,True,FULL),(True,False,SEQ)]:
    ids=np.where(((svals>0)==ss)&((fvals>0)==ff))[0]
    sc.scatter(svals[ids],fvals[ids],s=2.6**2,c=col,alpha=.85 if ss!=ff else .56,
               edgecolors='white',linewidths=.1,zorder=3 if ss!=ff else 2)
census=np.zeros((2,2),int)
for r,f,s in zip(cohort,fvals,svals):
    census[int(s<0),int(f<0)]+=1
    for arm,value in [('full',f),('sequence',s)]:rec('d','event',r['event_id'],arm,'source_preferred_minus_other_log_count',value)
    rec('d','event',r['event_id'],'full_sequence','Full_direction_agreement',int(f>0))
    rec('d','event',r['event_id'],'full_sequence','Sequence_direction_agreement',int(s>0))
assert census.tolist()==[[168,61],[43,124]]
sc.axhline(0,color='#A6ADB3',lw=.4);sc.axvline(0,color='#A6ADB3',lw=.4)
sc.set(xlim=(-.2,1.8),ylim=(-.2,1.8),xticks=[-.1,0,.1,1],yticks=[-.1,0,.1,1])
sc.set_xticklabels(['−.1','0','.1','1']);sc.set_yticklabels(['−.1','0','.1','1'])
sc.xaxis.set_minor_locator(NullLocator());sc.yaxis.set_minor_locator(NullLocator());clean(sc,True)
sc.set_xlabel('Sequence Δlog count',labelpad=.6);sc.set_ylabel('Full Δlog count',labelpad=.5)
for x,y,lab,col,ha in [(.025,.97,'43',FULL,'left'),(.97,.97,'168',INK,'right'),(.025,.025,'124',GRAY,'left'),(.97,.025,'61',SEQ,'right')]:
    sc.text(x,y,lab,transform=sc.transAxes,ha=ha,va='top' if y>.5 else 'bottom',fontsize=6.3,color=col,
            bbox={'fc':'white','ec':'none','pad':.2,'alpha':.88})
for i in range(2):
    for j in range(2):rec('d','census',str(i)+str(j),'full_sequence','event_count',census[i,j])
txt(163,72.4,'RBP-equal\naccuracy (%)',6.3,GRAY,ha='center',linespacing=1.05)
acc=ax(163,77.,14.,28.5)
for j,(arm,col) in enumerate(zip(ARMS,COLORS)):
    v,lo,hi=[float((summary if q=='accuracy' else two)[arm][q])*100 for q in ['accuracy','ci_low','ci_high']]
    check_values('d '+arm+' mean and two CI endpoints',[v,lo,hi],[40,78])
    acc.errorbar(v,j,xerr=[[v-lo],[hi-v]],fmt='o',ms=2.8,color=col,capsize=1.2,lw=.75)
    for q,val in [('accuracy_percent',v),('two_way_ci_low_percent',lo),('two_way_ci_high_percent',hi)]:rec('d','macro',arm,arm,q,val)
acc.axvline(50,color='#BDC4C9',lw=.4,zorder=0)
acc.set(xlim=(40,78),ylim=(4.5,-.5),yticks=range(5),yticklabels=['Full','Sequence','Flat Φ','Perm Φ','REF rule'],xticks=[40,60])
clean(acc);acc.tick_params(axis='y',length=0,pad=2)
acc.set_xlabel('95% CI',labelpad=.6)
for t in acc.get_yticklabels():t.set_color(INK)
txt(107,112.5,'85 RBPs · 287 sites · development-exposed',6.3,GRAY)

# e: true site groups define a single aligned identity-to-effect track atlas.
tag(1,117.5,'e','Disease-annotated variants · 8 sites / 13 events')
DS=HERE/'frozen_inputs/disease_overlap'
dr=rows(DS/'DISEASE_ALL13_EVENT_RESULTS.csv');vr=rows(DS/'DISEASE_ALL8_VARIANT_RESULTS.csv')
txt(5,123.4,'Site / gene',6.3,GRAY);txt(30.5,123.4,'REF>ALT',6.3,GRAY)
txt(45.5,123.4,r'$\Phi_0$',6.5,GRAY,ha='center')
txt(52,123.4,'RBP / cell',6.3,GRAY);txt(79,123.4,'ASB',6.3,GRAY,ha='center')
earms=['full','sequence','flat','perm']
ekeys=['A_native_ref_minus_alt','S_ref_minus_alt','A_flat_ref_minus_alt','A_permmean_ref_minus_alt']
eclrs=[FULL,SEQ,FLAT,PERM];elabs=['Full','Sequence','Flat Φ','Perm Φ'];ex=[87.,110.,133.,156.]
body_top,body_h,gap,pitch=127.,34.25,.25,2.5
eaxes=[]
agreement=[sum(float(r[k])>0 for r in dr) for k in ekeys]
assert agreement==[7,10,12,9]
for x,col,lab,key,n in zip(ex,eclrs,elabs,ekeys,agreement):
    ea=ax(x,body_top,20.,body_h);ea.set_xscale('asinh',linear_width=.02)
    ea.axvline(0,color='#9CA6AE',lw=.45,zorder=1)
    ea.set(xlim=(-.12,.88),ylim=(body_h,0),yticks=[],xticks=[-.1,0,.1,.8])
    ea.set_xticklabels(['−.1','0','.1','.8']);ea.xaxis.set_minor_locator(NullLocator());clean(ea)
    ea.tick_params(pad=.8,length=1.3);eaxes.append(ea)
    txt(x+10,122.0,lab,6.5,col,ha='center');txt(x+10,125.0,f'{n}/13',6.3,GRAY,ha='center')
    check_values('e '+key+' all 13 bar endpoints',[float(r[key]) for r in dr],[-.12,.88])
cursor=0.;ordered=[];e_negative=[];e_site_rows=[]
e_bar_physical_lengths=[];e_subprint_markers=[]
sign_threshold_mm,sign_offset_mm=.3,.8
for sid,v in enumerate(vr,1):
    er=[r for r in dr if r['variant_key']==v['variant_key']]
    block_h=len(er)*pitch
    cy=body_top+cursor+block_h/2
    if sid%2==0:
        fig.add_artist(Rectangle((4/W,1-(body_top+cursor+block_h)/H),173/W,block_h/H,
                                 transform=fig.transFigure,fc='#F1F4F6',ec='none',zorder=-1))
    for flag,x,marker,col in [('ClinVar_PLP',1.2,'o',GRAY),('TCGA',2.6,'^',FREQ)]:
        if int(v[flag]):
            fig.add_artist(Line2D([x/W],[1-cy/H],transform=fig.transFigure,marker=marker,ms=2.6,mfc=col,mec='none',linestyle='none'))
    txt(5,cy,str(sid),6.3,GRAY);txt(8.5,cy,v['gene'],6.3,weight='bold',style='italic')
    txt(30.5,cy,v['population_ref']+'>'+v['population_alt'],6.3)
    txt(45.5,cy,f'{float(v["reference_PhyloP"]):.2f}',6.3,ha='center')
    rec('e','variant',v['variant_key'],'reference','Phi_center',v['reference_PhyloP'])
    e_site_rows.append(dict(site=sid,gene=v['gene'],n_events=len(er),first_y_mm=body_top+cursor,last_y_mm=body_top+cursor+block_h))
    for k,r in enumerate(er):
        yy=cursor+(k+.5)*pitch
        txt(52,body_top+yy,r['RBP']+(' H' if r['cell']=='HepG2' else ' K'),6.3)
        assert r['observed_ASB_preference']=='REF'
        txt(79,body_top+yy,r['ref'],6.3,OBS,ha='center')
        for ea,key,arm,col in zip(eaxes,ekeys,earms,eclrs):
            v0=float(r[key]);agree=v0>0
            assert -.12<=v0<=.88
            ea.barh(yy,v0,height=1.52,color=col,edgecolor=col,linewidth=.65,zorder=3)
            # Same physical rule for every arm/event. A fixed offset marker
            # denotes sign only; the original bar and effect coordinate stay.
            x0_px=ea.transData.transform((0,yy))[0]
            xv_px=ea.transData.transform((v0,yy))[0]
            physical_length_mm=abs(xv_px-x0_px)/fig.dpi*25.4
            e_bar_physical_lengths.append(dict(event=r['event_id'],arm=arm,
                value=v0,length_mm=float(physical_length_mm)))
            if v0!=0 and physical_length_mm<sign_threshold_mm:
                direction=1 if v0>0 else -1
                shift=ScaledTranslation(direction*sign_offset_mm/25.4,0,fig.dpi_scale_trans)
                ea.plot(0,yy,transform=ea.transData+shift,marker='>' if direction>0 else '<',
                        ms=2.8,mfc='white',mec=col,mew=.65,ls='none',zorder=5)
                axis_left_px=ea.transAxes.transform((0,0))[0]
                marker_center_mm=(x0_px-axis_left_px)/fig.dpi*25.4+direction*sign_offset_mm
                marker_half_mm=(2.8+.65)*.5*25.4/72
                check_values('e '+r['event_id']+' '+arm+' sign-only symbol extent (mm)',
                             [marker_center_mm-marker_half_mm,marker_center_mm+marker_half_mm],[0,20])
                e_subprint_markers.append(dict(event=r['event_id'],variant=r['variant_key'],
                    arm=arm,value=v0,bar_length_mm=float(physical_length_mm),
                    marker='right_hollow_triangle' if direction>0 else 'left_hollow_triangle',
                    zero_relative_offset_mm=direction*sign_offset_mm,sign_only=True,
                    marker_center_from_axis_left_mm=float(marker_center_mm),marker_extent_inside_axis=True))
            rec('e','event',r['event_id'],arm,'REF_minus_ALT_log_count',v0)
            rec('e','event',r['event_id'],arm,'source_direction_agreement',int(agree))
            if not agree:e_negative.append(dict(event=r['event_id'],arm=arm,value=v0))
        ordered.append(r['event_id'])
    cursor+=block_h+(gap if sid<8 else 0)
assert len(ordered)==13 and len(set(ordered))==13 and np.isclose(cursor,body_h)
for x,marker,lab,col in [(126,'o','ClinVar P/LP',GRAY),(155,'^','TCGA',FREQ)]:
    fig.add_artist(Line2D([x/W],[1-117.5/H],transform=fig.transFigure,marker=marker,ms=2.7,mfc=col,mec='none',linestyle='none'))
    txt(x+2,117.5,lab,6.3,GRAY)
txt(5,166.0,'11 training / 2 validation',6.3,GRAY)
txt(49,166.0,'REF rule: 13/13',6.3)
txt(99,166.0,'REF − ALT Δlog count · asinh',6.3,GRAY)
for x,marker in [(139.5,'<'),(142.5,'>')]:
    fig.add_artist(Line2D([x/W],[1-166./H],transform=fig.transFigure,marker=marker,
                         ms=2.8,mfc='white',mec=GRAY,mew=.65,ls='none'))
txt(145.2,166.0,'Sign only <0.3 mm',6.3,GRAY)
assert len(e_bar_physical_lengths)==52
assert len(e_subprint_markers)==sum(r['value']!=0 and r['length_mm']<sign_threshold_mm for r in e_bar_physical_lengths)

# The baseline has the original recording order, enabling byte-level auditing
# in addition to exact key/value comparisons.
fig.canvas.draw();ren=fig.canvas.get_renderer();bad=[]
for t in fig.findobj(matplotlib.text.Text):
    if not t.get_visible() or not t.get_text():continue
    bb=t.get_window_extent(ren).transformed(fig.transFigure.inverted())
    if bb.x0<-.005 or bb.y0<-.005 or bb.x1>1.005 or bb.y1>1.005:bad.append(t.get_text())
assert not bad,bad
baseline=rows(baseline_path)
key=lambda r:(r['panel'],r['kind'],r['id'],r['arm'],r['quantity'])
old={key(r):float(r['value']) for r in baseline};new={key(r):float(r['value']) for r in records}
assert len(records)==len(new)==len(old)==2028
assert old.keys()==new.keys()
assert all(old[k]==new[k] for k in old)
assert [key(r) for r in baseline]==[key(r) for r in records]
writecsv(HERE/'PLOTTED_VALUES_V34.csv',records)
assert sha(HERE/'PLOTTED_VALUES_V34.csv')==baseline_sha
for ext in ['pdf','svg','png']:fig.savefig(HERE/f'fig5_v34.{ext}',dpi=400,facecolor='white')
pdf=fitz.open(HERE/'fig5_v34.pdf')
pdf_mm=[pdf[0].rect.width*25.4/72,pdf[0].rect.height*25.4/72]
assert np.allclose(pdf_mm,[W,H],atol=.001)
pdf[0].get_pixmap(dpi=160).save(HERE/'fig5_v34_pdf_160dpi.png')
pdf.close()
input_after=frozen_manifest()
assert input_before==input_after
assert sha(baseline_path)==baseline_sha
mapping=dict(version='v34',comparison_baseline='v32',native_mm=[W,H],pdf_measured_mm=pdf_mm,
    panels=['a','b','c','d','e'],plotted_records=len(records),
    exact_key_set_match=True,exact_value_match=True,same_record_order=True,
    coordinate_audit=coordinate_audit,all_original_points_and_CI_endpoints_within_axes=True,
    baseline_sha256=baseline_sha,output_value_sha256=sha(HERE/'PLOTTED_VALUES_V34.csv'),
    byte_identical_plotted_csv=True,out_of_page_text=bad,
    frozen_input_files=len(input_before),frozen_input_hashes_unchanged=input_before==input_after,
    frozen_input_manifest=input_after,
    a_modules=a_audit,b_nodes=30,b_edges=48,b_highlight_edges=highlight_edges,
    b_geometry='same x/y footprint and Z5 lift as v32; no data reorientation',
    b_response_scale=dict(palette=[FREQ,'#FFFFFF',FULL],linear_width=.02,limits=[-1.8,1.8]),
    c_model_control_responses=7,c_reporter_measurements=6,c_locus_split=cp['split'],
    c_endpoints=['K562 binding prediction','published HeLa normalized RNA'],
    c_readout_plot_mm=dict(prediction=[20,84,33,21.5],reporter=[77,84,23,21.5]),
    d_events=396,d_RBPs=85,d_sites=287,d_intervals=5,d_census=census.tolist(),
    d_population='development-exposed ASB-positive events',
    d_plot_mm=dict(scatter=[114,77,28.5,28.5],summary=[163,77,14,28.5]),
    e_sites=8,e_events=13,e_effects=52,e_agreement=agreement,constant_REF=13,
    e_event_pitch_mm=pitch,e_site_groups=e_site_rows,e_negative_effects=e_negative,
    e_scale=dict(linear_width=.02,limits=[-.12,.88],same_for_all_four_arms=True),
    e_bar_encoding='model colour for both signs; direction given by common zero and horizontal coordinates; no hatch',
    e_bar_physical_lengths=e_bar_physical_lengths,
    e_subprint_sign_rule=dict(threshold_mm=sign_threshold_mm,rule='nonzero and physical bar length below threshold, all 52 effects checked',
        marker_offset_from_zero_mm=sign_offset_mm,raw_bar_lengths_unchanged=True,
        markers=e_subprint_markers,legend='hollow left/right triangles: Sign only <0.3 mm'),
    e_splits={'training':11,'validation':2,'strict_test':0},
    visual_review_note='VISUAL_REVIEW_V34.md; actual full PNG and PDF raster inspected in two layout rounds',
    typography=dict(ordinary_pt=[6.3,6.5],panel_letter_pt=9,axis_lw_pt=.45,
                    baseline_lw_pt=[.4,.45],data_outline_lw_pt=[.65,.8],marker_diameter_pt=[2.6,3.2]),
    user_approval_claimed=False,new_analysis=False,new_experiments=False,
    manuscript_changed=False,Overleaf_changed=False,
    exported_artifacts={f'fig5_v34.{ext}':sha(HERE/f'fig5_v34.{ext}') for ext in ['pdf','svg','png']})
(HERE/'SOURCE_MAPPING_V34.json').write_text(json.dumps(mapping,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({k:mapping[k] for k in ['version','native_mm','plotted_records','byte_identical_plotted_csv','out_of_page_text','frozen_input_hashes_unchanged','e_agreement']},ensure_ascii=False))
