"""Reference-guided, vector/native 85-mm SNV figure extension. No invented tracks."""
import os
os.environ['OPENBLAS_NUM_THREADS']='4';os.environ['OMP_NUM_THREADS']='4'
import sys
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
from pathlib import Path
import json,shutil,runpy,xml.etree.ElementTree as ET,hashlib,gzip
import numpy as np,pandas as pd,fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parent;F=R/'figures';F.mkdir(exist_ok=True)
SOURCE=R/'b1_reference'
if not SOURCE.exists():SOURCE=Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/supervised_expansion_20261004/intervention_linked_20261005/direction_b')
P=R/'source_inputs'
if not P.exists():P=Path('M:/RBP/PaperSpine_ABCD_CARBON_20260922/reference_rebuild_20260924/expanded_functional_validation_20261003/06_ASB_expansion')
d=pd.read_csv(R/'EVENTS_WITH_CONSERVATION_AND_PREDICTIONS.csv');s=d[d.split=='strict_test']
selected=s.groupby('variant_key').agg(n_RBP=('RBP','nunique'),phi=('phi_center','first'),chrom=('chrom','first'),pos0=('pos0','first')).reset_index()
selected['chrom_numeric']=selected.chrom.str[3:].astype(int)
selected=selected.sort_values(['n_RBP','phi','chrom_numeric','pos0','variant_key'],ascending=[False,False,True,True,True])
key=selected.variant_key.iloc[0];assert key=='chr8:1844440:C:T'
ex=s[s.variant_key==key].sort_values('RBP');assert len(ex)==10
ex.to_csv(R/'FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv',index=False)
event_file=P/'EVENTS_INCLUDED.json'
events=json.loads(event_file.read_text('utf8') if event_file.exists() else gzip.decompress((P/'EVENTS_INCLUDED.json.gz').read_bytes()).decode('utf8'))
e=next(e for e in events if e['event_id']==ex.event_id.iloc[0])
FULL,SEQ,FLAT,PERM,REF='#E64B35','#4DBBD5','#00A087','#8491B4','#3C5488'
BASE={'A':REF,'C':FLAT,'G':'#7E6148','U':FULL,'T':FULL}
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.labelsize':6.5,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.linewidth':.55,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
def make(h):
    fig=plt.figure(figsize=(85/25.4,h/25.4),facecolor='white')
    def t(x,y,z,fs=6.5,**kw):return fig.text(x/85,1-y/h,z,fontsize=fs,va='center',**kw)
    def ax(x,y,w,hh):return fig.add_axes([x/85,1-(y+hh)/h,w/85,hh/h])
    return fig,t,ax
def save(fig,name):
    for ext in ['pdf','svg','png']:fig.savefig(F/(name+'.'+ext),dpi=600,facecolor='white')
    plt.close(fig)
fig,t,ax=make(49)
t(2.5,3,'a',fs=9,fontweight='bold');t(8,3,'ARHGEF10 · K562',fs=7.5,fontweight='bold')
t(8,6.3,'chr8:1,844,440  C>T')
seq=e['context']['reference_sequence_rna_direction_dna'][45:56].replace('T','U')
for j,b in enumerate(seq):t(8+2.05*j,10,b,fs=7,color=BASE[b],ha='center',fontweight='bold' if j==5 else 'normal')
t(8+5*2.05,12.9,'C→U',fs=6.5,ha='center')
ph=ax(45,4,37,7.5);phi=np.array(e['context']['phi_raw']);ph.plot(np.arange(1,102),phi,color=PERM,lw=.65)
ph.axvline(51,color='#B7BEC5',lw=.6,ls=(0,(2,2)));ph.scatter([51],[phi[50]],s=6,c=FULL,zorder=3)
ph.set_xlim(1,101);ph.set_xticks([]);ph.set_yticks([0]);ph.tick_params(length=1.5,pad=1,width=.4);ph.spines[['top','right','bottom']].set_visible(False)
t(64,2.0,'Φ fixed · center −3.50',ha='center',color=PERM)
t(8,17,'RBP');t(32,17,'BEAPR',ha='center')
t(53,17,'Full',color=FULL);t(69,17,'Sequence',color=SEQ)
chart=ax(45,19.5,37,22.5);chart.set_xlim(-.04,.04);chart.set_ylim(-.5,9.5);chart.set_yticks([])
chart.axvline(0,color='#9CA6AD',lw=.6,ls=(0,(2,2)));chart.spines[['top','right','left']].set_visible(False);chart.set_xticks([-.04,0,.04]);chart.tick_params(length=2,pad=1,width=.5)
for i,(_,row) in enumerate(ex.iterrows()):
    yc=9-i;ym=19.5+(9.5-yc)/10*22.5
    t(8,ym,row.RBP);t(32,ym,row.ref if row.ref_preferred else row.alt,ha='center',fontweight='bold')
    a,b=row.A_native_ref_minus_alt,row.S_ref_minus_alt
    chart.plot([a,b],[yc,yc],c='#BEC5C9',lw=.55,zorder=1)
    chart.scatter([a],[yc],s=8,c=FULL,zorder=3);chart.scatter([b],[yc],s=9,c=SEQ,marker='|',lw=1,zorder=4)
t(63.5,47.0,'Pred. Δlog count (REF−ALT)',ha='center')
save(fig,'FIG5_PANEL_a_NATURAL_SNV')

summ=pd.read_csv(R/'DIRECTION_SUMMARY.csv');summ=summ[summ.cohort=='strict_test'].set_index('arm')
ass=pd.read_csv(R/'CONSERVATION_ASSOCIATIONS.csv');ass=ass[(ass.cohort=='strict_test') & (ass.outcome=='credit')].set_index('arm')
fig,t,ax=make(42)
t(2.5,3,'e',fs=9,fontweight='bold')
colors=[FULL,SEQ,FLAT,PERM,REF];arms=['A_native','S','A_flat','A_permmean','ref'];labels=['Full','Sequence','Flat Φ','Perm. Φ','REF rule']
t(43,3,'396 events · 287 loci · 85 RBPs',ha='center')
left=ax(20,8,26,25);right=ax(57,8,25,25)
for p in [left,right]:
    p.set_ylim(-.5,4.5);p.set_yticks([]);p.spines[['top','right','left']].set_visible(False);p.tick_params(length=2,pad=1.5,width=.5)
left.set_xlim(40,80);left.set_xticks([40,60,80]);right.set_xlim(-.19,.19);right.set_xticks([-.15,0,.15]);right.axvline(0,color='#9CA6AD',lw=.6,ls=(0,(2,2)))
for i,(a,label,c) in enumerate(zip(arms,labels,colors)):
    y=4-i;t(8,8+(4.5-y)/5*25,label)
    r=summ.loc[a];left.plot([100*r.ci_low,100*r.ci_high],[y,y],c=c,lw=.8);left.scatter([100*r.accuracy],[y],s=12,c=c,zorder=3)
    if a!='ref':
        r=ass.loc[a];right.plot([r.ci_low,r.ci_high],[y,y],c=c,lw=.8);right.scatter([r.adjusted_standardized_rank_slope],[y],s=12,c=c,zorder=3)
t(33,38.2,'Direction accuracy (%)',ha='center');t(69.5,38.2,'Adjusted Φ slope',ha='center')
save(fig,'FIG5_PANEL_e_COHORT_CONSERVATION')

# Re-render only panel IDs of the accepted B1 in our OWN directory, retaining its sources/layout/data.
B=R/'b1_preserved';B.mkdir(exist_ok=True)
shutil.copyfile(SOURCE/'INTERVENTION_PACKET.npz',B/'INTERVENTION_PACKET.npz')
source=(SOURCE/'render.py').read_text('utf8')
source=source.replace("text(2.5,4,'a'","text(2.5,4,'b'").replace("text(2.5,64,'b'","text(2.5,64,'c'").replace("text(2.5,98,'c'","text(2.5,98,'d'")
(B/'render.py').write_text(source,'utf8');runpy.run_path(str(B/'render.py'))
def combine(name,parts):
    pt=72/25.4;h=sum(mm for _,mm in parts);doc=fitz.open();page=doc.new_page(width=85*pt,height=h*pt);y=0
    for p,mm in parts:
        src=fitz.open(p.with_suffix('.pdf'));page.show_pdf_page(fitz.Rect(0,y*pt,85*pt,(y+mm)*pt),src,0);src.close();y+=mm
    doc.save(F/(name+'.pdf'));page.get_pixmap(matrix=fitz.Matrix(600/72,600/72),alpha=False).save(F/(name+'.png'));doc.close()
    ns='http://www.w3.org/2000/svg';ET.register_namespace('',ns)
    root=ET.Element('{'+ns+'}svg',{'width':'85mm','height':str(h)+'mm','viewBox':f'0 0 {85*pt} {h*pt}'})
    y=0
    for index,(p,mm) in enumerate(parts):
        part=ET.parse(p.with_suffix('.svg')).getroot()
        ids={el.attrib['id']:f'p{index}_'+el.attrib['id'] for el in part.iter() if 'id' in el.attrib}
        for el in part.iter():
            if 'id' in el.attrib:el.attrib['id']=ids[el.attrib['id']]
            for attr,val in list(el.attrib.items()):
                if attr.endswith('href') and val.startswith('#'):el.attrib[attr]='#'+ids.get(val[1:],val[1:])
                elif 'url(#' in val:
                    for old,new in ids.items():val=val.replace('url(#'+old+')','url(#'+new+')')
                    el.attrib[attr]=val
        wrapper=ET.SubElement(root,'{'+ns+'}g',{'transform':f'translate(0,{y*pt})'})
        # Native coordinates are points; preserve vector labels and marks.
        for child in list(part):wrapper.append(child)
        y+=mm
    ET.ElementTree(root).write(F/(name+'.svg'),encoding='utf-8',xml_declaration=True)
combine('FIG5_EXPANDED_85mm',[(F/'FIG5_PANEL_a_NATURAL_SNV',49),(B/'figure',132),(F/'FIG5_PANEL_e_COHORT_CONSERVATION',42)])
combine('FIG5_ADDITIONAL_a_e_85mm',[(F/'FIG5_PANEL_a_NATURAL_SNV',49),(F/'FIG5_PANEL_e_COHORT_CONSERVATION',42)])
mapping=dict(width_mm=85,expanded_height_mm=223,additional_height_mm=91,natural_example=key,natural_selection='max distinct-RBP coverage; Phi and coordinate tie breakers; no model outcomes',
  manuscript_layout={'first_float':'original accepted B1 a/b/c, 85x132mm, byte-exact b1_reference/figure.pdf','continued_float':'FIG5_ADDITIONAL_d_e_85mm.pdf, natural variant d and cohort/conservation e, 85x91mm','expanded_composite':'engineering preview only, not manuscript'},
  preserve_B1='original 132-mm accepted panel geometry and all numeric data preserved; original a/b/c relabelled b/c/d only',
  reference_views=['PrismNet Figure4: linked sequence/perturbation/prediction/assay readout','PrismNet Figure7: natural SNV localized at shared coordinates with downstream RBP readout'],
  transferred='same-locus sequence/variant/Phi linked to paired model response and actual BEAPR direction; no structure assay, pull-down, odds ratio or causal conclusion inherited',
  source_files={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE/'render.py',SOURCE/'INTERVENTION_PACKET.npz',R/'FIG5_NATURAL_SNV_EXAMPLE_ALL10.csv',R/'DIRECTION_SUMMARY.csv',R/'CONSERVATION_ASSOCIATIONS.csv']})
(R/'FIGURE_SOURCE_MAPPING.json').write_text(json.dumps(mapping,indent=2)+'\n','utf8')
print('Rendered 85 x223 mm expanded; standalone 85x91 mm natural/conservation addition.')
