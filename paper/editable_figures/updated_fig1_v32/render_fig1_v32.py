"""Native 85 x 170 mm Figure 1. Architecture and interpretation schematics only.
All computation, dimensions and component roles follow the current Methods.
The two motif-logo exemplars are the original PPTX SVG vectors, cropped only
outside the letter geometry; no motif values or new model outputs are created.
"""
from pathlib import Path
import json, shutil, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle, FancyArrowPatch, PathPatch
from matplotlib.path import Path as MPath
from matplotlib.colors import to_rgb
import fitz

O=Path(__file__).resolve().parent
P=json.loads((O/'frozen_inputs/PALETTE_V3.json').read_text())
BASE=P['nucleotides']; NAVY=P['observed']; RED=P['models']['Full']; CYAN=P['models']['Sequence']
INK='#303A42'; MID='#74838D'; GRID='#D7DFE4'; BLUE='#CBD6E2'; PALE='#EDF2F5'; SALMON='#F5C4BC'; MINT='#C5E7E2'
W,H=85.,170.
plt.rcParams.update({'font.family':'Arial','font.size':5.5,'axes.linewidth':.45,'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','mathtext.fontset':'dejavusans'})
F=plt.figure(figsize=(W/25.4,H/25.4));A=F.add_axes([0,0,1,1]);A.set(xlim=(0,W),ylim=(H,0));A.axis('off')
TEXT=[];NODES={};EDGES=[]
def tx(x,y,s,fs=5.5,ha='center',va='center',color=INK,bold=False,**kw):
 a=A.text(x,y,s,fontsize=fs,ha=ha,va=va,color=color,fontweight='bold' if bold else 'normal',**kw);TEXT.append(a);return a

def box(id,x,y,w,h,label,fill=PALE,fs=5.5):
 A.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=.65',fc=fill,ec=MID,lw=.5,zorder=3));tx(x+w/2,y+h/2,label,fs=fs);NODES[id]={'bbox_mm':[x,y,w,h],'label':label};return (x,y,w,h)

def edge(id,pts,meaning='forward computation',head=True,lw=.6,col=INK):
 v=np.array(pts,dtype=float); codes=[MPath.MOVETO]+[MPath.LINETO]*(len(v)-1)
 patch=FancyArrowPatch(path=MPath(v,codes),arrowstyle='-|>' if head else '-',mutation_scale=4.5,color=col,lw=lw,zorder=2,capstyle='round',joinstyle='round');A.add_patch(patch);EDGES.append({'id':id,'points_mm':pts,'meaning':meaning});return patch

def dot(x,y):A.add_patch(Circle((x,y),.32,color=INK,zorder=4))
def op(id,x,y,s):
 A.add_patch(Circle((x,y),1.25,fc='white',ec=MID,lw=.5,zorder=4));tx(x,y,s,fs=6.3,zorder=6);NODES[id]={'center_mm':[x,y],'operation':s}
def panel(letter,y,title=None):
 tx(1,y,letter,fs=9,ha='left',bold=True)
 if title:tx(7,y,title,fs=6.0,ha='left',bold=True)
def rule(y):A.plot([1,84],[y,y],color=GRID,lw=.45,zorder=0)
def curve(x,y,w,h,values,col,fill=False,lw=.65):
 v=np.array(values,dtype=float);xx=np.linspace(x,x+w,len(v));yy=y+h-h*v
 if fill:A.fill_between(xx,y+h,yy,color=col,alpha=.32,lw=0,zorder=1)
 A.plot(xx,yy,color=col,lw=lw,zorder=2);A.plot([x,x+w],[y+h,y+h],color=MID,lw=.4,zorder=1)
def sequence(x,y,s,step=.98,highlight=None):
 for i,b in enumerate(s):
  if highlight==i:A.add_patch(Rectangle((x+i*step-step*.46,y-1.15),step*.92,2.3,fc=tuple(.22*v+.78 for v in to_rgb(BASE[b])),ec='none',zorder=0))
  tx(x+i*step,y,b,fs=5.5,color=BASE[b],fontfamily='DejaVu Sans Mono',bold=highlight==i)

# a: five-channel input -> MSME -> RCE -> profile/count factorisation.
panel('a',2.4)
tx(10,5,'101-nt input',bold=True)
seq='AACGUAAC'
for r,base in enumerate('ACGU'):
 for c,b in enumerate(seq):
  fc=BASE[base] if base==b else '#FFFFFF';fc=tuple(.33*v+.67 for v in to_rgb(fc))
  A.add_patch(Rectangle((3+c*1.75,8+r*1.8),1.7,1.75,fc=fc,ec='#9AA8B0',lw=.25))
tx(10,17.3,'Sequence · 4')
phy=[.05,.08,.04,.14,.55,.11,.13,.89,.21,.08,.2,.09,.12,.35,.07,.02]
curve(3,20.7,14,3.8,phy,'#8491B4',fill=True)
tx(10,26.8,'PhyloP · 1')
edge('a_sequence_to_input',[(17,11.6),(19.3,11.6),(19.3,17.5)],head=False)
edge('a_phylo_to_input',[(17,22.6),(19.3,22.6),(19.3,17.5)],head=False);dot(19.3,17.5)
box('a_MSME',24,13,14,9,'MSME\n5 → 96',BLUE)
box('a_RCE',43,13,14,9,'RCE\n96 → 128',PALE)
edge('a_input_MSME',[(19.3,17.5),(24,17.5)])
edge('a_MSME_RCE',[(38,17.5),(43,17.5)])
box('a_profile',62,6,20,7,'Profile · softmax',PALE)
box('a_count',62,19,20,7,'Count · exp(z)',MINT)
edge('a_RCE_profile',[(57,17.5),(59.8,17.5),(59.8,9.5),(62,9.5)])
edge('a_RCE_count',[(59.8,17.5),(59.8,22.5),(62,22.5)]);dot(59.8,17.5)
op('a_multiply',72,29,'×')
edge('a_profile_multiply',[(82,9.5),(84,9.5),(84,29),(73.25,29)])
edge('a_count_multiply',[(72,26),(72,27.75)])
edge('a_output',[(70.75,29),(61,29)])
curve(45,25.4,15,4.3,[0,.04,.08,.12,.25,.17,.1,.07,.2,.95,.35,.11,.03,.1,.3,.1,.02],RED,True)
tx(51.8,23.7,'Predicted profile',color=RED)
rule(31.5)

# b: unpooled V is gated, fused and added to the independent five-channel shortcut.
panel('b',34,'Multi-scale motif extraction')
for r in range(5):
 for c in range(4):A.add_patch(Rectangle((3+c*1.6,45+r*1.6),1.55,1.55,fc=BLUE if (r+c)%3 else '#9BAFC7',ec='white',lw=.2))
tx(6.1,55.4,'5 ch')
for k,yy in [(3,39),(5,48),(7,57)]:
 box(f'b_conv{k}',18,yy,13,5,f'k = {k} · 32',BLUE)
 edge(f'b_input_conv{k}',[(9.4,49),(13.2,49),(13.2,yy+2.5),(18,yy+2.5)])
 edge(f'b_conv{k}_concat',[(31,yy+2.5),(36,yy+2.5)])
box('b_concat',36,39,6,23,'',BLUE)
tx(39,50.5,'Concat · 96',rotation=90)
box('b_SE',46,38,34,8,'SE · global mean · 96\n24 ReLU → 96 sigmoid',SALMON)
edge('b_global_summary',[(42,50.5),(44,50.5),(44,42),(46,42)])
op('b_gate',61,50.5,'×')
edge('b_unpooled_to_gate',[(42,50.5),(59.75,50.5)])
tx(51.5,48.5,'Unpooled',color=MID)
edge('b_gate_weights',[(61,46),(61,49.25)])
box('b_fusion',66,47,15,7,'1×1 fusion\n96 → 96',BLUE)
edge('b_gate_fusion',[(62.25,50.5),(66,50.5)])
box('b_shortcut',37,62.4,25,5.4,'1×1 shortcut · 5 → 96',fill='#FFFFFF')
edge('b_input_shortcut',[(11.2,49),(11.2,65.1),(37,65.1)]);dot(11.2,49)
op('b_add',73.5,65.1,'+')
edge('b_fusion_add',[(73.5,54),(73.5,63.85)])
edge('b_shortcut_add',[(62,65.1),(72.25,65.1)])
edge('b_output',[(74.75,65.1),(83,65.1)])
tx(80.5,67.5,'F · 96')
rule(69.5)

# c: full recurrent residual, then separate profile and scalar count heads.
panel('c',71.5,'Context encoding and dual heads')
tx(5.7,77.3,'F · 96')
box('c_projection',13,76,14,7,'1×1\n96 → 64',BLUE)
box('c_BiGRU',31,75,16,10,'',PALE)
tx(39,77,'BiGRU',bold=True);tx(39,79.8,'64 each way')
edge('c_forward_direction',[(33,82.0),(45,82.0)],lw=.45).set_zorder(5)
edge('c_reverse_direction',[(45,83.8),(33,83.8)],lw=.45).set_zorder(5)
box('c_concat',51,76,12,7,'Concat\n128',BLUE)
edge('c_input_projection',[(7.5,79.5),(13,79.5)])
edge('c_projection_BiGRU',[(27,79.5),(31,79.5)])
edge('c_BiGRU_concat',[(47,79.5),(51,79.5)])
op('c_add',66,79.5,'+')
edge('c_concat_add',[(63,79.5),(64.75,79.5)])
box('c_shortcut',20,86,27,7,'1×1 shortcut\n96 → 128',fill='#FFFFFF')
edge('c_input_shortcut',[(8.3,79.5),(8.3,89.5),(20,89.5)]);dot(8.3,79.5)
edge('c_shortcut_add',[(47,89.5),(66,89.5),(66,80.75)])
box('c_LN_dropout',70,75.5,13,8,'LayerNorm\nDrop 0.2',PALE)
edge('c_add_LN',[(67.25,79.5),(70,79.5)])
edge('c_LN_to_heads',[(76.5,83.5),(76.5,95),(29.5,95),(29.5,97.5)])
edge('c_LN_to_count',[(76.5,95),(76.5,97.5)]);dot(76.5,95)
box('c_profile_head',10,97.5,39,6,'Profile · 1×1 conv\n128 → 64 → 32 → 1',BLUE)
box('c_count_head',54,97.5,29,6,'Count · global mean\n128 → 32 → 1',MINT)
edge('c_logits_softmax',[(29.5,103.5),(29.5,104.2)])
edge('c_z_exp',[(68.5,103.5),(68.5,104.2)])
tx(29.5,105.6,'Logits → softmax');tx(68.5,105.6,'z → exp(z)')
rule(107.7)

# d: joint objective; curves illustrate supervision, not additional measurements.
panel('d',109.6,'Joint objective')
obs=[0,.02,.01,.04,.12,.15,.7,1,.22,.07,.03,.02,.01,.02]
pred=[0,.02,.03,.08,.2,.72,.91,.35,.11,.09,.03,.01,.01,.02]
tx(12.5,112.0,'Observed');curve(5,113.5,15,4,obs,'#91D1C2',True)
tx(12.5,120.0,'Predicted');curve(5,121.5,15,4,pred,RED,True)
edge('d_observed_bus',[(20,115.5),(24,115.5)],head=False)
edge('d_predicted_bus',[(20,123.5),(24,123.5)],head=False)
edge('d_supervision_bus',[(24,114.5),(24,124.5)],head=False)
for name,y,label,weight in [('shape',112.5,'Profile shape','α'),('CDF',117.5,'Squared CDF','β'),('count',122.5,'Count','γ')]:
 box('d_'+name,29,y,22,4,label,PALE)
 edge('d_bus_'+name,[(24,y+2),(29,y+2)])
 edge('d_weight_'+name,[(51,y+2),(59,y+2),(59,119.5)],head=False)
 tx(54.8,y+.45,weight,fs=5.5)
box('d_total',65,115.2,18,8.6,'Weighted\nsum / s',MINT)
edge('d_joint_total',[(59,119.5),(65,119.5)])
rule(127.7)

# e: same original interpretation workflow and original motif glyph examples.
panel('e',129.5,'Motif recovery')
tx(10.5,132.5,'IG windows')
curve(4,134,13,7,[0,.05,.1,.12,.3,.8,.4,.12,.18,.32,.15,.08,.02],RED,True)
A.add_patch(Rectangle((9,134),3,7,fc='#8491B4',ec='none',alpha=.22,zorder=0))
edge('e_windows_segments',[(17,137.5),(21,137.5)])
tx(31.5,132.5,'Segments')
for y,seq in [(135.1,'GUAACAGCUGU'),(138.1,'GUAGCCGCUGU'),(141.1,'GUAACACCGGU')]:sequence(22,y,seq,step=1.8)
edge('e_segments_motifs',[(42,137.5),(45,137.5)])
box('e_TFMoDISco',45,133.5,17,8.5,'TF-MoDISco',PALE)
edge('e_TFMoDISco_logos',[(62,137.75),(65,137.75)])
tx(74,146.2,'CWM / PPM')
# Vector motif overlays occupy x65..83,y132.8..144.7, added after PDF export.
rule(148.0)

# f: same native/one-base-mutated sequence comparison; all output curves schematic.
panel('f',150.0,'Sequence substitution')
tx(11.0,153.1,'Observed');curve(5,154.5,12,4,[0,.01,.1,.8,1,.2,.12,.04,.02,.07,.01,.02],'#91D1C2',True)
tx(11.0,160.5,'Predicted');curve(5,162,12,4,[0,.02,.1,.9,.8,.2,.15,.08,.01,.03,.01,.03],RED,True)
edge('f_profiles_to_sequence',[(17,158.8),(21,158.8)])
sequence(23,155.7,'CAUGCACAGAU',step=1.65,highlight=5)
sequence(23,162.2,'CAUGCUCAGAU',step=1.65,highlight=5)
edge('f_mutation',[(31.25,157.1),(31.25,160.9)],col=MID,lw=.45)
tx(31.2,166.4,'A → U')
edge('f_sequence_to_prediction',[(40.5,158.8),(47.5,158.8)])
tx(65.5,153.1,'Re-predict')
A.add_patch(Rectangle((54.5,154.8),10,10.8,fc='#8491B4',ec='none',alpha=.18,zorder=0))
mut=[.02,.03,.05,.12,.7,.2,.14,.1,.15,.24,.13,.1,.06,.09,.04,.05,.03,.04,.02,.03]
wt=[.02,.03,.05,.13,1,.22,.15,.1,.15,.23,.14,.11,.06,.09,.04,.05,.03,.04,.02,.03]
curve(49,154.8,33,10.8,wt,RED,lw=.65);curve(49,154.8,33,10.8,mut,CYAN,lw=.65)
A.plot([55.95,55.95],[154.8,165.6],color=INK,lw=.4,ls=(0,(2,2)))
A.plot([53,56],[168.1,168.1],color=RED,lw=.7);tx(56.8,168.1,'WT',ha='left')
A.plot([65,68],[168.1,168.1],color=CYAN,lw=.7);tx(68.8,168.1,'Mutated',ha='left')

# Save editable native drawing; original accepted logo exemplars remain vectors.
F.savefig(O/'fig1_native.pdf',facecolor='white')
F.savefig(O/'fig1_native.svg',facecolor='white')
F.canvas.draw();renderer=F.canvas.get_renderer()
font_checks=[]
for t in TEXT:
 bb=t.get_window_extent(renderer);inv=A.transData.inverted();xy=inv.transform(bb.get_points());font_checks.append({'text':t.get_text(),'font_pt':t.get_fontsize(),'bounds_mm':[[float(x),float(y)] for x,y in xy]})
plt.close(F)
d=fitz.open(O/'fig1_native.pdf');p=d[0];logo=fitz.open(O/'frozen_inputs/motif_letters_only.svg');vlogo=fitz.open('pdf',logo.convert_to_pdf())
# The source is the accepted schematic's two motif examples. Remove outer axes only.
for dest,clip in [((65,132.7,83,138.8),(24,0,444.61,137)),((65,139.0,83,145.1),(24,173.14,444.61,311))]:
 p.show_pdf_page(fitz.Rect(*(v*72/25.4 for v in dest)),vlogo,0,clip=fitz.Rect(*clip),keep_proportion=True)
d.save(O/'fig1_v32.pdf')
(O/'fig1_v32.svg').write_text(p.get_svg_image(text_as_path=False),encoding='utf8')
for dpi in [96,300,600]:p.get_pixmap(dpi=dpi).save(O/f'fig1_v32_85mm_{dpi}dpi.png')
shutil.copy2(O/'fig1_v32_85mm_300dpi.png',O/'fig1_v32.png')
# Retain graph definitions to make scientific direction review concrete.
qa={'dimensions_mm':[85,170],'minimum_native_ordinary_font_pt':min(t['font_pt'] for t in font_checks),'font_objects':font_checks,'nodes':NODES,'edges':EDGES,'nucleotide_palette':BASE,'source_logo':'Original Figure1_PrismNet_nucleotides_v3.pptx ppt/media/image4.svg; all original 94 coloured glyph groups plus 2 baseline groups retained in motif_letters_only.svg; outer axes removed','no_new_training_inference_or_data':True,'role_map':{'a':'101-nt sequence+PhyloP -> MSME -> RCE -> profile softmax * exp(count log z)','b':'parallel k3/5/7,32 each; concat96; SE96->24ReLU->96sigmoid gates unpooled; fusion96 plus input shortcut5->96','c':'96->64; BiGRU64 each direction; concat128 plus shortcut96->128; LayerNorm/dropout.2; profile128->64->32->1 and countmean128->32->1','d':'weighted profile shape, squared CDF, count objective with global scale s','e':'IG windows -> sequence segments -> TF-MoDISco -> original CWM/PPM logo examples','f':'observed/prediction context -> one-base sequence substitution -> rerun comparison; schematic'}}
(O/'FIG1_NATIVE_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
print('Native figure1 85 x 170 mm; smallest native text',qa['minimum_native_ordinary_font_pt'],'pt; original motif vectors preserved')
