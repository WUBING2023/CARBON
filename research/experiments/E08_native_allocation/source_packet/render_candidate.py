from pathlib import Path
import os,sys,json,hashlib
sys.modules['numexpr']=None;sys.modules['bottleneck']=None
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parent;S=R/'nogru_matched';O=R/'figure';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':6.5,'axes.labelsize':6.5,'axes.titlesize':6.5,'xtick.labelsize':6,'ytick.labelsize':6,'axes.linewidth':.55,'lines.linewidth':.85,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
colors={'aligned':'#E64B35','nogru':'#3C5488','observed':'#333333','input':'#B0B0B0','phi':'#8491B4'}
fig=plt.figure(figsize=(85/25.4,119/25.4));sources=[]
selected=json.loads((O/'ILLUSTRATION_SELECTION.json').read_text())['examples'];examples=[v['selected'] for v in selected]
fig.text(.015,.983,'a',fontweight='bold',fontsize=8,va='top')
for ei,v in enumerate(examples):
 task=v['task'];ch=v['chromosome'].removeprefix('chr');src=R/'remote/data'/task if ch=='1' else R/'phase2_ambiguity/remote/data'/task/ch
 with np.load(src/'data.npz') as f:z={k:f[k] for k in f.files}
 meta=json.loads((src/'metadata.json').read_text());i=next(i for i,m in enumerate(meta) if m['index']==v['index'])
 with np.load(S/'remote/data'/task/ch/'data.npz') as f:nz={k:f[k] for k in f.files}
 assert int(nz['index'][i])==v['index']
 data={'observed':z['y'][i],'input':z['control'][i],'aligned':np.mean([z[f'aligned_seed{s}'][i].astype(float) for s in [11,23,47]],axis=0),'nogru':np.mean([nz[f'nogru_seed{s}'][i].astype(float) for s in [11,23,47]],axis=0),'phi':z['x'][i,:,4]}
 x0=.18+ei*.405;width=.34;top=.926
 fig.text(x0, .970,v['gene_name'],fontsize=6.5,fontstyle='italic');fig.text(x0,.954,v['cell'],fontsize=6);fig.text(x0,.939,'1 = 2: CTCT',fontsize=6)
 axes=[]
 for ti,key in enumerate(['observed','input','aligned','nogru','phi']):
  ax=fig.add_axes([x0,top-ti*.055-.045,width,.045]);axes.append(ax)
  for c in [v['c1'],v['c2']]:ax.axvspan(c-5,c+5,color='#EAEAEA',zorder=0)
  xx=np.arange(101);yy=data[key]
  if key in ['observed','input']:ax.fill_between(xx,0,yy,color=colors[key],lw=0)
  else:ax.plot(xx,yy,color=colors[key])
  ax.set_xlim(0,100);ax.spines[['top','right']].set_visible(False);ax.tick_params(length=2,pad=1)
  if key in ['observed','input']:ax.set_ylim(0,max(1,float(np.max(yy)))*1.05)
  if key in ['aligned','nogru']:ax.set_ylim(0,max(np.max(data['aligned']),np.max(data['nogru']))*1.10)
  max_tick=max(np.max(data['aligned']),np.max(data['nogru'])) if key in ['aligned','nogru'] else float(np.max(yy))
  ax.set_yticks([0,max_tick] if max_tick>0 else [0]);ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2g'))
  if ti<4:ax.set_xticks([])
  else:ax.set_xticks([0,50,100]);ax.set_xlabel('Window position (nt)',labelpad=1)
  if ei==0:fig.text(.126,top-ti*.055-.022,{'observed':'eCLIP','input':'Input','aligned':'Full','nogru':'NoGRU','phi':'PhyloP'}[key],ha='right',va='center',fontsize=6.5,color=colors[key])
 for num,c in enumerate([v['c1'],v['c2']],1):
  axes[0].text(c,axes[0].get_ylim()[1]*1.06,str(num),ha='center',fontsize=6,clip_on=False)
 shares={key:float(data[key][v['c1']-5:v['c1']+6].sum()/(data[key][v['c1']-5:v['c1']+6].sum()+data[key][v['c2']-5:v['c2']+6].sum())) for key in ['aligned','nogru']}
 Y=float(data['observed'].sum());C=float(data['input'].sum());m=[float(data['observed'][v[c]-5:v[c]+6].sum()-Y*data['input'][v[c]-5:v[c]+6].sum()/C) for c in ['c1','c2']];shares['observed']=m[0]/sum(m)
 assert abs(shares['aligned']-v['aligned_ensemble_share'])<1e-10 and abs(shares['nogru']-v['nogru_ensemble_share'])<1e-10 and abs(shares['observed']-v['q_adjusted'])<1e-10
 sources.append(dict(**v,native_array_row=i,observed_signal=data['observed'].tolist(),input_signal=data['input'].tolist(),full_profile=data['aligned'].tolist(),nogru_profile=data['nogru'].tolist(),phi=data['phi'].tolist(),plotted_array_recomputed_shares=shares))
fig.text(.015,.619,'b',fontweight='bold',fontsize=8,va='top')
df=pd.read_csv(S/'PAIR_RESULTS_ALL.csv');sd=df[df.adj_good&df.unique_gene]
for ei,v in enumerate(examples):
 row=sd[(sd.task==v['task'])&(sd.chromosome==v['chromosome'])&(sd['index']==v['index'])].iloc[0]
 ax=fig.add_axes([.18+ei*.405,.505,.34,.091]);vals=[row.q_adjusted,row.aligned_ensemble_share,row.nogru_ensemble_share]
 for j,(val,col) in enumerate(zip(vals,[colors['observed'],colors['aligned'],colors['nogru']])):
  ax.plot([0,1],[j,j],color='#DDDDDD',lw=3);ax.plot([0,val],[j,j],color=col,lw=3);ax.plot(val,j,'o',color=col,ms=2.3)
  ax.text(.985,j-.14,f'{val:.2f}',ha='right',va='bottom',fontsize=6)
 ax.set_xlim(0,1);ax.set_ylim(2.5,-.5);ax.set_xticks([0,.5,1]);ax.set_yticks([0,1,2],['eCLIP','Full','NoGRU'] if ei==0 else [])
 ax.spines[['left','top','right']].set_visible(False);ax.tick_params(length=2,pad=2);ax.set_xlabel('Site 1 share',labelpad=1)
fig.text(.015,.455,'c',fontweight='bold',fontsize=8,va='top')
ax=fig.add_axes([.19,.135,.745,.27]);plot=[]
for ci,(cell,col) in enumerate([('K562','#3C5488'),('HepG2','#F39B7F')]):
 sub=sd[sd.cell==cell];g=sub.groupby('gene_id')[['aligned_ensemble_adjusted_loss','nogru_ensemble_adjusted_loss']].mean();delta=g.iloc[:,0]-g.iloc[:,1]
 jit=np.array([int(hashlib.sha256(n.encode()).hexdigest()[:8],16)/2**32*.36-.18 for n in g.index])
 ax.scatter(delta,ci+jit,s=3,color=col,alpha=.38,linewidths=0,rasterized=True)
 ax.plot(float(delta.mean()),ci,'D',ms=3.5,color=col,markeredgecolor='white',markeredgewidth=.4)
 for gene,val in delta.items():plot.append({'cell':cell,'gene_id':gene,'loss_delta':val})
 ax.text(.27,ci+0.31,f'n = {len(g)} genes',ha='right',fontsize=6)
ax.axvline(0,color='#666666',ls='--',lw=.6);ax.set_xlim(-.35,.28);ax.set_ylim(1.48,-.48);ax.set_yticks([0,1],['K562','HepG2']);ax.set_xticks([-.3,-.15,0,.15]);ax.spines[['top','right']].set_visible(False);ax.set_xlabel('Full − NoGRU paired-share error',labelpad=3);ax.tick_params(length=2,pad=2)
con=pd.read_csv(S/'CONTRASTS.csv');rr=con[(con.scene=='all_core')&(con.cohort=='adjusted')&(con.group=='PTBP1')&(con.seed=='ensemble')&(con.control=='nogru')].iloc[0]
fig.text(.19,.058,f"Δ = {rr.gene_equal_delta:.4f}  (95% CI {rr.bonf1_low:.4f}, {rr.bonf1_high:.4f})",fontsize=6.5)
for ext in ['png','pdf','svg']:fig.savefig(O/('CONTEXT_CANDIDATE.'+ext),dpi=450,facecolor='white')
(O/'FIGURE_SOURCE.json').write_text(json.dumps({'examples':sources,'all_cohort':{'pairs':len(sd),'genes':sd.gene_id.nunique(),'task_gene_units':len(plot),'gene_equal_delta':rr.gene_equal_delta,'ci95':[rr.bonf1_low,rr.bonf1_high]},'reference':'RBPNet Fig4 coordinate-linked profile logic, original CARBON layout','dimensions_mm':[85,119],'minimum_font_pt':6,'posthoc_example_selection':True,'position_capacity_control':'training not complete; absent from figure'},indent=2),encoding='utf-8')
(O/'FIGURE_QA.json').write_text(json.dumps({'all_array_shares_match_csv_within':1e-10,'same_probability_ylim_within_each_example':True,'inclusive_11position_supports':True,'examples_selection':'post hoc, median of preexisting two transition groups, fixed before new capacity-control outcomes','minimum_font_pt':6,'required_visual_reinspection_after_this_export':True},indent=2))
pd.DataFrame(plot).to_csv(O/'PLOTTED_GENES.csv',index=False,float_format='%.12g')
print(json.dumps({'png':str(O/'CONTEXT_CANDIDATE.png'),'pdf':str(O/'CONTEXT_CANDIDATE.pdf'),'size_mm':[85,119]}))
