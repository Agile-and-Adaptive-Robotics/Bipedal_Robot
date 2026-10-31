from pathlib import Path
import csv,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Polygon,FancyArrowPatch
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'AGENTS.md').exists())
rows=list(csv.DictReader((HERE/'s3k_human92_edges.csv').open()))
EDGES={(r['source'],r['target']):float(r['reversal']) for r in rows}
plt.rcParams.update({'font.family':'Arial','font.size':12.5,'pdf.fonttype':42})
fig,axs=plt.subplots(2,1,figsize=(7.2,8.8),gridspec_kw={'height_ratios':[1,1.08]})
fig.subplots_adjust(left=.035,right=.965,top=.97,bottom=.06,hspace=.19)
BLUE='#4477AA';ORANGE='#CC6677';INK='#252c32';shown=[]
for ax in axs:ax.set(xlim=(0,10),ylim=(0,6));ax.set_aspect('equal');ax.axis('off')
def node(ax,xy,label,kind='neuron'):
 x,y=xy;r=.45
 ax.add_patch(Circle((x,y),r,fc='white' if kind=='neuron' else '#e3eee4',ec=INK,lw=1.5,zorder=4))
 ax.text(x,y+.09,label,ha='center',va='center',fontsize=11.3,fontweight='bold',zorder=5)
 xx=np.linspace(-.22,.22,24);yy=.04*np.sin(np.linspace(0,2*np.pi,24))
 ax.plot(x+xx,y-.18+yy,color=INK,lw=1.1,zorder=5)
def edge(ax,positions,a,b,rad=0,ls='-'):
 assert (a,b) in EDGES,(a,b)
 p=np.array(positions[a]);q=np.array(positions[b]);v=q-p;u=v/np.linalg.norm(v);start=p+.48*u;end=q-.49*u
 excit=EDGES[(a,b)]>=0
 ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-' if not excit else '-|>',connectionstyle=f'arc3,rad={rad}',mutation_scale=14,lw=1.5,color=INK,linestyle=ls,zorder=2,fill=False))
 if not excit:ax.add_patch(Circle(end,.065,fc=INK,ec=INK,zorder=3))
 shown.append({'source':a,'target':b,'reversal_mV':EDGES[(a,b)]})
# Selected, verified source-network feedforward connections.
a=axs[0];a.text(.1,5.9,'(A) Rhythm and pattern inputs to right-knee motor pools',fontsize=13.2,fontweight='bold',va='top')
pos={'RG_E_r':(1.2,4.7),'RG_F_r':(8.8,4.7),'InE_r':(4.1,5.05),'InF_r':(5.9,4.05),'PF_E1_r':(1.2,2.7),'PF_E2_r':(3.8,2.7),'PF_F1_r':(8.8,2.7),'MN_vas_lat_r':(2.5,.9),'MN_bifemsh_r':(8.3,.9)}
labels={'RG_E_r':'RG E','RG_F_r':'RG F','InE_r':'IN E','InF_r':'IN F','PF_E1_r':'PF E1','PF_E2_r':'PF E2','PF_F1_r':'PF F1','MN_vas_lat_r':'MN','MN_bifemsh_r':'MN'}
for name,xy in pos.items():node(a,xy,labels[name])
for x,y,rad in [('RG_E_r','InE_r',-.03),('InE_r','RG_F_r',-.06),('RG_F_r','InF_r',0),('InF_r','RG_E_r',.03),('RG_E_r','PF_E1_r',0),('RG_E_r','PF_E2_r',0),('RG_F_r','PF_F1_r',0),('PF_E1_r','MN_vas_lat_r',0),('PF_E2_r','MN_vas_lat_r',0),('PF_F1_r','MN_bifemsh_r',0)]:edge(a,pos,x,y,rad)
a.text(2.5,.18,'Vastus lateralis',ha='center',fontsize=11.5);a.text(8.3,.18,'Biceps femoris short head',ha='center',fontsize=11.5)
# Selected source-network Ia/II excitation and Ib inhibition on those same MNs.
b=axs[1];b.text(.1,5.9,'(B) Proprioceptive inputs to the same motor pools',fontsize=13.2,fontweight='bold',va='top')
pos={'Ia_vas_lat_r':(.8,4.7),'II_vas_lat_r':(.8,3.1),'Ib_vas_lat_r':(.8,1.5),'Ia_bifemsh_r':(9.2,4.7),'II_bifemsh_r':(9.2,3.1),'Ib_bifemsh_r':(9.2,1.5),'MN_vas_lat_r':(3.1,3.1),'MN_bifemsh_r':(6.9,3.1)}
for name,xy in pos.items():node(b,xy,name.split('_')[0],kind='neuron' if name.startswith('MN') else 'afferent')
for side,m in [('vas_lat_r','MN_vas_lat_r'),('bifemsh_r','MN_bifemsh_r')]:
 for s in ['Ia','II','Ib']:edge(b,pos,s+'_'+side,m)
# Reciprocal Ia inhibition is separately shown, not inferred from anatomy.
edge(b,pos,'Ia_vas_lat_r','MN_bifemsh_r',-.25)
edge(b,pos,'Ia_bifemsh_r','MN_vas_lat_r',.25)
def bio(ax,x,y,label,col):
 xx=np.linspace(-.58,.58,60);rr=.21*(1-(xx/.58)**2)**.72
 ax.add_patch(Polygon(np.c_[np.r_[x+xx,x+xx[::-1]],np.r_[y+rr,y-rr[::-1]]],fc=col,ec=INK,lw=1,zorder=4))
 for q in np.linspace(-.40,.25,5):ax.plot([x+q,x+q+.2],[y-.11,y+.11],color=INK,lw=.65,zorder=5)
 ax.plot([x-.88,x-.58],[y,y],color=INK,lw=2);ax.plot([x+.58,x+.88],[y,y],color=INK,lw=2)
 ax.text(x,y-.52,label,ha='center',fontsize=11.5)
for x,label,col in [(3.1,'Knee extensor',ORANGE),(6.9,'Knee flexor',BLUE)]:
 b.add_patch(FancyArrowPatch((x,2.65),(x,1.2),arrowstyle='-|>',mutation_scale=13,lw=1.5,color=INK))
 bio(b,x,.92,label,col)
fig.text(.08,.025,'Open arrowhead: excitation     Filled terminal: inhibition     Selected connections only',fontsize=11.3)
out=ROOT/'output/pdf/Gait2392_Neural_Pathways_Preview_20261030.pdf'
fig.savefig(out);fig.savefig(HERE/'preview.png',dpi=170)
(HERE/'PREVIEW_MANIFEST.json').write_text(json.dumps({'source':'reports_20261028_family_audit/s3k_human92_edges.csv','source_sha256':hashlib.sha256((HERE/'s3k_human92_edges.csv').read_bytes()).hexdigest(),'controller':'s3k, human92,410neurons,1186edges','glyph_source':'Code/Matlab/SNS_Simscape/sns_build_library.m and sns_icon_geom.m','selected_edges':shown,'omissions':'Other drive,interleg,adaptation,load-sharing,andmuscleinputs omitted; this is not a full connectome.','active_include':False},indent=2))
print(out)
