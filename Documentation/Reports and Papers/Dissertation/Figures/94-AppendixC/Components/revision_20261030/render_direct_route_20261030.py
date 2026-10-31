"""Render unchanged saved direct-route arrays; zoom to contact region, without bones."""
from pathlib import Path
import re,json,hashlib
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[8]
if R.name!='Bipedal_Robot':
 R=next(p for p in Path(__file__).resolve().parents if p.name=='Bipedal_Robot')
F=R/'Documentation/Reports and Papers/Dissertation/Figures/94-AppendixC'
c=re.findall(r"c\{\d+\}\s*=\s*'(#[A-Fa-f0-9]{6})'",(R/'Code/Matlab/Colors.m').read_text())
source=R/'Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_20260925.mat'
S=loadmat(source,simplify_cells=True,variable_names=['ctx','predBest']);ctx=S['ctx'];ri=S['predBest']['routeInfo'];geo=ctx['geo'];a=np.asarray(ri['active'],bool)
idx=np.unique(np.r_[0,np.flatnonzero(np.any(a[:,:-1]&~a[:,1:],axis=0)),a.shape[1]-1]);assert len(idx)==8
plt.rcParams.update({'font.family':'Arial','font.size':11.5,'axes.labelsize':11.5,'axes.titlesize':11.5,'xtick.labelsize':11.5,'ytick.labelsize':11.5,'pdf.fonttype':42})
fig,axs=plt.subplots(2,4,figsize=(8.4,4.2),sharex=True,sharey=True);th=np.linspace(0,2*np.pi,200)
def trans(P,T):return np.asarray(P)@T[:3,:3].T+T[:3,3]
def to_fem(P,i):return trans(trans(P,ctx['T_ICR_t1'][:,:,i]),ctx['T_Pam'][:,:,i])
for k,(ax,i) in enumerate(zip(axs.flat,idx)):
 P=ri['raw'][:,:,i].copy();P[5:]=to_fem(P[5:],i)
 C=np.asarray(geo['femurCylCenter']);rad=geo['femurCylClearRadius'];ax.plot(C[0]+rad*np.cos(th),C[1]+rad*np.sin(th),color=c[0],lw=1.5)
 ax.plot([geo['femurLineX']]*2,geo['femurLineY'],color=c[1],lw=1.5)
 Q=np.vstack((geo['femurOffsetBoundary'],geo['femurOffsetBoundary'][0]));ax.plot(Q[:,0],Q[:,1],color=c[2],lw=1.5)
 if 'femurCondyleClipY' in geo:ax.plot([geo['femurCondyleClipX']]*2,geo['femurCondyleClipY'],color=c[2],lw=1.5)
 for j,key in enumerate(['Lower','Upper']):
  C=np.asarray(geo['tibia'+key+'Center']);rad=geo['tibia'+key+'ClearRadius'];pts=to_fem(np.c_[C[0]+rad*np.cos(th),C[1]+rad*np.sin(th),np.zeros(len(th))],i);ax.plot(pts[:,0],pts[:,1],color=c[3+j],lw=1.5)
 if a[1,i]:ax.plot([geo['femurCylCenter'][0],P[1,0]],[geo['femurCylCenter'][1],P[1,1]],'--',color=c[5],lw=1.8)
 if a[7,i]:
  tc=to_fem(np.r_[geo['tibiaLowerCenter'],0][None,:],i)[0];ax.plot([tc[0],P[7,0]],[tc[1],P[7,1]],'--',color=c[5],lw=1.8)
 ax.plot(P[:,0],P[:,1],'-o',color=c[4],lw=1.8,ms=3,mec='none');ax.scatter(P[-1,0],P[-1,1],s=35,c=c[1],marker='D',edgecolors='none',zorder=5)
 nxt=np.flatnonzero(a[:,i]&~a[:,i+1]) if i<a.shape[1]-1 else []
 event=', '.join(f'p{j+1}' for j in nxt)
 title=f'({chr(65+k)}) '+(f'Before {event} release\n' if event else '\n')+f'{ctx["phiD"][i]:.1f} deg'
 ax.set_title(title,pad=8,fontsize=11.5)
 for j in nxt:
  ax.annotate(f'p{j+1}',P[j,:2],xytext=(16,-15),textcoords='offset points',fontsize=11.5,arrowprops={'arrowstyle':'-','lw':.6,'color':'.35'})
 if k==0:ax.annotate('p2',P[1,:2],xytext=(9,5),textcoords='offset points',fontsize=11.5)
 ax.set_xlim(-.24,.18);ax.set_ylim(-.67,-.20);ax.set_aspect('equal',adjustable='box');ax.spines[['top','right']].set_visible(False);ax.set_xticks([-.2,0,.15]);ax.set_yticks([-.6,-.4,-.2]);ax.tick_params(length=3)
 if k%4==0:ax.set_ylabel('Femur-frame y (m)')
 if k>=4:ax.set_xlabel('Femur-frame x (m)')
handles=[Line2D([],[],color=c[j],lw=1.6,label=lab) for j,lab in enumerate(['Femur cylinder','Femur wall','Condyle envelope','Tibia lower','Tibia upper','Bend radius','']) if j<6]
handles[5].set_linestyle('--');handles.append(Line2D([],[],color=c[4],lw=1.8,marker='o',ms=3,label='BPA/tendon route'));handles.append(Line2D([],[],color=c[1],marker='D',ls='none',label='Insertion p9'))
fig.legend(handles=handles,loc='lower center',ncol=4,fontsize=11.5,frameon=False,bbox_to_anchor=(.5,.005),handlelength=2.0,columnspacing=1.1)
fig.subplots_adjust(left=.065,right=.99,top=.87,bottom=.24,hspace=.72,wspace=.35)
fig.savefig(F/'direct_route_contact_detail_20261030.pdf');fig.savefig(F/'direct_route_contact_detail_20261030.png',dpi=220);plt.close(fig)
(F/'direct_route_contact_detail_20261030.alt.txt').write_text('Eight unchanged saved direct-route poses, cropped to the knee contact region with no bone overlay. Panels show full flexion, each pre-release frame, and full extension. All saved route vertices and all clearance curves are retained; the proximal origin is outside this close-up. A preceding anatomical-context figure shows the full path. A magenta polyline shows the BPA/tendon route; purple dashes are local bend radii, and panel titles identify the next released contact. Labels are limited to that contact and the p2 landmark to avoid collisions.\n')
print(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'indices_matlab':(idx+1).tolist(),'angles':ctx['phiD'][idx].tolist(),'saved_route_values_changed':False}))
