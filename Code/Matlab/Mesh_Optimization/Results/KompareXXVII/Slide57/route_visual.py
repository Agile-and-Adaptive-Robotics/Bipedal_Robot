from pathlib import Path
from scipy.io import loadmat
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
W=Path(__file__).parent
K=W.parent
with (K/'route_configurations.csv').open(encoding='utf-8-sig',newline='') as f:
    selected=[q for q in csv.DictReader(f) if q['side']=='selected']
assert len(selected)==27 and all(q['points_identical_to_baseline']=='true' for q in selected)
M=loadmat(W/'route_projection.mat',squeeze_me=True,struct_as_record=False)
routes=M['Pel']; entries=M['sel']
def group(q):
    if max(int(q['tackleLineParts_c1']),int(q['tackleLineParts_c2']))>1:return 'tackle'
    return 'bundle' if int(q['nBPA'])>1 else 'single'
colors={'single':'#157349','bundle':'#1958AA','tackle':'#BD520F'}
styles={'single':'-','bundle':'-','tackle':'--'}
def project(x):
    x=np.atleast_2d(x)
    return np.c_[.80*x[:,0]-.60*x[:,2],x[:,1]]
def line(ax,p,color,style='-',width=2.0,alpha=1):
    p=project(p)
    ax.plot(p[:,0],p[:,1],linestyle=style,color=color,lw=width,alpha=alpha)
plt.rcParams.update({'font.family':'Arial'})
fig=plt.figure(figsize=(5.8,5.5),facecolor='white')
ax=fig.add_axes([.0,.02,.73,.96])
hip=M['Thip'][:3,3]; knee=M['Ttib'][:3,3]
ank=M['Ttib'][:3,:3]@M['ankP']+knee
sub=M['Ttib'][:3,:3]@M['subP']+knee
toe=M['Ttib'][:3,:3]@M['toeP']+knee
line(ax,[[0,.36,0],[0,0,0],hip,knee,ank,sub,toe],'#d2d7da',width=8)
for i,q in enumerate(selected):
    assert entries[i].site==q['site']
    if q['site']=='psoas_r':continue
    g=group(q);p=np.asarray(routes[i].P)
    line(ax,p,colors[g],styles[g],2.1,.89)
    pp=project(p)
    ax.scatter(pp[[0,-1],0],pp[[0,-1],1],s=12,color=colors[g],zorder=3)
ax.set_aspect('equal');ax.axis('off')
inset=fig.add_axes([.72,.42,.26,.28])
q=next(q for q in selected if q['site']=='psoas_r')
i=int(q['actuator_id'])-1
line(inset,entries[i].Pw,colors[group(q)],styles[group(q)],3)
inset.set_aspect('equal');inset.axis('off')
inset.set_title('Psoas\nlocal frame',fontsize=12,color='#2b3a42',pad=2)
fig.savefig(W/'route_schematic.png',dpi=220);plt.close(fig)

