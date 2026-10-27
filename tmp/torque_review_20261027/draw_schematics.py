from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Arc, FancyArrowPatch, PathPatch
from matplotlib.path import Path as MPath

ROOT=Path(r'D:\GitHub\Bipedal_Robot')
OUT=ROOT/'tmp/torque_review_20261027'
DEST=[ROOT/'Documentation/Reports and Papers/Dissertation/Figures/30-torque',
      ROOT/'Documentation/Reports and Papers/Knee_Torque_Test/Figures']
plt.rcParams.update({'font.family':'Arial','font.size':11,'mathtext.fontset':'dejavusans',
                     'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.8})
INK='#20252b'; BLUE='#4477AA'; MAGENTA='#AA3377'; CYAN='#66CCEE'; GRAY='#e3e6e8'; GOLD='#997700'

def arrow(ax,a,b,color=INK,lw=1.4,style='-|>',ms=10):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=ms,lw=lw,color=color,zorder=8))

def save(fig,name,alt):
    for ext in ['pdf','svg','png']:
        f=OUT/f'{name}.{ext}'
        fig.savefig(f,dpi=240,facecolor='white',bbox_inches='tight',pad_inches=.08)
        for dest in DEST: shutil.copy2(f,dest/f.name)
    for dest in DEST: (dest/f'{name}.alt.txt').write_text(alt,encoding='utf-8')
    plt.close(fig)

# Projection: a compact force-path view plus a true geometric projection inset.
fig=plt.figure(figsize=(7.2,3.9))
ax=fig.add_axes([.01,.45,.98,.53]);ax.set(xlim=(-.4,7.6),ylim=(-.65,1.8));ax.axis('off')
for x,side in [(0,1),(7,-1)]:
    ax.plot([x,x],[-.35,1.2],color=INK,lw=2)
    for y in np.arange(-.3,1.2,.22):ax.plot([x,x-.18*side],[y,y-.14],color=INK,lw=.7)
    ax.plot([x,x+.65*side],[.4,.4],color=INK,lw=4)
ax.plot([.65,6.35],[.4,.4],color='#7b8085',ls=(0,(4,3)),lw=1)
ax.add_patch(Polygon([[2,.21],[2.2,.08],[4.8,.08],[5,.21],[5,.59],[4.8,.72],[2.2,.72],[2,.59]],closed=True,fc='#eef4fa',ec=BLUE,lw=1.5))
for x in np.arange(2.1,4.8,.25):
    ax.plot([x,x+.2],[.12,.67],color=BLUE,lw=.65)
    ax.plot([x,x+.2],[.67,.12],color=BLUE,lw=.65)
ax.plot([.65,2],[.4,.4],color=INK,lw=1.7)
ax.plot([5,6.35],[.4,.4],color=INK,lw=1.7)
ax.text(3.5,.95,'BPA',ha='center',color=BLUE)
ax.text(5.7,.03,'Tendon',ha='center')
arrow(ax,(.65,.9),(1.45,.9));ax.text(1.02,1.17,r'$F\hat{\mathbf{u}}_1$',ha='center')
arrow(ax,(6.35,.9),(5.55,.9));ax.text(5.96,1.17,r'$F\hat{\mathbf{u}}_2$',ha='center')
ax.text(0,-.51,'Bracket 1',ha='center',fontsize=10);ax.text(7,-.51,'Bracket 2',ha='center',fontsize=10)
for x0,x1,lab in [(0.65,1.25,r'$\delta_{br,1}$'),(5.35,6.15,r'$\delta_{ten}$'),(6.35,5.75,r'$\delta_{br,2}$')]:
    y=-.15 if 'ten' not in lab else -.62
    arrow(ax,(x0,y),(x1,y),color=MAGENTA,style='<->',ms=8)
    ax.text((x0+x1)/2,y-.07,lab,color=MAGENTA,ha='center',va='top',fontsize=10)

ax=fig.add_axes([.04,.02,.42,.4]);ax.set_aspect('equal');ax.set(xlim=(-.3,2.4),ylim=(-.3,1.65));ax.axis('off')
o=np.array([0.,0.]);eta=np.array([1.5,1.1]);proj=np.array([1.5,0.])
arrow(ax,o,(2.1,0));ax.text(2.13,0,r'$\hat{\mathbf{u}}$',va='center')
arrow(ax,o,eta,MAGENTA);ax.text(.38,.85,r'$\boldsymbol{\eta}_{br}$',color=MAGENTA)
ax.plot([eta[0],proj[0]],[eta[1],proj[1]],color='#777777',ls='--',lw=1)
ax.plot([1.34,1.34,1.5],[0,.16,.16],color='#777777',lw=.8)
arrow(ax,(0,-.15),(1.5,-.15),BLUE,style='<->',ms=8)
ax.text(.75,-.22,r'$\delta_{br}$',color=BLUE,ha='center',va='top')
ax.text(0,1.48,'Local projection',fontsize=10)
fig.text(.49,.25,r'$\delta_{br,j}=\hat{\mathbf{u}}_j^{T}\boldsymbol{\eta}_{br,j}$',fontsize=12)
fig.text(.49,.10,r'$\delta=\delta_{br,1}+\delta_{ten}+\delta_{br,2}$',fontsize=11)
save(fig,'xiProjection_review','Two compliant brackets connected by a BPA and tendon. Opposed local tension directions act toward the actuator. An inset projects an off-axis bracket displacement onto its force direction. The two projected displacements and tendon elongation sum to the scalar actuator shortening. All displacements are schematic and exaggerated.')

# Sagittal knee: a continuous, braided actuator passes down the anterior
# femur, around a circular contact, then onto a distal artificial tendon.
fig,ax=plt.subplots(figsize=(7.2,4.8))
ax.set_aspect('equal');ax.axis('off');ax.set(xlim=(-2.85,4.55),ylim=(-2.6,3.15))
femur=[(-.46,2.9),(.35,2.9),(.3,1.2),(.55,.55),(.7,.2),(.56,-.48),(-.1,-.7),(-.65,-.45),(-.8,.1),(-.45,.7)]
ax.add_patch(Polygon(femur,closed=True,fc=GRAY,ec='#828b93',lw=1.1))
v=np.array([-.9063,-.4226]);n=np.array([.4226,-.9063])
tibia0=np.array([-.1,-.43]);tibia1=np.array([-2.4,-1.9])
ax.add_patch(Polygon([tibia0+np.array([-.43,.16]),tibia0+np.array([.25,-.32]),tibia1+np.array([.22,-.18]),tibia1+np.array([-.18,.2])],closed=True,fc=GRAY,ec='#828b93',lw=1.1))
ax.add_patch(Circle((0,0),.87,fc='#edf0f2',ec='#828b93',lw=1.1,zorder=2))
ax.add_patch(Circle((0,0),.09,fc='white',ec=INK,lw=1,zorder=3))
ax.text(-1.0,2.2,'Femur',ha='right')
ax.text(-2.3,-1.12,'Tibia',ha='center')
ax.annotate('Knee axis',xy=(0,0),xytext=(-1.3,.62),arrowprops=dict(arrowstyle='-',lw=.7),ha='right',fontsize=10)

R=1.08;theta=np.linspace(0,np.deg2rad(-65),160)
top=np.column_stack([np.full(80,R),np.linspace(2.62,0,80)])
arc=np.column_stack([R*np.cos(theta),R*np.sin(theta)])
end=arc[-1];tail=np.array([end+t*v for t in np.linspace(0,1.2,65)])
center=np.vstack([top,arc[1:],tail[1:]])
dx=np.gradient(center,axis=0);normals=np.column_stack([-dx[:,1],dx[:,0]])
normals/=np.linalg.norm(normals,axis=1)[:,None]
width=.125;poly=np.vstack([center+width*normals,(center-width*normals)[::-1]])
body=Polygon(poly,closed=True,fc='#d8e8f5',ec=BLUE,lw=1.1,zorder=5);ax.add_patch(body)
# Alternating diagonal braid traces clipped to the actuator envelope.
dist=np.r_[0,np.cumsum(np.linalg.norm(np.diff(center,axis=0),axis=1))]
for ss in np.arange(0,dist[-1]-.15,.18):
    ix=np.searchsorted(dist,ss);jx=np.searchsorted(dist,ss+.15)
    for sign in [-1,1]:
        q1=center[ix]+sign*width*normals[ix];q2=center[jx]-sign*width*normals[jx]
        line,=ax.plot([q1[0],q2[0]],[q1[1],q2[1]],color=BLUE,lw=.65,zorder=6)
        line.set_clip_path(body)

# Fittings and anchor connections.
for pt,tangent in [(center[0],np.array([0.,-1.])),(center[-1],v)]:
    normal=np.array([-tangent[1],tangent[0]])
    coords=[pt+u*tangent+w*normal for u,w in [(-.12,-.16),(.12,-.16),(.12,.16),(-.12,.16)]]
    ax.add_patch(Polygon(coords,fc='#b9c2c9',ec=INK,lw=.8,zorder=7))
anchor=center[-1]+.66*v
ax.plot([center[-1,0],anchor[0]],[center[-1,1],anchor[1]],color=INK,lw=1.6,zorder=4)
ax.plot(anchor[0],anchor[1],marker='D',ms=5,color=INK,zorder=8)
ax.plot([.12,R],[2.9,2.9],color=INK,lw=2)
ax.plot([R,R],[2.9,2.65],color=INK,lw=1.5)
ax.text(1.48,2.82,'Proximal attachment',fontsize=10)
ax.annotate('Braided BPA',xy=(R,1.7),xytext=(2.02,1.7),arrowprops=dict(arrowstyle='-',lw=.7),va='center')
ax.annotate('Distal tendon',xy=tuple(center[-1]+.3*v),xytext=(.55,-2.2),arrowprops=dict(arrowstyle='-',lw=.7),fontsize=10)

for th in [0,-65]:
    r=np.deg2rad(th);q=R*np.array([np.cos(r),np.sin(r)])
    ax.plot([0,q[0]],[0,q[1]],color=INK,lw=.9,ls=(0,(4,3)),zorder=8)
arrow(ax,(.08,.08),(R-.08,.08),GOLD,style='<->',ms=8)
ax.text(.58,.19,r'$R$',ha='center',color=GOLD,fontsize=12)
ax.add_patch(Arc((0,0),.88,.88,theta1=-65,theta2=0,color=INK,lw=1.25,zorder=9))
ax.text(.49,-.32,r'$\alpha$',fontsize=14,zorder=9)
ax.text(2.0,.25,'Contact geometry',fontsize=11,weight='bold')
ax.text(2.0,-.15,r'$\alpha=|\theta_{\mathrm{wrap}}-\theta_k|$',fontsize=12)
ax.text(2.0,-.57,r'$\Delta l=\chi_3 R\alpha(1-\varepsilon^*)^2$',fontsize=12)
ax.text(2.0,-1.0,'Fitted usable-length loss',fontsize=10)
ax.text(-2.6,-2.45,'Sagittal view  |  schematic, not to scale',fontsize=10,color='#555555')
fig.subplots_adjust(left=.01,right=.99,top=.99,bottom=.01)
save(fig,'bpa_knee_wrapping_review','Sagittal robot knee schematic with a proximal femur, flexed tibia, and a braided extensor BPA running down the anterior femur and around the knee into a distal tendon. Two dashed contact radii enclose the full bend angle alpha. R labels the radius of the actuator path. The fitted usable-length loss is chi3 times R times alpha times the squared complement of normalized contraction. The figure uses a single circular contact approximation and is not to scale.')
print('Wrote two vector schematics and PNG previews to both figure trees.')
