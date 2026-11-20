"""Render the sagittal wrapping schematic using the original joint and actuator geometry."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Arc, FancyArrowPatch

plt.rcParams.update({'font.family':'Arial','font.size':11,'mathtext.fontset':'custom',
                     'mathtext.it':'Arial:style=italic','mathtext.rm':'Arial',
                     'pdf.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'knee-wrap'})
INK='#20252b'; BLUE='#4477AA'; GRAY='#e3e6e8'; GOLD='#997700'
leader={'arrowstyle':'-|>', 'lw':.85, 'color':INK, 'mutation_scale':9.5,
        'shrinkA':0, 'shrinkB':0, 'joinstyle':'round'}

def arrow(ax,a,b,color=INK,lw=1.4,style='-|>',ms=10):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=ms,lw=lw,color=color,zorder=8))

def draw(folder):
    fig,ax=plt.subplots(figsize=(5.33,4.38))
    ax.set_aspect('equal');ax.axis('off');ax.set(xlim=(-2.76,3.5),ylim=(-2.25,3.03))
    femur=[(-.46,2.9),(.35,2.9),(.3,1.2),(.55,.55),(.7,.2),(.56,-.48),(-.1,-.7),(-.65,-.45),(-.8,.1),(-.45,.7)]
    ax.add_patch(Polygon(femur,closed=True,fc=GRAY,ec='#828b93',lw=1.1))
    v=np.array([-.9063,-.4226]);n=np.array([.4226,-.9063])
    tibia0=np.array([-.1,-.43]);tibia1=np.array([-2.4,-1.9])
    ax.add_patch(Polygon([tibia0+np.array([-.43,.16]),tibia0+np.array([.25,-.32]),tibia1+np.array([.22,-.18]),tibia1+np.array([-.18,.2])],closed=True,fc=GRAY,ec='#828b93',lw=1.1))
    ax.add_patch(Circle((0,0),.87,fc='#edf0f2',ec='#828b93',lw=1.1,zorder=2))
    ax.add_patch(Circle((0,0),.09,fc='white',ec=INK,lw=1,zorder=3))
    ax.text(-1.0,2.2,'Femur',ha='right')
    ax.text(-2.3,-1.12,'Tibia',ha='center')
    ax.annotate('Knee axis', xy=(-.19,.075), xytext=(-1.30,.65), ha='right', va='center',
                arrowprops={**leader, 'shrinkA':3}, fontsize=11, zorder=10)

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
    ax.annotate('Braided BPA', xy=(R+.23,1.70), xytext=(1.93,1.70), va='center',
                arrowprops={**leader, 'shrinkA':3}, fontsize=11, zorder=10)
    # Place the shaft below the tendon. Its tip ends on the normal to the
    # middle of the tendon, leaving a visible clearance in both print sizes.
    target = center[-1]+.34*v
    tip = target + .125*np.array([.4226,-.9063])
    elbow = (-.38,-2.00)
    ax.text(.24,-2.00,'Distal tendon',va='center',ha='left',fontsize=11)
    ax.plot([.18,elbow[0]],[-2.00,elbow[1]],color=INK,lw=.85,solid_capstyle='butt',zorder=10)
    ax.annotate('', xy=tip, xytext=elbow, arrowprops=leader, zorder=10)

    for th in [0,-65]:
        r=np.deg2rad(th);q=R*np.array([np.cos(r),np.sin(r)])
        ax.plot([0,q[0]],[0,q[1]],color=INK,lw=.9,ls=(0,(4,3)),zorder=8)
    arrow(ax,(.08,.08),(R-.08,.08),GOLD,style='<->',ms=8)
    ax.text(.58,.19,r'$R$',ha='center',color=GOLD,fontsize=12)
    ax.add_patch(Arc((0,0),.88,.88,theta1=-65,theta2=0,color=INK,lw=1.25,zorder=9))
    ax.text(.49,-.32,r'$\alpha$',fontsize=14,zorder=9)
    fig.subplots_adjust(left=.012,right=.988,bottom=.018,top=.98)
    for extension in ('pdf','svg','png'):
        metadata={'CreationDate':None} if extension=='pdf' else {'Date':None} if extension=='svg' else None
        fig.savefig(folder/('bpa_knee_wrapping_review.'+extension),facecolor='white',dpi=220,metadata=metadata)
    plt.close(fig)

if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('directory',type=Path)
    folder=cli.parse_args().directory
    folder.mkdir(parents=True,exist_ok=True)
    draw(folder)
