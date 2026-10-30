from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle,Polygon,FancyBboxPatch,Arc,FancyArrowPatch
plt.rcParams.update({'font.family':'Arial','font.size':11,'pdf.fonttype':42})
f,ax=plt.subplots(figsize=(7.2,5.8)); f.subplots_adjust(left=.025,right=.99,bottom=.04,top=.98);ax.set(xlim=(0,10),ylim=(0,8));ax.set_aspect('equal');ax.axis('off')
red='#CC6677';blue='#0072B2';grey='#555555'
# Stylized proximal and distal bones; linkage indicates a migrating-center knee, not a simple hinge.
ax.add_patch(Polygon([[.6,3.8],[.6,4.5],[5.9,4.6],[6.7,4.4],[7.1,3.8],[6.7,3.1],[6.0,3.5]],closed=True,fc='#f0ece4',ec='#7f7b73',lw=1.6))
ax.add_patch(Polygon([[6.9,3.2],[7.45,3.4],[8.75,.65],[8.15,.4],[6.7,2.9]],closed=True,fc='#f0ece4',ec='#7f7b73',lw=1.6))
for x,y in [(6.25,3.85),(6.65,4.05),(6.75,3.2),(7.05,3.35)]:ax.plot(x,y,'o',mfc='white',mec=grey,ms=5,zorder=6)
ax.plot([6.25,7.05],[3.85,3.35],color=grey,lw=3);ax.plot([6.65,6.75],[4.05,3.2],color=grey,lw=3)
ax.text(2.1,3.93,'Femur',ha='center');ax.text(8.85,1.5,'Tibia',ha='center');
# Two recognizable braided actuators in parallel, attached to a common moving crosshead.
ax.plot([.85,.85],[4.5,6.9],color=grey,lw=4)
for y in [5.55,6.45]:
 ax.plot([.85,1.15],[y,y],color=grey,lw=2)
 ax.add_patch(FancyBboxPatch((1.15,y-.23),2.65,.46,boxstyle='round,pad=0,rounding_size=.22',fc='#f5dfe5',ec=red,lw=1.5))
 for x in np.arange(1.35,3.5,.23):
  ax.plot([x,x+.25],[y-.16,y+.16],color=red,lw=.6,alpha=.7);ax.plot([x,x+.25],[y+.16,y-.16],color=red,lw=.6,alpha=.7)
 ax.plot([3.8,4.08],[y,y],color=grey,lw=2)
ax.plot([4.08,4.08],[5.55,6.45],color=grey,lw=3);ax.plot([4.08,4.8],[6,6],color=grey,lw=2)
ax.add_patch(Circle((4.8,6),.45,ec=grey,fc='white',lw=1.5))
# Cable: anchor on proximal mount, around moving block, then fixed redirect to distal insertion.
ax.plot([4.8,7.35],[6.45,6.45],color=blue,lw=2);ax.add_patch(Arc((4.8,6),.9,.9,theta1=90,theta2=270,color=blue,lw=2));ax.plot([4.8,7.35],[5.55,5.55],color=blue,lw=2)
ax.plot([7.35,7.35],[4.45,6.65],color=grey,lw=3,zorder=0);ax.plot([7.35,6.45],[4.45,4.45],color=grey,lw=3,zorder=0)
ax.plot(7.35,6.45,'o',color=blue,ms=5)
ax.add_patch(Circle((7.35,5.1),.45,ec=grey,fc='white',lw=1.5));ax.add_patch(Arc((7.35,5.1),.9,.9,theta1=0,theta2=90,color=blue,lw=2));ax.plot([7.8,7.8],[5.1,2.35],color=blue,lw=2);ax.plot(7.8,2.35,'D',color=blue,ms=6)
ax.text(2.5,7.15,'Two parallel BPAs',color=red,ha='center',weight='bold');ax.text(.85,7.58,'Origin mount',ha='center');ax.text(.85,7.32,r'$F_{R,\mathrm{origin}}=F_{\mathrm{BPA}}$',ha='center',fontsize=10);ax.text(4.8,7.15,'Moving\nsheave',ha='center');ax.plot([4.8,4.8],[6.85,6.5],color=grey,lw=.8)
ax.text(8.0,6.6,'Fixed tackle\nmount on femur',va='center');ax.text(8.1,4.45,'Artificial\ntendon',color=blue);ax.text(8.15,2.35,'Insertion',color=blue,va='center')
def arrow(a,b,c):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=13,lw=1.5,color=c))
arrow((4.65,5),(3.65,5),red);ax.text(3.55,4.8,r'$\Delta\ell_{\mathrm{BPA}}$',color=red,ha='right',va='center')
arrow((9.25,2.6),(9.25,3.7),blue);ax.text(9.45,3.15,r'Tendon travel $\Delta s_t=n\Delta\ell_{\mathrm{BPA}}$',rotation=90,ha='center',va='center',color=blue)
arrow((7.8,3.8),(7.8,4.3),blue);ax.text(8.0,4.0,r'$F_t=F_{\mathrm{BPA}}/n$',color=blue,fontsize=10)

ax.text(2.5,6.92,r'$F_{\mathrm{BPA}}=2F_{20}$',color=red,ha='center'); ax.text(5.95,6.8,r'$n=2$ cable parts',ha='center');arrow((4.75,6),(4.2,6),red)
f.savefig(Path(__file__).resolve().parents[2]/'reverse_pulley_biomimetic.pdf');f.savefig(Path(__file__).resolve().parent/'reverse_pulley_biomimetic_preview.png',dpi=180)
