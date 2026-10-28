from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle,Arc,FancyArrowPatch
R=Path('D:/Github/Bipedal_Robot');D=R/'Documentation/Reports and Papers/Dissertation';A=D/'Archive/appendix_c_corrections_20261028'
plt.rcParams.update({'font.family':'Arial','font.size':13,'pdf.fonttype':42})
fig=plt.figure(figsize=(7.2,5.3));ax=fig.add_axes([.02,.13,.61,.80]);ax.set(xlim=(0,10),ylim=(.5,8));ax.axis('off');ax.set_aspect('equal',adjustable='box');blue='#0072B2';red='#CC6677';gray='.35'
ax.text(.05,7.65,'(A) Reverse tackle: two line parts',fontweight='bold',fontsize=14)
ax.plot([.6,.6],[3.0,6.9],color=gray,lw=4)
for yy in [3.6,6.3]:
 ax.add_patch(Rectangle((.8,yy-.18),2.1,.36,fc='#f7e5e9',ec=red,lw=1.5));ax.plot([.6,.8],[yy,yy],color=gray,lw=2);ax.plot([2.9,3.2],[yy,yy],color=gray,lw=2)
ax.plot([3.2,3.2],[3.6,6.3],color=gray,lw=2);ax.plot([3.2,4],[5,5],color=gray,lw=2)
ax.add_patch(Rectangle((6.6,3.1),1.4,3.0,ec=gray,fc='#eeeeee',lw=1.5,zorder=0))
ax.add_patch(Circle((4,5),.6,ec=blue,fc='white',lw=1.5));ax.plot([4,7],[5.6,5.6],color=blue,lw=2);ax.add_patch(Arc((4,5),1.2,1.2,theta1=90,theta2=270,color=blue,lw=2));ax.plot([4,7],[4.4,4.4],color=blue,lw=2)
ax.add_patch(Circle((7,3.8),.6,ec=gray,fc='white',lw=1.5));ax.add_patch(Arc((7,3.8),1.2,1.2,theta1=0,theta2=90,color=blue,lw=2));ax.plot([7.6,7.6],[3.8,1.2],color=blue,lw=2);ax.plot(7,5.6,'o',color=gray,ms=6);ax.plot(7.6,1.2,'D',color=blue,ms=7)
def arrow(a,b,c=gray):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=14,color=c,lw=1.5))
arrow((4,2.7),(3.0,2.7),red);ax.text(3.5,2.15,r'BPA contraction $\Delta\ell_b$',ha='center',fontsize=12)
arrow((8.8,1.2),(8.8,3.4),blue);ax.text(9.15,2.25,r'$\Delta s_t=2\Delta\ell_b$',rotation=90,va='center',fontsize=12)
ax.text(1.85,6.85,r'$n$ parallel BPAs',ha='center',color=red);ax.text(.45,2.05,'BPA origin\nmount',ha='center',fontsize=12);ax.text(4.3,6.6,'Moving\nsheave',ha='center',fontsize=12,va='bottom');ax.text(7.3,6.55,'Fixed tackle\nmount',ha='center',fontsize=12);ax.text(7.45,.68,'Tendon insertion',ha='center',fontsize=12)
arrow((4.1,5.0),(3.4,5.0),red);ax.text(3.1,4.6,r'$nF_b$',ha='right',color=red,fontsize=12)
for yy in [5.6,4.4]:arrow((4.8,yy),(5.7,yy),blue)
ax.text(5.3,5.83,r'$F_t$',color=blue);ax.text(5.3,4.03,r'$F_t$',color=blue)
ax2=fig.add_axes([.67,.12,.31,.81]);ax2.axis('off');ax2.text(0,.98,'(B) Ideal relations',weight='bold',fontsize=14,va='top')
ax2.text(0,.84,'Travel multiplier',fontsize=12)
ax2.text(0,.76,r'$G=\Delta s_t/\Delta\ell_b$',fontsize=15)
ax2.text(0,.63,'Force multiplier',fontsize=12)
ax2.text(0,.55,r'$F_t/F_b=n/G$',fontsize=15)
ax2.text(0,.42,'Fixed tackle cable load',fontsize=12)
ax2.text(0,.34,r'$\mathbf{R}_c=nF_b\hat{\mathbf{u}}_b+F_t\hat{\mathbf{u}}_t$',fontsize=13)
ax2.text(0,.24,r'$\hat{\mathbf{u}}_b$: left; $\hat{\mathbf{u}}_t$: down',fontsize=12)
ax2.text(0,.12,r'Support reaction: $-\mathbf{R}_c$',fontsize=12)
ax2.text(0,.02,'Separate BPA origin load',fontsize=12)
fig.text(.03,.035,'Ideal cable and rigid blocks shown. Compliance and pulley losses are omitted.',fontsize=12)
fig.savefig(A/'reverse_pulley_principle.pdf',bbox_inches='tight');fig.savefig(A/'reverse_pulley_principle.png',dpi=160,bbox_inches='tight')
