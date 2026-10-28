from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
f=Path(__file__).resolve().parents[3]
a=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'Arial','font.size':11,'pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(6.5,6.3));fig.subplots_adjust(.015,.015,.985,.985);ax.set(xlim=(0,6.5),ylim=(0,6.3));ax.axis('off')
def box(x,y,w,h,t,fc='#f4f4f4'):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.04',fc=fc,ec='#333333',lw=1));ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=11,linespacing=1.3)
def arrow(x,y,xx,yy,both=False):ax.add_patch(FancyArrowPatch((x,y),(xx,yy),arrowstyle='<->' if both else '-|>',mutation_scale=12,lw=1.1,color='#333333'))
box(.15,5.52,2.8,.6,'SOLIDWORKS assembly\ngeometry, joints, and inertias')
box(3.55,5.52,2.8,.6,'OpenSim human / robot models\ngeometry and muscle paths')
box(.15,4.26,2.8,.82,'CAD export and model assembly\nmeshes, joint definitions,\nand anatomical frames')
box(3.55,4.26,2.8,.82,'MyoConverter and model repairs\npaths, joint directions,\nforce scales, and contacts')
arrow(1.55,5.48,1.55,5.12);arrow(4.95,5.48,4.95,5.12)
box(.15,2.93,2.8,.72,'Simscape Multibody model\nSNS-Simscape actuator blocks')
box(3.55,2.93,2.8,.72,'MuJoCo model (MJCF)\nbiological-muscle or BPA forces')
arrow(1.1,4.22,1.1,3.69);arrow(4.95,4.22,4.95,3.69)
# Separate CAD reconstruction route to MuJoCo, routed through the central gap.
ax.plot([2.5,2.5,3.25,3.25],[4.22,3.98,3.98,3.3],color='#333333',lw=1.1);arrow(3.25,3.3,3.5,3.3)
box(.15,1.0,6.2,.78,'SNS-Simulink block editor\nSNS-Simscape neuron, synapse, afferent, and activation blocks',fc='#eaf2f8')
arrow(1.55,1.82,1.55,2.89,True);ax.text(.2,2.34,'commands /\nsensory return',fontsize=10,va='center')
arrow(4.95,1.82,4.95,2.89,True);ax.text(3.38,2.34,'MuJoCo–Simulink\nbridge',fontsize=10,va='center')
ax.text(3.25,.55,'Alternative scripted controller: SNS-Toolbox (Python) ↔ MuJoCo',ha='center',fontsize=10)
ax.text(3.25,.12,'Component tests → coupled-model checks → controlled experiments',ha='center',fontsize=11,weight='bold')
fig.savefig(f/'60-sns/simulation_workflow.pdf');fig.savefig(a/'simulation_workflow.png',dpi=150);plt.close(fig)
