"""Plot actual neural/force traces; continued neural bursts are not foot steps."""
from pathlib import Path
import numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent;root=H/'runs/fourcell_hunt_20261029_v1/ground_isolation';O=H/'network_figures'
n=np.load(root/'full_feedback.neurons.npz');z=np.load(root/'full_feedback.npz');f=np.load(root/'full_feedback.assistance.npz')['samples']
plt.rcParams.update({'font.family':'Arial','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axes=plt.subplots(2,2,figsize=(12,7.2),sharex=True);colors=['#0077BB','#EE7733','#0077BB','#EE7733']
names=n['cells'].tolist();keep=(n['t']>=5)&(n['t']<13)
for ax,side,title in zip(axes[0],('r','l'),('A  Right paired PF','B  Left paired PF')):
    cells=['PF_BI_JOINT_RECT_'+side,'PF_BI_JOINT_HAM_'+side,'PF_IN_BI_PAIR_RECT_'+side,'PF_IN_BI_PAIR_HAM_'+side]
    for j,(cell,label) in enumerate(zip(cells,('RF PF','HAM PF','RF inhibitory relay','HAM inhibitory relay'))):ax.plot(n['t'][keep],n['voltage'][keep,names.index(cell)],color=colors[j],ls='-' if j<2 else '--',lw=1.5,label=label)
    ax.set(title=title,ylabel='Model voltage (mV)');ax.legend(fontsize=10,ncol=2,loc='upper right');ax.grid(alpha=.18)
keep=(z['t']>=5)&(z['t']<13);names=z['neuro_names'].tolist();ax=axes[1,0]
for side,color,style in [('r','#0077BB','-'),('l','#CC3311','--')]:
    contrast=z['neuro'][:,names.index('RG_E_'+side)]-z['neuro'][:,names.index('RG_F_'+side)]
    ax.plot(z['t'][keep],contrast[keep],color=color,ls=style,lw=1.2,label='Right' if side=='r' else 'Left')
ax.set(title='C  Rhythm generators continue bursting',ylabel='RG E − F voltage (mV)',xlabel='Time (s)');ax.legend(fontsize=10);ax.grid(alpha=.18)
keep=(f[:,0]>=5)&(f[:,0]<13);bw=78.4646*9.81;ax=axes[1,1]
for col,color,style,label in [(7,'#0077BB','-','Right foot'),(8,'#CC3311','--','Left foot'),(1,'#777777','-.','Vertical harness')]:ax.plot(f[keep,0],f[keep,col]/bw,color=color,ls=style,lw=1.4,label=label)
ax.axhline(.05,color='#333333',lw=.7,ls=':',label='Contact-on threshold');ax.set(title='D  Both feet remain loaded',ylabel='Force / body weight',xlabel='Time (s)',ylim=(-.02,.9));ax.legend(fontsize=10,ncol=2);ax.grid(alpha=.18)
fig.suptitle('Four-cell RF/hamstring PF: grounded diagnostic, not successful walking',fontsize=14)
fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(O/'fourcell_ground_diagnosis.png',dpi=180);fig.savefig(O/'fourcell_ground_diagnosis.pdf')
(O/'fourcell_ground_diagnosis.alt.txt').write_text('Human Gait2392 single-cadence four-cell controller, grounded full-feedback diagnostic. Right RF and left hamstring PF activity are biased despite continuing bilateral RG bursts. Once-counted GRFs remain positive on both feet, so these neural cycles are not alternating ground steps. Source: runs/fourcell_hunt_20261029_v1/ground_isolation/full_feedback files. The harness remains active.')
print(O/'fourcell_ground_diagnosis.png')
