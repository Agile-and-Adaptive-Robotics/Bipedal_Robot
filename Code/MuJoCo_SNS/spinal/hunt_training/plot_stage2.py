import os,sys,json
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent;O=H/'runs/stage2_diagnosis_v1'
plt.rcParams.update({'font.family':'Arial','font.size':10})
fig,axes=plt.subplots(2,2,figsize=(11,6),sharex=True,sharey=True,layout='constrained')
for ax,(name,title) in zip(axes.flat,[('both_off','Both feedback blocks off (diagnostic control)'),('timing_only','Timing feedback only'),('reflex_only','Motor reflex feedback only'),('both_on','Both feedback blocks on')]):
    z=np.load(O/(name+'.npz'),allow_pickle=True);keep=(z['t']>=5)&(z['t']<13)
    for side,color,style in [('r','#4477AA','-'),('l','#CC6677','--')]:
        q=z['q'][:,list(z['key_joints']).index('knee_angle_'+side)]
        ax.plot(z['t'][keep],q[keep],color=color,ls=style,label=side.upper()+' knee')
    ax.set_title(title);ax.grid(alpha=.2);ax.set_ylim(-100,20)
axes[0,0].legend();fig.supxlabel('Simulation time (s)');fig.supylabel('Knee coordinate (degrees; flexion negative)')
fig.savefig(H/'STAGE2_ABLATIONS.png',dpi=160)
(H/'STAGE2_ABLATIONS.alt.txt').write_text('Four physical knee-motion panels isolate feedback blocks. With both off, both knees cycle. Timing-only gives bilateral motion with poor coordination. Motor-reflex-only leaves one knee near extension. Both on leaves both knees near extension. Deafferented success does not qualify stage 2.')
