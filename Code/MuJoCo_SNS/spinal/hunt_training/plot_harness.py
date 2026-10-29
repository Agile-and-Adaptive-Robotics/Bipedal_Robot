"""Plot measured restraint force with actual ground loading and knee motion."""
import os,sys
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent;stem=H/'runs/corrected_curriculum_v3/s4_height-0.06'
z=np.load(str(stem)+'.npz',allow_pickle=True);a=np.load(str(stem)+'.assistance.npz',allow_pickle=True)
x=a['samples'];t=x[:,0];bw=769.737726;keep=(t>=5)&(t<13)
plt.rcParams.update({'font.family':'Arial','font.size':10})
fig,axs=plt.subplots(3,1,figsize=(10,7),sharex=True,layout='constrained')
axs[0].plot(t[keep],-x[keep,10]/bw,color='#4477AA',label='Harness force, anatomical right')
axs[0].plot(t[keep],x[keep,11]/bw,color='#228833',ls='--',label='Harness force, superior')
axs[0].set_ylabel('Force / body weight');axs[0].legend(loc='upper right',fontsize=9)
for side,col,color,ls in [('Right',7,'#4477AA','-'),('Left',8,'#CC6677','--')]:axs[1].plot(t[keep],x[keep,col]/bw,color=color,ls=ls,label=side+' vertical GRF')
axs[1].set_ylabel('Foot force / body weight');axs[1].legend(loc='upper right',fontsize=9)
for side,color,ls in [('r','#4477AA','-'),('l','#CC6677','--')]:
    q=z['q'][:,list(z['key_joints']).index('knee_angle_'+side)];k=(z['t']>=5)&(z['t']<13)
    axs[2].plot(z['t'][k],q[k],color=color,ls=ls,label=side.upper()+' knee')
axs[2].set_ylabel('Knee coordinate (degrees)');axs[2].set_xlabel('Simulation time (s)');axs[2].legend(loc='upper right',fontsize=9)
for ax in axs:ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False)
fig.suptitle('Diagnostic supported stepping: bilateral motion, failed gait gates',fontsize=12)
fig.savefig(H/'HARNESS_DIAGNOSTICS.png',dpi=180)
(H/'HARNESS_DIAGNOSTICS.alt.txt').write_text('Three panels show anatomical-right and superior harness force, individual foot vertical forces, and knee coordinates over 5–13 seconds. Both knees cycle, but periods differ and foot-loading alternation fails. This is an assisted diagnostic, not accepted gait.')
