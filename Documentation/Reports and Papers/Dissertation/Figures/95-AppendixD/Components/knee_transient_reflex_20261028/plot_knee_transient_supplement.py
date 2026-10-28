from pathlib import Path
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Path('Documentation/Reports and Papers/Dissertation');A=D/'Archive/knee_and_robot_update_20261028';M=loadmat(A/'knee_transient_reflex_20261028/replay_summary.mat',squeeze_me=True);t=M['t'];P=M['P'];C=M['C'];ix=(t>=29)&(t<=40);x=t[ix]-30
plt.rcParams.update({'font.family':'Arial','font.size':14,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
f,ax=plt.subplots(3,2,figsize=(7.2,8.0),sharex=True)
ys=[P[ix,8],-np.column_stack([P[ix,0],C[ix,0]])*180/np.pi,-P[ix,1],P[ix,2:4],P[ix,4:6],P[ix,6:8]];labs=['Ia input, nA','Knee angle, deg','Velocity, rad/s','MN voltage, mV','Activation','Actuator force, N']
for k,(a,y,lab) in enumerate(zip(ax.flat,ys,labs)):
 a.plot(x,y if y.ndim==1 else y[:,0],color='#0072B2',lw=1.5);a.set_ylabel(lab);a.set_xlim(-1,10);a.grid(alpha=.2)
 if y.ndim>1:a.plot(x,y[:,1],'--',color='black',lw=1.3);a.legend(['Stimulus','No stimulus'] if k==1 else ['Extensor','Flexor'],fontsize=10,loc='upper right')
for a in ax[-1]:a.set_xlabel('Time from pulse onset, s')
f.tight_layout();f.savefig(D/'Figures/95-AppendixD/knee_transient_response.pdf');f.savefig(A/'knee_transient_response_final.png',dpi=140)
