from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
d=np.load(HERE/'dirusso_pso_warning_checked.npz')
plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42})
fig,axes=plt.subplots(3,2,figsize=(12,9),layout='constrained',sharex=True)
for j,name in enumerate(('Hip','Knee','Ankle')):
    ax=axes[j,0]
    ax.plot(d['t'],np.degrees(d['q'][:,j]),color='#4477AA',label='Right')
    ax.plot(d['t'],np.degrees(d['q'][:,j+3]),color='#EE6677',ls='--',label='Left')
    ax.set_ylabel(name+' angle (deg)');ax.grid(alpha=.2)
axes[0,0].legend(ncol=2)
axes[0,1].plot(d['t'],d['com'][:,2],color='#228833')
h0=d['com'][0,2]
axes[0,1].axhline(.8*h0,color='#AA3377',ls='--',label='Fall threshold')
axes[0,1].set_ylabel('COM height (m)');axes[0,1].legend()
axes[1,1].plot(d['t'],d['heel'][:,0],color='#EE6677',ls='--',label='Left')
axes[1,1].plot(d['t'],d['heel'][:,1],color='#4477AA',label='Right')
axes[1,1].set_ylabel('Heel force (N)');axes[1,1].legend(ncol=2)
axes[2,1].plot(d['t'],np.degrees(d['tilt']),color='#AA3377')
axes[2,1].set_ylabel('Pelvis pitch (deg)')
for ax in axes[-1]:
    ax.set_xlabel('Time (s)')
fig.suptitle('Rejected Di Russo-inspired gait2392 attempt\n90% feedforward support; excessive body rotation, then a fall at 6.982 s',fontsize=12)
for suffix in ('png','pdf'):
    fig.savefig(HERE/('dirusso_attempt.'+suffix),dpi=180,bbox_inches='tight')
(HERE/'dirusso_attempt.alt.txt').write_text('Six panels show hip, knee, and ankle joint angles, center-of-mass height with the 80 percent fall threshold, left and right heel forces, and pelvis pitch. The converted gait2392 model undergoes excessive hip flexion and body rotation beyond 180 degrees before falling at 6.982 seconds with 90 percent feedforward weight support. Each foot has only two merged heel-contact events. The survival-only search objective is rejected because it rewards this behavior. This is not stable walking.')
print('q ranges deg',np.degrees(d['q']).min(0),np.degrees(d['q']).max(0))
print('pitch range deg',np.degrees(d['tilt']).min(),np.degrees(d['tilt']).max())
