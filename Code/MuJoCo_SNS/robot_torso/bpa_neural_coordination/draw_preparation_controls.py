"""Plot matched neural-state pressure-initialization controls, not synthetic data."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');FOLDER=ROOT/'pinned_preparation'

def main():
 plt.rcParams.update({'font.family':'Arial','font.size':17,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':240})
 fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True,constrained_layout=True)
 for col,(suffix,title) in enumerate([('live','(A) Live neural feedback'),('replay','(B) Nominal sensory replay')]):
  for case,label,color,style,marker in [('harsh_'+suffix,'Prepared pressure','#0077BB','-','o'),('harsh_zero_'+suffix,'Pressure reset','#EE7733','--','^'),('harsh_cold_'+suffix,'Pressure + history reset','#555555',':','s')]:
   prefix=FOLDER/('syn6_prepared_'+case+'_v1');a=np.load(prefix.with_suffix('.apparatus.npz'));p=np.load(prefix.with_suffix('.pressure.npz'));r=json.loads(prefix.with_suffix('.result.json').read_text())
   t=a['samples'][:,0];idx=p['names'].tolist().index('vas_med_l');ch=p['channels'].tolist().index('actual_pressure_kpa')
   axes[0,col].plot(t,np.degrees(a['samples'][:,2]),label=label,color=color,ls=style,marker=marker,markevery=20,ms=5)
   axes[1,col].plot(p['time'],p['samples'][:,ch,idx],color=color,ls=style,marker=marker,markevery=20,ms=5)
   if r['failure']:axes[0,col].plot(t[-1],np.degrees(a['samples'][-1,2]),'kx',ms=8,mew=1.5)
  axes[0,col].set_title(title,loc='left');axes[0,col].set_ylabel('Left knee angle (°)');axes[1,col].set_ylabel('Pressure (kPa)');axes[1,col].set_xlabel('Physics time after release (s)')
  for row in range(2):axes[row,col].set_xlim(0,.8);axes[row,col].grid(alpha=.18)
 axes[0,0].legend(loc='best',fontsize=17)
 out=ROOT/'figures'
 for ext in ('png','pdf','svg'):fig.savefig(out/('virtual_bpa_preparation_controls.'+ext))
 plt.close(fig)
 (out/'virtual_bpa_preparation_controls.alt.txt').write_text('Matched neural-state controls under the same harsh virtual actuator changes. Prepared pressure maintains supported standing in both live-feedback and nominal sensory-replay conditions. Resetting pressure, with or without resetting queued commands, causes the left knee to flex until the characterized pressure-supported equilibrium limit is reached. Black crosses mark the last recorded steps, not observed falls. The lower plots show the difference in pressure reserve at release. All controls have identical neural states before the intentional actuator-state reset.',encoding='utf-8')
 (out/'PREPARATION_FIGURE_NOTES.md').write_text('# Prepared-state control figure\n\nMachine: easteregg2. Source files are the six `syn6_prepared_harsh*_v1` apparatus and pressure traces under `'+str(FOLDER)+'`.\n\nThe neural states and mechanical pose are matched before the deliberate pressure/command-history resets; each live/replay pair also has exactly matched initialized states after its reset. These are controlled initial-state manipulations, not identified physical venting procedures. The actuator uncertainty is the same fixed hypothetical slow/weak/direction-memory case. The rig, knee stops, and contact model remain present. Black crosses denote the recorded endpoint before an unsupported force-equilibrium update, not an observed fall.\n\nThe result is initialization sensitivity in this virtual Gait2392 pressure-actuated standing model, not a calibrated CAD result or a claim that feedback is unnecessary generally. Exact state audit: `'+str(FOLDER/'STATE_MATCH_AUDIT.json')+'`.\n',encoding='utf-8')
 print(out/'virtual_bpa_preparation_controls.png')

if __name__=='__main__':main()
