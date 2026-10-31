"""Matched saved ground-trial neural and mechanical signals, no new simulation."""
from pathlib import Path
import os,json,hashlib,shutil
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=next(p for p in Path(__file__).resolve().parents if (p/'AGENTS.md').exists());D=R/'Documentation/Reports and Papers/Dissertation';N=Path(__file__).parent;O=N/'previews'
plt.rcParams.update({'font.family':'Arial','font.size':17,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axes=plt.subplots(6,2,figsize=(7.2,8.5),sharex=True,layout='constrained');provenance=[]
colors=('#D55E00','#0072B2');styles=('-','--')
for col,(file,title) in enumerate([('spinal_run_spkbase_s3k.npz','Graded'),('spinal_run_spkbase_spiking.npz','Spiking (untuned)')]):
 p=R/'Code/MuJoCo_SNS/spinal'/file;z=np.load(p);t=z['t'];w=(t>=10)&(t<=12.5)
 nn=z['neuro_names'].tolist();an=z['key_acts'].tolist();jn=z['key_joints'].tolist()
 rows=[(z['neuro'],[nn.index('RG_E_r'),nn.index('RG_F_r')],'RG\n(mV)'),(z['neuro'],[nn.index('PF_E1_r'),nn.index('PF_F1_r')],'PF\n(mV)'),(z['act'],[an.index('vas_lat_r'),an.index('semimem_r')],'Muscle\nactivation')]
 rows += [(z['q'],[jn.index(j+'_r'),jn.index(j+'_l')],label+'\n(deg)') for j,label in [('hip_flexion','Hip'),('knee_angle','Knee'),('ankle_angle','Ankle')]]
 for row,(array,indices,label) in enumerate(rows):
  ax=axes[row,col]
  for index,color,style in zip(indices,colors,styles):ax.plot(t[w],array[w,index],color=color,ls=style,lw=1.1)
  ax.set_xlim(10,12.5);ax.grid(alpha=.15);ax.set_xticks([10,11,12])
  if col==0:ax.set_ylabel(label)
 axes[0,col].set_title(title,fontsize=17,weight='bold');axes[-1,col].set_xlabel('Time (s)')
 provenance.append({'source':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'window_s':[10,12.5],'neural':['RG_E_r','RG_F_r','PF_E1_r','PF_F1_r'],'activation':['vas_lat_r','semimem_r'],'q_units':'degrees, as stored; confirmed by original generator','neural_scope':'archived output/readout channels, not individual spiking membrane voltage','no_smoothing_or_time_shift':True})
for row in range(6):
 limits=[axes[row,c].get_ylim() for c in range(2)];bounds=(min(a[0] for a in limits),max(a[1] for a in limits))
 for c in range(2):axes[row,c].set_ylim(*bounds)
from matplotlib.lines import Line2D
fig.legend([Line2D([],[],color=colors[i],ls=styles[i]) for i in range(2)],['Extensor / Right','Flexor / Left'],loc='outside upper center',ncol=2,fontsize=15,frameon=False)
fig.savefig(O/'spiking_ground_neuromechanics_short_20261031.png',dpi=220);fig.savefig(O/'spiking_ground_neuromechanics_short_20261031.pdf');plt.close(fig)
target=D/'Figures/95-AppendixD/spiking_ground_neuromechanics_short_20261031.pdf';shutil.copy2(O/target.name,target)
(O/'SPIKING_GROUND_COMPARISON_PROVENANCE.json').write_text(json.dumps(provenance,indent=2))
print('Created matched neural/activation/kinematic figure from saved 16s ground trials; plotted10–12.5s.')
