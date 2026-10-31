"""Shorten a saved W2L feedback comparison without changing any simulation."""
from pathlib import Path
import os,hashlib,json,shutil
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=next(p for p in Path(__file__).resolve().parents if (p/'AGENTS.md').exists());D=R/'Documentation/Reports and Papers/Dissertation';N=Path(__file__).parent;O=N/'previews'
S=R/'Code/MuJoCo_SNS/spinal/reports_20261027_codex'
plt.rcParams.update({'font.family':'Arial','font.size':15.5,'legend.fontsize':15.5,'axes.titlesize':15.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axes=plt.subplots(5,2,figsize=(6.5,7.5),sharex=True,layout='constrained');manifest=[]
for col,(name,title) in enumerate([('panel_control','Feedback intact'),('panel_del_all_feedback','Feedback disabled')]):
 p=S/(name+'.npz');z=np.load(p);t=z['t'];w=(t>=12)&(t<=15)
 for row,(kind,indices,label) in enumerate([('rgE',(0,1),'RG extensor\n(mV)'),('q',(0,3),'Hip angle\n(deg)'),('q',(1,4),'Knee angle\n(deg)'),('q',(2,5),'Ankle angle\n(deg)'),('heelF',(0,1),'Heel normal\nforce (N)')]):
  ax=axes[row,col]
  for index,color,style,leg in zip(indices,('#0072B2','#D55E00'),('-','--'),('Left','Right')):ax.plot(t[w],z[kind][w,index],color=color,ls=style,lw=1,label=leg)
  ax.grid(alpha=.17);ax.set_xlim(12,15)
  if col==0:ax.set_ylabel(label)
 axes[0,col].set_title(title);axes[-1,col].set_xlabel('Time (s)')
 manifest.append(dict(source=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),window_s=[12,15],channels=['rgE','q','heelF'],note='No PF or MN data in this saved record; no inferred traces. No smoothing or realignment.'))
for row in range(5):
 limits=[axes[row,col].get_ylim() for col in range(2)];lo=min(v[0] for v in limits);hi=max(v[1] for v in limits)
 for col in range(2):axes[row,col].set_ylim(lo,hi)
fig.legend(*axes[0,0].get_legend_handles_labels(),ncol=2,loc='outside upper center',frameon=False)
fig.savefig(O/'w2l_feedback_short.png',dpi=220);fig.savefig(O/'w2l_feedback_short.pdf');plt.close(fig)
target=D/'Figures/60-sns/w2l_feedback_short_20261031.pdf';backup=D/'Archive/session1_sns_20261031/Figures/60-sns/w2l_feedback_verified.pdf'
if not backup.exists():
 backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(D/'Figures/60-sns/w2l_feedback_verified.pdf',backup)
shutil.copy2(O/'w2l_feedback_short.pdf',target)
(O/'W2L_FEEDBACK_SHORT_PROVENANCE.json').write_text(json.dumps(manifest,indent=2))
print('Saved 3-second figure from existing records; old figure archived.')
