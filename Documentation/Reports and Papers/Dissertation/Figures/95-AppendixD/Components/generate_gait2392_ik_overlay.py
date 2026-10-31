from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'Testing_Data').exists())
source=ROOT/'Code/MuJoCo_SNS/spinal/reports_20261028_publication_flow/human_ik/gait2392_ik_overlay.csv'
a=np.genfromtxt(source,delimiter=',',names=True)
out=Path(__file__).resolve().parent.parent
plt.rcParams.update({'font.family':'Arial','font.size':11.5,'axes.titlesize':12,'axes.labelsize':11.5,'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':11.5,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
f,axs=plt.subplots(4,2,figsize=(6.5,5.95),sharex=True,gridspec_kw={'height_ratios':[1,1,1,.45]})
x=a['cycle_percent']
for col,side in enumerate(['r','l']):
 for row,j in enumerate(['hip','knee','ankle']):
  ax=axs[row,col];mu=a[f'{side}_{j}_model_mean_deg'];sd=a[f'{side}_{j}_model_sd_deg'];ref=a[f'{side}_{j}_reference_deg']
  ax.plot(x,mu,color='#0072B2',lw=1.5,label='Model mean');ax.fill_between(x,mu-sd,mu+sd,color='#0072B2',alpha=.16,edgecolor='none');ax.plot(x,ref,color='black',ls='--',lw=1.4,label='Reference IK');ax.set_ylabel(j.capitalize()+' (deg)');ax.grid(alpha=.18);ax.set_xlim(0,100)
  if row==0:ax.set_title('Right leg' if side=='r' else 'Left leg')
 ax=axs[3,col];ax.plot(x,a[side+'_model_stance_probability'],color='#0072B2',lw=1.5);ax.plot(x,a[side+'_reference_stance'],color='black',ls='--',lw=1.4);ax.set_ylim(-.05,1.05);ax.set_yticks([0,1]);ax.set_ylabel('Stance\nfraction');ax.set_xlabel('Gait cycle (%)');ax.set_xlim(0,100);ax.set_xticks([0,25,50,75,100]);ax.grid(alpha=.18)
f.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(0.55,1.0),ncol=2,frameon=False)
f.subplots_adjust(left=.13,right=.96,bottom=.10,top=.88,hspace=.18,wspace=.36)
f.savefig(out/'gait2392_ik_overlay.pdf');f.savefig(Path(__file__).resolve().parent/'gait2392_ik_overlay_preview.png',dpi=150)
(Path(__file__).resolve().parent/'gait2392_ik_overlay_provenance.json').write_text(json.dumps({'source':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'operation':'Replot saved mean, standard deviation, reference, and stance columns only; no simulation or alignment changes.'},indent=2))
