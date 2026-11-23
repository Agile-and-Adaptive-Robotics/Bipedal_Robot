"""Vector active-phase stimulus panels from native saved aligned arrays.

Read only: no resampling, smoothing, amplitude scaling, decimation, or dynamics.
Run with D:\Anaconda\python.exe.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
ROOT=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104')
paths={'Baseline':ROOT/'hunt_layer_stimulus_2nA_v2/Baseline/capture/aligned_trace.npz'}
paths.update({k:ROOT/f'hunt_active_phase_stimulus_v3/{k}/capture/aligned_trace.npz'for k in ('RG','PF','MN')})
data={key:np.load(p,allow_pickle=False)for key,p in paths.items()}
styles={'Baseline':('#222222','-'),'RG':('#CC6677','--'),'PF':('#4477AA','-.'),'MN':('#228833',':')}
names={'Baseline':'Baseline','RG':'RG pulse','PF':'PF pulse','MN':'MN pulse'}
window=(7.5,11.05)
plt.rcParams.update({'font.family':'Arial','font.size':10.2,'axes.labelsize':10.2,'xtick.labelsize':10.2,'ytick.labelsize':10.2,'legend.fontsize':10.2,'axes.linewidth':1.0,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
proof={'generator':str(Path(__file__).resolve()),'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'display_window_s':list(window),'response_scoring_window_s':[7.65,11.05],'response_end_excluded':True,'stimulus_input_window_s':[7.65,7.75],'stimulus_targets':{'RG':'RG_F_r','PF':'PF_HIP-F_r','MN':'MN_iliacus_r'},'source_files':{},'operations':['native sample selection on each signal clock','radians to degrees for joint coordinates'],'not_applied':['smoothing','interpolation','phase shift','amplitude rescaling','decimation','new dynamics'],'panels':{}}

def finish(fig,axes,stem,labels):
    for i,(ax,label) in enumerate(zip(axes,labels)):
        ax.set_ylabel(label)
        ax.grid(axis='both',alpha=.22,lw=.5)
        ax.axvspan(7.65,7.75,color='#DDCC77',alpha=.26,lw=0)
        ax.set_xlim(*window)
        ax.spines[['right','top']].set_visible(False)
        ax.text(-.205,1.01,f'({chr(65+i)})',transform=ax.transAxes,ha='left',va='bottom',weight='bold',fontsize=10.2)
        ax.tick_params(width=1,length=3)
    axes[-1].set_xlabel('Recorded time (s)')
    axes[-1].set_xticks([7.5,8.,9.,10.,11.])
    for suffix in ('pdf','svg','png'):
        fig.savefig(OUT/f'{stem}.{suffix}',dpi=220)
    plt.close(fig)

fig,axes=plt.subplots(4,1,figsize=(14/2.54,15/2.54),sharex=True)
fig.subplots_adjust(left=.225,right=.985,bottom=.095,top=.905,hspace=.31)
cells=['RG_F_r','PF_HIP-F_r','MN_iliacus_r']
panel_counts={}
for label,z in data.items():
    color,ls=styles[label]
    tn=z['t_neural'];ti=z['t_neural_input'];keep=(tn>=window[0])&(tn<=window[1]);ik=(ti>=window[0])&(ti<=window[1]);nn=z['neural_names'].tolist()
    if label!='Baseline':
        stimulus=z['stimulus_nA'][:,nn.index(proof['stimulus_targets'][label])]
        axes[0].plot(ti[ik],stimulus[ik],color=color,ls=ls,lw=1.2)
        actual=np.where(np.nan_to_num(stimulus)!=0.)[0]
        assert len(actual)==50 and np.all(stimulus[actual]==2.)
    for row,cell in enumerate(cells,1):
        axes[row].plot(tn[keep],z['neural_voltage'][keep,nn.index(cell)],color=color,ls=ls,lw=1.1,label=names[label])
    panel_counts[label]={'neural_samples':int(keep.sum()),'input_samples':int(ik.sum())}
axes[0].set_ylim(-.15,2.3);axes[0].set_yticks([0,1,2])
handles,legend_labels=axes[1].get_legend_handles_labels()
fig.legend(handles,legend_labels,loc='upper center',bbox_to_anchor=(.53,.994),ncol=4,frameon=False,columnspacing=.85,handlelength=1.8,handletextpad=.35)
finish(fig,axes,'AIR_stimulus_neural',['Added current\n(nA)','Right RG-F\n(mV)','Right hip-F PF\n(mV)','Right iliacus MN\n(mV)'])
proof['panels']['neural']={'dimensions_cm':[14,15],'minimum_font_pt':10.2,'counts':panel_counts,'clocks':['t_neural_input for current','t_neural for voltage'],'cells':cells}

fig,axes=plt.subplots(3,1,figsize=(14/2.54,12.5/2.54),sharex=True)
fig.subplots_adjust(left=.225,right=.985,bottom=.12,top=.89,hspace=.31)
joints=['hip_flexion_r','knee_angle_r','ankle_angle_r'];counts={}
for label,z in data.items():
    color,ls=styles[label];t=z['t_state'];keep=(t>=window[0])&(t<=window[1]);jn=z['joint_names'].tolist()
    for row,joint in enumerate(joints):
        qcol=z['joint_qpos_addresses'][jn.index(joint)]
        values=np.degrees(z['qpos'][keep,qcol])
        axes[row].plot(t[keep],values,color=color,ls=ls,lw=1.1,label=names[label])
    counts[label]=int(keep.sum())
handles,legend_labels=axes[0].get_legend_handles_labels()
fig.legend(handles,legend_labels,loc='upper center',bbox_to_anchor=(.53,.994),ncol=4,frameon=False,columnspacing=.85,handlelength=1.8,handletextpad=.35)
finish(fig,axes,'AIR_stimulus_mechanics',['Right hip\n(deg)','Right knee\n(deg)','Right ankle\n(deg)'])
proof['panels']['mechanics']={'dimensions_cm':[14,12.5],'minimum_font_pt':10.2,'counts':counts,'clock':'t_state','joint_names':joints,'knee_sign':'Native flexion-negative; no sign reversal'}
for key,p in paths.items():
    proof['source_files'][key]={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
for stem,alt in [('AIR_stimulus_neural','Separate +2 nA,100 ms current pulses at7.65 s to the right RG flexor, right hip-flexor PF, or right iliacus MN in the suspended BiW2Lvar+ model. Baseline and three stimulated trajectories retain all samples on their native clocks. Panels show applied current and right RG-F, hip-F PF, and iliacus MN membrane voltages.'),('AIR_stimulus_mechanics','Right hip, native flexion-negative knee, and ankle coordinates over7.5–11.05 s in matched baseline and separate RG/PF/MN stimulus trials. The shaded7.65–7.75 s interval marks the actual current pulse; curves retain all recorded samples with no smoothing or phase shift.')]:
    (OUT/f'{stem}.alt.txt').write_text(alt+'\n',encoding='utf-8')
    for suffix in ['pdf','svg','png','.alt.txt']:
        p=OUT/f'{stem}.{suffix}' if suffix!='.alt.txt' else OUT/f'{stem}.alt.txt'
        proof.setdefault('outputs',{})[p.name]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
(OUT/'AIR_stimulus_provenance.json').write_text(json.dumps(proof,indent=2),encoding='utf-8')
print(json.dumps({'panels':proof['panels'],'outputs':list(proof['outputs'])},indent=2))
