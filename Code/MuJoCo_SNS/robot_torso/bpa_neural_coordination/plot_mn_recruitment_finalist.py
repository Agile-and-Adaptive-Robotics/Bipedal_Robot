"""Selected right-leg neural/actuator trace from byte-verified fresh challenge."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\source_domain_recruitment')

def main():
 audit=json.loads((ROOT/'STRONG_FRESH_AUDIT.json').read_text())
 if not(audit['accepted'] and audit['all_training_array_hashes_exact']):raise ValueError('Fresh strong proof required')
 stem=ROOT/'pressure/fresh_training_strong/w2lvar_case';cfg=json.loads((ROOT/'PREREGISTRATION.json').read_text());spec=next(s for s in cfg['training'] if s['name']=='train_strong');act=json.loads(Path(spec['path']).read_text())['actuators']['pect_r']
 with np.load(stem.with_suffix('.final_neural_state.npz')) as z:r=z['recruit_logs']
 with np.load(stem.with_suffix('.pressure.npz')) as z:p=z['samples'];pt=z['time'];i=z['names'].tolist().index('pect_r')
 with np.load(stem.with_suffix('.motor_state.npz')) as z:mt=z['time'];u=z['excitation'][:,z['actuator_names'].tolist().index('pect_r')]
 nt=(np.arange(len(r))-300+1)*.002
 plt.rcParams.update({'font.family':'Arial','font.size':11,'axes.labelsize':11,'axes.titlesize':11,'legend.fontsize':10,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig,axes=plt.subplots(4,1,figsize=(8.5,8.2),sharex=True,constrained_layout=True);blue='#4477AA';orange='#EE7733';purple='#AA3377';grey='#777777'
 axes[0].plot(nt,10*r[:,0,4],color=purple,lw=1.4,label='Recruitment-neuron input');axes[0].set_ylabel('Input (nA)');axes[0].set_title('(A) Predicted extension / measured-tension error');axes[0].legend(loc='upper left')
 axes[1].plot(nt,r[:,0,5],color=orange,lw=1.4,label='Recruitment relay');axes[1].set_ylabel('Relay (mV)');right=axes[1].twinx();right.plot(mt,u,color=blue,lw=1.3,ls='--',label='MN command');right.set_ylim(0,.75);right.set_ylabel('MN command (0–1)');lines=axes[1].lines+right.lines;axes[1].legend(lines,[l.get_label() for l in lines],loc='upper left');axes[1].set_title('(B) Neural relay and right-pectineus motor command')
 axes[2].plot(pt,p[:,2,i],color=blue,lw=1.5,label='Actual tension');axes[2].plot(mt,u*act['force_command_scale_N'],color=orange,lw=1.3,ls='--',label='Neural force reference');axes[2].plot(nt,r[:,0,1]+5.,color=purple,lw=1.1,ls=':',label='Predicted minimum + headroom');axes[2].set_ylabel('Tension (N)');axes[2].set_title('(C) Actual force and neural / sensory references');axes[2].legend(loc='upper left')
 axes[3].plot(pt,p[:,1,i],color=blue,lw=1.4,label='Actual pressure');axes[3].plot(pt,p[:,0,i],color=orange,lw=1.2,ls='--',label='Final pressure reference');axes[3].set_ylabel('Pressure (kPa)');axes[3].set_xlabel('Released-physics time (s)');axes[3].set_title('(D) Unchanged pressure-response dynamics');axes[3].legend(loc='upper left')
 for ax in axes:ax.set_xlim(17.15,17.85);ax.grid(True,alpha=.2);ax.axvline(17.294,color=grey,lw=1.,ls='-.');ax.spines['top'].set_visible(False)
 axes[0].set_ylim(0,13.);axes[0].legend(loc='lower left');axes[0].text(17.31,12.4,'Valve-only control ends',fontsize=10,color=grey,va='top')
 mask=(pt>=17.15)&(pt<=17.85);mmask=(mt>=17.15)&(mt<=17.85);nmask=(nt>=17.15)&(nt<=17.85);axes[2].set_ylim(0,1.25*max(float(np.max(p[mask,2,i])),float(np.max(u[mmask]*act['force_command_scale_N'])),float(np.max(r[nmask,0,1]+5.))));axes[3].set_ylim(0,1.3*float(np.max(p[mask,0:2,i])));fig.suptitle('Guided Gait2392 virtual-BPA standing: right pectineus\nFrozen 5 nA recruitment; known combined uncertainty; fresh exact repeat',fontsize=12)
 for suffix in ('png','pdf','svg'):fig.savefig(ROOT/('MN_RECRUITMENT_STRONG.'+suffix),dpi=220)
 plt.close(fig)
 caption='Right pectineus only; no bilateral rhythm or locomotion claim. The selected controller has frozen shared5nA extra MN current and10kPa valve correction. Known training challenge: force-map multiplier0.8, per-route fill72–80ms and delay14–16ms. This fresh20s study (19.4s released) is byte-identical to its training trajectories/states. Neural ticks align to2ms physics-tick ends; sensory transport adds6ms and samples have native engine latency. Predicted minimum uses delayedL/velocity, so it is not a current-pose force floor. MN command is logged normalized excitation, not native activation or raw membrane voltage. Neural force reference is excitation times declared pooled force scale; actual tension is measured physical output. Dashed/dotted styles distinguish lines independently of hue. The vertical marker denotes the separate matched MN-off/valve-only experiment ending at17.294s; the active trajectory is not assumed identical to that control. Strong lateral/rotational apparatus guidance, virtual actuator assumptions, and small foot-mesh penetration remain. Frozen independent validation3/4 passes; weaker0.78map case fails. No generalrobustness/hardware/gait claim.'
 (ROOT/'MN_RECRUITMENT_STRONG_caption.txt').write_text(caption,encoding='utf-8');print(str(ROOT/'MN_RECRUITMENT_STRONG.png'))

if __name__=='__main__':main()
