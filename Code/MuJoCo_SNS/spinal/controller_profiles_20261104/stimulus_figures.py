"""Few-cycle plots and paired metrics from actually logged layer stimuli."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(folder):
    z=np.load(folder/'capture/aligned_trace.npz');return dict(z=z,t=z['t_state'],tn=z['t_neural'],ti=z['t_neural_input'],names=z['neural_names'].tolist(),jn=z['joint_names'].tolist())

def main():
    p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path,required=True);p.add_argument('--trials',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    base=read(a.baseline);sources=dict(RG=a.trials/'RG',PF=a.trials/('PF_logged'if(a.trials/'PF_logged').exists()else'PF'),MN=a.trials/'MN');runs={label:read(path)for label,path in sources.items()}
    pf='PF_HIP-F_r'if all('PF_HIP-F_r'in r['names']for r in[base,*runs.values()])else'PF_BI_JOINT_RECT_r'
    cells=('RG_F_r',pf,'MN_iliacus_r');joints=('hip_flexion_r','knee_angle_r','ankle_angle_r');rows=[];proof=[]
    for label,r in runs.items():
        z=r['z'];stim=z['stimulus_nA'];active=np.where(np.any(np.nan_to_num(stim)!=0.,axis=1))[0]
        if not len(active):raise ValueError('No actually recorded pulse for '+label)
        currents=stim[active];currents=currents[np.isfinite(currents)]
        observed=dict(layer=label,source=str(sources[label]),pulse_samples=len(active),start_s=float(r['ti'][active[0]]),last_input_s=float(r['ti'][active[-1]]),currents_nA=sorted({float(value)for value in currents if value!=0.}))
        onset=observed['start_s'];before=(r['tn']>=onset-1.1)&(r['tn']<onset-.1);after=(r['tn']>=onset)&(r['tn']<onset+3.4)
        error={};prior=0.
        for cell in cells:
            y=z['neural_voltage'][:,r['names'].index(cell)];reference=np.interp(r['tn'],base['tn'],base['z']['neural_voltage'][:,base['names'].index(cell)]);delta=y-reference
            prior=max(prior,float(abs(delta[before]).max()));error[cell]=dict(maximum_delta_mV=float(abs(delta[after]).max()),rms_delta_mV=float(np.sqrt(np.mean(delta[after]**2))))
        observed.update(neural_response=error,maximum_pre_pulse_voltage_difference_mV=prior)
        if prior>1e-9:raise ValueError('Unmatched pre-pulse trajectory: '+label+' '+str(prior))
        proof.append(observed)
    if max(row['start_s']for row in proof)-min(row['start_s']for row in proof)>1e-9:raise ValueError('Use a common pulse onset for the layer comparison')
    start=proof[0]['start_s']-.6;end=proof[0]['start_s']+3.4
    plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none'})
    fig,axes=plt.subplots(7,1,figsize=(8.3,10),sharex=True,constrained_layout=True)
    styles=dict(Baseline=('#222222','-'),RG=('#CC6677','--'),PF=('#4477AA','-.'),MN=('#228833',':'))
    for label,r in [('Baseline',base),*runs.items()]:
        color,style=styles[label];mask=(r['tn']>=start)&(r['tn']<=end)
        if label!='Baseline':
            stimulus=r['z']['stimulus_nA'];total=np.nansum(stimulus,axis=1);axes[0].plot(r['ti'][mask],total[mask],style,color=color,label=label+' stimulation')
        for i,cell in enumerate(cells,1):axes[i].plot(r['tn'][mask],r['z']['neural_voltage'][mask,r['names'].index(cell)],style,color=color,lw=1.1,label=label)
        for i,joint in enumerate(joints,4):
            column=int(r['z']['joint_qpos_addresses'][r['jn'].index(joint)]);y=np.degrees(r['z']['qpos'][:,column]);y=-y if joint.startswith('knee')else y
            axes[i].plot(r['t'][mask],y[mask],style,color=color,lw=1.1)
    for ax in axes:
        ax.grid(alpha=.2);ax.set_xlim(start,end);ax.axvspan(proof[0]['start_s'],proof[0]['last_input_s']+.002,color='#DDCC77',alpha=.25)
    for ax,label in zip(axes,('Added current\n(nA)','RG-F\n(mV)',('Hip-F PF'if pf=='PF_HIP-F_r'else'RF PF')+'\n(mV)','Iliacus MN\n(mV)','Hip coordinate\n(deg)','Knee flexion\n(deg)','Ankle coordinate\n(deg)')):ax.set_ylabel(label)
    axes[1].legend(loc='upper right',ncol=4,fontsize=10);axes[-1].set_xlabel('Recorded time (s)')
    amplitudes={value for row in proof for value in row['currents_nA']}
    amplitude_label=', '.join(format(value,'g')for value in sorted(amplitudes))
    fig.suptitle('Gait2392: separately applied RG, PF, and MN current pulses\nGraded BiW2Lvar+; Bolen curriculum stage 2 (air); '+amplitude_label+' nA for 100 ms',fontsize=12)
    for ext in('png','pdf','svg'):fig.savefig(a.output/('HUNT_LAYER_STIMULUS.'+ext),dpi=180)
    plt.close(fig)
    (a.output/'PAIRED_STIMULUS_METRICS.json').write_text(json.dumps(dict(baseline=str(a.baseline),trials=proof,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Actual current pulses and recorded neural/joint responses; no inferred stimulus and no new dynamics during plotting'),indent=2))
if __name__=='__main__':main()
