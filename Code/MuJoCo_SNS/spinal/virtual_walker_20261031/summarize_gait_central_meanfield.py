"""Mechanism/phase diagnostics from frozen studies; no fitting or body integration."""
import sys,json,numpy as np
from pathlib import Path
from scipy.signal import find_peaks
import gait_central_meanfield as M
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
R=M.ROOT;C=M.C
def load(path):
    with np.load(path) as z:return {k:z[k].copy() for k in z.files}
study=json.loads((R/'rustuna24/SUMMARY.json').read_text());assert C.X.digest(M.__file__)==study['provenance']['adapter_sha256']
trials=[json.loads(l) for l in (R/'rustuna24/TRIALS.jsonl').read_text().splitlines()];assert len(trials)==24
passing=[t for t in trials if t['diagnostics']['passed_offline_gate']];assert not passing,'Do not use failed-candidate report for a passing study'
frozen=load(R/'MEANFIELD_FROZEN.npz');soma=load(R/'SOMA_CONTROL.npz');reference=load(M.SOURCE/'GRADED_REFERENCE.npz')
best=study['best'];selected={};fresh={}
if '--plot-only' not in sys.argv:
    protocol=(R/'PROTOCOL.json').read_bytes();M.prepare();(R/'PROTOCOL.json').write_bytes(protocol)
    for key,u in M.STATE['inputs'].items():
        a,_=M.replay(C.clone(M.STATE['source']),M.STATE['names'],u,best['params'],fixture=key=='deaff');metrics=C.metrics(a,M.STATE['reference'][key],key)
        if key=='standing':metrics['RG_h_mean_absolute_error']=np.mean(abs(a['h_mean'][500:,[0,1,8,9]]-M.STATE['graded_h'][key][500:,[0,1,8,9]]),axis=0).tolist()
        else:
            hs=M.h_profile(M.STATE['reference'][key]['trace'],M.STATE['graded_h'][key]);hp=M.h_profile(a['trace'],a['h_mean']);metrics['h_peak_cycle_phase_error']=(np.full(4,.5) if hs is None or hp is None else abs((np.argmax(hp,axis=1)/64-np.argmax(hs,axis=1)/64+.5)%1-.5)).tolist()
        assert metrics==best['diagnostics']['fixtures'][key],'Fresh selected failed metrics differ'
        selected.update({key+'_'+f:a[f] for f in ('trace','events','somas','h_mean')});fresh[key]=metrics
    np.savez_compressed(R/'RUSTUNA_BEST_FAILED.npz',**selected,configuration_json=json.dumps(best['params']),dt_s=.002)
else:
    assert json.loads((R/'CENTRAL_CHECKPOINT.json').read_text())['fresh_best_metrics_exact']
    selected=load(R/'RUSTUNA_BEST_FAILED.npz');assert json.loads(str(selected['configuration_json']))==best['params']
def profile(trace):
    v=np.clip(trace[:,:16]/5.,0,1);p,_=find_peaks(v[1000:,1],prominence=.3,distance=100);p=p+1000
    cycles=[np.array([np.interp(np.linspace(a,b,64,endpoint=False),np.arange(len(v)),v[:,i]) for i in range(16)]) for a,b in zip(p[:-1],p[1:])]
    return np.mean(cycles,axis=0) if cycles else None
sourceprofile=profile(reference['deaff']);diagnostic={}
for label,z in (('soma_frozen',soma),('meanfield_frozen',frozen),('meanfield_fitted_failed',selected)):
    p=profile(z['deaff_trace']);diagnostic[label]=dict(phase_normalized_chemical_fraction_RMSE=None if p is None else float(np.sqrt(np.mean((p-sourceprofile)**2))),interpretation='Diagnostic waveform shape after each ownRG-Fcycle registration; does NOT replace raw-time acceptance gate or rescue failed candidates.')
checkpoint=dict(machine='easteregg2',physical_trials=0,completed_trials=24,offline_passes=len(passing),frozen_soma_control_exact_v1=True,source_adapter_sha256=study['provenance']['adapter_sha256'],frozen_mechanism=json.loads((R/'MEANFIELD_FROZEN.json').read_text()),best_failed=best,fresh_best_metrics_exact=True,phase_normalized_diagnostics=diagnostic,study_wall_s=study['this_process_wall_s'],offline_neural_seconds_per_wall_second=24*13/study['this_process_wall_s'],scope='Spike-envelopeNaP availability restores oscillation/within-cyclehphase but no prototype passes fullsource waveform+standingh gates. Zero bodytrials, no standing/walkingpromotion.')
C.X.write(R/'CENTRAL_CHECKPOINT.json',checkpoint)
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'svg.fonttype':'none','pdf.fonttype':42})
fig,ax=plt.subplots(6,2,figsize=(6.5,9.3),sharex='col');handles=[]
for j,key in enumerate(('standing','deaff')):
    t=np.arange(1,len(reference[key])+1)*.002
    for label,z,color,ls in (('Soma h',soma,'#d55e00','--'),('Envelope h',frozen,'#0072b2','-')):
        ax[0,j].plot(t,z[key+'_events'][:,0],color=color,ls=ls,lw=.8)
        for row,column in ((1,0),(2,1),(4,2),(5,20)):
            line,=ax[row,j].plot(t,z[key+'_trace'][:,column],color=color,ls=ls,lw=1.,label=label)
            if row==1 and j==0:handles.append(line)
        ax[3,j].plot(t,z[key+'_h_mean'][:,0],color=color,ls=ls,lw=1.)
    for row,column in ((1,0),(2,1),(4,2),(5,20)):
        line,=ax[row,j].plot(t,reference[key][:,column],color='#555555',ls=':',lw=1.2,label='Graded source')
        if row==1 and j==0:handles.append(line)
    ax[3,j].plot(t,frozen[key+'_graded_h'][:,0],color='#555555',ls=':',lw=1.2)
    ax[0,j].set_title('Recorded standing\nFrozen parameters' if key=='standing' else 'Deafferented drive\nFrozen parameters',fontsize=10);ax[0,j].set_xlim((1.,3.) if key=='standing' else (2.,9.8))
    ax[5,j].set_xlabel('Neural time (s)')
    for a in ax[:,j]:a.grid(alpha=.18)
for row,label in enumerate(['RG-E spikes\n(per2ms)','RG-E readout\n(mV)','RG-F readout\n(mV)','RG-E NaP h\n(availability)','PF S1 readout\n(mV)','Soleus MN\n(graded mV)']):ax[row,0].set_ylabel(label)
fig.subplots_adjust(left=.15,right=.98,top=.92,bottom=.23,wspace=.32,hspace=.20)
fig.legend(handles,['Soma-driven h','Spike-envelope h','Graded source'],ncol=2,frameon=False,loc='lower center',bbox_to_anchor=(.5,.12))
fig.text(.035,.025,'Only NaP inactivation changes; central chemical readouts derive from spikes.\nThe spike envelope is a declared mean-field signal, not a single-cell voltage.\nIt restores oscillations and h timing but fails the complete fidelity gates.\n256 central spiking cells; remaining778 cells graded. No body trial.',fontsize=10)
for suffix in ('pdf','svg','png'):fig.savefig(R/('FROZEN_H_MECHANISM.'+suffix),dpi=220,bbox_inches='tight',facecolor='white')
plt.close(fig)
text='# Mean-field NaP mechanism comparison\n\nMachine: easteregg2. **No candidate passed the complete offline gates; zero physical trials.** This is separate from preserved v1.\n\n'
text+='The frozen comparison changes only h_inf’s voltage signal from actual threshold/reset soma to the actual spike-filtered population envelope. Fast activation still uses actual soma voltage. This is a declared phenomenological population mechanism, not literal single-cell physiology. No original graded presynaptic voltage is used. All source weights, shared parameters, 200 µs clock, resets, adaptation, decoder, and fixtures are unchanged.\n\n'
text+='The soma-h control exactly reproduces every saved v1 trace/event/soma array. It latches in the driven fixture. Envelope-h restores bilateral cycles (2.482–2.633 s versus source2.132 s) and h peak phase (0–0.03125-cycle error), but standing h and waveform fidelity fail. This supports a causal role for the inactivation signal, not complete model equivalence.\n\n'
text+='The preregistered 24-trial calibration completed in '+f"{study['this_process_wall_s']:.2f}"+' s with no passer. Its best failed candidate has 2.388/2.482 s RG periods, zero bilateral phase error relative to the in-phase source, but driven chemical waveform RMSE0.62088 and standing extensor-h error0.25426. A fresh fixed-candidate replay exactly reproduces its metrics. No body or pulse trial was launched.\n\n'
text+='Phase-normalized waveform diagnostics (not replacement gates):\n\n'
for label,v in diagnostic.items():text+='- '+label+': '+str(v['phase_normalized_chemical_fraction_RMSE'])+'\n'
text+='\nElectrical checks:64-step disabled/default equality,64-step discarded-graded-voltage poison equality,32-step complete-state continuation equality, and exact recheck of all6500×16 graded reference central samples pass. Mean-field policy is included in state snapshots; importing a soma-policy snapshot into it is rejected. PF→MN keeps its nonzero fitted0.013344897908126326 value; no zero-output shortcut.\n\n'
text+='Actual architecture:256 central spiking cells (fourRG+12PF pools×16),778 graded cells includingMN/reflex pathways. The standing source is tonic; the deafferented source is nearly in-phase, not alternating gait. Measured study throughput is '+f"{24*13/study['this_process_wall_s']:.3f}"+' simulated neural-seconds per wall second including setup/driver overhead, not coupled-physics speedup.\n\n'
text+='Evidence/full paths:\n\n'
for name in ('PROTOCOL.json','STATE_AUDIT.json','SOMA_CONTROL.json','MEANFIELD_FROZEN.json','CENTRAL_CHECKPOINT.json','rustuna24/TRIALS.jsonl','rustuna24/ADAPTER_SOURCE.py','FROZEN_H_MECHANISM.pdf'):
    text+='- `'+str(R/name)+'`\n'
text+='\nReplay/figure program: `'+str(Path(__file__))+'`; invoke with `D:\\Anaconda\\envs\\myo\\python.exe`. Default repeats only the fixed failed candidate; `--plot-only` uses saved traces with no neural replay. No manuscript, canonical model, original controller, or v1 file changed. No owned worker remains.\n'
(R/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8');print(json.dumps(dict(offline_passes=0,physical_trials=0,phase_normalized_diagnostics=diagnostic,record=str(R/'EXPERIMENT_RECORD.md'))))
