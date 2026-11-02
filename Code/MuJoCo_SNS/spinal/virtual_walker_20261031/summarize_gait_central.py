"""Saved-study audit and fixed-candidate replay; never starts fitting or mechanics."""
from pathlib import Path
import json,hashlib,copy,sys,numpy as np
import gait_central_spiking as C
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
ROOT=C.ROOT
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
study=json.loads((ROOT/'rustuna24/SUMMARY.json').read_text());rows=[json.loads(s) for s in (ROOT/'rustuna24/TRIALS.jsonl').read_text().splitlines()]
assert digest(C.__file__)==study['provenance']['adapter_sha256'],'Replay must use exact numerical adapter tested in the completed study'
assert len(rows)==24 and not any(r['diagnostics']['passed_offline_gate'] for r in rows)
best=study['best'];config=best['params'];fresh={};data={}
if '--plot-only' not in sys.argv:
    C.prepare()
    for key,u in C.STATE['inputs'].items():
        a,controller=C.run_replay(C.clone(C.STATE['source']),C.STATE['names'],u,config,fixture=key=='deaff');metrics=C.metrics(a,C.STATE['reference'][key],key)
        assert metrics==best['diagnostics']['fixtures'][key],'Fixed-candidate fresh numerical metrics must reproduce exactly'
        fresh[key]=metrics;data[key]=a
    np.savez_compressed(ROOT/'RUSTUNA_BEST_FAILED_TRACE.npz',**{k+'_'+f:a[f] for k,a in data.items() for f in ('trace','events','somas')},names=np.array(C.STATE['names'])[data['standing']['selected']],configuration_json=json.dumps(config),dt_s=.002)
else:
    assert json.loads((ROOT/'CENTRAL_CHECKPOINT.json').read_text())['fresh_best_metrics_exact']
    with np.load(ROOT/'GRADED_REFERENCE.npz') as z:C.STATE['reference']={k:dict(trace=z[k]) for k in ('standing','deaff')}
    with np.load(ROOT/'RUSTUNA_BEST_FAILED_TRACE.npz') as z:data={k:{f:z[k+'_'+f] for f in ('trace','events','somas')} for k in ('standing','deaff')}
checkpoint=dict(machine='easteregg2',status='No central prototype passed BOTH preregistered source replay fixtures; physical gate remains closed.',physical_simulations=0,completed_trials=24,offline_gate_passes=0,study_wall_s=study['this_process_wall_s'],best_objective_candidate=best,fresh_best_metrics_exact=True,nominal=json.loads((ROOT/'PROTOTYPE_RESULT.json').read_text()),source_configuration=str(C.X.BEST),source_configuration_sha256=digest(C.X.BEST),source_graph=str(C.X.DIAGRAM),source_graph_sha256=digest(C.X.DIAGRAM),logical_source_cells=794,central_logical_populations=16,actual_central_spiking_cells=256,remaining_graded_cells=778,actual_MN_projection_gain=C.X.settings()['pf_mn_gain'],graded_MN_scope='Six sampled rightMN traces in waveform scorer; complete native MN/reflex graph remains graded and unaltered.',state_audit=str(ROOT/'ELECTRICAL_STATE_AUDIT.json'),candidate_replay_scope='Exact fresh deterministic neural-only replay of objective-ranked failed candidate; no robust walking/standing claim.')
C.X.write(ROOT/'CENTRAL_CHECKPOINT.json',checkpoint)
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'svg.fonttype':'none','pdf.fonttype':42})
def figure(data,label,configuration):
    fig,axes=plt.subplots(4,2,figsize=(6.5,7.6),sharex='col');colors=['#0072b2','#d55e00'];handles=[]
    for j,key in enumerate(('standing','deaff')):
        a=data[key];ref=C.STATE['reference'][key];t=np.arange(1,len(a['trace'])+1)*.002
        for k,color in enumerate(colors):
            axes[0,j].plot(t,a['events'][:,k],color=color,lw=.8,label=('E' if k==0 else 'F')+' spikes /2ms')
            for row,i in ((1,k),(2,2+k),(3,20+k)):
                line,=axes[row,j].plot(t,a['trace'][:,i],color=color,lw=1.,label=('E' if k==0 else 'F')+' spiking')
                baseline,=axes[row,j].plot(t,ref['trace'][:,i],color=color,ls=':',lw=1.2,label=('E' if k==0 else 'F')+' source')
                if row==1 and j==0:handles.extend([line,baseline])
        axes[0,j].set_title('Recorded standing inputs\nTonic source' if key=='standing' else 'Deafferented drive\nIn-phase source rhythm',fontsize=10)
        axes[3,j].set_xlabel('Neural time (s)');axes[0,j].set_xlim((1.,3.) if key=='standing' else (4.,10.))
        for ax in axes[:,j]:ax.grid(alpha=.18)
    for row,labely in enumerate(['RG E/F counts\n(actual events)','RG E/F readout\n(mV)','PF S1/S2 readout\n(mV)','Soleus / TA MN\n(graded mV)']):axes[row,0].set_ylabel(labely)
    fig.subplots_adjust(left=.15,right=.98,top=.89,bottom=.30,wspace=.34,hspace=.22)
    fig.legend(handles,['E spiking','E graded source','F spiking','F graded source'],loc='lower center',bbox_to_anchor=(.5,.16),ncol=2,frameon=False)
    fig.text(.04,.035,'Offline '+label+'. Sixteen central pools ×16 reset/refractory cells; all MN/reflex cells\nremain graded. Continuous lines use actual spike-filter presynaptic readouts.\nPF panel: S1 blue, S2 orange; MN panel: soleus blue, tibialis anterior orange.\nNo candidate passed both gates; no physical trial, alternating gait, or weaning claim.',fontsize=10)
    for suffix in ('pdf','svg','png'):fig.savefig(ROOT/(label+'_OFFLINE_DIAGNOSTIC.'+suffix),dpi=220,bbox_inches='tight',facecolor='white')
    plt.close(fig)
with np.load(ROOT/'PROTOTYPE_TRACE.npz') as z:nominal={k:{f:z[k+'_'+f] for f in ('trace','events','somas')} for k in ('standing','deaff')}
figure(nominal,'NOMINAL',C.DEFAULT);figure(data,'RUSTUNA_BEST_FAILED',config)
text='# Central-spiking syn6 prototype\n\nMachine: easteregg2. **No prototype passed both offline gates; no physical trial was launched.** Twenty-four mirrored-parameter Rustuna trials completed in '+f"{study['this_process_wall_s']:.2f}"+' s. Failed probes and all trial metrics remain recorded.\n\n'
text+='The accepted syn6 standing network has tonic bilateral RG-E and suppressed RG-F. Its archived source trace cannot establish a cadence/phase target. A separate deafferented source fixture uses DRIVE 2.5 nA plus a named 1 nA/100 ms left-RG-F startup pulse, and oscillates at 2.132 s almost bilaterally in-phase. This is a neural test fixture, not alternating gait.\n\n'
text+='Four NaP RG cells and twelve graded synergy PF cells are replaced by 256 actual reset/refractory LIF cells with source-derived conductance inputs and spike-filter readouts. The remaining 778 source cells, including MNs and reflex pathways, stay graded. W2L contributes the population/reset/adaptation implementation pattern, not its calcium parameters or successful controller. The syn6 PF→MN multiplier remains its fitted nonzero '+str(C.X.settings()['pf_mn_gain'])+'.\n\n'
text+='The first common-threshold/linear-readout probe suppressed a mixed PF whose source steady voltage is below its spike threshold. That failure is preserved as INITIAL_LINEAR_READOUT_*. The corrected PF decoder analytically converts measured firing rate and actual conductance to a passive-LIF voltage estimate; no graded central target supplies an output. Its nominal standing chemical-fraction RMSE is 0.03507 and six sampled MN RMSE is 0.000534 mV, but the driven prototype latches.\n\n'
text+='The objective-ranked failed candidate emits real bilateral RG/PF spikes, but driven RG periods are about 0.77–0.79 s versus 2.132 s, chemical waveform error remains excessive, and standing fidelity is lost. A fixed-candidate fresh neural replay exactly reproduces all its reported metrics. This is not a neural-training success.\n\n'
text+='Disabled/default64-step output parity, discarded-shadow-voltage poisoning64-step parity, and full local electrical-state32-step continuation all pass exactly. The latter includes soma membranes, NaP slow gates, refractory/adaptation/filter states, counts, graded/sparse backend continuation, and counters, without RNG. Physiological interpretation remains limited: membrane/readout voltages are in the source relative-rest frame.\n\n'
text+='A plausible failure mechanism is the changed NaP inactivation operating point: spike resets hold somas near 2.5–4 mV, below the graded burst envelope, so slow-channel availability differs. This remains a hypothesis; no extra parameter campaign or physical workaround was started.\n\n'
text+='Files (full paths):\n\n'
for filename in ('PROTOCOL.json','SOURCE_AUDIT.json','ELECTRICAL_STATE_AUDIT.json','CENTRAL_CHECKPOINT.json','RUSTUNA_MANIFEST.json','rustuna24/PROVENANCE.json','rustuna24/TRIALS.jsonl','rustuna24/ADAPTER_SOURCE.py','NOMINAL_OFFLINE_DIAGNOSTIC.pdf','RUSTUNA_BEST_FAILED_OFFLINE_DIAGNOSTIC.pdf'):
    text+='- `'+str(ROOT/filename)+'`\n'
text+='\nSource candidate: `'+str(C.X.BEST)+'`. Numerical adapter: `'+str(Path(C.__file__))+'`. Replay/figure program: `'+str(Path(__file__))+'`. No manuscript, canonical source model, or original controller file changed. No owned numerical process remains after this bounded replay.\n'
(ROOT/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8');print(json.dumps(dict(completed_trials=24,offline_passes=0,fresh_failed_candidate_exact=True,physical_trials=0,record=str(ROOT/'EXPERIMENT_RECORD.md'))))
