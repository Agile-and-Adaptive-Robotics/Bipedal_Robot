"""One source-supported SFA ablation and saved-data dynamic-observer summary."""
import json,sys,time,inspect,hashlib,zipfile,numpy as np
from pathlib import Path
from scipy.signal import find_peaks
import gait_dynamic_readout as D
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
C=D.C;M=D.M;R=D.ROOT

class ZeroSFA(D.DynamicCentral):
    def __init__(self,inner,names,configuration=None,**kwargs):
        actual=dict(C.DEFAULT if configuration is None else configuration);actual['rg_adapt_jump_nA']=0.
        warm=actual.copy();warm['rg_adapt_jump_nA']=C.DEFAULT['rg_adapt_jump_nA']
        super().__init__(inner,names,warm,**kwargs)
        self.configuration=actual;self.config[0]=0. # Warm-up above uses copies only; all actual adaptation stays exactlyzero.

class TimingZeroSFA(ZeroSFA):
    def __init__(self,inner,names,configuration=None,**kwargs):
        actual=dict(C.DEFAULT if configuration is None else configuration);scale=float(actual.pop('rg_h_tau_scale',2.133/2.246));assert .85<=scale<=1.05
        assert actual.get('rg_adapt_jump_nA',0)==0
        super().__init__(inner,names,actual,**kwargs);self.params[self.is_rg,11]*=scale;self.configuration['rg_h_tau_scale']=scale

def timing_calibration():
    scale=2.133/2.246;cfg=dict(C.DEFAULT,rg_adapt_jump_nA=0.,rg_h_tau_scale=scale)
    protocol=dict(machine='easteregg2',hypothesis='One analyticHunt-like timing calibration, notsearch: observedsourcehperiod2.133/uncorrectedactualhperiod2.246 sets sharedNaPh tau multiplier. AllfourRGpools shareit; fixed200usclock unchanged.',configuration=cfg,source_tau_s=.35,effective_tau_s=.35*scale,source_h_period_s=2.133,measured_candidate_h_period_s=2.246,change='Only h timeconstant versus savedSFA0 condition; no clock/weight/decoder/geometry/current/support change.',original_gates_retained=True,physical_trials=0,optimization_trials=0)
    C.X.write(R/'TIMING_PROTOCOL.json',protocol);prior=(R/'NETWORK_PROTOCOL.json').read_bytes();D.prepare_network();(R/'NETWORK_PROTOCOL.json').write_bytes(prior)
    c=TimingZeroSFA(C.clone(M.STATE['source']),M.STATE['names'],cfg)
    for u in M.STATE['inputs']['standing'][:64]:c.forward(u)
    initial=c.export_state()
    for u in M.STATE['inputs']['standing'][64:96]:c.forward(u)
    first=c.export_state();c.import_state(initial)
    for u in M.STATE['inputs']['standing'][64:96]:c.forward(u)
    second=c.export_state();assert all(np.array_equal(v,second[k],equal_nan=True) if v.dtype.kind in 'fc' else np.array_equal(v,second[k]) for k,v in first.items())
    rows={};data={};old=M.MeanFieldCentral;M.MeanFieldCentral=TimingZeroSFA;t0=time.perf_counter()
    try:
        for key,u in M.STATE['inputs'].items():
            a,_=M.replay(C.clone(M.STATE['source']),M.STATE['names'],u,cfg,fixture=key=='deaff');r=C.metrics(a,M.STATE['reference'][key],key);assert a['peak_adapt_nA']==0
            if key=='standing':r['RG_h_mean_absolute_error']=np.mean(abs(a['h_mean'][500:,[0,1,8,9]]-M.STATE['graded_h'][key][500:,[0,1,8,9]]),axis=0).tolist()
            else:
                hp=M.h_profile(a['trace'],a['h_mean']);hs=M.h_profile(M.STATE['reference'][key]['trace'],M.STATE['graded_h'][key]);r['h_peak_cycle_phase_error']=(np.full(4,.5) if hp is None else abs((np.argmax(hp,axis=1)/64-np.argmax(hs,axis=1)/64+.5)%1-.5)).tolist()
            rows[key]=r;data[key]=a
    finally:M.MeanFieldCentral=old
    passed=all(r['passed'] for r in rows.values()) and max(rows['standing']['RG_h_mean_absolute_error'])<=.15 and max(rows['deaff']['h_peak_cycle_phase_error'])<=.15
    result=dict(configuration=cfg,passed_original_offline_gate=passed,complete_state32_exact=True,fixtures=rows,wall_s=time.perf_counter()-t0,physical_trials=0,optimization_trials=0)
    C.X.write(R/'TIMING.json',result);np.savez_compressed(R/'TIMING.npz',**{k+'_'+f:a[f] for k,a in data.items() for f in ('trace','events','somas','h_mean')},**{k+'_graded_h':v for k,v in M.STATE['graded_h'].items()},configuration_json=json.dumps(cfg),dt_s=.002);print(json.dumps(result),flush=True)

def fresh_timing():
    cfg=json.loads((R/'TIMING.json').read_text())['configuration'];prior=(R/'NETWORK_PROTOCOL.json').read_bytes();D.prepare_network();(R/'NETWORK_PROTOCOL.json').write_bytes(prior)
    old=M.MeanFieldCentral;M.MeanFieldCentral=TimingZeroSFA;checks={}
    try:
        with np.load(R/'TIMING.npz') as saved:
            for key,u in M.STATE['inputs'].items():
                a,_=M.replay(C.clone(M.STATE['source']),M.STATE['names'],u,cfg,fixture=key=='deaff');checks[key]={f:bool(np.array_equal(a[f],saved[key+'_'+f])) for f in ('trace','events','somas','h_mean')};assert all(checks[key].values())
    finally:M.MeanFieldCentral=old
    C.X.write(R/'TIMING_FRESH.json',dict(fresh_process_all_recorded_arrays_exact=True,configuration=cfg,checks=checks,physical_trials=0));print(json.dumps(dict(fresh_timing_exact=True,checks=checks)),flush=True)

def zero_ablation():
    protocol=dict(machine='easteregg2',hypothesis='Original graded sourceNaPRG has no extra spike-triggered adaptation. SFA0.3 was added during the earlier soma-reset workaround; test exactlyzero while sourceNaPh and the independently validated dynamicdecoder remain active.',change='Only shared RG adaptationjump0.3→0; allweights,hkinetics,observer,memory,reset,refractory,clock,initialization,PFdecoder unchanged. Exactzero, notnearzero.',source_default=str(R/'DYNAMIC_FROZEN.json'),configuration=dict(C.DEFAULT,rg_adapt_jump_nA=0.),physical_trials=0,optimization_trials=0,acceptance='Originalbothoffline+h gates unchanged. No search, support, MNoverride, or bodytrial.')
    C.X.write(R/'ZERO_SFA_PROTOCOL.json',protocol);original_protocol=(R/'NETWORK_PROTOCOL.json').read_bytes();D.prepare_network();(R/'NETWORK_PROTOCOL.json').write_bytes(original_protocol)
    a=ZeroSFA(C.clone(M.STATE['source']),M.STATE['names']);b=ZeroSFA(C.clone(M.STATE['source']),M.STATE['names'],poison_shadow=True)
    for u in M.STATE['inputs']['standing'][:64]:assert np.array_equal(a.forward(u),b.forward(u))
    start=a.export_state()
    for u in M.STATE['inputs']['standing'][64:96]:a.forward(u)
    first=a.export_state();a.import_state(start)
    for u in M.STATE['inputs']['standing'][64:96]:a.forward(u)
    second=a.export_state();checks={k:bool(np.array_equal(v,second[k],equal_nan=True)) if v.dtype.kind in 'fc' else bool(np.array_equal(v,second[k])) for k,v in first.items()};assert all(checks.values()) and np.all(a.adapt==0)
    C.X.write(R/'ZERO_SFA_STATE_AUDIT.json',dict(complete_continuation32_exact=True,shadow_poison64_exact=True,adaptation_exactly_zero=True,arrays=checks))
    rows={};data={};old=M.MeanFieldCentral;M.MeanFieldCentral=ZeroSFA;t0=time.perf_counter()
    try:
        for key,u in M.STATE['inputs'].items():
            a,_=M.replay(C.clone(M.STATE['source']),M.STATE['names'],u,protocol['configuration'],fixture=key=='deaff');r=C.metrics(a,M.STATE['reference'][key],key);assert a['peak_adapt_nA']==0
            if key=='standing':r['RG_h_mean_absolute_error']=np.mean(abs(a['h_mean'][500:,[0,1,8,9]]-M.STATE['graded_h'][key][500:,[0,1,8,9]]),axis=0).tolist()
            else:
                hp=M.h_profile(a['trace'],a['h_mean']);hs=M.h_profile(M.STATE['reference'][key]['trace'],M.STATE['graded_h'][key]);r['h_peak_cycle_phase_error']=(np.full(4,.5) if hp is None else abs((np.argmax(hp,axis=1)/64-np.argmax(hs,axis=1)/64+.5)%1-.5)).tolist()
            rows[key]=r;data[key]=a
    finally:M.MeanFieldCentral=old
    passed=all(r['passed'] for r in rows.values()) and max(rows['standing']['RG_h_mean_absolute_error'])<=.15 and max(rows['deaff']['h_peak_cycle_phase_error'])<=.15
    report=dict(configuration=protocol['configuration'],passed_offline_gate=passed,fixtures=rows,wall_s=time.perf_counter()-t0,physical_trials=0,optimization_trials=0)
    C.X.write(R/'ZERO_SFA.json',report);np.savez_compressed(R/'ZERO_SFA.npz',**{k+'_'+f:a[f] for k,a in data.items() for f in ('trace','events','somas','h_mean')},**{k+'_graded_h':v for k,v in M.STATE['graded_h'].items()},configuration_json=json.dumps(protocol['configuration']),dt_s=.002)
    print(json.dumps(report),flush=True)

def load(path):
    with np.load(path) as z:return {k:z[k].copy() for k in z.files}

def cycle_diagnostics(trace,h,events=None):
    peaks,_=find_peaks(h[1000:,0],prominence=.15,distance=250);peaks=peaks+1000;v=np.clip(trace[:,:16]/5,0,1);profiles=[]
    for a,b in zip(peaks[:-1],peaks[1:]):profiles.append(np.array([np.interp(np.linspace(a,b,64,endpoint=False),np.arange(len(v)),v[:,i]) for i in range(16)]))
    out=dict(h_peaks_s=(peaks*.002).tolist(),h_period_s=float(np.median(np.diff(peaks))*.002) if len(peaks)>2 else None)
    if events is not None:
        density=np.convolve(events[:,0]/(.002*16),np.ones(50)/50,mode='same');p,_=find_peaks(density[1000:],prominence=10,distance=300);out['actual_event_burst_period_s']=float(np.median(np.diff(p))*.002) if len(p)>2 else None
    return out,np.mean(profiles,axis=0) if profiles else None

def burst_onsets(v,dt=.002):
    if not np.isfinite(v).all():raise ValueError('Nonfinite signal cannot be a neural burst')
    q=np.convolve(np.clip(v/5,0,1),np.ones(round(.1/dt))/round(.1/dt),mode='same');active=False;on=0;off=0;times=[];dwell=round(.05/dt)
    for k,value in enumerate(q):
        if not active:
            on=on+1 if value>=.7 else 0
            if on>=dwell:
                active=True;off=0;t=(k-dwell+1)*dt
                if not times or t-times[-1]>=.6:times.append(t)
        else:
            off=off+1 if value<=.3 else 0
            if off>=dwell:active=False;on=0;off=0
    return np.array([t for t in times if t>=2.])

def audit_burst_scorer():
    dt=.002;t=np.arange(6000)*dt;base=np.where((t%2.132)/2.132<.75,5.,0.);ripple=np.where(base>0,np.clip(base+1.5*np.sin(2*np.pi*t/.24),0,5),.2*np.maximum(0,np.sin(2*np.pi*t/.24)))
    clean=burst_onsets(base);noisy=burst_onsets(ripple);assert len(clean)==len(noisy) and np.max(abs(clean-noisy))<.05
    assert len(burst_onsets(np.zeros_like(t)))==0 and len(burst_onsets(np.ones_like(t)*5))==0
    rejected=False
    try:burst_onsets(np.full_like(t,np.nan))
    except ValueError:rejected=True
    assert rejected
    unit=dict(clean_burst_count=len(clean),ripple_burst_count=len(noisy),max_onset_shift_s=float(np.max(abs(clean-noisy))),flat_dead_and_tonic_rejected=True,nonfinite_rejected=True,scope='Independent synthetic ripple/flat/NaN check. Scorer audit after observedfailure; not retroactivepreregistration.')
    reference=load(M.SOURCE/'GRADED_REFERENCE.npz');candidate=load(R/'TIMING.npz');original=json.loads((R/'TIMING.json').read_text());names=reference['names'][:16].tolist()
    measurements={};alltimes={}
    for label,trace in (('source',reference['deaff']),('candidate',candidate['deaff_trace'])):
        rows={};times={}
        for i in (0,1,8,9):
            a=burst_onsets(trace[:,i]);times[i]=a;rows[names[i]]=dict(onsets_s=a.tolist(),cycles=len(a),period_s=float(np.median(np.diff(a))) if len(a)>2 else None)
        alltimes[label]=times;measurements[label]=rows
    period_pass=all(measurements['candidate'][name]['cycles']>=3 and measurements['candidate'][name]['period_s'] is not None and abs(measurements['candidate'][name]['period_s']/measurements['source'][name]['period_s']-1)<=.2 for name in measurements['source'])
    def phase(times):
        rows=[]
        for a,b in zip(times[1][:-1],times[1][1:]):
            other=times[9][(times[9]>=a)&(times[9]<b)]
            if len(other):rows.append((other[0]-a)/(b-a))
        return float(np.angle(np.mean(np.exp(2j*np.pi*np.array(rows))))/(2*np.pi)%1) if rows else None
    ps=phase(alltimes['source']);pc=phase(alltimes['candidate']);phaseerror=abs((pc-ps+.5)%1-.5) if ps is not None and pc is not None else .5
    def hp(h,times):
        rows=[]
        for a,b in zip(times[1][:-1],times[1][1:]):rows.append(np.array([np.interp(np.linspace(a,b,64,endpoint=False),np.arange(len(h))*.002,h[:,i]) for i in (0,1,8,9)]))
        return np.mean(rows,axis=0)
    sourceh=hp(candidate['deaff_graded_h'],alltimes['source']);newh=hp(candidate['deaff_h_mean'],alltimes['candidate']);he=abs((np.argmax(newh,axis=1)/64-np.argmax(sourceh,axis=1)/64+.5)%1-.5)
    waveform=original['fixtures']['deaff']['chemical_fraction_RMSE']<.2;standing=original['fixtures']['standing']['passed'] and max(original['fixtures']['standing']['RG_h_mean_absolute_error'])<=.15
    result=dict(status='SCORER_REVIEW_PREVIEW_ONLY_NOT_ACCEPTED',rule='Same source/candidate100ms averagedchemicalsignal;0.7/0.3on/offhysteresis;50msdwell;0.6sminimumcycle;onsetsafter2s. No alteredcontroller/physics/curve. hphase uses these coherentburstboundaries rather than localripplemaxima.',independent_unit=unit,measurements=measurements,period_pass=period_pass,bilateral_phase_error=phaseerror,h_cycle_peak_errors=he.tolist(),original_raw_waveform_gate_pass=waveform,original_standing_gate_pass=standing,preview_all_gates_pass=bool(period_pass and phaseerror<=.1 and max(he)<=.15 and waveform and standing),original_gate_pass=original['passed_original_offline_gate'],physical_trials=0,scope='Originalnegativecheckpointandrawmetricsretained. Rootreviewrequired before replacingburstcounter; no bodypromotion here.')
    C.X.write(R/'BURST_SCORER_REVIEW.json',result);print(json.dumps(result),flush=True)

def render_timing():
    assert json.loads((R/'TIMING_FRESH.json').read_text())['fresh_process_all_recorded_arrays_exact']
    z=load(R/'TIMING.npz');ref=load(M.SOURCE/'GRADED_REFERENCE.npz');result=json.loads((R/'TIMING.json').read_text());review=json.loads((R/'BURST_SCORER_REVIEW.json').read_text())
    plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(5,2,figsize=(6.5,8.2),sharex='col');handles=[];colors=['#0072b2','#d55e00']
    for j,key in enumerate(('standing','deaff')):
        t=np.arange(1,len(z[key+'_trace'])+1)*.002
        for k,color in enumerate(colors):
            axes[0,j].plot(t,z[key+'_events'][:,8+k],color=color,lw=.8)
            for row,col in ((1,8+k),(2,11-k),(3,20+k)):
                actual,=axes[row,j].plot(t,z[key+'_trace'][:,col],color=color,lw=1.)
                target,=axes[row,j].plot(t,ref[key][:,col],color=color,ls=':',lw=1.2)
                if row==1 and j==0:handles.extend([actual,target])
            axes[4,j].plot(t,z[key+'_h_mean'][:,8+k],color=color,lw=1.);axes[4,j].plot(t,z[key+'_graded_h'][:,8+k],color=color,ls=':',lw=1.2)
        axes[0,j].set_title('Standing inputs\nNeural replay' if key=='standing' else 'Deafferented drive\nNeural replay',fontsize=10);axes[0,j].set_xlim((1.,3.) if key=='standing' else (2.,9.8));axes[4,j].set_xlabel('Neural time (s)')
        for ax in axes[:,j]:ax.grid(alpha=.18)
    for row,label in enumerate(['Right RG events\n(per 2 ms)','Right RG E/F\n(readout mV)','PF S2/S1\n(readout mV)','Soleus / TA MN\n(graded mV)','Right RG E/F h\n(availability)']):axes[row,0].set_ylabel(label)
    fig.subplots_adjust(left=.15,right=.98,top=.91,bottom=.26,hspace=.22,wspace=.32)
    fig.legend(handles,['E: spiking readout','E: graded source','F: spiking readout','F: graded source'],loc='lower center',bbox_to_anchor=(.5,.125),ncol=2,frameon=False)
    fig.text(.035,.025,'Timing-calibrated central-spiking candidate; actual AP counts shown separately.\nPF: S2 blue/S1 orange. MN: soleus blue/TA orange; all right-side.\nRaw driven chemical RMSE 0.14977; original peak scorer fails (under review).\n256 spiking cells; 778 graded. No mechanical simulation or gait claim.',fontsize=10)
    for suffix in ('pdf','svg','png'):fig.savefig(R/('TIMING_CANDIDATE_NEURAL.'+suffix),dpi=220,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    params=result['configuration'];core=hashlib.sha256((inspect.getsource(ZeroSFA)+inspect.getsource(TimingZeroSFA)).encode()).hexdigest()
    dependencies=[dict(path=str(Path(p)),sha256=C.X.digest(p)) for p in (D.__file__,D.Q.__file__,M.__file__,C.__file__,C.X.__file__)]
    cfg=dict(machine='easteregg2',source_family='syn6',source_standing_configuration=str(C.X.BEST),source_standing_configuration_sha256=C.X.digest(C.X.BEST),source_standing_parameters=C.X.settings(),source_graph=str(C.X.DIAGRAM),source_graph_sha256=C.X.digest(C.X.DIAGRAM),source_model=str(C.F.R.MODEL),source_model_sha256=C.X.digest(C.F.R.MODEL),controller_file=str(Path(__file__)),controller_class='TimingZeroSFA',controller_numerical_core_sha256=core,parameters=params,observer_constants=dict(cells_per_population=16,neural_step_s=.002,spiking_step_s=.0002,spike_history_tau_s=.025,support_min_events=2),actual_spiking_cells=256,remaining_graded_cells=778,dependencies=dependencies,fresh_neural_arrays_exact=True,original_scorer_gate_pass=result['passed_original_offline_gate'],corrected_scorer_preview_pass=review['preview_all_gates_pass'],status='Neural candidate; scorer review pending; physical_trials0. No automatic promotion.')
    C.X.write(R/'TIMING_CANDIDATE.json',cfg);print(json.dumps(dict(configuration=str(R/'TIMING_CANDIDATE.json'),figure=str(R/'TIMING_CANDIDATE_NEURAL.pdf'),physical_trials=0)),flush=True)

def finalize_record():
    timing=json.loads((R/'TIMING.json').read_text());review=json.loads((R/'BURST_SCORER_REVIEW.json').read_text());fresh=json.loads((R/'TIMING_FRESH.json').read_text());p=json.loads((R/'TIMING_PROTOCOL.json').read_text())
    checkpoint=json.loads((R/'CENTRAL_CHECKPOINT.json').read_text());checkpoint.update(timing_candidate=timing,timing_protocol=p,fresh_timing=fresh,burst_scorer_review=review,physical_trials=0,optimization_trials=0,status='Original scorer remains negative. Corrected burst-scorer preview passes; parent review pending. No body trial or automatic promotion.');C.X.write(R/'CENTRAL_CHECKPOINT.json',checkpoint)
    text='''# Causal dynamic spike readout: neural checkpoint

Machine: **easteregg2**. One timing-calibrated candidate now matches standing and rhythmic source waveforms. Its original local-maximum cadence scorer still fails; a separately labelled burst-scorer review passes. **No mechanical trial or automatic promotion.** There have been zero optimizer trials and one analytically derived timing calibration in this study.

## Mechanism and independent validation

Actual AP events update firing-rate, total-conductance, and NaP-inactivation histories with the same 25 ms kernel. The inverse firing-rate readout uses event-conditioned conductance and h, rather than pairing old rate with instantaneous conductance. A declared two-event uncertainty regularizer shrinks poorly supported estimates toward silence. Its Poisson-derived confidence model is an engineering assumption; the simulated cells remain deterministic. No driving-current target or original graded presynaptic voltage enters the decoder, and no spikes are fabricated.

Seven observable held-out conductance/current/channel step and ramp fixtures pass the preregistered transient criteria: worst chemical-activation RMSE 0.13310, h-infinity RMSE 0.11618, and inhibitory off lag 54–80 ms. Four factor controls use exactly the same actual spikes. Instantaneous inverse: 3/7; event matching alone: 3/7; confidence regularization alone: 7/7; combined: 7/7. These tests therefore do not establish that event matching alone was necessary. Zero observed events produce exact silence. Positive subthreshold voltage without any AP remains explicitly unidentifiable.

## Full-network comparisons and one timing calibration

| Neural condition | Standing chemical RMSE | Driven raw chemical RMSE | Status |
|---|---:|---:|---|
| Frozen dynamic observer, added spike SFA 0.3 nA/event | 0.01879 | 0.56685 | Original gates fail |
| Same observer, extra spike SFA exactly zero | 0.01050 | 0.41458 | Original gates fail |
| Zero SFA plus shared NaP-h timing calibration | 0.01073 | 0.14977 | Waveform gates pass; original cadence scorer fails |

The original graded NaP source has no separate spike-triggered adaptation current. Removing the added SFA is one source-supported mechanism contrast; its actual current is exactly zero. No synapse, geometry, muscle law, support condition, or motor override changes in these comparisons.

The single timing calibration uses source h-cycle period 2.133 s divided by the observed zero-SFA candidate period 2.246 s: multiplier 0.9496883348174533. All four RG populations share this multiplier. Their h time constant changes from 0.35 s to 0.3323909171861087 s. The spiking clock remains 200 microseconds; the chemical/network step remains 2 ms. This is an explicit intrinsic-kinetics calibration, not an exact-source claim or clock change. Source PF-to-MN multiplier remains nonzero at 0.013344897908126326.

The calibrated candidate gives standing sampled-MN RMSE 0.0000605680 mV and driven sampled-MN RMSE 0.0228328 mV. Every recorded trace, actual-event-count, soma, and h array repeats bit-exactly in a fresh process for both fixtures. Complete-state continuation is exact. Earlier default, shadow-poison, scalar/JIT, and zero-event audits remain saved.

## Scorer review, not retroactive acceptance

The old local-peak counter splits spike-readout ripple into extra bursts. Its negative result is preserved. The review uses the same 100 ms averaged chemical signal, 0.7/0.3 hysteresis, 50 ms dwell, and 0.6 s minimum spacing for both source and candidate. It reports E periods 2.132 versus 2.130 s and F periods 2.132 versus 2.149 s, four E/three F onsets per side, bilateral phase error zero, and h-cycle phase errors 0/0.03125. The source rhythmic fixture is bilateral in-phase; this is not alternating locomotion.

An independent ripple check preserves all five true bursts with maximum 22 ms onset shift; dead, tonic, and nonfinite signals are rejected. The preview passes the original numerical waveform and period/phase thresholds with coherent burst boundaries. Its status remains `SCORER_REVIEW_PREVIEW_ONLY_NOT_ACCEPTED`. No gate has silently replaced the original scorer, and no body run follows without review.

## Scope, timing, and reproducibility

This is neural replay: 256 actual central spiking cells represent 16 RG/PF logical pools; 778 remaining cells, including MN/reflex pathways, remain graded. It is not an all-spiking controller or physical standing/walking result. Source muscle policy remains 92 MN channels/86 force-active routes; no physical model was integrated here.

Warm-JIT replay and assessment of 3 s standing inputs plus 10 s driven inputs took 6.4413 s wall: 2.02 times neural duration. Compilation and source construction were outside that timer. This is not a full-coupled MuJoCo speed comparison.

All files below are under `D:\\Users\\Ben Bolen\\NeuroCampaign_20261031\\spiking\\Gait_dynamic_readout_v4\\`:

- `TIMING_CANDIDATE.json`: all four consumed parameters, source/dependency hashes, class and full source paths, and explicit negative/pending gate status.
- `TIMING_CANDIDATE_NEURAL.pdf`, `.svg`, and `.png`: actual right-side RG spike counts, decoded RG/PF, graded MN, and NaP-h versus source; no motor replay or joint-motion claim.
- `TIMING_PROTOCOL.json`, `TIMING.json`, `TIMING.npz`, and `TIMING_FRESH.json`: frozen protocol, original negative scorer, all neural data, and exact fresh-array checks.
- `BURST_SCORER_REVIEW.json`: separate scorer preview and independent checks.
- `TRANSIENT_PROTOCOL.json`, `TRANSIENT_VALIDATION.json`, and `TRANSIENT_TRACES.npz`: preregistered component tests and factor controls.
- `DYNAMIC_FROZEN.json` and `ZERO_SFA.json`: prior frozen negatives remain unchanged; v1/v2/v3 records also remain intact.
- `CENTRAL_CHECKPOINT.json` and `PDF_LAYOUT_QA.json`: combined status and figure QA.

Numerical program: `D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal\\virtual_walker_20261031\\gait_dynamic_readout.py`.
Analysis/calibration program: `D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal\\virtual_walker_20261031\\summarize_gait_dynamic.py`.
Interpreter: `D:\\Anaconda\\envs\\myo\\python.exe`.

Reproduce the candidate with `summarize_gait_dynamic.py --timing`; verify saved arrays in a new process with `--fresh-timing`. Audit the scorer with `--audit-scorer`. Render saved data only with `--render-timing`. None of these flags starts a physical simulation or optimizer. No owned worker remains.
'''
    text+='\n`SOURCE_SNAPSHOT.zip` freezes this numerical/analysis version and its five direct controller dependencies; the candidate JSON records their exact source paths and hashes. Raw traces remain outside the repository.\n'
    (R/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8')
    with zipfile.ZipFile(R/'SOURCE_SNAPSHOT.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
        for path in (Path(__file__),Path(D.__file__),Path(D.Q.__file__),Path(M.__file__),Path(C.__file__),Path(C.X.__file__),R/'TIMING_CANDIDATE.json',R/'TIMING_PROTOCOL.json'):archive.write(path,path.name)
    checkpoint['source_snapshot']=dict(path=str(R/'SOURCE_SNAPSHOT.zip'),sha256=C.X.digest(R/'SOURCE_SNAPSHOT.zip'));C.X.write(R/'CENTRAL_CHECKPOINT.json',checkpoint)
    C.X.write(R/'FILE_MANIFEST.json',dict(files=[dict(path=str(q),bytes=q.stat().st_size,sha256=C.X.digest(q)) for q in sorted(R.iterdir()) if q.is_file() and q.name!='FILE_MANIFEST.json']))
    print(json.dumps(dict(record=str(R/'EXPERIMENT_RECORD.md'),status=checkpoint['status'],physical_trials=0,optimization_trials=0)),flush=True)

def summarize():
    reference=load(M.SOURCE/'GRADED_REFERENCE.npz');v3=load(R/'V3_CONTROL.npz');default=load(R/'DYNAMIC_FROZEN.npz');zero=load(R/'ZERO_SFA.npz');fixtures=json.loads((R/'TRANSIENT_VALIDATION.json').read_text())
    counts={str(mode):sum(c['modes'][str(mode)]['passed'] for c in fixtures['cases'].values() if c['observable']) for mode in range(4)}
    sr,sp=cycle_diagnostics(reference['deaff'],default['deaff_graded_h']);diag={'source':sr}
    for label,z in (('v3',v3),('dynamic_default',default),('dynamic_zero_SFA',zero)):
        r,p=cycle_diagnostics(z['deaff_trace'],z['deaff_h_mean'],z['deaff_events']);r['h_registered_chemical_RMSE']=float(np.sqrt(np.mean((p-sp)**2))) if p is not None else None;diag[label]=r
    checkpoint=dict(machine='easteregg2',independent_transients=fixtures,transient_factor_pass_counts=counts,unchanged_v3_control_exact=True,frozen=json.loads((R/'DYNAMIC_FROZEN.json').read_text()),source_supported_zero_SFA=json.loads((R/'ZERO_SFA.json').read_text()),posthoc_cycle_diagnostics=diag,posthoc_scope='Availability/event-baseddiagnostics only; originalreadoutpeak/rawtimegates retained. No candidatepromoted.',physical_trials=0,optimization_trials=0,numerical_process_remaining=False)
    C.X.write(R/'CENTRAL_CHECKPOINT.json',checkpoint)
    plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,ax=plt.subplots(6,2,figsize=(6.5,9.3),sharex='col');handles=[]
    for j,key in enumerate(('standing','deaff')):
        t=np.arange(1,len(default[key+'_trace'])+1)*.002
        for label,z,color,style in (('V3',v3,'#d55e00','--'),('Matched+support',default,'#0072b2','-'),('Matched+support,SFA0',zero,'#009e73','-.')):
            ax[0,j].plot(t,z[key+'_events'][:,0],color=color,ls=style,lw=.7)
            for row,col in ((1,0),(2,1),(4,2),(5,20)):
                line,=ax[row,j].plot(t,z[key+'_trace'][:,col],color=color,ls=style,lw=1.,label=label)
                if row==1 and j==0:handles.append(line)
            ax[3,j].plot(t,z[key+'_h_mean'][:,0],color=color,ls=style,lw=1.)
        for row,col in ((1,0),(2,1),(4,2),(5,20)):
            line,=ax[row,j].plot(t,reference[key][:,col],color='#555555',ls=':',lw=1.2)
            if row==1 and j==0:handles.append(line)
        ax[3,j].plot(t,default[key+'_graded_h'][:,0],color='#555555',ls=':',lw=1.2)
        ax[0,j].set_title('Recorded standing\nNo fitting' if key=='standing' else 'Deafferented drive\nNo fitting',fontsize=10);ax[0,j].set_xlim((1.,3.) if key=='standing' else (2.,9.8));ax[5,j].set_xlabel('Neural time (s)')
        for a in ax[:,j]:a.grid(alpha=.18)
    for row,label in enumerate(['RG-E spikes\n(per2ms)','RG-E readout\n(mV)','RG-F readout\n(mV)','RG-E NaP h\n(availability)','PF S1 readout\n(mV)','Soleus MN\n(graded mV)']):ax[row,0].set_ylabel(label)
    fig.subplots_adjust(left=.15,right=.98,top=.92,bottom=.24,wspace=.32,hspace=.21);fig.legend(handles,['V3','Matched+support','Matched+support,SFA0','Graded source'],loc='lower center',bbox_to_anchor=(.5,.125),ncol=2,frameon=False)
    fig.text(.035,.025,'Spike/history observer passes independent transients; full source gates still fail.\nSFA0 removes only the added spike adaptation absent from the graded NaP source.\nActual events preserved; current/teacher voltage never supply an output carrier.\n256 central spiking cells; MN/reflex layers graded. No body trial or gain search.',fontsize=10)
    for suffix in ('pdf','svg','png'):fig.savefig(R/('DYNAMIC_OBSERVER_NETWORK.'+suffix),dpi=220,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    text='# Causal dynamic spike readout\n\nMachine: easteregg2. **Independent transient gates pass; original fullnetwork gates still fail. Zero fitting/body trials.**\n\n'
    text+='Actual AP events drive matched rate, leak/synaptic-conductance, and NaP-h histories with the same25ms kernel. Low-event inverse estimates shrink toward silence using a declared two-event support regularizer (effective history2tau). This Poisson-derived uncertainty model is not the actual deterministic spiking law. No current/teacherVm input enters the decoder, no fabricated spikes, and no original graded presynaptic voltage.\n\n'
    text+='All seven observable held-out step/ramp/current/channel fixtures pass: worst chemicalRMSE0.13310, h∞RMSE0.11618, inhibitory off lag54–80ms. Actual200µs clock and shared16-cell parameters unchanged. Four-mode factor pass counts (instantaneousinverse,event-matchedonly,supportonly,combined): '+str(counts)+'. All modes use identical actual spikes. Initializedzero-observed-event controls remain exactlysilent despite known driving current. Stable positive subthreshold voltage without spikes is explicitly unidentifiable, not voltage-fidelity success.\n\n'
    text+='Frozen fullnetwork: standing chemicalRMSE0.01879,herror0.04537; RG-E driven readout period2.236s versus source2.132s, but RG-F readout splitting and rawtimewaveformRMSE0.56685 fail. V3control exactly reproduces all stored traces/events/somas/h; complete-state/default/shadow/scalar checks pass.\n\n'
    text+='One source-supported ablation removes only the extra0.3nA/spike adaptation, because the original graded source RG has NaP inactivation but no separate spike-SFA mechanism. Actual adaptation remains exactly0. Standing h error9.53e−5, chemicalRMSE0.01050,and sampledMN-RMSE6.25e−5mV improve, but driven rawwaveformRMSE0.41458 and readout-peak gates still fail. This was a single mirroredmechanism contrast, not optimization. No candidate is promoted.\n\n'
    text+='Separate posthoc burst/cycle diagnostics (no acceptance substitution):\n\n'
    for label,r in diag.items():text+='- '+label+': '+str(r)+'\n'
    text+='\nPhysical/architecture scope:256 actual central spikingcells and778 gradedcells, includingMN/reflex pathways. Originalsyn6 source/weights, nonzeroPF→MN multiplier0.013344897908126326, geometry, and support unchanged. No body was integrated in this study.\n\nEvidence/full paths:\n\n'
    for name in ('TRANSIENT_PROTOCOL.json','TRANSIENT_VALIDATION.json','NETWORK_PROTOCOL.json','STATE_AUDIT.json','DYNAMIC_FROZEN.json','ZERO_SFA_PROTOCOL.json','ZERO_SFA_STATE_AUDIT.json','ZERO_SFA.json','CENTRAL_CHECKPOINT.json','DYNAMIC_OBSERVER_NETWORK.pdf'):
        text+='- `'+str(R/name)+'`\n'
    text+='\nNumericalprogram: `'+str(Path(D.__file__))+'`; use `D:\\Anaconda\\envs\\myo\\python.exe` with `--validate-transients --compare-network`. This is boundedvalidation/frozenreplay, not optimization. Source-SFAablation/analysisprogram: `'+str(Path(__file__))+'`; `--zero-sfa` repeats that singleablation; `--plot-only` uses saveddata only. V1/v2/v3 preserved. No ownedworker remains.\n'
    (R/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8');print(json.dumps(dict(physical_trials=0,optimization_trials=0,transient_factor_pass_counts=counts,posthoc_cycle_diagnostics=diag)),flush=True)
    if (R/'TIMING_FRESH.json').exists():finalize_record()

if __name__=='__main__':
    if '--zero-sfa' in sys.argv:zero_ablation()
    if '--timing' in sys.argv:timing_calibration()
    if '--fresh-timing' in sys.argv:fresh_timing()
    if '--audit-scorer' in sys.argv:audit_burst_scorer()
    if '--render-timing' in sys.argv:render_timing()
    if '--finalize' in sys.argv:finalize_record()
    if '--plot-only' in sys.argv:summarize()
