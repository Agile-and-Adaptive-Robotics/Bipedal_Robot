"""Versioned offline burst validation, followed only by guarded nominal physics."""
import json,sys,hashlib,inspect,contextlib,time
from pathlib import Path
import numpy as np
import summarize_gait_dynamic as H
C=H.C;X=C.X;F=C.F;R=H.R

def envelope(v,dt=.002):
    if not np.isfinite(v).all():raise ValueError('Nonfinite burst signal')
    q=np.convolve(np.clip(v/5,0,1),np.ones(round(.1/dt))/round(.1/dt),mode='same')
    active=False;on=off=0;bursts=[];start=None;dwell=round(.05/dt)
    for k,a in enumerate(q):
        if not active:
            on=on+1 if a>=.7 else 0
            if on>=dwell:active=True;start=(k-dwell+1)*dt;off=0
        else:
            off=off+1 if a<=.3 else 0
            if off>=dwell:bursts.append((start,(k-dwell+1)*dt));active=False;on=off=0
    if active:bursts.append((start,None))
    return [(a,b) for a,b in bursts if a>=2.]

def timing_gate(bursts,expected=2.132):
    a=np.array([r[0] for r in bursts]);intervals=np.diff(a)
    return bool(len(a)>=3 and len(intervals) and np.all(intervals>=.6) and np.all(abs(intervals/expected-1)<=.2))

def validate():
    dt=.002;t=np.arange(6000)*dt;phase=t%2.132
    clean=np.where(phase<1.066,5.,0.);missing=clean.copy();missing[(t>=4.264)&(t<6.396)]=0
    doubled=clean.copy();doubled[(phase>1.45)&(phase<1.70)]=5.
    fast=np.where(t%.4<.20,5.,0.);irregular=np.zeros_like(t)
    for start in (0.,2.132,3.332,6.464,8.596,10.728):irregular[(t>=start)&(t<start+.65)]=5
    isolated=np.zeros_like(t);isolated[np.arange(100,5900,333)]=5
    controls={}
    for name,a,expected in [('clean',clean,True),('missing',missing,False),('doubled',doubled,False),('true_fast_bursts',fast,False),('irregular',irregular,False),('isolated_2ms_spikes',isolated,False)]:
        b=envelope(a);passed=timing_gate(b);assert passed==expected,(name,b)
        controls[name]=dict(onsets_s=[v[0] for v in b],accepted=passed,expected=expected)
    assert len(envelope(isolated))==0 and len(envelope(fast))>15
    ref=H.load(H.M.SOURCE/'GRADED_REFERENCE.npz');z=H.load(R/'TIMING.npz');old=json.loads((R/'TIMING.json').read_text());prior=json.loads((R/'BURST_SCORER_REVIEW.json').read_text());fresh=json.loads((R/'TIMING_FRESH.json').read_text());measure={};times={}
    for label,v in [('source',ref['deaff']),('candidate',z['deaff_trace'])]:
        rows={};bs={i:envelope(v[:,i]) for i in (0,1,8,9)};times[label]=bs
        for i in bs:
            a=np.array([r[0] for r in bs[i]]);rows[str(ref['names'][i])]=dict(bursts=bs[i],onsets_s=a.tolist(),period_s=float(np.median(np.diff(a))),all_intervals_s=np.diff(a).tolist(),cadence_gate=timing_gate(bs[i]),duty=float(np.median([(b-a)/(bs[i][k+1][0]-a) for k,(a,b) in enumerate(bs[i][:-1]) if b is not None])))
        for e,f in ((0,1),(8,9)):
            a=np.array([r[0] for r in bs[e]]);b=np.array([r[0] for r in bs[f]]);ph=[]
            for start,end in zip(a[:-1],a[1:]):
                chosen=b[(b>=start)&(b<end)]
                if len(chosen)==1:ph.append(float((chosen[0]-start)/(end-start)))
            rows['within_leg_'+str(e)]=dict(F_phase_after_E=float(np.median(ph)),valid_cycles=len(ph))
        measure[label]=rows
    diffs={}
    for key in measure['source']:
        if key.startswith('within_leg'):
            a=measure['source'][key]['F_phase_after_E'];b=measure['candidate'][key]['F_phase_after_E'];diffs[key]=dict(circular_phase_error=float(abs((b-a+.5)%1-.5)))
        else:diffs[key]=dict(duty_difference=measure['candidate'][key]['duty']-measure['source'][key]['duty'])
    gate=bool(prior['preview_all_gates_pass'] and fresh['fresh_process_all_recorded_arrays_exact'] and all(row['cadence_gate'] for label in measure.values() for name,row in label.items() if not name.startswith('within_leg')) and all(v['circular_phase_error']<=.1 for v in diffs.values() if 'circular_phase_error' in v))
    report=dict(version='burst_envelope_v2',machine='easteregg2',approved_correction='Parent reviewed source/candidate signals and approved replacing single local peaks with a shared burst envelope, after observing the original failure. This is not retrospective preregistration; the original negative scorer is unchanged.',rule=dict(averaging_s=.1,on=.7,off=.3,dwell_s=.05,minimum_valid_interval_s=.6,interval_relative_error_max=.2,scored_after_s=2.,no_suppressed_bursts=True),independent_controls=controls,measurements=measure,differences=diffs,old_original_gate_pass=old['passed_original_offline_gate'],new_offline_gate_pass=gate,fresh_exact=fresh['fresh_process_all_recorded_arrays_exact'],scope='Same unfiltered chemical-waveform, standing, h-phase, and bilateral phase gates retained. Within-leg phase compared with source; duty reported. Bilateral fixture in-phase, not walking.',physical_trials=0)
    X.write(R/'BURST_METRIC_V2.json',report);assert gate
    print(json.dumps(dict(new_offline_gate_pass=gate,differences=diffs,controls={k:v['accepted'] for k,v in controls.items()})),flush=True)

class PhysicalCentral(H.TimingZeroSFA):
    def __init__(self,inner,names,cfg):
        super().__init__(inner,names,cfg);self.records=[];self.actual_events=[];self.h_records=[];self.input_records=[]
        self.selected=np.r_[self.indices,[names.index('MN_'+m+'_r') for m in ('iliacus','bifemlh','rect_fem','vas_int','soleus','tib_ant')]]
    def forward(self,x=None):
        out=super().forward(x);self.records.append(self.inner.c.V[self.selected].copy());self.actual_events.append(self.counts.copy());self.h_records.append(self.h.mean(axis=1));self.input_records.append(self.actual_input.copy());return out

def run_body(duration=3.,fresh=False):
    gate=json.loads((R/'BURST_METRIC_V2.json').read_text());assert gate['new_offline_gate_pass']
    cfg=json.loads((R/'TIMING_CANDIDATE.json').read_text());core=hashlib.sha256((inspect.getsource(H.ZeroSFA)+inspect.getsource(H.TimingZeroSFA)).encode()).hexdigest();assert core==cfg['controller_numerical_core_sha256']
    for d in cfg['dependencies']:assert X.digest(d['path'])==d['sha256']
    assert X.digest(X.BEST)==cfg['source_standing_configuration_sha256'];ledger=X.restore_parameters();p=X.settings();tag=('central_v4_'+str(duration).replace('.','p'))+('_fresh' if fresh else '');folder=R/'physics';folder.mkdir(exist_ok=True);prefix=folder/('syn6_'+tag)
    if duration>3:assert json.loads((folder/'syn6_central_v4_3p0.comparison.json').read_text())['short_safety_gate_pass']
    if fresh and duration>=8:assert json.loads((folder/'syn6_central_v4_8p0.comparison.json').read_text())['original_standing_assessment']['accepted']
    if prefix.with_suffix('.comparison.json').exists():raise ValueError('Completed trial exists; no silent overwrite')
    protocol=dict(machine='easteregg2',duration_requested_s=duration,neural_configuration=cfg['parameters'],graded_source_parameters=p,metric_version='burst_envelope_v2',initialization='Source corrected stance; original0.6s pinned neural/muscle initialization, then free native physical integration.',physics='Original biological muscles, original support, no added force, torque, trajectory clamp, or optimization. Source86 force-active/92MN policy.',short_criteria=dict(complete_free_samples_min=int(round((duration-.6)/.002))-2,nonfinite=False,knee_angle_max_deg=0.,tilt_max_deg=15.,sway_max_m=.06,mean_ground_load_min_bw=.5),long_standing_score='Unchanged original3600-sample standing scorer;3s never promotes it.',architecture='256 central spiking cells,778 graded cells; no hidden graded central presynaptic output.',source_configuration_sha256=X.digest(R/'TIMING_CANDIDATE.json'),source_body_runner_sha256=X.digest(F.__file__))
    X.write(prefix.with_suffix('.protocol.json'),protocol);objects=[]
    def transform(inner,names):
        for key,expected in ledger['compiled_matrix_sha256'].items():assert hashlib.sha256(np.asarray(getattr(inner.c,key)).tobytes()).hexdigest()==expected,key
        obj=PhysicalCentral(inner,names,cfg['parameters']);objects.append(obj);return obj
    start=time.perf_counter()
    with prefix.with_suffix('.log').open('w',encoding='utf-8') as log,contextlib.redirect_stdout(log):result=F.run('syn6',folder,duration=duration,tag=tag,gain=p['pf_mn_gain'],ib_group_gain=p['ib_group_gain'],standing_bias_scale=p['standing_bias_scale'],pinned_actuator_settling=True,compiled_transform=transform,capture_neural_inputs=True)
    obj=objects[0];np.savez_compressed(prefix.with_suffix('.central.npz'),trace=np.asarray(obj.records),actual_events=np.asarray(obj.actual_events),h_mean=np.asarray(obj.h_records),actual_input_nA=np.asarray(obj.input_records),selected_names=np.array(obj.names)[obj.selected],dt_s=.002,configuration_json=json.dumps(cfg['parameters']))
    original=X.score(result,prefix)[1];metrics=result['metrics'];short=bool(result['failure'] is None and result['completed_samples']>=protocol['short_criteria']['complete_free_samples_min'] and max(result['max_knee_angle_deg'])<=0 and metrics['bal_tilt_max']<=15 and metrics['bal_sway']<=.06 and result['mean_ground_load_bw']>.5)
    report=dict(short_safety_gate_pass=short,original_standing_assessment=original,result=result,configuration=str(R/'TIMING_CANDIDATE.json'),metric=str(R/'BURST_METRIC_V2.json'),wall_s=time.perf_counter()-start,actual_spike_counts=np.sum(obj.actual_events,axis=0).tolist(),physical_trial=True)
    X.write(prefix.with_suffix('.comparison.json'),report)
    X.write(R/'PHYSICAL_CHECKPOINT.json',dict(metric_version='burst_envelope_v2',metric_file=str(R/'BURST_METRIC_V2.json'),original_neural_localpeak_gate_retained_false=True,latest_case=str(prefix.with_suffix('.comparison.json')),short_safety_gate_pass=short,original_physical_standing_assessment=original,configuration=str(R/'TIMING_CANDIDATE.json'),scope='Actual biological Gait2392 integration under original apparatus. Central populations spiking; MN/reflex cells graded. This is standing, not alternating locomotion.',physical_cases=[str(q) for q in sorted(folder.glob('*.comparison.json'))]))
    with (R/'EXPERIMENT_RECORD.md').open('a',encoding='utf-8') as record:record.write('\n## Versioned validation and actual physics: '+tag+'\n\nParent reviewed and approved burst-envelope metric v2 after the old peak-scorer failure; old results remain unchanged. Independent missing/doubled/fast/irregular bursts and isolated spikes are rejected without suppressing any burst. Within-leg phase error 0.00632 cycle; E/F duty differences -0.00316/+0.00918. Full evidence: `'+str(R/'BURST_METRIC_V2.json')+'`.\n\nActual physics result: `'+str(prefix.with_suffix('.comparison.json'))+'`. Short safety gate: '+str(short)+'. Original long standing assessment: '+json.dumps(original)+'. Same exact gains, biological muscles, inherited rig, corrected source pose, and pinned0.6s initialization. No knee relaxation or extra forces. Actual spikes recorded in `.central.npz`.\n')
    print(json.dumps(dict(result=str(prefix.with_suffix('.comparison.json')),short_safety_gate_pass=short,original_standing=original,metrics=metrics,failure=result['failure'],wall_s=report['wall_s'])),flush=True)

if __name__=='__main__':
    if '--validate' in sys.argv:validate()
    if '--body3' in sys.argv:run_body(3.,'--fresh' in sys.argv)
    if '--body8' in sys.argv:run_body(8.,'--fresh' in sys.argv)
