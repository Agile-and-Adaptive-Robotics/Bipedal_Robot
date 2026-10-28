"""Build result tables and numerical replay checks from saved evidence."""
import json,hashlib
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;OUT=HERE/'reports_20261028_family_audit';H=OUT/'human_trials';R=OUT/'robot_trials'
def compare(left,right):
    a,b=np.load(left,allow_pickle=True),np.load(right,allow_pickle=True)
    result={}
    for key in a.files:
        if key not in b.files:result[key]={'same':False};continue
        x,y=a[key],b[key]
        same=x.shape==y.shape and np.array_equal(x,y)
        delta=float(np.max(abs(x-y))) if x.shape==y.shape and x.size and x.dtype.kind in 'fiu' else None
        result[key]=dict(same=bool(same),max_difference=delta)
    return dict(bit_exact=all(r['same'] for r in result.values()) and set(a.files)==set(b.files),arrays=result)
rows=[];human=[];parity={};fast=[];robot=[]
for family in ('s3k','w2lvar','syn6'):
    for path in sorted(H.glob(family+'_*.json')):
        if path.stem.endswith(('_reference','_preallocated')):continue
        s=json.loads(path.read_text());m=s['metrics'];k=m.get('kine') or {}
        human.append(dict(family=family,case=s['case'],kine_score=m['kine_score'],cycles_r=k.get('n_cycles_r',0),cycles_l=k.get('n_cycles_l',0),fell=m['bal_fell'],nan=m['nan']))
    parity[family]={}
    for label in ('reference','preallocated'):
        path=H/(family+'_baseline_'+label+'.npz')
        if path.exists():parity[family][label]=compare(H/(family+'_baseline.npz'),path)
for path in sorted(OUT.glob('*_preallocated.json')):fast.append(json.loads(path.read_text()))
for path in sorted(R.glob('*authority*.json')):
    if path.stem.endswith(('_assessment','_abduction_assessment','_geometry_audit')):continue
    s=json.loads(path.read_text());robot.append(dict(run=path.stem,completed=s['failed'] is None,t_end=s['duration_simulated_s'],failed=s['failed']))
causality={}
control=R/'w2lvar_authority0p1_ramp.npz';ablated=R/'w2lvar_authority0p1_ramp_pf_deleted.npz'
if control.exists() and ablated.exists():
    a,b=np.load(control),np.load(ablated);n=min(len(a['time']),len(b['time']));time=a['time'][:n]
    before=time<15.;after=time>=15.
    checks={key:float(np.max(abs(a[key][:n][before]-b[key][:n][before]))) for key in ('qpos','qvel','activation','control','rg')}
    changes={key:float(np.max(abs(a[key][:n][after]-b[key][:n][after]))) for key in ('qpos','qvel','activation','control','rg')}
    early=(time>=15.)&(time<15.1)
    early_changes={key:float(np.max(abs(a[key][:n][early]-b[key][:n][early]))) for key in ('activation','control','rg')}
    causality=dict(intervention='PF-to-MN conductances deleted at 15 s',
        identical_before_intervention=all(v==0. for v in checks.values()),before_max_differences=checks,
        after_max_differences=changes,first_100ms_max_differences=early_changes,motor_effect=changes['control']>1e-6,
        scope='Local causal response within a 90% engineering MN clamp; not proof of a physiological gait controller')
negative_control={}
left=R/'w2lvar_authority0p0_ramp.npz';right=R/'w2lvar_authority0p0_ramp_pf_deleted.npz'
if left.exists() and right.exists():
    a,b=np.load(left),np.load(right);t=a['time'];before=t<15.;after=t>=15.;early=(t>=15)&(t<15.1)
    negative_control=dict(before_max_differences={k:float(np.max(abs(a[k][before]-b[k][before]))) for k in ('qpos','qvel','activation','control','rg')},
        after_max_differences={k:float(np.max(abs(a[k][after]-b[k][after]))) for k in ('qpos','qvel','activation','control','rg')},
        first_100ms_max_differences={k:float(np.max(abs(a[k][early]-b[k][early]))) for k in ('activation','control','rg')})
robot_replay={}
replay=R/'w2lvar_authority0p1_ramp_preallocated.npz'
if replay.exists():robot_replay=compare(control,replay)
result=dict(human=human,full_trajectory_parity=parity,backend_timings=fast,robot=robot,pf_causality=causality,pf_negative_control=negative_control,partial_neural_robot_replay=robot_replay)
(OUT/'RESULTS.json').write_text(json.dumps(result,indent=2))
lines=['### Human-model screen','', '| Family | Configuration | Kinematic score | R/L cycles | Fell |','|---|---|---:|---:|---|']
for s in human:lines.append('| {family} | {case} | {kine_score:.3f} | {cycles_r}/{cycles_l} | {fell} |'.format(**s))
lines+=['','Higher score is better. These short supported tests do not establish walking. Every targeted change scored below its own frozen baseline. No gait winner was replaced.','','### Full-trajectory numerical checks','']
for family,entries in parity.items():
    for label,s in entries.items():lines.append('- '+family+' versus '+label+': '+('all saved arrays bit-exact.' if s['bit_exact'] else '**Mismatch; reject this backend result.**'))
lines+=['','### Preallocated neural evaluator','', '| Family/model | Reference / preallocated speed ratio | Previous fast / preallocated ratio | State parity |','|---|---:|---:|---|']
for s in fast:lines.append('| %s/%s | %.2f | %.2f | %s |'%(s['family'],s['body'],s['reference_speedup'],s['previous_fast_speedup'],s['bit_exact']))
lines+=['','Ratios above one indicate faster preallocated stepping. These are warmed neural-step timings on a shared machine, not whole-trial speedups or controlled throughput benchmarks. For s3k, retain the earlier exact-fast evaluator if it is faster.','','### CAD neural-authority probes','', '| Trial | Simulated duration (s) | Outcome |','|---|---:|---|']
for s in robot:lines.append('| %s | %.3f | %s |'%(s['run'],s['t_end'],'Completed requested duration' if s['completed'] else s['failed']))
if causality:
    lines+=['','Delayed PF deletion: identical history before 15 s = **%s**; maximum subsequent control difference = **%.6g**. This is a causal effect of PF transmission inside a largely engineering-controlled system, not an accepted gait result.'%(causality['identical_before_intervention'],causality['after_max_differences']['control'])]
    lines+=['','Within the first 100 ms, the exposed-current case changes control by **%.6g**.'%causality['first_100ms_max_differences']['control']]
if negative_control:lines+=['The fully clamped negative control changes control by **%.6g** in that same first 100 ms (maximum over the full post-deletion interval: **%.6g**).'%(negative_control['first_100ms_max_differences']['control'],negative_control['after_max_differences']['control'])]
if robot_replay:lines+=['','Fresh preallocated replay of the 25 s partial-neural CAD run: **%s**.'%('all saved arrays bit-exact' if robot_replay['bit_exact'] else 'mismatch; do not accept numerical replay')]
report=OUT/'AUDIT_REPORT.md';text=report.read_text();marker='<!-- GENERATED_RESULTS -->';start=text.index(marker)+len(marker);end=text.index('\n## Recommended next changes',start)
report.write_text(text[:start]+'\n\n'+'\n'.join(lines)+'\n'+text[end:])
print(json.dumps(dict(human_trials=len(human),parity={f:{b:r['bit_exact'] for b,r in e.items()} for f,e in parity.items()},robot_trials=len(robot),pf_causality=causality,pf_negative_control=negative_control,partial_neural_robot_replay=robot_replay.get('bit_exact')),indent=2))
