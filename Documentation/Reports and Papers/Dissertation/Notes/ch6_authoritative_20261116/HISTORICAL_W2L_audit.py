"""Independent read-only scoring of retained historical W2L recordings."""
from pathlib import Path
import json
import hashlib
import numpy as np

OUT=Path(__file__).resolve().parent
ROOT=Path(r'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal')
R=ROOT/'reports_20261027_codex'
sources={}
def bind(p):
    p=Path(p);sources[str(p)]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size};return p
def js(p):return json.loads(bind(p).read_text(encoding='utf-8-sig'))
def nz(p):return np.load(bind(p),allow_pickle=True)

def primary(f,t):
    on=f>20.;rise=np.flatnonzero(on[1:]&~on[:-1])+1
    indices=[]
    for k in rise:
        if not indices or t[k]-t[indices[-1]]>.45:indices.append(int(k))
    return indices

def score(z):
    keep=z['t']>=4.5;t=z['t'][keep];heel=z['heelF'][keep];toe=z['toeF'][keep];com=z['com'][keep]
    L=primary(heel[:,0],t);Rr=primary(heel[:,1],t)
    phases=[]
    for k in range(len(Rr)-1):
        for j in L:
            if t[Rr[k]]<=t[j]<=t[Rr[k+1]]:phases.append(float((t[j]-t[Rr[k]])/(t[Rr[k+1]]-t[Rr[k]])))
    finite=bool(z['finite']);fall=z['fall_t'].item() is not None
    phase=float(np.median(phases))if phases else None
    fwd=float(com[-1,0]-com[0,0]);v=t>=8
    return {'cycL':len(L),'cycR':len(Rr),'phase_med':phase,'phase_n':len(phases),'fwd':fwd,'slope':float(np.polyfit(t[v],com[v,0],1)[0]),'harness_percent_legacy411N':float(100*z['rigFz'][keep].mean()/411.),'harness_mean_N':float(z['rigFz'][keep].mean()),'dutycL':float(((heel[:,0]>5)|(toe[:,0]>5)).mean()),'dutycR':float(((heel[:,1]>5)|(toe[:,1]>5)).mean()),'finite':finite,'fell':fall,'survive':bool(finite and not fall and min(len(L),len(Rr))>=6 and phase is not None and .35<=phase<=.65 and fwd>.05),'contact_window_s':[float(t[0]),float(t[-1])],'speed_window_s':[float(t[v][0]),float(t[v][-1])],'heel_threshold_N':20,'primary_merge_s':.45,'duty_threshold_N':5,'phase_definition':'left primary strike phase within consecutive right primary strikes; median across retained events','primary_times_L':t[L].tolist(),'primary_times_R':t[Rr].tolist()}

summary={'scorer':'Reimplemented exact retained _primary/_phase and metric formulas; no evaluator import or simulation','runs':{},'repeat_comparisons':{},'perturbation_conditions':{},'sources':sources}
for p in [ROOT/'reports_20261016/w2l_wean/_w2l_wean_eval.py',R/'recover_w2l.py',R/'perturb_w2l.py',R/'full_basin_w2l.py',R/'verification.json',R/'w2l_verified_config.json',ROOT/'w2l_mujoco/test_w2l_ground.py',ROOT/'w2l_mujoco/w2l_ground.xml']:
    bind(p)
for name in ['trial458_d30','trial458_d30_ws0','panel_control','panel_del_all_feedback','panel_del_Renshaw','panel_del_commissural']:
    z=nz(R/f'{name}.npz');record=js(R/f'{name}.json');s=score(z)
    s['saved_metrics']=record['metrics'];s['job']=record['job'];s['saved_metric_max_difference']=max(abs(s[k]-record['metrics'][k]) for k in ['cycL','cycR','phase_med','fwd','slope','dutycL','dutycR'])
    summary['runs'][name]=s
base=nz(R/'trial458_d30.npz')
for name in ['repeat_0','repeat_1','repeat_2','canonical_repeat']:
    z=nz(R/f'{name}.npz');keys=[k for k in base.files if np.issubdtype(base[k].dtype,np.number)]
    summary['repeat_comparisons'][name]={'numeric_arrays':len(keys),'arrays_equal':{k:bool(np.array_equal(base[k],z[k],equal_nan=True))for k in keys}}
results=js(R/'full_basin_results.json');summary['full_basin']={'count':len(results),'saved_pass_count':sum(bool(r['metrics']['survive'])for r in results),'recomputed':{}}
for row in results:
    name=row['job']['name'];s=score(nz(R/f'{name}.npz'));summary['full_basin']['recomputed'][name]={k:s[k]for k in ['survive','cycL','cycR','phase_med','fwd','finite','fell']}
pert=js(R/'perturbation_results.json')
for row in pert:
    name=row['job']['name'];s=score(nz(R/f'{name}.npz'));summary['perturbation_conditions'][name]={k:s[k]for k in ['survive','cycL','cycR','phase_med','fwd','finite','fell']}
control=js(R/'panel_control.json')['job']['knobs'];removed=js(R/'panel_del_all_feedback.json')['job']['knobs']
summary['all_feedback_config_changes']={k:{'control':control[k],'removed':removed[k]}for k in control if control[k]!=removed[k]}
push=nz(R/'panel_push_fwd.npz');intact=nz(R/'panel_control.npz')
summary['forward_push_equals_control']={k:bool(np.array_equal(push[k],intact[k],equal_nan=True))for k in intact.files if np.issubdtype(intact[k].dtype,np.number)}
summary['forward_push_source_definition']={'runner_defaults':{'pushN':40.,'pushT':8.,'pushAxis':'+x'},'pulse_duration_s':.15,'control_job_omits_push_keys':not any(k in control for k in ['pushN','pushT','pushAxis']),'limitation':'No separate recorded horizontal applied-force array; inference uses retained runner/defaults and exact trajectory identity.'}
(OUT/'HISTORICAL_W2L_RESULTS.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps({'runs':{k:{x:v[x]for x in ['cycL','cycR','phase_med','fwd','slope','harness_percent_legacy411N','harness_mean_N','saved_metric_max_difference']}for k,v in summary['runs'].items()},'repeat_all_numeric_equal':{k:all(v['arrays_equal'].values())for k,v in summary['repeat_comparisons'].items()},'basin_pass':sum(v['survive']for v in summary['full_basin']['recomputed'].values()),'perturbations':summary['perturbation_conditions']},indent=2))
