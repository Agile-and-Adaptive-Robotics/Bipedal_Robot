"""Isolated staged Gait2392 evaluations of functional PF variants."""
import os,sys,json,argparse,time,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;S=HERE.parent;sys.path[:0]=[str(S),str(HERE),str(S/'w2l_cpg')];os.chdir(S)
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent);os.environ['AARL_NET']='w2lvar'
OUT=S/'reports_20261028_functional_pf';(OUT/'trials').mkdir(exist_ok=True)
p=argparse.ArgumentParser();p.add_argument('variant',choices=['baseline','axes','bi','combined','static_axes','primitive','s3k_gate']);p.add_argument('--stage',choices=['air_deaff','air_aff','stand','air_pattern','ground'],default='ground');p.add_argument('--axis-gain',type=float,default=.1);p.add_argument('--bi-gain',type=float,default=1.);p.add_argument('--feedback',type=float,default=1.);p.add_argument('--kick',type=float,default=0.);p.add_argument('--curriculum-seed',action='store_true');p.add_argument('--tau',type=float);p.add_argument('--primitive-gain',type=float,default=1.);p.add_argument('--label');a=p.parse_args()
import numpy as np
import params as P
import runner as R
import build_network as BN
from preallocated_exact_eval import PreallocatedExactForward
import network as F
if a.variant=='s3k_gate':
    os.environ['AARL_NET']='';import _curriculum as C
    base=dict(json.loads((S/'best_walk_params_v10.json').read_text())['multipliers']);base.update(json.loads((S/'reports_20260923/s3k_trial34_full_params.json').read_text())['params']);base.update(renshaw=.5,full_rules=0.);C.set_stage(4,base)
else:
    import _curriculum_w2lvar as C
    seed_stage={'air_deaff':1,'air_aff':2,'stand':3,'air_pattern':4,'ground':5}[a.stage] if a.curriculum_seed else 5
    base=dict(json.loads((S/'best_walk_params_v10.json').read_text())['multipliers']);base.update(json.loads((S/('curriculum_w2lvar_stage'+str(seed_stage)+'.json')).read_text())['params']);C.set_stage(seed_stage,base)
if a.tau is not None:P.TAU['rg_nap_h']=a.tau;base['rg_nap_h']=a.tau
opts=dict(axes=a.variant in ('axes','combined','static_axes','primitive'),biarticular=a.variant in ('bi','combined','primitive'),
          pose_dependent=a.variant!='static_axes',axis_gain=a.axis_gain,trunk_gain=a.axis_gain,bi_gain=a.bi_gain,axis_feedback=a.feedback,
          primitive=a.variant=='primitive',primitive_gain=a.primitive_gain)
label=a.label or a.variant+'_'+a.stage+'_a'+str(a.axis_gain).replace('.','p')+'_b'+str(a.bi_gain).replace('.','p')
path=OUT/'trials'/label;os.environ['AARL_NPZ']=str(path.with_suffix('.npz'))
sources=[Path(__file__),HERE/'network.py',HERE/'mechanics.py',S/'runner.py',S/'build_network_w2lvar.py',S/'w2l_cpg/preallocated_exact_eval.py']
source_hashes={str(p.relative_to(S)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
snapshot=OUT/'source_snapshots';snapshot.mkdir(exist_ok=True)
for source in sources:
    target=snapshot/(hashlib.sha256(source.read_bytes()).hexdigest()+source.suffix)
    if not target.exists():target.write_bytes(source.read_bytes())
capture={};original=BN.build;logs=[];netref=[]
def build(acts,dt,interleg):
    net=original(acts,dt=dt,interleg=interleg) if a.variant=='s3k_gate' else F.build(acts,dt=dt,interleg=interleg,options=opts)
    if a.kick:
        net.net.add_input('RG_E_l');net.inputs.append('INIT_RG_E_l');net.compile(dt)
    net.compiled=PreallocatedExactForward(net.compiled);netref.append(net)
    capture.update(neurons=int(net.compiled.V.size),inputs=len(net.inputs),gains=dict(P.G),tau=dict(P.TAU))
    if hasattr(net,'extra_pairs') or a.kick:
        pf_names=[cell for cells in getattr(net,'extra_pairs',{}).values() for cell in cells];capture['functional_pf_cells']=pf_names
        pf_idx=np.array([net.idx[name] for name in pf_names],dtype=int);step=net.step;counter=0
        def tracked(u):
            nonlocal counter
            if a.kick:u[net.inputs.index('INIT_RG_E_l')]=a.kick if 2.<=counter*dt<2.5 else 0.
            value=step(u)
            if counter%5==0:logs.append(np.r_[counter*dt,net.compiled.V[pf_idx]])
            counter+=1;return value
        net.step=tracked
    return net
BN.build=build
command=['--eval','--drive',str(base['drive'])]
if a.stage.startswith('air'):command+=['--no-ground']
if a.stage=='air_deaff':command+=['--no-afferents']
if a.stage=='stand':command+=['--stand-eval','8']
start=time.perf_counter();metrics=R.main(command)
result=dict(machine='easteregg2',model=str(R.MODEL),variant=a.variant,stage=a.stage,options=opts,
    initialization_pulse=dict(amplitude_nA=a.kick,start_s=2.,end_s=2.5),curriculum_seed=a.curriculum_seed,seed_parameters=base,
    metrics=metrics,compiled=capture,wall_s=time.perf_counter()-start,source_hashes=source_hashes,
    force_and_physics='Unchanged existing Gait2392 runner, contact, muscle model, and support; fixed 2 ms physics/neural step',
    controller='Model-aware PF recruitment and native SNS motor outputs; no CAD MN target-current clamp',
    stage_scope='Staged diagnostics seeded from the recorded w2lvar configuration; not a completed five-stage optimization')
if netref and hasattr(netref[0],'change_log'):result['routing_updates']=netref[0].change_log
if logs:np.savez_compressed(path.with_name(path.name+'_functional.npz'),trace=np.array(logs),pf_names=capture['functional_pf_cells'])
path.with_suffix('.json').write_text(json.dumps(result,indent=2))
print(json.dumps(dict(label=label,wall_s=result['wall_s'],metrics=metrics),indent=2),flush=True)
