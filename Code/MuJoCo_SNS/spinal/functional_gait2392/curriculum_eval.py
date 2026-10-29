"""Fresh-process evaluation for the functional PF five-stage curriculum."""
import os,sys,json,time,hashlib
from pathlib import Path
H=Path(__file__).resolve().parent;S=H.parent
sys.path[:0]=[str(H),str(S),str(S/'w2l_cpg')];os.chdir(S)
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent);os.environ['AARL_NET']='w2lvar'
configpath=Path(sys.argv[1]).resolve();cfg=json.loads(configpath.read_text());out=configpath.with_suffix('')
if out.with_suffix('.result.json').exists():raise FileExistsError(str(out))
os.environ['AARL_NPZ']=str(out.with_suffix('.npz'))
import numpy as np
import params as P
import runner as R
import build_network as BN
import build_network_w2lvar as W
import optuna_walk_v10 as OW
import network as F
from preallocated_exact_eval import PreallocatedExactForward
from sparse_pairwise_eval import SparsePairwiseForward
from curriculum_metrics import assess
from frontal_metrics import assess_frontal
p=cfg['parameters'];stage=cfg['stage']
OW.load_fitted_baseline();OW.set_params(p)
P.TAU['rg_nap_h']=p['rg_nap_h']
# Full explicit load, including keys omitted or clobbered by the older curriculum.
for key in ('rg_mutual_inh','rg_weak_exc','ia_to_mn','ia_to_antagonist','ii_to_mn','ib_to_mn_inh',
            'contact_onset','contra_swing','ib_group_exc','ib_exc_to_mn','posture_to_mn','coact','vest_prop'):
    P.G[key]=float(p[key])
for key in W.GAINS:
    if key in p:W.GAINS[key]=float(p[key])
opts=dict(axes=True,biarticular=True,primitive=False,axis_gain=p['axis_gain'],trunk_gain=p['trunk_gain'],
          axis_feedback=0. if stage==1 else p['axis_feedback'],bi_gain=p['bi_gain'],bi_afferent=0. if stage==1 else p['bi_afferent'],
          bi_stance={f:p['stance_'+f] for f in F.BI_FAMILIES},
          mono_gains={g:p['mono_'+g] for g in ('hip_ext','hip_flex','knee_ext','knee_flex','ankle_pf','ankle_df')})
capture={}
def build(acts,dt,interleg):
    net=F.build(acts,dt=dt,interleg=interleg,options=opts)
    net.fast_io=cfg.get('fast_io',True)
    # Only initial finite perturbation. No continuing phase machine or motor override.
    net.net.add_input('RG_E_l');net.inputs.append('INIT_RG_E_l');net.compile(dt)
    backend=SparsePairwiseForward if cfg.get('backend','sparse_pairwise')=='sparse_pairwise' else PreallocatedExactForward
    net.compiled=backend(net.compiled)
    capture.update(neurons=len(net.compiled.V),inputs=len(net.inputs),gains=dict(P.G),variant_gains=dict(W.GAINS),tau=dict(P.TAU))
    step=net.step;ports={name:i for i,name in enumerate(net.inputs)};counter=0
    postports=np.array([i for name,i in ports.items() if name.startswith('POST_')],dtype=int)
    def advance(u):
        nonlocal counter
        u[ports['INIT_RG_E_l']]=p['kick'] if 2.<=counter*dt<2.5 and stage!=3 else 0.
        u[postports]*=p['preset_scale']
        value=step(u);counter+=1;return value
    net.step=advance
    return net
BN.build=build
args=['--eval','--drive',str(p['drive'])]
if stage in (1,2,4):args+=['--no-ground']
if stage==1:args+=['--no-afferents','--no-interleg']
if stage==3:args=['--stand-eval','8','--rig-scale',str(p['rig_scale'])]
if stage in (3,5):
    os.environ['AARL_PELVIS_TY']=str(p['pelvis_ty'])
    if stage==5:args+=['--rig-scale',str(p['rig_scale'])]
start=time.perf_counter();metrics=R.main(args)
assessment=assess(out.with_suffix('.npz'),metrics,stage)
assessment['pre_frontal_screen_accepted']=assessment['accepted']
assessment.update(assess_frontal(out.with_suffix('.npz')))
assessment['accepted']=bool(assessment['accepted'] and assessment['frontal_gate'])
result=dict(stage=stage,parameters=p,options=opts,metrics=metrics,assessment=assessment,
            compiled=capture,wall_s=time.perf_counter()-start,machine='easteregg2',
            physics='MuJoCo 2.3.7 / implicitfast / 2 ms; unchanged muscle forces; native neural outputs',
            backend=cfg.get('backend','sparse_pairwise'),
            source_hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (Path(__file__),H/'curriculum_metrics.py',H/'frontal_metrics.py',H/'network.py',S/'w2l_cpg/sparse_pairwise_eval.py')})
out.with_suffix('.result.json').write_text(json.dumps(result,indent=2))
print(json.dumps(dict(score=assessment['score'],accepted=assessment['accepted'],wall_s=result['wall_s'])),flush=True)
