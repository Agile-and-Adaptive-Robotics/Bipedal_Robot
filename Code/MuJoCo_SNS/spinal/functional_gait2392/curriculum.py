"""Resumable, process-isolated, five-stage physical-motion optimization.

All studies and full configurations are campaign-local. Each trial starts
fresh; fresh replay gates control promotion. Diagnostic continuation is explicit.
"""
import os,sys,json,time,argparse,subprocess,hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;S=H.parent
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import optuna
import numpy as np
parser=argparse.ArgumentParser();parser.add_argument('--trials',type=int,default=100);parser.add_argument('--workers',type=int,default=4)
parser.add_argument('--start',type=int,default=1);parser.add_argument('--end',type=int,default=5);parser.add_argument('--campaign',default='curriculum_20261028_v2')
parser.add_argument('--diagnostic-continue',action='store_true',help='Explicitly explore later environments after failed gates; does not certify progression.')
a=parser.parse_args();O=S/'reports_20261028_functional_pf'/a.campaign;O.mkdir(exist_ok=True)
(O/'trials').mkdir(exist_ok=True)
optuna.logging.set_verbosity(optuna.logging.WARNING)
base=json.loads((S/'best_walk_params_v10.json').read_text())['multipliers']
base.update(json.loads((S/'curriculum_w2lvar_stage5.json').read_text())['params'])
base.update(rg_weak_exc=0.,posture_to_mn=1.,coact=0.,vest_prop=0.,rig_scale=1.,pelvis_ty=.92,
            axis_gain=.1,trunk_gain=.1,axis_feedback=1.,bi_gain=1.,bi_afferent=.03,kick=2.,preset_scale=1.,
            w2l_rg_to_pf=2.4,w2l_pf_to_mn=2.,w2l_pf_recip=2.749,
            w2l_ia_central=.53,w2l_ii_central=.53,w2l_ib_central=.5,
            w2l_v3e_inE=1.,w2l_ini_inh=2.,w2l_v0d_inh=1.2,w2l_alpha=2.)
for f,v in dict(HAM=.25,RECT=.8,SAR=.1,GRAC=.2,TFL=.3,GAS=.9).items():base['stance_'+f]=v
for g in ('hip_ext','hip_flex','knee_ext','knee_flex','ankle_pf','ankle_df'):base['mono_'+g]=1.
RANGES1=dict(drive=(1.5,4.),rg_nap_h=(.22,.7),desc_e=(.7,2.),desc_f=(.7,2.),
    w2l_rg_to_pf=(1.5,4.),w2l_pf_to_mn=(.5,4.),w2l_pf_recip=(1.5,5.),
    bi_gain=(.3,2.),posture_to_mn=(0.,1.),preset_scale=(0.,1.),
    mono_hip_ext=(.2,2.),mono_hip_flex=(.3,2.5),mono_knee_ext=(.1,2.),mono_knee_flex=(.3,3.),
    mono_ankle_pf=(.1,1.5),mono_ankle_df=(.5,3.),
    stance_HAM=(0.,.7),stance_RECT=(.3,1.),stance_SAR=(0.,.5),stance_GRAC=(0.,.6),stance_TFL=(0.,.7),stance_GAS=(.3,1.))
RANGES2=dict(ia_to_mn=(0.,.9),ia_to_antagonist=(0.,.6),ii_to_mn=(0.,.6),ib_to_mn_inh=(0.,.5),
    w2l_ia_central=(0.,.55),w2l_ii_central=(0.,.55),w2l_ib_central=(0.,.5),
    w2l_v3e_inE=(0.,1.2),w2l_ini_inh=(.2,3.),w2l_v0d_inh=(.1,2.),w2l_alpha=(.1,2.5),
    axis_gain=(.01,.6),trunk_gain=(.03,.8),axis_feedback=(.1,2.),bi_afferent=(0.,.15))
RANGES3=dict(rig_scale=(.2,1.),axis_gain=(.01,1.),trunk_gain=(.1,2.),axis_feedback=(.3,3.),
    posture_to_mn=(.2,2.),preset_scale=(.3,1.5),coact=(0.,.15),vest_prop=(0.,.8),pelvis_ty=(.84,.95))
RANGES4={**RANGES1,
         'axis_gain':(.01,.6),'axis_feedback':(.1,2.)}
RANGES5=dict(contact_onset=(0.,.5),contra_swing=(0.,1.),ib_group_exc=(0.,.8),ib_exc_to_mn=(0.,.9),
             w2l_ib_central=(0.,.5),w2l_v3e_inE=(0.,1.),pelvis_ty=(.84,.95),rig_scale=(.5,1.),
             posture_to_mn=(0.,1.5),preset_scale=(0.,1.),axis_gain=(.02,.8),axis_feedback=(.2,2.5))
RANGES={1:RANGES1,2:RANGES2,3:RANGES3,4:RANGES4,5:RANGES5}
sources=[*H.glob('*.py'),S/'runner.py',S/'build_network_w2lvar.py',S/'w2l_cpg/preallocated_exact_eval.py',S/'w2l_cpg/sparse_pairwise_eval.py',S/'reports_20261028_functional_pf/routing_audit.json']
manifest={str(p.relative_to(S)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
mp=O/'source_manifest.json'
if mp.exists():
    if json.loads(mp.read_text())!=manifest:raise RuntimeError('Campaign sources changed. Use a new campaign name.')
else:
    mp.write_text(json.dumps(manifest,indent=2));(O/'sources').mkdir(exist_ok=True)
    for p in sources:(O/'sources'/(hashlib.sha256(p.read_bytes()).hexdigest()+p.suffix)).write_bytes(p.read_bytes())

def evaluate(stage,p,label):
    path=O/'trials'/label;config=path.with_suffix('.json');result=path.with_suffix('.result.json')
    if result.exists():return json.loads(result.read_text())
    config.write_text(json.dumps(dict(stage=stage,parameters=p),indent=2))
    with path.with_suffix('.log').open('w') as log:
        run=subprocess.run([sys.executable,str(H/'curriculum_eval.py'),str(config)],stdout=log,stderr=subprocess.STDOUT)
    if run.returncode or not result.exists():
        return dict(assessment=dict(score=-2000.,accepted=False,error='Evaluation failed; inspect '+str(path.with_suffix('.log'))))
    return json.loads(result.read_text())

def report():
    stages=[]
    for n in range(1,6):
        p=O/f'stage{n}_best.json'
        if p.exists():stages.append(json.loads(p.read_text()))
    text=['# Functional PF five-stage curriculum — Session 1','',
          'Machine: easteregg2. Model: articulated human Gait2392. Native SNS outputs; preallocated Numba/NumPy backend.','',
          'Fresh replay acceptance controls stage promotion. Explicit diagnostic continuation is not graduation. Failed gates are retained.','',
          '| Stage | Completed trials | Best objective | Fresh replay accepted |','|---|---:|---:|---|']
    for s in stages:text.append(f"| {s['stage']} | {s['completed_trials']} | {s['score']:.3f} | {s['replay']['assessment']['accepted']} |")
    text+=['','Acceptance uses actual bilateral knee/hip excursions and cycle timing, not RG burst counts. Standing has separate tilt/sway/load gates. Ground adds contact-cycle and duty gates. Forward progression and reduced support are not assumed.','',
           'The stage-one trials remove muscle afferents, added joint feedback, and interleg wiring. Later stages enable them explicitly. The finite left RG initialization pulse ends at 2.5 s.','',
           'Full parameters, metrics, raw trajectories, and logs are stored beside this report. No production winner or database is replaced.']
    (O/'REPORT.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    # Update only Session 1's owned block through its byte-preserving helper.
    own=S/'reports_20261027_robot_control/SESSION_1_UPDATE_20261027.md'
    value=own.read_text(encoding='utf-8');start='<!-- FUNCTIONAL_PF_CURRICULUM_START -->';end='<!-- FUNCTIONAL_PF_CURRICULUM_END -->'
    block=start+'\n### Five-stage functional PF tuning campaign\n\nEvidence: `Code\\MuJoCo_SNS\\spinal\\reports_20261028_functional_pf\\'+a.campaign+'\\REPORT.md`.\n\n'+'\n'.join(text[6:6+len(stages)+2])+'\n\nThis table records optimization and fresh-process gate outcomes, not a blanket walking claim.\n'+end
    if start in value:value=value[:value.index(start)]+block+value[value.index(end)+len(end):]
    else:value=value.replace('### Scope and next steps',block+'\n\n### Scope and next steps')
    own.write_text(value,encoding='utf-8')
    helper=S/'reports_20261027_robot_control/update_session1_report.ps1'
    subprocess.run(['pwsh','-NoProfile','-File',str(helper)],stdout=subprocess.DEVNULL,check=True)

with ThreadPoolExecutor(max_workers=a.workers) as pool:
    for stage in range(a.start,a.end+1):
        previous=O/f'stage{stage-1}_best.json'
        if stage>1 and not a.diagnostic_continue:
            if not previous.exists() or not json.loads(previous.read_text())['replay']['assessment']['accepted']:
                raise RuntimeError('Prior stage has not passed its fresh replay. Stop for module training; diagnostic continuation requires --diagnostic-continue.')
            if stage==5 and not json.loads(previous.read_text())['replay']['assessment'].get('frontal_gate',False):
                raise RuntimeError('Ground-stage promotion requires an explicit frontal-plane gate, not the old sagittal-only screen.')
        if previous.exists():base.update(json.loads(previous.read_text())['parameters'])
        # Standing tone/support tuning must not overwrite the successful air drive.
        ranges=RANGES[stage]
        study=optuna.create_study(direction='maximize',storage='sqlite:///'+str(O/'curriculum.db').replace('\\','/'),
            study_name=f'functional_s{stage}',load_if_exists=True,sampler=optuna.samplers.TPESampler(seed=2810+stage,n_startup_trials=16,multivariate=True))
        if not study.trials:
            study.enqueue_trial({k:float(np.clip(base[k],*v)) for k,v in ranges.items()})
            pilot=S/'reports_20261028_functional_pf/curriculum_20261028_pilot/stage1_best.json'
            if stage==1 and pilot.exists():
                candidate=json.loads(pilot.read_text())['parameters']
                for tau in (candidate['rg_nap_h'],.65):
                    candidate['rg_nap_h']=tau
                    study.enqueue_trial({k:float(np.clip(candidate[k],*v)) for k,v in ranges.items()})
        best=study.best_trial if any(t.state==optuna.trial.TrialState.COMPLETE for t in study.trials) else None
        p={**base,**best.params} if best is not None else dict(base)
        state={}
        while len([t for t in study.trials if t.state==optuna.trial.TrialState.COMPLETE])<a.trials:
            remaining=a.trials-len([t for t in study.trials if t.state==optuna.trial.TrialState.COMPLETE])
            batch=[]
            for _ in range(min(a.workers,remaining)):
                trial=study.ask();p=dict(base)
                for k,(lo,hi) in ranges.items():p[k]=trial.suggest_float(k,lo,hi)
                batch.append((trial,p,pool.submit(evaluate,stage,p,f's{stage}_t{trial.number:04}')))
            for trial,p,job in batch:
                result=job.result();assessment=result['assessment'];trial.set_user_attr('physical_gate',assessment)
                study.tell(trial,float(assessment['score']))
                print(json.dumps(dict(stage=stage,trial=trial.number,score=assessment['score'],accepted=assessment['accepted'],
                                     knees=assessment.get('knee_cycles'),phase_error=assessment.get('cycle_phase_error'))),flush=True)
            best=study.best_trial;p={**base,**best.params}
            state=dict(stage=stage,trial=best.number,score=best.value,parameters=p,
                       completed_trials=len([t for t in study.trials if t.state==optuna.trial.TrialState.COMPLETE]),
                       replay=dict(assessment=dict(accepted=False,pending=True)))
            (O/f'stage{stage}_progress.json').write_text(json.dumps(state,indent=2))
        feasible=[t for t in study.trials if t.state==optuna.trial.TrialState.COMPLETE and
                  t.user_attrs.get('physical_gate',{}).get('accepted',False)]
        if feasible:
            best=max(feasible,key=lambda t:t.value);p={**base,**best.params}
        replay=evaluate(stage,p,f's{stage}_best_t{best.number:04}_fresh')
        state=dict(stage=stage,trial=best.number,score=best.value,parameters=p,
                   completed_trials=len([t for t in study.trials if t.state==optuna.trial.TrialState.COMPLETE]))
        state['replay']=replay
        state['diagnostic_continuation']=a.diagnostic_continue
        (O/f'stage{stage}_best.json').write_text(json.dumps(state,indent=2));base=p;report()
        print('STAGE FINISHED '+str(stage)+' accepted='+str(replay['assessment']['accepted']),flush=True)
        if not replay['assessment']['accepted'] and not a.diagnostic_continue:
            print('STOPPED: failed stage gate. Later stages were not started.',flush=True)
            break
