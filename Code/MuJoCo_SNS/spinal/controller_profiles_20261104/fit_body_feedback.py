"""Bounded actual-body Type-II/expression fit; fresh-process curriculum gates."""
import argparse,hashlib,json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'virtual_walker_20261031'))
sys.path.insert(0,r'D:\Users\Ben Bolen\ControllerPreparation_20261104\rustuna_runtime_py311')
from rustuna_transfer import Worker
import rustuna

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--passed-air',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--trials',type=int,default=12);a=p.parse_args()
    if not 1<=a.trials<=24:raise ValueError('Bounded body fit required')
    from pf_design_gate import require_pf_design_ready
    source_C=Path(json.loads((a.fit/'MANIFEST.json').read_text())['fit'])
    require_pf_design_ready(json.loads((source_C/'MANIFEST.json').read_text())['configuration']['family'])
    a.output.mkdir(parents=True,exist_ok=False)
    rawmodel=HERE.parents[3]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
    files=[Path(__file__),HERE/'body_fit_adapter.py',HERE/'body_probe.py',HERE/'dominant_cycles.py',HERE/'profile.py',HERE.parent/'runner.py',HERE.parent/'hunt_training/physical_eval.py',rawmodel]
    manifest=dict(reflex_fit=str(a.fit),passed_air=str(a.passed_air),output=str(a.output),seed=61204,space=dict(expression=[.1,.3],local_II_scale=[.05,.35]),trials=a.trials,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in files},scope='Actual-body air-stage feedback calibration; unchanged mechanics and matched left/right gains',target_period_s=1.233,target_period_tolerance_fraction=.15,objective='5*relative_period_error+5*phase_error+0.02*sum(abs(knee_ROM-35)); failed gates add100, ROM deficits add2/degree, fast-spectrum rejection adds500+50*Hz, guarded/incomplete run10000')
    from kinematic_response import reference,self_check
    ref=reference();self_check()
    manifest.update(objective_profile='stage_specific_kinematics_v2',reference=ref['identity'],reference_source_sha256=ref['source_sha256'],objective='Air: 5*bilateral_phase_error + period_cv + normalized_cycle_variability; existing motion/safety failures add1000; invalid cycles10000. Ground IK amplitude/phase/cadence remain diagnostic in air.',target_period_s=None,target_period_tolerance_fraction=None)
    manifest['source_sha256'].update(ref['source_sha256'])
    for f in(HERE/'kinematic_response.py',HERE.parent/'kine_ref.py'):
        manifest['source_sha256'][str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
    (a.output/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    study=rustuna.create_study(direction='minimize',sampler=rustuna.samplers.TPESampler(seed=61204,n_startup_trials=6))
    for expression,gain in((.15,.1),(.2,.1),(.25,.1),(.15,.2),(.2,.2),(.25,.2)):study.enqueue_trial(dict(expression=expression,local_II_scale=gain))
    best=None;w=Worker(HERE/'body_fit_adapter.py',a.output);deadline=time.time()+1800
    try:
        ready=w.request(dict(op='preflight',manifest=manifest),180.);(a.output/'PREFLIGHT.json').write_text(json.dumps(ready,indent=2))
        if not ready['ready_for_optimization']:raise ValueError('Prerequisite gate failed')
        for i in range(a.trials):
            if time.time()>=deadline:break
            trial=study.ask();parameters={key:trial.suggest_float(key,*bounds)for key,bounds in manifest['space'].items()};result=w.request(dict(op='evaluate',trial=i,params=parameters),240.);study.tell(trial.number,result['value']);row=dict(trial=i,parameters=parameters,**result)
            with(a.output/'TRIALS.jsonl').open('a',encoding='utf-8')as f:f.write(json.dumps(row)+'\n')
            if best is None or(result['diagnostics']['accepted']and not best['diagnostics']['accepted'])or(result['diagnostics']['accepted']==best['diagnostics']['accepted']and result['value']<best['value']):best=row;(a.output/'BEST.json').write_text(json.dumps(best,indent=2))
            print(json.dumps(dict(trial=i,value=result['value'],accepted=result['diagnostics']['accepted'],parameters=parameters)),flush=True)
    finally:w.close()
    accepted=False
    if best['diagnostics']['accepted']:
        with(a.output/'fresh_validation.log').open('a',encoding='utf-8')as log:subprocess.run([r'D:\Anaconda\envs\myo\python.exe',str(HERE/'air_validation.py'),'--fit',str(a.fit),'--output',str(a.output/'fresh_validation'),'--expression',str(best['parameters']['expression']),'--local-ii-scale',str(best['parameters']['local_II_scale'])],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=900,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        rows=json.loads((a.output/'fresh_validation/COMPLETE.json').read_text());accepted=len(rows)==2 and all(r['accepted']for r in rows)
    (a.output/'FINAL.json').write_text(json.dumps(dict(accepted=bool(accepted),best=best,scope=manifest['scope']),indent=2))
if __name__=='__main__':main()
