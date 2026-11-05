"""Bounded robust postural tuning, fresh repeat, and unseen-seed verification."""
import argparse,json,hashlib,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'virtual_walker_20261031'))
sys.path.insert(0,r'D:\Users\Ben Bolen\ControllerPreparation_20261104\rustuna_runtime_py311')
from rustuna_transfer import Worker
import rustuna

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--trials',type=int,default=24);p.add_argument('--training-noise',type=float,default=.025);p.add_argument('--fit-neural-tau',action='store_true');a=p.parse_args()
    if not 1<=a.trials<=48:raise ValueError('Bounded trial count required')
    a.output.mkdir(parents=True,exist_ok=False)
    space=dict(kp=(35.,100.),kd=(5.,25.),reciprocal_scale=(.3,2.),Ia_scale=(.3,2.),II_scale=(.3,2.),Ib_scale=(.3,2.))
    if not 0<=a.training_noise<=.1:raise ValueError('Use a declared bounded noise amplitude')
    if a.fit_neural_tau:space['tau_scale']=(.3,1.5)
    manifest=dict(seed=61200,space=space,trials=a.trials,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in(HERE/'model.py',HERE/'run.py',HERE/'robust_adapter.py',Path(__file__))},training='Mirrored quiet tests and +/-16 Nm, 30/70 ms added-delay tests with 0.025 nA noise, seeds61200/61201',held_out='0/30/70 ms added delay, 0.05 nA noise, +/-16 Nm, +/-0.08rad, new seeds61210/61211/61212',objective='Mean cost + 0.25 worst cost; each cost=1000*failed_gate+RMS_deg+0.2*max_deg+0.2*settle_s+0.1*mean_squared_excitation',scope='PD task gains and mirrored neural gains only; no tendon, damping, actuator, or walker-model changes')
    manifest['training_noise_nA']=a.training_noise;manifest['training']='Mirrored quiet tests and +/-16 Nm, 30/70ms added delay, declared training_noise_nA, seeds61200/61201'
    (a.output/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    study=rustuna.create_study(direction='minimize',sampler=rustuna.samplers.TPESampler(seed=61200,n_startup_trials=8))
    for kp,kd in((80.,10.),(45.,10.),(55.,15.),(40.,7.),(65.,12.)):
        candidate=dict(kp=kp,kd=kd,reciprocal_scale=1.,Ia_scale=1.,II_scale=1.,Ib_scale=1.)
        if a.fit_neural_tau:candidate['tau_scale']=.5
        study.enqueue_trial(candidate)
    preflight=dict(op='preflight',training_noise_nA=a.training_noise)
    best=None;w=Worker(HERE/'robust_adapter.py',a.output);deadline=time.time()+1800
    try:
        ready=w.request(preflight,180.);(a.output/'PREFLIGHT.json').write_text(json.dumps(ready,indent=2))
        if not ready['ready_for_optimization']:raise ValueError('Nonfinite baseline')
        for i in range(a.trials):
            if time.time()>=deadline:break
            trial=study.ask();parameters={key:trial.suggest_float(key,*bounds)for key,bounds in space.items()}
            result=w.request(dict(op='evaluate',params=parameters),180.);study.tell(trial.number,result['value']);row=dict(trial=i,parameters=parameters,**result)
            with(a.output/'TRIALS.jsonl').open('a',encoding='utf-8')as f:f.write(json.dumps(row)+'\n')
            if best is None or row['value']<best['value']:best=row;(a.output/'BEST.json').write_text(json.dumps(best,indent=2))
            print(json.dumps(dict(trial=i,value=result['value'],passed=result['diagnostics']['passed'],total=result['diagnostics']['total'])),flush=True)
    finally:w.close()
    fresh=Worker(HERE/'robust_adapter.py',a.output)
    try:
        fresh.request(preflight,180.)
        repeat=fresh.request(dict(op='evaluate',params=best['parameters']),180.)
        held=fresh.request(dict(op='held_out',params=best['parameters']),180.)
    finally:fresh.close()
    matched=abs(repeat['value']-best['value'])<1e-10
    (a.output/'FRESH_REPEAT.json').write_text(json.dumps(dict(matched=matched,**repeat),indent=2));(a.output/'HELD_OUT.json').write_text(json.dumps(held,indent=2))
    (a.output/'FINAL.json').write_text(json.dumps(dict(accepted=bool(best['diagnostics']['accepted']and repeat['diagnostics']['accepted']and held['diagnostics']['accepted']and matched),parameters=best['parameters'],best=best,fresh_process_matched=matched,held_out_passed=held['diagnostics']['passed'],held_out_total=held['diagnostics']['total'],scope=manifest['scope']),indent=2))
if __name__=='__main__':main()
