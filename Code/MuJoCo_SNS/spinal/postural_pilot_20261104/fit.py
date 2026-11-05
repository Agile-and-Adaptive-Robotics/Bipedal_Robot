"""Rustuna fit of shared neural scales, followed by held-out physical tests."""
import argparse,json,os,sys,time,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'virtual_walker_20261031'))
sys.path.insert(0,r'D:\Users\Ben Bolen\ControllerPreparation_20261104\rustuna_runtime_py311')
from rustuna_transfer import Worker
import rustuna

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--trials',type=int,default=20);a=p.parse_args()
    if not 1<=a.trials<=48:raise ValueError('Bounded fit required')
    a.output.mkdir(parents=True,exist_ok=False)
    space=dict(reciprocal_scale=(.2,2.),Ia_scale=(.1,3.),II_scale=(.1,3.),Ib_scale=(.1,2.))
    manifest=dict(machine=os.environ.get('COMPUTERNAME'),seed=61104,space=space,trial_count=a.trials,training_conditions=[dict(initial_rad=.04,push_Nm=8.),dict(initial_rad=-.04,push_Nm=-8.)],task_gains=dict(kp=80.,kd=10.),source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in HERE.glob('*.py')},objective='mean(1000*failure+final_rms_deg+0.2*max_angle_deg+0.2*settling_s+0.1*mean_squared_excitation)',scope='Single-ankle bench neural-scale fit with shared antagonistic weights')
    (a.output/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    study=rustuna.create_study(direction='minimize',sampler=rustuna.samplers.TPESampler(seed=61104,n_startup_trials=8))
    study.enqueue_trial({k:1. for k in space});best=None
    w=Worker(HERE/'adapter.py',a.output)
    try:
        ready=w.request(dict(op='preflight'),180.);(a.output/'PREFLIGHT.json').write_text(json.dumps(ready,indent=2))
        if not ready['ready_for_optimization']:raise RuntimeError('Premotor baseline failed')
        for i in range(a.trials):
            t=study.ask();params={k:t.suggest_float(k,*v)for k,v in space.items()}
            result=w.request(dict(op='evaluate',params=params),180.)
            study.tell(t.number,result['value']);row=dict(trial=i,parameters=params,**result)
            with(a.output/'TRIALS.jsonl').open('a',encoding='utf-8')as f:f.write(json.dumps(row)+'\n')
            if best is None or row['value']<best['value']:best=row;(a.output/'BEST.json').write_text(json.dumps(best,indent=2))
            print(json.dumps(dict(trial=i,value=result['value'],accepted=result['diagnostics']['accepted'])),flush=True)
    finally:w.close()
    (a.output/'FIT_COMPLETE.json').write_text(json.dumps(best,indent=2))
if __name__=='__main__':main()
