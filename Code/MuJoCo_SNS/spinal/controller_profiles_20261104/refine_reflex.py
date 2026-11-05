"""Small Rustuna Ia/Ib refinement after frozen mirrored per-muscle peak fitting."""
import argparse,hashlib,json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'virtual_walker_20261031'))
sys.path.insert(0,r'D:\Users\Ben Bolen\ControllerPreparation_20261104\rustuna_runtime_py311')
from rustuna_transfer import Worker
import rustuna

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--trials',type=int,default=24);p.add_argument('--include-II',action='store_true');p.add_argument('--preserve-peaks',action='store_true');a=p.parse_args()
    if not 1<=a.trials<=64:raise ValueError('Bounded fit required')
    from pf_design_gate import require_pf_design_ready
    require_pf_design_ready(json.loads((a.fit/'MANIFEST.json').read_text())['configuration']['family'])
    a.output.mkdir(parents=True,exist_ok=False)
    preflight=dict(op='preflight',fit=str(a.fit),fit_sha256=hashlib.sha256((a.fit/'FINAL.json').read_bytes()).hexdigest())
    manifest=dict(**preflight,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in (Path(__file__),HERE/'motor_adapter.py',HERE/'replay_adapter.py',HERE/'profile.py')},space=dict(local_Ia=[.05,3.],local_Ib=[.05,3.]),seed=61104,scope='Offline Hunt D; peak solution frozen, Ia/Ib refined, not a body trial')
    manifest['source_sha256'][str(HERE/'pose_replay.py')]=hashlib.sha256((HERE/'pose_replay.py').read_bytes()).hexdigest()
    manifest['scope']='Ordered D waveform refinement; fixed timing and feedforward weights, mirrored local sensory gains; no body integration'
    if a.include_II:manifest['space']['local_II']=[.01,1.]
    manifest['preserve_peak_gate']=a.preserve_peaks
    def apply_peak_guard(result):
        if a.preserve_peaks:
            d=result['diagnostics'];d['waveform_gate_passed']=d['accepted'];d['peak_preservation_gate']=d['mn_peak_rmse']<=.25;d['accepted']=bool(d['accepted']and d['peak_preservation_gate'])
        return result
    (a.output/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    study=rustuna.create_study(direction='minimize',sampler=rustuna.samplers.TPESampler(seed=61104,n_startup_trials=8));study.enqueue_trial(dict(local_Ia=1.,local_Ib=1.,**({'local_II':1.}if a.include_II else{})));best=None
    if a.include_II:
        for ia,ib,ii in((1.,.1,.1),(1.,.1,.01),(.1,.1,.1),(.1,.1,.01),(1.,1.,.1)):
            study.enqueue_trial(dict(local_Ia=ia,local_Ib=ib,local_II=ii))
    w=Worker(HERE/'motor_adapter.py',a.output)
    try:
        ready=w.request(preflight,180.);(a.output/'PREFLIGHT.json').write_text(json.dumps(ready,indent=2))
        if not ready['ready_for_optimization']:raise ValueError('Fresh-process C check failed')
        deadline=time.time()+900.
        for i in range(a.trials):
            if time.time()>=deadline:break
            trial=study.ask();parameters={key:trial.suggest_float(key,*bounds,log=True)for key,bounds in manifest['space'].items()}
            result=apply_peak_guard(w.request(dict(op='evaluate',params=parameters),180.));study.tell(trial.number,result['value']);row=dict(trial=i,parameters=parameters,**result)
            with(a.output/'TRIALS.jsonl').open('a',encoding='utf-8')as f:f.write(json.dumps(row)+'\n')
            if best is None or(row['diagnostics']['accepted']and not best['diagnostics']['accepted'])or(row['diagnostics']['accepted']==best['diagnostics']['accepted']and row['value']<best['value']):best=row;(a.output/'BEST.json').write_text(json.dumps(best,indent=2))
            print(json.dumps(dict(trial=i,rmse=result['diagnostics'].get('mn_command_rmse'),accepted=result['diagnostics']['accepted'])),flush=True)
    finally:w.close()
    fresh=Worker(HERE/'motor_adapter.py',a.output)
    try:
        fresh.request(preflight,180.)
        repeat=apply_peak_guard(fresh.request(dict(op='capture',params=best['parameters'],path=str(a.output/'BEST_TRACE.npz')),180.))
    finally:fresh.close()
    matched=abs(repeat['value']-best['value'])<1e-10
    (a.output/'FRESH_REPEAT.json').write_text(json.dumps(dict(matched=matched,**repeat),indent=2))
    (a.output/'FINAL.json').write_text(json.dumps(dict(accepted=bool(best['diagnostics']['accepted']and repeat['diagnostics']['accepted']and matched),best=best,parameters=best['parameters'],fresh_process_matched=matched,scope=manifest['scope']),indent=2))
if __name__=='__main__':main()
