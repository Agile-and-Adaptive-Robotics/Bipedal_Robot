"""Fresh-process peak-fit validation and frozen-feedforward reflex refinement."""
import hashlib,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
from replay_adapter import Replay
from muscle_map import classify
from readout_keys import multiplier
INSTANCE=None;FIXED=None

def handle_request(request):
    global INSTANCE,FIXED
    if request['op']=='preflight':
        fit=Path(request['fit']);raw=(fit/'FINAL.json').read_bytes()
        if hashlib.sha256(raw).hexdigest()!=request['fit_sha256']:raise ValueError('Peak checkpoint changed')
        final=json.loads(raw);saved=json.loads((fit/'MANIFEST.json').read_text())
        if not final['accepted']:raise ValueError('Peak fit failed')
        cfg=saved['configuration'];FIXED=json.loads((Path(saved['source_timing_block'])/'FINAL.json').read_text())['parameters']
        if cfg.get('pose_recruitment_replay',False)or cfg.get('feedforward_only',False):
            from pose_replay import ContractReplay
            INSTANCE=ContractReplay(cfg)
        else:INSTANCE=Replay(cfg)
        r=INSTANCE
        for i,(pre,post)in enumerate(zip(r.c.pre,r.c.post)):
            source,target=r.names[pre],r.names[post]
            if target.startswith('MN_')and source.startswith('PF_')and not source.startswith('PF_IN'):
                r.base_g[i]*=multiplier(source,target,final)
        check=r.evaluate(FIXED,time.time()+180.)
        matched=abs(check['value']-final['best']['value'])<1e-10
        r.phase='D'
        if hasattr(r,'feedforward_only'):r.feedforward_only=False
        return dict(ok=True,ready_for_optimization=bool(check['diagnostics']['accepted']and matched),fresh_process_peak_check=check,peak_objective_abs_difference=abs(check['value']-final['best']['value']))
    if INSTANCE is None or request['op']not in ('evaluate','capture'):raise ValueError('Preflight required')
    result=INSTANCE.evaluate({**FIXED,**request['params']},time.time()+180.)
    if request['op']=='capture':
        path=Path(request['path']).resolve()
        if not path.is_relative_to(Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104')):raise ValueError('Unowned output')
        r=INSTANCE;np.savez_compressed(path,t=r.t,rg=r.rg_trace,pf=r.pf_trace,mn=r.output,desired=r.desired,acts=np.array(r.acts),rg_names=np.array([r.names[i]for i in r.rg]),pf_names=np.array([r.names[i]for i in r.pf]))
    return result
