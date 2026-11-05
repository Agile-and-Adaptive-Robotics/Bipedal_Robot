"""Mirrored per-muscle calibration after a coarse six-group fit fails.

Use a conductance transfer estimate for simultaneous bounded gain updates,
and score each proposal by integrating the actual compiled circuit. No target
voltage is clamped. Original target files and saved controllers are unchanged.
"""
import argparse,json,os,sys,time,hashlib
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
from replay_adapter import Replay
from muscle_map import classify

def main():
    p=argparse.ArgumentParser();p.add_argument('--timing-block',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--iterations',type=int,default=16);p.add_argument('--pose-aware',action='store_true');p.add_argument('--feedforward-only',action='store_true');a=p.parse_args()
    if not 1<=a.iterations<=40:raise ValueError('Bounded fit required')
    a.output.mkdir(parents=True,exist_ok=False)
    prior=json.loads((a.timing_block/'FINAL.json').read_text())
    if not prior['accepted']:raise ValueError('Timing block must pass')
    cfg=json.loads((a.timing_block/'MANIFEST.json').read_text());cfg['phase']='C';fixed=prior['parameters']
    from pf_design_gate import require_pf_design_ready
    require_pf_design_ready(cfg['family'])
    replay=Replay
    if a.pose_aware or a.feedforward_only:
        from pose_replay import ContractReplay
        replay=ContractReplay;cfg['pose_recruitment_replay']=a.pose_aware;cfg['feedforward_only']=a.feedforward_only
    r=replay(cfg);original=r.base_g.copy();keep=r.t>=2*r.T
    bases=sorted({classify(act).base for act,selected in zip(r.acts,r.selected)if selected})
    act_ids={base:np.array([i for i,act in enumerate(r.acts)if classify(act).base==base],dtype=int)for base in bases}
    edge_ids={base:np.array([i for i,(pre,post)in enumerate(zip(r.c.pre,r.c.post))if r.names[post].startswith('MN_')and classify(r.names[post][3:]).base==base and r.names[pre].startswith('PF_')and not r.names[pre].startswith('PF_IN')],dtype=int)for base in bases}
    # Report incompatible target amplitudes instead of silently changing data.
    desired={base:float(np.mean(r.desired[keep][:,ids].max(axis=0)))for base,ids in act_ids.items()}
    target_audit={base:dict(actuators=[r.acts[i]for i in ids],target_peaks=r.desired[keep][:,ids].max(axis=0).tolist(),target_peak_spread=float(np.ptp(r.desired[keep][:,ids].max(axis=0))))for base,ids in act_ids.items()}
    (a.output/'TARGET_AUDIT.json').write_text(json.dumps(target_audit,indent=2))
    manifest=dict(source_timing_block=str(a.timing_block),configuration=cfg,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in (Path(__file__),Path(__file__).parent/'replay_adapter.py',Path(__file__).parent/'profile.py')},bounds=[.02,20.],update='Log-damped inverse DC transfer ratio, proposal scored with the actual SNS dynamics',scope='Shared per-muscle PF->MN gain refit, not independent left/right weights or full PF-basis refit')
    if a.pose_aware or a.feedforward_only:
        helper=Path(__file__).parent/'pose_replay.py';manifest['source_sha256'][str(helper)]=hashlib.sha256(helper.read_bytes()).hexdigest()
    (a.output/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    gains={base:1. for base in bases};best=None;rows=[]
    def install(values):
        r.base_g[:]=original
        for base,ids in edge_ids.items():r.base_g[ids]=original[ids]*values[base]
    for step in range(a.iterations):
        install(gains);result=r.evaluate(fixed,time.time()+180.)
        if not result['diagnostics']['finite']:
            if best is None:raise RuntimeError('Baseline is nonfinite')
            gains={base:np.sqrt(gains[base]*best['gains'][base])for base in bases}
            rows.append(dict(iteration=step,gains=gains,**result));continue
        observed={base:float(np.mean(r.output[keep][:,ids].max(axis=0)))for base,ids in act_ids.items()}
        row=dict(iteration=step,gains=gains.copy(),**result,observed=observed)
        rows.append(row)
        if best is None or result['value']<best['value']:
            best=row;(a.output/'BEST.json').write_text(json.dumps(best,indent=2))
            np.savez_compressed(a.output/'BEST_TRACE.npz',t=r.t,rg=r.rg_trace,pf=r.pf_trace,mn=r.output,desired=r.desired,acts=np.array(r.acts),rg_names=np.array([r.names[i]for i in r.rg]),pf_names=np.array([r.names[i]for i in r.pf]))
        (a.output/'ITERATIONS.json').write_text(json.dumps(rows,indent=2))
        print(json.dumps(dict(iteration=step,value=result['value'],peak_rmse=result['diagnostics']['mn_peak_rmse'],accepted=result['diagnostics']['accepted'])),flush=True)
        for base in bases:
            target=np.clip(desired[base],.005,1.);actual=np.clip(observed[base],.005,1.)
            ratio=(target/(1.6-target))/(actual/(1.6-actual))
            gains[base]=float(np.clip(gains[base]*np.clip(ratio,.3,3.)**.6,.02,20.))
    # Independent reconstruction, not just resetting the fitted in-memory object.
    repeat=replay(cfg);repeat_base=repeat.base_g.copy()
    for base,ids in edge_ids.items():repeat.base_g[ids]=repeat_base[ids]*best['gains'][base]
    fresh=repeat.evaluate(fixed,time.time()+180.)
    matched=abs(fresh['value']-best['value'])<1e-10
    (a.output/'FINAL.json').write_text(json.dumps(dict(accepted=bool(best['diagnostics']['accepted']and fresh['diagnostics']['accepted']and matched),fresh_reconstruction_matched=matched,best=best,repeat=fresh,scope=manifest['scope']),indent=2))
    # A real fresh-process check remains separate from reconstruction above.
if __name__=='__main__':main()
