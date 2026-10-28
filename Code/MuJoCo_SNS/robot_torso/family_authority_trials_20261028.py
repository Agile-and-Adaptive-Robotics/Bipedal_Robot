"""Expose natural PF/reflex MN currents in isolated supported CAD trials."""
import os,sys,json,argparse,hashlib
from pathlib import Path
from types import SimpleNamespace
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'spinal'))
import run_neural_robot as R
OUT=HERE.parent/'spinal/reports_20261028_family_audit/robot_trials';OUT.mkdir(parents=True,exist_ok=True)
def main(family,authority,delete_pf,duration,ramp,backend):
    adapter_bytes=Path(__file__).read_bytes();adapter_sha=hashlib.sha256(adapter_bytes).hexdigest()
    snapshots=OUT/'source_snapshots';snapshots.mkdir(exist_ok=True)
    snapshot=snapshots/(adapter_sha+'.py')
    if not snapshot.exists():snapshot.write_bytes(adapter_bytes)
    if backend=='preallocated':
        from preallocated_exact_eval import PreallocatedExactForward
        R.ExactFastForward=PreallocatedExactForward
    original_build=R.bn.build;original_clamp=R.motor_target_current;capture={};indices={}
    def build(acts,dt,interleg):
        if family=='s3k':net=original_build(acts,dt=dt,interleg=interleg)
        elif family=='w2lvar':
            import build_network_w2lvar as B
            net=B.build(acts,dt=dt,interleg=interleg)
        else:
            import build_network_syn6 as B
            net=B.build(acts,dt=dt,interleg=interleg)
        names=[p['name'] for p in net.net.populations]
        pf=np.array([i for i,n in enumerate(names) if n.startswith('PF_') and not n.startswith('PF_IN')],dtype=int)
        mn=np.array([net.idx['MN_'+act] for act in acts]);matrix=net.compiled.g_max_non
        before=matrix[np.ix_(mn,pf)].copy()
        indices.update(mn=mn,pf=pf)
        capture.update(true_neurons=int(net.compiled.V.size),pf_to_mn_edges_before=int(np.count_nonzero(before)),
            pf_to_mn_edges_after=int(np.count_nonzero(matrix[np.ix_(mn,pf)])),effective_gains=dict(R.params.G))
        return net
    call_count=0;deleted=False
    def current(compiled,mn,desired,tau):
        nonlocal call_count,deleted
        call_count+=1
        # The inherited warm-up calls this function 100 times before the
        # first physical step. One call per 2 ms physical step follows.
        t=max(0,call_count-101)*.002
        if delete_pf and not deleted and t>=15.:
            # Delayed deletion preserves the entire pre-intervention history.
            target,source=indices['mn'],indices['pf']
            compiled.g_max_non[np.ix_(target,source)]=0.
            mask=np.isin(compiled.post,target)&np.isin(compiled.pre,source)
            compiled.gmax[mask]=0.
            compiled.g[np.ix_(target,source)]=0.;compiled.ge[np.ix_(target,source)]=0.
            capture['pf_to_mn_edges_after']=0;capture['pf_deletion_time_s']=t;deleted=True
        fraction=authority*np.clip((t-5.)/10.,0.,1.) if ramp else authority
        # Equilibrium postural bias, capped at desired activation in [0,1].
        # At authority=1, no instantaneous PF/reflex synaptic current is canceled.
        bias=compiled.g_m[mn]*(R.params.E_HI*np.clip(desired,0,1)-compiled.V_rest[mn])-compiled.i_b[mn]
        if fraction==1.:return bias
        return (1-fraction)*original_clamp(compiled,mn,desired,tau)+fraction*bias
    R.bn.build=build;R.motor_target_current=current
    settings=json.loads((HERE/'runs/abd_defense_exact.json').read_text())['run_options']
    label=family+'_authority'+str(authority).replace('.','p')+('_ramp' if ramp else '')+('_pf_deleted' if delete_pf else '')+('_preallocated' if backend=='preallocated' else '')
    settings.update(duration=duration,output=OUT/label);settings['model']=Path(settings['model'])
    summary=R.run(SimpleNamespace(**settings))
    summary.update(family=family,neural_authority=authority,authority_ramp_5_to_15_s=ramp,pf_to_mn_deleted=delete_pf,
        controller='Engineering allocated posture bias with uncanceled native PF/reflex MN synaptic currents' if authority==1. else 'Partial MN-clamp release; engineered posture target remains',
        adapter_sha256=adapter_sha,adapter_snapshot=str(snapshot),experimental_neural_backend=backend,compiled_audit=capture,
        scope='Exploratory supported CAD neural-authority test. Same original mount and 80% post-ramp support; not a calibrated BPA model.')
    settings['output'].with_suffix('.json').write_text(json.dumps(summary,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('family',choices=['s3k','w2lvar','syn6']);p.add_argument('--authority',type=float,default=1.);p.add_argument('--delete-pf',action='store_true');p.add_argument('--duration',type=float,default=25.);p.add_argument('--ramp',action='store_true');p.add_argument('--backend',choices=['exact_fast','preallocated'],default='exact_fast');a=p.parse_args()
    assert 0<=a.authority<=1
    main(a.family,a.authority,a.delete_pf,a.duration,a.ramp,a.backend)
