"""Prepare a side-corrected Di Russo-inspired imitation fit; opt-in execution.

No simulation is launched. This fits an engineering activation target using the
existing wrapped primitive basis, not the paper's full objective or reflex fit.
"""
from pathlib import Path
import os,sys,json,hashlib,argparse
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import controller as C
DATA=HERE.parent/'bsolve_out.npz'
GRF=HERE.parents[3]/'Solid_Models/OpenSim/Gait2392_Robotbody/subject01_walk1_grf.mot'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def measured_onsets(path):
    names=None;rows=[];ready=False
    for line in path.read_text(encoding='utf-8',errors='replace').splitlines():
        if line.startswith('endheader'):ready=True;continue
        if not ready or not line.strip():continue
        if line.startswith('time'):names=line.split();continue
        rows.append([float(x) for x in line.split()])
    values=np.array(rows);assert names and np.all(np.diff(values[:,0])>0.)
    result={}
    for side,key in [('r','ground_force_vy'),('l','1_ground_force_vy')]:
        on=values[:,names.index(key)]>20.;raw=values[1:,0][on[1:]&~on[:-1]];keep=[]
        for t in raw:
            if not keep or t-keep[-1]>.45:keep.append(float(t))
        result[side]=keep
    return result

def prepare(out,execute=False):
    out.mkdir(parents=True,exist_ok=True)
    if execute and (out/'FIT_RESULT.json').exists():raise FileExistsError('Preserve previous fit; choose a new output directory.')
    z=np.load(DATA,allow_pickle=True);t=z['t'];acts=z['acts'];names=z['act_names'].tolist()
    assert np.all(np.isfinite(acts)) and np.all(np.diff(t)>0.) and acts.min()>=0. and acts.max()<=1.+1e-12
    onset=measured_onsets(GRF);cycles={};targets={};binding={}
    for side in ('r','l'):
        targets[side]=np.empty((len(t),9));binding[side]={}
        for mi,group in enumerate(C.MUSCLES):
            channels=[name.format(s=side) for name in C.MAP[group]]
            assert all(name.endswith('_'+side) for name in channels)
            columns=[names.index(name) for name in channels]
            targets[side][:,mi]=acts[:,columns].mean(axis=1)
            binding[side][group]=dict(actuator_names=channels,columns=columns,aggregation='arithmetic mean of reference activation; unchanged legacy target definition')
        cycles[side]=[]
        for a,b in zip(onset[side][:-1],onset[side][1:]):
            if a<t[0] or b>t[-1]:continue
            indices=np.flatnonzero((t>=a)&(t<b));phase=(t[indices]-a)/(b-a)
            if len(indices)<10:continue
            cycles[side].append(dict(start_s=a,end_s=b,indices=indices,phase=phase,P=np.stack([C.primitives(p) for p in phase])))
    assert cycles['r'],'No complete canonical-right cycle in activation data.'
    # Entire right-cycle holdout where available; otherwise transparently limited
    # interleaved sample checks. Left never contributes to fitted parameters.
    training=cycles['r'][:-1] if len(cycles['r'])>=2 else cycles['r']
    full_holdout=cycles['r'][-1:] if len(cycles['r'])>=2 else []
    X=[];Y=[];right_test=[]
    for cycle in training:
        test=np.arange(len(cycle['indices']))%5==0 if not full_holdout else np.zeros(len(cycle['indices']),bool)
        train=~test;weight=1./np.sqrt(int(train.sum()))
        X.append(cycle['P'][train]*weight);Y.append(targets['r'][cycle['indices'][train]]*weight)
        if test.any():right_test.append((cycle['P'][test],targets['r'][cycle['indices'][test]]))
    for cycle in full_holdout:right_test.append((cycle['P'],targets['r'][cycle['indices']]))
    left_test=[(c['P'],targets['l'][c['indices']]) for c in cycles['l']]
    def desc(c):return dict(start_s=c['start_s'],end_s=c['end_s'],samples=len(c['indices']),phase_min=float(c['phase'].min()),phase_max=float(c['phase'].max()))
    plan=dict(format='side_correct_bounded_mirrored_primitive_fit_v1',machine=os.environ.get('COMPUTERNAME'),training_executed=execute,canonical_side='right',basis='existing periodic/wrapped five raised cosines; not a new NMF basis',weights=dict(shape=[18,5],bounds=[-1.,1.],shared_parameters=45,left_equals_right=True,ordering=C.MUSCLES),sources={str(p):sha(p) for p in (DATA,GRF,HERE/'controller.py',HERE/'dirusso_W_stage1.npy',Path(__file__))},dataset_time_s=[float(t[0]),float(t[-1])],dataset_frames=len(t),heel_cue='Measured TOTAL vertical foot GRF crossing20N upward, merged only within450ms; not a separate regional heel sensor.',GRF_columns={'r':'ground_force_vy','l':'1_ground_force_vy'},onsets_s=onset,complete_cycles={s:[desc(c) for c in cycles[s]] for s in ('r','l')},binding=binding,training_samples=sum(len(y) for y in Y),validation_right_samples=sum(len(y) for _,y in right_test),validation_left_samples=sum(len(y) for _,y in left_test),independent_right_cycle_holdout=bool(full_holdout),validation_limitation='With one canonical cycle, every fifth right sample tests interpolation only; entire left cycle is held out for mirrored functional comparison, not independent-subject/cycle generalization.',target='Back-solved model activation, not measured EMG or published Di Russo training data.',fixed_preregistered_reconstruction_checks=dict(per_muscle_RMSE_max=.20,correlation_min_if_target_std_at_least_0p05=.50),body_trial_allowed=False,reflexes_fitted=False,paper_eq13_objective_implemented=False,interface='paper_checked_variant.py --weights ABSOLUTE_W_MIRRORED_NPY; existing default W remains unchanged')
    (out/'FIT_PLAN.json').write_text(json.dumps(plan,indent=2))
    if not execute:print(json.dumps({k:plan[k] for k in ('canonical_side','complete_cycles','training_samples','validation_right_samples','validation_left_samples','independent_right_cycle_holdout','training_executed')},indent=2));return plan
    from scipy.optimize import lsq_linear
    x=np.vstack(X);y=np.vstack(Y);W=np.empty((9,5));solver=[]
    for mi in range(9):
        fitted=lsq_linear(x,y[:,mi],bounds=(-1.,1.),method='bvls',tol=1e-12,max_iter=100)
        assert fitted.success and np.all(np.isfinite(fitted.x));W[mi]=fitted.x
        solver.append(dict(muscle=C.MUSCLES[mi],status=int(fitted.status),optimality=float(fitted.optimality)))
    full=np.vstack((W,W));assert np.array_equal(full[:9],full[9:]) and abs(full).max()<=1.+1e-12
    path=out/'W_mirrored.npy';np.save(path,full)
    def assess(data):
        if not data:return dict(available=False,pass_reconstruction=False)
        p=np.vstack([a for a,b in data]);target=np.vstack([b for a,b in data]);prediction=p@W.T
        rmse=np.sqrt(np.mean((prediction-target)**2,axis=0));corr=[]
        for i in range(9):corr.append(float(np.corrcoef(prediction[:,i],target[:,i])[0,1]) if np.std(target[:,i])>=.05 and np.std(prediction[:,i])>1e-12 else None)
        passed=bool(np.all(rmse<=.2) and all(c is None or c>=.5 for c in corr))
        return dict(available=True,samples=len(target),muscle_order=C.MUSCLES,RMSE=rmse.tolist(),correlation=corr,pass_reconstruction=passed)
    result=dict(training_executed=True,weights_path=str(path),weights_sha256=sha(path),solver=solver,right_validation=assess(right_test),left_mirror_validation=assess(left_test),body_trial_allowed=False,physical_success=False,paper_replication=False,plan=str(out/'FIT_PLAN.json'))
    (out/'FIT_RESULT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--execute',action='store_true');a=p.parse_args();prepare(a.out,a.execute)
