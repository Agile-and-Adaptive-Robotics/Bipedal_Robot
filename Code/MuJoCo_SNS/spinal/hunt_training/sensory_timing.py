"""Hunt step B: fit six afferent projections with intrinsic parameters frozen."""
import os,sys,json,time,hashlib,argparse
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parent;S=H.parent;sys.path.insert(0,str(S))
import numpy as np,optuna
from scipy.signal import find_peaks
from intrinsic import build,arrays,simulate,DT
import bsolve_ik as B
import kine_ref as KR
def phase_inputs(dataset,side,cycles=12):
    z=np.load(H/'runs/targets/targets.npz',allow_pickle=True)
    tg,names,grf=B.read_mot(B.GRF_MOT)
    vy=grf[:,names.index(('' if side==0 else '1_')+'ground_force_vy')]
    onsets=KR._loading_onsets(tg,vy);begin,end=KR._first_pair_in_ik(onsets,z['t']);period=end-begin
    t=np.arange(int(cycles*period/DT))*DT;phi=t%period
    features=np.column_stack([np.interp(begin+phi,z['t'],z['features'][:,side,k]) for k in range(6)])
    stance=np.interp(begin+phi,tg,vy)>50.
    u=np.zeros((len(t),10));u[:50,0]=2.;u[:,4:]=5.*np.clip(features,0.,1.5)
    return u,stance,period
def metrics(v,stance,period):
    keep=np.arange(len(v))*DT>4*period;contrast=v[:,0]-v[:,1]
    prediction=1./(1.+np.exp(-np.clip(contrast,-30,30)))
    mismatch=float(np.mean((prediction[keep]-stance[keep])**2))
    transitions=np.flatnonzero((contrast[1:]>0)&(contrast[:-1]<=0))+1
    transitions=transitions[transitions*DT>4*period]
    intervals=np.diff(transitions)*DT;T=float(np.median(intervals)) if len(intervals) else 0.
    phase=((transitions*DT)%period)/period;phase_error=float(np.mean(np.minimum(phase,1.-phase))) if len(phase) else 1.
    loss=mismatch+abs(T/period-1.)+phase_error
    return dict(loss=loss,stance_probability_mse=mismatch,period_s=T,target_period_s=period,
        transition_phase_error=phase_error,cycles=len(transitions),finite=bool(np.all(np.isfinite(v))))
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--side',type=int,choices=(0,1),default=0);args=parser.parse_args()
    out=H/('runs/step_B' if args.side==0 else 'runs/step_B_left');out.mkdir(parents=True,exist_ok=True)
    if (out/'winner.json').exists():
        prior=(out/'winner.json').read_bytes();(out/('prior_'+hashlib.sha256(prior).hexdigest()+'.json')).write_bytes(prior)
    a=json.loads((H/'runs/step_A/winner.json').read_text());assert a['accepted']
    u,stance,period=phase_inputs(None,args.side);ul,sl,pl=phase_inputs(None,1-args.side)
    c=build(a['tau_h'],a['mutual_inhibition'],a['nap'],np.ones(6));p=arrays(c)
    initial_v=c.V.copy();initial_h=c.b_gate[0].copy()
    destinations=np.array([1,0,1,0,0,0]);sources=np.arange(4,10)
    def evaluate(gains,inputs=u,desired=stance,T=period):
        matrix=p[0].copy();matrix[destinations,sources]=gains
        v,h=simulate((matrix,*p[1:]),initial_v.copy(),initial_h.copy(),inputs)
        return metrics(v,desired,T),v
    study=optuna.create_study(direction='minimize',storage='sqlite:///'+str(out/'timing.db').replace('\\','/'),
        study_name='frozen_intrinsic_sensory',load_if_exists=True,sampler=optuna.samplers.TPESampler(seed=281028))
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    right=H/'runs/step_B/winner.json'
    if args.side==1 and not study.trials and right.exists():
        seed=json.loads(right.read_text())['afferent_gains'];study.enqueue_trial({'afferent_'+str(i):seed[i] for i in range(6)})
    def objective(trial):
        gains=[trial.suggest_float('afferent_'+str(i),.001,20.,log=True) for i in range(6)]
        result,_=evaluate(gains);trial.set_user_attr('metrics',result)
        if trial.number%32==0:print(json.dumps(dict(trial=trial.number,**result)),flush=True)
        return result['loss']
    remaining=max(0,512-len(study.trials))
    study.optimize(objective,n_trials=remaining)
    gains=np.array([study.best_params['afferent_'+str(i)] for i in range(6)])
    train,v=evaluate(gains);heldout,vl=evaluate(gains,ul,sl,pl)
    c=build(a['tau_h'],a['mutual_inhibition'],a['nap'],gains);native=np.empty_like(v)
    for i in range(len(u)):c.forward(np.r_[u[i,:2],u[i,4:]]);native[i]=c.V
    native_metric=metrics(native,stance,period);error=float(np.max(abs(native-v)))
    rng=np.random.default_rng(281028);robust=[]
    for k in range(8):
        up=u.copy();up[:,4:]*=rng.uniform(.99,1.01,(1,6));metric,_=evaluate(gains,up,stance,period);robust.append(metric)
    accepted=all(r['finite'] and abs(r['period_s']/r['target_period_s']-1)<.05 and r['transition_phase_error']<.08 and r['stance_probability_mse']<.15 for r in (train,native_metric,*robust))
    result=dict(step='B_sensory_timing',side='r' if args.side==0 else 'l',afferent_gains=gains.tolist(),intrinsic=a,train=train,opposite_leg_transfer=heldout,
        validation_scope='Eight independent +/-1% afferent-amplitude perturbations; opposite leg is separately calibrated after the shared-weight transfer failure, not called an independent subject holdout',robustness=robust,native_replay=native_metric,
        max_voltage_disagreement=error,accepted=accepted,trials=len(study.trials),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'winner.json').write_text(json.dumps(result,indent=2));np.savez_compressed(out/'validation.npz',v=v,left_v=vl,native_v=native,stance=stance)
    print('STEP_B',json.dumps({k:v for k,v in result.items() if k!='intrinsic'}),flush=True)
if __name__=='__main__':main()
