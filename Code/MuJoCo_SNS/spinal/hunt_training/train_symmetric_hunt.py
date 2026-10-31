"""Fit six shared timing parameters against both legs; refit tied C/D.

No independent left-side parameter vector. Loaded offset cases are stress
probes, not loaded-task ground-truth targets or accepted walking evidence.
"""
import os,sys,json,time,shutil,runpy,subprocess
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parent)]
import numpy as np,optuna
from scipy.ndimage import gaussian_filter1d
import intrinsic as A,sensory_timing as B
from bilateral import canonical,assert_tied
name=os.environ.get('HUNT_SYMMETRIC_CAMPAIGN','symmetric_hunt_20261030_s1_v1');out=H/'runs'/name;W=out/'weights'
if out.exists():raise SystemExit('Existing symmetric fit retained')
W.mkdir(parents=True);start=time.perf_counter()
source=H/'runs/tonic_abduction_20261030_s1_v1/0p02/weights'
z=np.load(H/'runs/swing_release_20261030_v1/feasible_targets/fast.teacher_targets.npz',allow_pickle=True)
T=1.233;grid=np.arange(256)/256;phase=((z['t']-2)%T)/T;cyc=np.floor((z['t']-2)/T).astype(int);features=[]
for side in (0,1):
 total=np.zeros((256,6));count=0
 for cycle in np.unique(cyc[z['t']>=8]):
  ix=np.flatnonzero((z['t']>=8)&(cyc==cycle))
  if len(ix)<.85*T/.01:continue
  ph=phase[ix];values=z['features'][ix,side];total+=np.column_stack([np.interp(grid,np.r_[ph[-1]-1,ph,ph[0]+1],np.r_[values[-1,k],values[:,k],values[0,k]]) for k in range(6)]);count+=1
 features.append(gaussian_filter1d(total/count,1,axis=0,mode='wrap'))
a=json.loads((H/'runs/step_A/winner.json').read_text());c=A.build(a['tau_h'],a['mutual_inhibition'],a['nap'],np.ones(6));p=A.arrays(c);v0=c.V.copy();h0=c.b_gate[0].copy()
t=np.arange(int(12*T/.002))*.002;ph=t/T%1;stance=ph<.61;cases=[]
for side in (0,1):
 u=np.zeros((len(t),10));u[:50,0]=2.;q=(ph+.5*side)%1
 for k in range(6):u[:,4+k]=5*np.interp(q,np.r_[grid,1],np.r_[features[side][:,k],features[side][0,k]])
 cases.append(u)
labels=['r','l']
if os.environ.get('HUNT_TIMING_STRESS')=='1':
 for side,u in zip(('r','l'),cases.copy()):
  shifted=u.copy();shifted[:,5]=7.5;cases.append(shifted);labels.append('loaded_offset_'+side)
def evaluate(g,u):
 matrix=p[0].copy();matrix[np.array([1,0,1,0,0,0]),np.arange(4,10)]=g
 v,_=A.simulate((matrix,*p[1:]),v0.copy(),h0.copy(),u);return B.metrics(v,stance,T),v
optuna.logging.set_verbosity(optuna.logging.WARNING)
study=optuna.create_study(sampler=optuna.samplers.TPESampler(seed=301031))
old=json.loads((source/'timing_weights.json').read_text())['gains']['r'];study.enqueue_trial({str(k):min(x,2.) if k<2 else x for k,x in enumerate(old)})
def objective(trial):
 g=[trial.suggest_float(str(k),.0001,2. if k<2 else 20.,log=True) for k in range(6)]
 losses=[evaluate(g,u)[0]['loss'] for u in cases]
 return max(losses)+.25*np.mean(losses)
study.optimize(objective,n_trials=2048)
g=[study.best_params[str(k)] for k in range(6)];metrics={}
for side,u in zip(labels,cases):
 r,v=evaluate(g,u);native=A.build(a['tau_h'],a['mutual_inhibition'],a['nap'],g);nv=np.empty_like(v)
 for i,row in enumerate(u):native.forward(np.r_[row[:2],row[4:]]);nv[i]=native.V
 r['native_max_error']=float(abs(v-nv).max());r['accepted']=bool(abs(r['period_s']/T-1)<.05 and r['transition_phase_error']<.08 and r['stance_probability_mse']<.15 and r['native_max_error']<1e-10);metrics[side]=r
accepted=all(r['accepted'] for r in metrics.values())
timing=dict(gains=dict(r=g,l=g),metrics=metrics,accepted=accepted,bilateral_tied=True,free_parameters=6,source='Joint fit of both expected-air afferent sets; intrinsic A frozen')
(W/'timing_weights.json').write_text(json.dumps(timing,indent=2));print('SHARED_B',json.dumps(timing),flush=True)
if not accepted:raise SystemExit('Shared timing failed gates; no motor or physical promotion')
pf=json.loads((source/'pf_weights.json').read_text());values={}
for edge,value in pf['edges'].items():values.setdefault(canonical(edge),value)
for edge in pf['edges']:pf['edges'][edge]=pf['edges'].get(canonical(edge),values[canonical(edge)])
assert_tied(pf['edges']);pf['scope']='Explicit right-to-left PF symmetry projection before tied C/D refit; not independently fitted left PF.'
(W/'pf_weights.json').write_text(json.dumps(pf,indent=2));(W/'circuit_options.json').write_text(json.dumps(dict(require_bilateral_symmetry=True),indent=2))
os.environ.update(HUNT_FEASIBLE_CAMPAIGN=name,HUNT_FEASIBLE_TARGET='mono_command',HUNT_FEASIBLE_MODES='afferented',HUNT_BILATERAL_TIED='1')
runpy.run_path(str(H/'train_feasible_hunt.py'),run_name='__main__')
for p in W.glob('*weights.json'):
 if 'edges' in (obj:=json.loads(p.read_text())):assert_tied(obj['edges'])
subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',str(W),'--output-directory',str(out/'air_pair')],check=True)
print('SYMMETRIC_HUNT_COMPLETE',time.perf_counter()-start,flush=True)
