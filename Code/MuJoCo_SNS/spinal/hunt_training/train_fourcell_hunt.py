"""Versioned Hunt A/B revalidation, fresh four-cell C fitting, then D fitting."""
import os,sys,json,hashlib,shutil,time
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parent;O=H/'runs/fourcell_hunt_20261029_v1';O.mkdir(exist_ok=True);fits=O/'weights';fits.mkdir(exist_ok=True)
if (O/'TRAINING_COMPLETE.json').exists():print('Existing completed training retained:',O);raise SystemExit(0)
import numpy as np,mujoco
from scipy.signal import find_peaks
import motor_fit as M,controller as C,intrinsic as A,sensory_timing as B
start=time.perf_counter()
sources=O/'sources';sources.mkdir(exist_ok=True)
for p in (Path(__file__),H/'motor_fit.py',H/'controller.py',H/'physical_eval.py'):
    target=sources/p.name
    if target.exists() and target.read_bytes()!=p.read_bytes():raise RuntimeError('Changed training sources require a new run directory: '+p.name)
    if not target.exists():shutil.copyfile(p,target)
a=json.loads((H/'runs/step_A/winner.json').read_text());u=np.zeros((20000,4));u[:50,0]=2.
c=A.build(a['tau_h'],a['mutual_inhibition'],a['nap']);v,h=A.simulate(A.arrays(c),c.V.copy(),c.b_gate[0].copy(),u)
av=A.measure(v);av['accepted']=bool(av['finite'] and abs(av['period_s']/a['target_undriven_period_s']-1)<.02 and av['cycles']>=8)
(O/'step_A.json').write_text(json.dumps(av,indent=2));print('A',json.dumps(av),flush=True)
if not av['accepted']:raise RuntimeError('Intrinsic gate failed')
b_results={}
for side in (0,1):
    winner=json.loads((H/('runs/step_B' if side==0 else 'runs/step_B_left')/'winner.json').read_text());u,stance,period=B.phase_inputs(None,side)
    c=A.build(a['tau_h'],a['mutual_inhibition'],a['nap'],winner['afferent_gains']);v,h=A.simulate(A.arrays(c),c.V.copy(),c.b_gate[0].copy(),u)
    metric=B.metrics(v,stance,period);metric['accepted']=bool(metric['finite'] and abs(metric['period_s']/period-1)<.05 and metric['transition_phase_error']<.08)
    b_results[str(side)]=metric
(O/'step_B_isolated.json').write_text(json.dumps(b_results,indent=2));print('B isolated',json.dumps(b_results),flush=True)
if not all(r['accepted'] for r in b_results.values()):raise RuntimeError('Sensory timing gate failed')
real=C.build
def build(*args,**kwargs):
    options=kwargs.pop('options',{}) or {}
    options.update(weights_directory=str(fits),biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True)
    return real(*args,options=options,**kwargs)
C.build=build;M.O=fits
# Revalidate timing in the actual paired/coupled network before fitting outputs.
z,m,d,n,t,q,vel,act,pftarget=M.setup();u=n.make_inputs();u[n.ports['DRIVE']]=2.5
trace=np.empty((len(t),4))
for k in range(len(t)):
    d.qpos[:]=q[k];n.arm_model.project(d);d.qvel[:]=vel[k];d.act[:]=act[k];mujoco.mj_forward(m,d);M.afferents(n,u,d)
    u[n.ports['HUNT_INIT_E_r']]=2. if k<50 else 0.;u[n.ports['HUNT_INIT_F_l']]=2. if k<50 else 0.
    n.step(u);trace[k]=[n.compiled.V[n.idx['RG_'+phase+'_'+side]] for side in ('r','l') for phase in ('E','F')]
coupled={}
for side,col in [('r',0),('l',2)]:
    peaks=find_peaks(trace[:,col]-trace[:,col+1],prominence=2.,distance=100)[0];peaks=peaks[t[peaks]>4*1.233]
    period=float(np.median(np.diff(t[peaks]))) if len(peaks)>1 else 0.
    coupled[side]=dict(period_s=period,target_s=1.233,cycles=len(peaks),accepted=bool(abs(period/1.233-1)<.05 and len(peaks)>=3))
(O/'step_B_coupled.json').write_text(json.dumps(coupled,indent=2));print('B coupled',json.dumps(coupled),flush=True)
if not all(r['accepted'] for r in coupled.values()):raise RuntimeError('Coupled expected-afferent timing failed before C')
for kind in ('PF','BI','MN'):
    filename={'PF':'pf_weights.json','BI':'bi_weights.json','MN':'mn_weights.json'}[kind]
    if not (fits/filename).exists():M.replay(kind)
M.validate(.8,'step_C_native')
if not (fits/'reflex_weights.json').exists():M.replay('REFLEX')
M.validate()
manifest=dict(pair_gain=.3,commissural_gain=.3,learned_reciprocal=True,sartorius_gracilis_separate=True,
    source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources.iterdir()},wall_s=time.perf_counter()-start,
    weights={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in fits.glob('*weights.json')},scope='A/B revalidated and frozen; C PF, four-cell integration, MN output, then D reflex fitted sequentially. Physical curriculum still required.')
(O/'TRAINING_COMPLETE.json').write_text(json.dumps(manifest,indent=2));print('TRAINING_COMPLETE',json.dumps(manifest),flush=True)
