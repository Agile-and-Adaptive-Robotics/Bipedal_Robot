"""Fit the same frozen A/B circuit's motor mapping across slow air and loaded timing."""
import os,sys,json,time,hashlib,shutil
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
H=Path(__file__).resolve().parent;campaign=os.environ.get('HUNT_CAMPAIGN','fourcell_multicondition_20261029_v1')
if Path(campaign).name!=campaign:raise ValueError('Campaign must be a directory name under hunt_training/runs')
O=H/'runs'/campaign;O.mkdir(exist_ok=True);fits=O/'weights';fits.mkdir(exist_ok=True)
if (O/'TRAINING_COMPLETE.json').exists():print('Existing completed training retained:',O);raise SystemExit(0)
assert (H/'runs/fourcell_hunt_20261029_v1/TRAINING_COMPLETE.json').exists()
import numpy as np
import motor_fit as M,controller as C
start=time.perf_counter();real_build=C.build;real_setup=M.setup;real_aff=M.afferents;real_fit=M.lsq_linear;MODE='afferented';captured=[]
pf_target_mode=os.environ.get('HUNT_PF_TARGET','gross_muscle_torque')
if pf_target_mode not in ('gross_muscle_torque','net_active_torque'):raise ValueError(pf_target_mode)
sources=O/'sources';sources.mkdir(exist_ok=True)
for p in (Path(__file__),H/'motor_fit.py',H/'controller.py'):
    dst=sources/p.name
    if dst.exists() and dst.read_bytes()!=p.read_bytes():raise RuntimeError('Changed sources require a new run directory')
    if not dst.exists():shutil.copyfile(p,dst)
def build(*args,**kwargs):
    n=real_build(*args,options=dict(weights_directory=str(fits),biarticular_pair_gain=.3,commissural_gain=0. if MODE=='deafferented' else .3,learn_reciprocal=True,use_timing_afferents=MODE!='deafferented',pf_target_mode=pf_target_mode),**kwargs)
    if MODE=='deafferented':n.reflex_gain[:]=0.
    return n
def setup(cycles=8):
    z,m,d,n,t,q,v,a,pf=real_setup(cycles)
    if MODE!='deafferented':return z,m,d,n,t,q,v,a,pf
    tt=np.arange(len(t)*2)*M.DT;query=tt/2
    def stretch(x):return np.column_stack([np.interp(query,t,x[:,k]) for k in range(x.shape[1])])
    return z,m,d,n,tt,stretch(q),stretch(v)/2,stretch(a),{key:np.interp(query,t,value) for key,value in pf.items()}
def aff(n,u,d):
    if MODE!='deafferented':real_aff(n,u,d)
def capture(A,b,**kwargs):
    captured.append((A.copy(),b.copy(),kwargs));return real_fit(A,b,**kwargs)
C.build=build;M.setup=setup;M.afferents=aff
for kind,filename in [('PF','pf_weights.json'),('BI','bi_weights.json'),('MN','mn_weights.json')]:
    if (fits/filename).exists():continue
    blocks=[];files=[]
    for condition in ('afferented','deafferented'):
        MODE=condition;M.O=O/(kind.lower()+'_'+condition);M.O.mkdir(exist_ok=True);captured=[];M.lsq_linear=capture
        M.replay(kind,2.466 if MODE=='deafferented' else 1.233)
        blocks.append(captured);files.append(json.loads((M.O/filename).read_text()))
    assert len(blocks[0])==len(blocks[1]);keys=list(files[0]['edges']);post_names=list(files[0]['regression']);weights={};diagnostics={};offset=0
    for j,pair in enumerate(zip(*blocks)):
        matrices=[];targets=[]
        for A,b,kw in pair:
            count=A.shape[0]-A.shape[1];matrices.append(A/np.sqrt(count));targets.append(b/np.sqrt(count))
        solution=real_fit(np.vstack(matrices),np.concatenate(targets),**pair[0][2]);size=pair[0][0].shape[1]
        for key,value in zip(keys[offset:offset+size],solution.x):weights[key]=max(0.,float(value))
        offset+=size;diagnostics[post_names[j]]={condition:float(np.sqrt(np.mean((A[:-size]@solution.x-b[:-size])**2))) for condition,(A,b,kw) in zip(('afferented','deafferented'),pair)}
    assert offset==len(keys)
    (fits/filename).write_text(json.dumps(dict(block=kind,edges=weights,regression=diagnostics,timing_frozen=True,scope='Equal-condition joint fit of expected-afferent 1.233 s and deafferented 2.466 s targets; no runtime teacher clamp.'),indent=2));print('JOINT',kind,'edges',len(weights),flush=True)
MODE='afferented';M.lsq_linear=real_fit;M.O=fits;M.validate(.8,'step_C_afferented')
MODE='deafferented';M.validate(.8,'step_C_deafferented')
MODE='afferented'
if not (fits/'reflex_weights.json').exists():M.replay('REFLEX')
M.validate()
(O/'TRAINING_COMPLETE.json').write_text(json.dumps(dict(wall_s=time.perf_counter()-start,pf_target_mode=pf_target_mode,conditions=['Expected afferents at 1.233 s','Deafferented at 2.466 s'],upstream_validation=str(H/'runs/fourcell_hunt_20261029_v1'),
    source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources.iterdir()},weights={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in fits.glob('*weights.json')},physical_curriculum_passed=False),indent=2));print('TRAINING_COMPLETE',flush=True)
