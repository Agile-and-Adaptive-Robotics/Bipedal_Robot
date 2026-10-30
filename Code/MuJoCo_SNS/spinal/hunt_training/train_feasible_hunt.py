"""Hunt A/B/C/D using commands and expected afferents from feasible air teachers.

Clock and inverse control are confined to teacher-data generation. The final
network receives only the original neural inputs and native afferents.
"""
import os,sys,json,time,hashlib,shutil
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parent;O=H/'runs'/os.environ.get('HUNT_FEASIBLE_CAMPAIGN','feasible_hunt_20261030_v1');O.mkdir(exist_ok=True);W=O/'weights';W.mkdir(exist_ok=True)
inactive_voltage=float(os.environ.get('HUNT_PF_INACTIVE','0'))
fit_modes=os.environ.get('HUNT_FEASIBLE_MODES','afferented,deafferented').split(',')
assert fit_modes in (['afferented'],['afferented','deafferented'])
if os.environ.get('HUNT_JOINT_PF_INHIBITION'):
    (W/'circuit_options.json').write_text(json.dumps(dict(joint_pf_inhibition_scale=float(os.environ['HUNT_JOINT_PF_INHIBITION'])),indent=2))
if (O/'TRAINING_COMPLETE.json').exists():print('Existing completed fit retained');raise SystemExit(0)
import numpy as np,optuna,mujoco
from scipy.ndimage import gaussian_filter1d
import motor_fit as M,controller as C,intrinsic as A,sensory_timing as B
from pf_targets import mono_command_targets
start=time.perf_counter();datasets={};source=H/'runs/swing_release_20261030_v1/feasible_targets'
def folded(z,key,T):
    phase=((z['t']-2.)%T)/T;valid=z['t']>=8.;grid=np.arange(256)/256
    x=z[key];shape=x.shape[1:];flat=x.reshape(len(x),-1);result=np.zeros((256,flat.shape[1]));count=np.zeros(256)
    # Average each successful physical cycle, not every phase sample from startup.
    cyc=np.floor((z['t']-2.)/T).astype(int)
    for cycle in np.unique(cyc[valid]):
        ids=np.flatnonzero(valid&(cyc==cycle))
        if len(ids)<int(.85*T/.01):continue
        ph=phase[ids];order=np.argsort(ph);ph=ph[order];vals=flat[ids][order]
        xx=np.r_[ph[-1]-1,ph,ph[0]+1];yy=np.vstack([vals[-1],vals,vals[0]])
        result+=np.column_stack([np.interp(grid,xx,yy[:,j]) for j in range(flat.shape[1])]);count+=1
    assert np.all(count>0),(key,T)
    return gaussian_filter1d((result/count[:,None]).reshape((256,)+shape),1.,axis=0,mode='wrap')
for condition,filename,T in [('afferented','fast',1.233),('deafferented','slow',2.466)]:
    path=source/(filename+'.teacher_targets.npz');z=np.load(path,allow_pickle=True)
    proof=json.loads((source/(filename+'.result.json')).read_text());assert proof['assessment']['teacher_motion_gate']
    datasets[condition]={key:folded(z,key,T) for key in ('qpos','qvel','activation','command','rho','passive','capacity','features')}
    datasets[condition]['period']=T
a=json.loads((H/'runs/step_A/winner.json').read_text());assert a['accepted']
timing_path=W/'timing_weights.json'
if not timing_path.exists() and os.environ.get('HUNT_TIMING_SOURCE'):
    shutil.copyfile(Path(os.environ['HUNT_TIMING_SOURCE']),timing_path)
if not timing_path.exists():
    optuna.logging.set_verbosity(optuna.logging.WARNING);gains={};metrics={}
    for side,col in [('r',0),('l',1)]:
        T=1.233;t=np.arange(int(12*T/.002))*.002;phase=(t/T)%1;shift=0. if side=='r' else .5
        query=(phase+shift)%1;features=datasets['afferented']['features'][:,col,:]
        inputs=np.zeros((len(t),10));inputs[:50,0]=2.
        for f in range(6):inputs[:,f+4]=5.*np.interp(query,np.r_[np.arange(256)/256,1.],np.r_[features[:,f],features[0,f]])
        stance=phase<.61;c=A.build(a['tau_h'],a['mutual_inhibition'],a['nap'],np.ones(6));p=A.arrays(c)
        initial_v=c.V.copy();initial_h=c.b_gate[0].copy();dest=np.array([1,0,1,0,0,0]);src=np.arange(4,10)
        def evaluate(values):
            matrix=p[0].copy();matrix[dest,src]=values;v,h=A.simulate((matrix,*p[1:]),initial_v.copy(),initial_h.copy(),inputs)
            return B.metrics(v,stance,T),v
        study=optuna.create_study(direction='minimize',sampler=optuna.samplers.TPESampler(seed=301030+col))
        old=json.loads((H/('runs/step_B' if side=='r' else 'runs/step_B_left')/'winner.json').read_text())['afferent_gains'];study.enqueue_trial({str(i):v for i,v in enumerate(old)})
        def objective(trial):
            values=[trial.suggest_float(str(i),.001,20.,log=True) for i in range(6)];r,_=evaluate(values);return r['loss']
        study.optimize(objective,n_trials=768)
        values=[study.best_params[str(i)] for i in range(6)];r,v=evaluate(values)
        native=A.build(a['tau_h'],a['mutual_inhibition'],a['nap'],values);nv=np.empty_like(v)
        for k,u in enumerate(inputs):native.forward(np.r_[u[:2],u[4:]]);nv[k]=native.V
        r['max_native_error']=float(np.max(abs(nv-v)));r['accepted']=bool(abs(r['period_s']/T-1)<.05 and r['transition_phase_error']<.08 and r['stance_probability_mse']<.15 and r['max_native_error']<1e-10)
        gains[side]=values;metrics[side]=r;print('B',side,json.dumps(r),flush=True)
    timing=dict(gains=gains,metrics=metrics,accepted=all(r['accepted'] for r in metrics.values()),source='Expected afferents from feasible air teacher; stance proxy is reference phase, not physical air contact')
    timing_path.write_text(json.dumps(timing,indent=2))
timing=json.loads(timing_path.read_text());assert timing['accepted'],timing
MODE='afferented';real_build=C.build;real_setup=M.setup;real_aff=M.afferents;real_fit=M.lsq_linear;captured=[]
def build(*args,**kwargs):
    n=real_build(*args,options=dict(weights_directory=str(W),biarticular_pair_gain=.3,commissural_gain=0. if MODE=='deafferented' else .3,learn_reciprocal=True,use_timing_afferents=MODE!='deafferented'),**kwargs)
    if MODE=='deafferented':n.reflex_gain[:]=0.
    return n
def setup(cycles=8):
    z,m,d,n,*_=real_setup(cycles);record=datasets[MODE];T=record['period'];t=np.arange(int(cycles*T/.002))*.002;phase=(t/T)%1
    def replay(key):
        v=record[key];flat=v.reshape(256,-1);grid=np.r_[np.arange(256)/256,1.]
        return np.column_stack([np.interp(phase,grid,np.r_[flat[:,j],flat[0,j]]) for j in range(flat.shape[1])]).reshape((len(t),)+v.shape[1:])
    q=replay('qpos');v=replay('qvel');command=replay('command');n.expected_activation=replay('activation')
    strength=record['rho']*(record['capacity']*record['command'])[:,None,:];pf={}
    for joint,pair in n.hunt_pairs.items():
        row=C.JOINTS.index(joint)
        for sign,cell in zip((1.,-1.),pair):
            tau=np.maximum(sign*strength[:,row].sum(axis=1),0.);activity=np.clip(tau/max(float(np.percentile(tau,95)),1.),0,1)
            wave=2.75*activity+inactive_voltage*(1.-activity)
            pf[cell]=np.interp(phase,np.r_[np.arange(256)/256,1.],np.r_[wave,wave[0]])
    if os.environ.get('HUNT_FEASIBLE_TARGET')=='mono_command':pf.update(mono_command_targets(command,z['act_names'],n.hunt_pairs))
    return z,m,d,n,t,q,v,command,pf
def aff(n,u,d):
    if MODE!='deafferented':real_aff(n,u,d)
def capture(X,y,**kw):captured.append((X.copy(),y.copy(),kw));return real_fit(X,y,**kw)
C.build=build;M.setup=setup;M.afferents=aff
sources=O/'sources';sources.mkdir(exist_ok=True)
for path in (Path(__file__),H/'motor_fit.py',H/'controller.py',H/'motor_teacher.py'):shutil.copyfile(path,sources/path.name)
for block,filename in [('PF','pf_weights.json'),('BI','bi_weights.json'),('MN','mn_weights.json')]:
    if (W/filename).exists():continue
    blocks=[];files=[]
    for mode in fit_modes:
        MODE=mode;M.O=O/(block.lower()+'_'+mode);M.O.mkdir(exist_ok=True);captured=[];M.lsq_linear=capture;M.replay(block,datasets[MODE]['period']);blocks.append(captured);files.append(json.loads((M.O/filename).read_text()))
    keys=list(files[0]['edges']);assert all(keys==list(f['edges']) for f in files);weights={};offset=0
    if len(files)==1:weights=files[0]['edges'];offset=len(keys)
    else:
        for left,right in zip(*blocks):
            xx=[];yy=[]
            for X,y,kw in (left,right):factor=np.sqrt(X.shape[0]-X.shape[1]);xx.append(X/factor);yy.append(y/factor)
            solution=real_fit(np.vstack(xx),np.r_[yy[0],yy[1]],**left[2]);width=left[0].shape[1]
            weights.update({key:max(0.,float(value)) for key,value in zip(keys[offset:offset+width],solution.x)});offset+=width
    assert offset==len(keys);(W/filename).write_text(json.dumps(dict(edges=weights,scope='Feasible-air teacher motor commands; fixed timing',fit_modes=fit_modes,target_source=os.environ.get('HUNT_FEASIBLE_TARGET','net_active_torque'),inactive_pf_voltage=inactive_voltage),indent=2));print('JOINT',block,len(weights),flush=True)
    if block=='PF' and os.environ.get('HUNT_FEASIBLE_PF_ONLY')=='1':
        waveform={}
        for mode in ('afferented','deafferented'):
            MODE=mode;M.O=W;M.validate(.8,'pf_only_'+mode)
            measured=np.load(W/('pf_only_'+mode+'.npz'));teacher=np.load(O/('pf_'+mode)/'pf_teacher.npz');keep=measured['t']>=4*datasets[mode]['period'];cells=teacher['pf_cells'].tolist()
            rows={}
            for i,cell in enumerate(cells):
                if not cell.startswith(('PF_HIP','PF_KNEE')):continue
                actual=np.maximum(measured['pf'][keep,i],0.);target=np.maximum(teacher['pf_target'][keep,i],0.)
                corr=float(np.corrcoef(actual,target)[0,1]) if np.std(actual)>1e-8 and np.std(target)>1e-8 else 0.
                duty=float(np.mean(actual>.1));peak=float(actual.max())
                rows[cell]=dict(correlation=corr,positive_fraction=duty,peak=peak,accepted=bool(corr>.5 and .1<duty<.9 and peak>.5))
            waveform[mode]=rows
        (O/'PF_COMPLETE.json').write_text(json.dumps(dict(inactive_pf_voltage=inactive_voltage,native_waveform=waveform,accepted=all(r['accepted'] for rows in waveform.values() for r in rows.values()),scope='Preliminary free PF waveform check only; no motor fit or physical curriculum result'),indent=2))
        print('PF_ONLY_COMPLETE',flush=True);raise SystemExit(0)
MODE='afferented';M.O=W;M.lsq_linear=real_fit;M.validate(.8,'step_C_afferented');MODE='deafferented';M.validate(.8,'step_C_deafferented');MODE='afferented'
if not (W/'reflex_weights.json').exists():M.replay('REFLEX')
M.validate();(O/'TRAINING_COMPLETE.json').write_text(json.dumps(dict(wall_s=time.perf_counter()-start,teacher_driven_motion_is_neural_result=False,physical_curriculum_passed=False,source=str(source),timing=timing,fit_modes=fit_modes,target_source=os.environ.get('HUNT_FEASIBLE_TARGET','net_active_torque')),indent=2));print('FEASIBLE_HUNT_COMPLETE',flush=True)
