"""Read-only shadow-network audit of PT input→MN→source-force transfer.

The live standing trajectory is unchanged. Independent neural copies receive
fixed inputs, and independent MjData objects receive source steady forces.
All temporary model force-law arrays are restored before the live step resumes.
This diagnoses a local transfer, not recovery or walking performance.
"""
import os,sys,json,copy,hashlib,time,argparse,runpy
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
HERE=Path(__file__).resolve().parent;S=HERE.parent;REPO=S.parents[2]
sys.path[:0]=[str(HERE),str(S/'hunt_training'),str(S/'functional_gait2392'),str(S),str(S/'w2l_cpg')]
import numpy as np,mujoco
import controller as C
from thelen_dynamic import DynamicThelenBank,parameters_from_osim
SOURCE=REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody.osim'

def probe(n,u,out):
    tic=time.perf_counter();m,d=n.live;c=n.compiled;names=n.names;acts=list(n.muscles);mn=np.array([n.idx[n.mn_names[a]] for a in acts]);ports=n.balance.ports
    fields=('actuator_gaintype','actuator_gainprm','actuator_biastype','actuator_biasprm');saved={key:getattr(m,key).copy() for key in fields};data_fields=('qpos','qvel','qacc','act','ctrl','actuator_force','xfrc_applied','qacc_warmstart');actual={key:getattr(d,key).copy() for key in data_fields};voltage=c.V.copy();gmax=c.gmax.copy()
    p=parameters_from_osim(SOURCE,[m.actuator(i).name for i in range(m.nu)],0.);bank=DynamicThelenBank(p,float(m.opt.timestep));q=mujoco.MjData(m);jac=np.zeros((3,m.nv));rows=[]
    def physical(activation):
        q.qpos[:]=d.qpos;q.qvel[:]=0.;q.act[:]=activation;q.ctrl[:]=activation;q.xfrc_applied[:]=d.xfrc_applied;q.qacc_warmstart[:]=0.
        mujoco.mj_forward(m,q);bank.initialize(q.actuator_length,q.act)
        if np.any(bank.bounds):raise RuntimeError('Source fiber bound in static PT transfer probe')
        m.actuator_gaintype[:]=mujoco.mjtGain.mjGAIN_FIXED;m.actuator_gainprm[:]=0.;m.actuator_biastype[:]=mujoco.mjtBias.mjBIAS_AFFINE;m.actuator_biasprm[:]=0.;m.actuator_biasprm[:,0]=-bank.force;q.qacc_warmstart[:]=0.
        mujoco.mj_forward(m,q);mujoco.mj_jacSubtreeCom(m,q,jac,0);return (jac@q.qacc).copy()
    def neural(scale,label=None,amplitude=0.):
        # Avoid delegated __getattr__ during copy reconstruction of the
        # sparse wrapper; deep-copy its ordinary state dictionary instead.
        shadow=object.__new__(type(c));shadow.__dict__.update(copy.deepcopy(c.__dict__))
        if np.shares_memory(c.V,shadow.V) or np.shares_memory(c.gmax,shadow.gmax):raise RuntimeError('Neural shadow shares mutable state')
        selected=np.array([i for i in n.balance.edges if names[c.pre[i]].startswith('PT_BAL_')],dtype=int);shadow.gmax[selected]*=scale;shadow.c.g_max_non[shadow.post,shadow.pre]=shadow.gmax
        inputs=u.copy();inputs[ports]=0.
        if label is not None:inputs[ports[label]]=amplitude
        for _ in range(500):shadow.forward(inputs)
        if not np.all(np.isfinite(shadow.V)):raise RuntimeError('Nonfinite shadow neuron voltage')
        return np.clip(shadow.V[mn]/C.P.E_HI,0.,1.)
    try:
        alpha0=neural(1.);a0=physical(alpha0)
        for scale in (.5,1.,2.,4.,8.):
            for amplitude in (.02,.05):
                for label,axis,sign in [(0,0,1.),(1,0,-1.),(2,1,1.),(3,1,-1.)]:
                    alpha=neural(scale,label,amplitude);acceleration=(physical(alpha)-a0)/amplitude;target=np.zeros(3);target[axis]=sign
                    rows.append(dict(excitatory_gain_multiplier=scale,command=amplitude,label=('XP','XN','YP','YN')[label],actual_static_acceleration_per_command=acceleration.tolist(),desired_acceleration_per_command=target.tolist(),unit_response_squared_error=float(np.sum((acceleration-target)**2)),maximum_activation_change=float(np.max(abs(alpha-alpha0))),neuronal_activation_delta=(alpha-alpha0).tolist()))
    finally:
        for key,value in saved.items():getattr(m,key)[:]=value
    model_restored=all(np.array_equal(getattr(m,key),value) for key,value in saved.items());data_unchanged=all(np.array_equal(getattr(d,key),value) for key,value in actual.items());neural_unchanged=np.array_equal(c.V,voltage) and np.array_equal(c.gmax,gmax)
    if not(model_restored and data_unchanged and neural_unchanged):raise RuntimeError('Shadow transfer audit modified the live state')
    errors={str(scale):float(np.mean([row['unit_response_squared_error'] for row in rows if row['excitatory_gain_multiplier']==scale])) for scale in (.5,1.,2.,4.,8.)}
    result=dict(scope='Local shadow-neural and source-static forward transfer at actual standing pose. No production change, no dynamic perturbation/recovery validation; fixed sensory inputs and frozen contacts/pose.',probe_module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),physical_time_s=float(d.time),native_frame='x anterior,y left,z superior; PT_YP is anatomically left.',source='Bilateral original Thelen2003; steady compliant tendon equilibrium,beta0.',unit_convention='PT command currently enters the SNS as input current directly. No presumed unit-gain map from desired acceleration to actual MN activation.',source_basis_initial_calibration=n.balance.calibration,live_pose_source_static_baseline_acceleration_mps2=a0.tolist(),common_excitation_multiplier_mean_squared_errors=errors,lowest_static_error_multiplier=float(min(errors,key=errors.get)),all_model_arrays_restored=model_restored,live_data_unchanged=data_unchanged,live_neural_arrays_unchanged=neural_unchanged,wall_s=time.perf_counter()-tic,actuators=acts,rows=rows)
    (out/'PT_CURRENT_CHAIN.json').write_text(json.dumps(result,indent=2));return result

def install(out):
    original=C.build
    def build(*args,**kwargs):
        n=original(*args,**kwargs);native=n.step;done=False
        def step(u):
            nonlocal done
            result=native(u)
            if not done and n.step_count>=3500:
                if not hasattr(n,'balance'):raise RuntimeError('PT pathway required')
                probe(n,u,out);done=True
            return result
        n.step=step;return n
    C.build=build

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    if (args.out/'BACKEND_RESULT.json').exists():raise FileExistsError('Preserve completed audit')
    install(args.out);sys.argv=[str(HERE/'live_tendon_trial.py'),'--protocol','hunt-stand','--physics','compliant','--source-task-basis','--config',str(args.config),'--out',str(args.out)];runpy.run_path(str(HERE/'live_tendon_trial.py'),run_name='__main__')
