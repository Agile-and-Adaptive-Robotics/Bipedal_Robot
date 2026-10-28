"""Measure Gait2392 routing across IK and hip flexion before PF assignment."""
import os,sys,json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;SPINAL=HERE.parent;sys.path[:0]=[str(SPINAL),str(HERE)]
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
import bsolve_ik as B
from mechanics import CoupledArms,JOINTS,BI_FAMILIES
OUT=SPINAL/'reports_20261028_functional_pf';OUT.mkdir(exist_ok=True)
model,data=B.build_air_model(id_mode=False);arms=CoupledArms(model)
source=np.load(SPINAL/'bsolve_out.npz',allow_pickle=True)
assert source['qpos'].shape[1]==model.nq
names=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(model.nu)]
assert names==source['act_names'].tolist()
samples=[];errors=[]
for index in range(0,len(source['t']),4):
    data.qpos[:]=source['qpos'][index];arms.project(data);mujoco.mj_forward(model,data)
    measured=arms.update(data).copy();samples.append(measured)
    if index in (0,40,80):errors.append(float(np.max(abs(measured-arms.finite_difference(data)))))
assert max(errors)<2e-6,errors
samples=np.array(samples)
tg,gn,gv=B.read_mot(B.GRF_MOT)
phase_weights={}
temporal_primitives={};phase_metadata={}
for side,prefix in (('r',''),('l','1_')):
    force,_,_=B.grf_at(source['t'],tg,gn,gv,pref=prefix)
    stance=force[:,2]>50.
    onsets=source['t'][np.flatnonzero(np.diff(stance.astype(int))==1)+1]
    period=float(np.median(np.diff(onsets))) if len(onsets)>=2 else float(np.median(np.diff(source['hs_r'])))
    origin=float(onsets[0]) if len(onsets) else float(source['t'][0])
    phase=((source['t']-origin)/period)%1.;order=np.argsort(phase)
    phase_metadata[side]=dict(period_s=period,duty=float(stance.mean()),heel_strike_origin_s=origin)
    from scipy.ndimage import gaussian_filter1d
    for family,members in BI_FAMILIES.items():
        ids=[i for i,n in enumerate(names) if n.rsplit('_',1)[0] in members and n.endswith('_'+side)]
        target=source['acts'][:,ids].mean(axis=1)
        phi=np.r_[phase[order]-1,phase[order],phase[order]+1];values=np.tile(target[order],3)
        wave=np.interp(np.linspace(0,1,101),phi,values)
        wave=gaussian_filter1d(wave,2.,mode='wrap');wave[-1]=wave[0]
        temporal_primitives[family+'_'+side]=wave.tolist()
    for ai,name in enumerate(names):
        if not name.endswith('_'+side):continue
        e=float(source['acts'][stance,ai].mean());f=float(source['acts'][~stance,ai].mean())
        phase_weights[name]=dict(stance=e/(e+f) if e+f>1e-6 else .5,
                                swing=f/(e+f) if e+f>1e-6 else .5,
                                fallback=bool(e+f<=1e-6))
sweep=[]
base=source['qpos'][0].copy()
angles=np.arange(-20,101,10)
for angle in angles:
    data.qpos[:]=base
    for side in ('r','l'):
        data.qpos[arms.addresses[JOINTS.index('hip_flexion_'+side)]]=np.deg2rad(angle)
        data.qpos[arms.addresses[JOINTS.index('hip_adduction_'+side)]]=0.
        data.qpos[arms.addresses[JOINTS.index('hip_rotation_'+side)]]=0.
        data.qpos[arms.addresses[JOINTS.index('knee_angle_'+side)]]=np.deg2rad(-20.)
    arms.project(data);mujoco.mj_forward(model,data);sweep.append(arms.update(data).copy())
sweep=np.array(sweep);mapping={};crossovers=[]
for ai,name in enumerate(names):
    side=name[-1];base_name=name.rsplit('_',1)[0];leg=0 if side=='r' else 7
    maxarm=np.max(abs(samples[:,:,ai]),axis=0)
    hip=np.linalg.norm(maxarm[leg:leg+3]);knee=maxarm[leg+3];ankle=np.linalg.norm(maxarm[leg+4:leg+7])
    family=next((key for key,members in BI_FAMILIES.items() if base_name in members),None)
    confirmed=bool(family is not None and ((hip>.002 and knee>.002) or (knee>.002 and ankle>.002)))
    mapping[name]=dict(family=family if confirmed else None,hip_arm_norm_m=float(hip),knee_arm_m=float(knee),distal_arm_norm_m=float(ankle),confirmed_biarticular=confirmed)
    values=sweep[:,leg+2,ai]
    if values.min()<-.002 and values.max()>.002:crossovers.append(dict(actuator=name,min_rotation_arm_m=float(values.min()),max_rotation_arm_m=float(values.max())))
np.savez_compressed(OUT/'routing_samples.npz',actuator_names=names,joint_names=JOINTS,
    ik_moment_arms=samples,hip_flexion_deg=angles,sweep_moment_arms=sweep)
report=dict(machine='easteregg2',model=str(B.runner.MODEL),actuator_names=names,joint_names=JOINTS,phase_weights=phase_weights,
    temporal_primitives=temporal_primitives,phase_metadata=phase_metadata,
    finite_difference_max_errors_m=errors,hip_rotation_sign_crossovers=crossovers,actuator_mapping=mapping,
    frame='OpenSim anatomical joint coordinates; legacy MuJoCo vector adapter (x,y,z)_OS -> (x,-z,y)_MJ',
    force_model='Existing converted Gait2392 muscle model; no CAD or reference-tension robot substitution',
    scope='Measured geometric action; sweep -20 to 100 degrees hip flexion with knee -20 degrees, hip adduction/rotation zero. Not an anatomical validation of every route.',
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'mechanics.py',SPINAL/'bsolve_out.npz',B.runner.MODEL)})
(OUT/'routing_audit.json').write_text(json.dumps(report,indent=2))
print(json.dumps(dict(confirmed_biarticular=[n for n,v in mapping.items() if v['confirmed_biarticular']],rotation_crossovers=crossovers,max_fd_error=max(errors)),indent=2))
