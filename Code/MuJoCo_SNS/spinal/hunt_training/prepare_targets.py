"""Model-consistent supported-in-place targets, with explicit IK sign adapter.

Retargets measured joint motion to the actual converted geometry. GRF shape
and foot-local CoP come from the source; total mean loading is scaled to a
declared 40% body-weight support condition. This is an engineering training
reference, not a replay of the experimental subject or a new experimental result.
"""
import os,sys,json,hashlib
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parent;S=H.parent;sys.path[:0]=[str(S),str(S/'functional_gait2392')]
import numpy as np,mujoco,opensim as osim
from scipy.signal import savgol_filter
from scipy.optimize import lsq_linear
from mechanics import CoupledArms,JOINTS
import bsolve_ik as B
from motion_io import read_motion

def main():
    out=H/'runs/targets';out.mkdir(parents=True,exist_ok=True)
    t,names,values,unit=read_motion(B.IK_MOT);assert unit=='yes'
    dt=float(np.median(np.diff(t)));model,data=B.build_air_model(False);arms=CoupledArms(model)
    signs={r['coordinate']:r['sign'] for r in json.loads((H/'runs/ik_audit/fk_signs.json').read_text())}
    q=[]
    for row in values:
        for j,adr in zip(JOINTS,arms.addresses):data.qpos[adr]=np.deg2rad(row[names.index(j)])*signs[j]
        arms.project(data);q.append(data.qpos.copy())
    q=np.asarray(q);velocity=savgol_filter(q,11,3,deriv=1,delta=dt,axis=0);acceleration=savgol_filter(q,11,3,deriv=2,delta=dt,axis=0)
    assert model.nq==model.nv
    tg,gn,gv=B.read_mot(B.GRF_MOT)
    foot=[B.grf_at(t,tg,gn,gv,pref=p) for p in ('','1_')]
    mass=float(model.body_mass.sum());bw=mass*9.81
    load_scale=.4*bw/np.mean(foot[0][0][:,2]+foot[1][0][:,2])
    om=osim.Model(str(B.OSIM_DIR/B.OSIM_MODEL));state=om.initSystem();cs=om.getCoordinateSet()
    acts=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(model.nu)]
    Fmax=np.maximum(model.actuator_gainprm[:,2],5.)
    length=[];speed=[];force=[];activation=[];torques=[];residual=[];rhos=[];passive=[];capacity=[]
    jp=np.zeros((3,model.nv));jr=np.zeros_like(jp)
    for k,row in enumerate(values):
        data.qpos[:]=q[k];data.qvel[:]=velocity[k];data.qacc[:]=acceleration[k];data.act[:]=0.
        mujoco.mj_inverse(model,data);tau=data.qfrc_inverse.copy()
        for ci in range(cs.getSize()):
            c=cs.get(ci);name=c.getName()
            if name in names and not c.getLocked(state):c.setValue(state,float(row[names.index(name)]) if name in B.TRANS_DRIVERS else np.deg2rad(row[names.index(name)]))
        om.realizePosition(state)
        for side,(F,P,T) in zip(('r','l'),foot):
            body='calcn_'+side;bid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,body)
            transform=om.getBodySet().get(body).getTransformInGround(state);p=transform.p();r=transform.R().asMat33()
            p_os=np.array([p.get(i) for i in range(3)]);R_os=np.array([[r.get(i,j) for j in range(3)] for i in range(3)])
            local_cop=R_os.T@(np.array([P[k,0],P[k,2],-P[k,1]])-p_os)
            point=data.xpos[bid]+data.xmat[bid].reshape(3,3)@local_cop
            mujoco.mj_jac(model,data,jp,jr,point,bid)
            tau-=load_scale*(jp.T@F[k]+jr.T@T[k])
        rho=arms.update(data).copy();target=arms.E.T@tau
        data.act[:]=0.;mujoco.mj_forward(model,data);f0=-data.actuator_force.copy()
        L=data.actuator_length.copy();V=data.actuator_velocity.copy()
        data.act[:]=1.;mujoco.mj_forward(model,data);unit_force=-data.actuator_force-f0
        A=rho*unit_force[None,:];rhs=target-rho@f0
        norm=np.linalg.norm(A,axis=0);lam=.025*np.median(norm[norm>1e-6])
        sol=lsq_linear(np.vstack([A,lam*np.eye(model.nu)]),np.r_[rhs,np.zeros(model.nu)],bounds=(0.,.98),method='bvls',tol=1e-7,max_iter=120)
        a=sol.x;activation.append(a);length.append(L);speed.append(V);force.append(f0+unit_force*a);torques.append(target);rhos.append(rho)
        passive.append(f0);capacity.append(unit_force)
        residual.append(float(np.linalg.norm(A@a-rhs)/max(np.linalg.norm(rhs),1e-6)))
    activation=np.asarray(activation);length=np.asarray(length);speed=np.asarray(speed);force=np.asarray(force);rhos=np.asarray(rhos)
    # Sensor definitions are fitted to these mechanics and saved for runtime.
    mid=length.mean(axis=0);half=np.maximum(np.ptp(length,axis=0)/2,.003)
    features=[];feature_members={}
    for side,offset in (('r',0),('l',7)):
        source_rho=rhos[:,offset,:]*signs['hip_flexion_'+side]
        masks=[(source_rho.mean(axis=0)>.005),(source_rho.mean(axis=0)<-.005)]
        apf=rhos[:,offset+4,:].mean(axis=0)<-.005;kext=rhos[:,offset+3,:].mean(axis=0)>.005
        indices=[np.flatnonzero(mask) for mask in (masks[0],masks[1],masks[0],masks[1],apf,kext)]
        feature_members[side]=[i.tolist() for i in indices]
        channels=[]
        for f,ids in enumerate(indices):
            raw=(length-mid)/half if f<2 else speed/.6 if f<4 else force/Fmax
            channels.append(np.maximum(raw[:,ids],0.).mean(axis=1))
        features.append(np.stack(channels,axis=1))
    np.savez_compressed(out/'targets.npz',t=t,qpos=q,qvel=velocity,qacc=acceleration,activation=activation,
        length=length,velocity=speed,force=force,Fmax=Fmax,torque=np.asarray(torques),rho=rhos,
        features=np.stack(features,axis=1),act_names=acts,joints=JOINTS,residual=residual,passive_force=passive,unit_force=capacity)
    meta=dict(kind='Supported-in-place engineering reference retargeted from experimental IK',
        source_ik=str(B.IK_MOT),source_ik_sha256=hashlib.sha256(B.IK_MOT.read_bytes()).hexdigest(),
        source_grf=str(B.GRF_MOT),native_signs=signs,mass_kg=mass,mean_contact_fraction_bw=.4,
        grf_scale=load_scale,kinematics='All 17 coordinates retained; pelvis fixed to native initial support pose; measured foot-local CoP transferred without anthropometric scaling',
        inverse='Actual forward model, including its equality followers, armature, and force-velocity effects',
        residual_mean=float(np.mean(residual)),residual_p95=float(np.percentile(residual,95)),
        residual_scope='Projected 17-coordinate inverse torque reconstruction with bounded native muscle activation',
        sensor_mid=mid.tolist(),sensor_half=half.tolist(),feature_members=feature_members,
        feature_names=['hip_flexor_II','hip_extensor_II','hip_flexor_Ia','hip_extensor_Ia','ankle_PF_Ib','knee_EXT_Ib'],
        foot_target_limit='The experimental source holds subtalar/MTP essentially constant; it does not identify their periodic PF patterns')
    (out/'metadata.json').write_text(json.dumps(meta,indent=2));print(json.dumps({k:v for k,v in meta.items() if k not in ('sensor_mid','sensor_half','feature_members')},indent=2))
if __name__=='__main__':main()
