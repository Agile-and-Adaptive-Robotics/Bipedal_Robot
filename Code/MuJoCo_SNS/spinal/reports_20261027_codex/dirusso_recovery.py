"""Corrected gait2392 experiment with compiled neural updates.

Converted-model adaptation, not a verified reproduction of the SCONE study.
No linear-motor replacement, replay forcing, or reference kinematics imposed
after release. All controls and support forces are logged.
"""
from __future__ import annotations
import os
import sys
import json
from pathlib import Path
import time
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
SPINAL = HERE.parent
os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)
for var in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[var] = '1'
import numpy as np
import mujoco
sys.path.insert(0,str(SPINAL/'dirusso_gait2392'))
import controller as C
import plant as P
import dirusso_run as D
from controller_numba import CompiledDiRusso

MODEL = None
REF = None
IMITATION = None


def model_and_reference():
    global MODEL, REF
    if MODEL is None:
        sys.path.insert(0,str(SPINAL))
        import runner
        raw=Path(D.XML).read_text(encoding='utf-8')
        original=mujoco.MjModel.from_xml_path(D.XML)
        initial=mujoco.MjData(original)
        mujoco.mj_resetDataKeyframe(original,initial,0)
        repaired=runner.patch_xml(raw,original,initial,kxy=0.,kz=0.,ky=0.,
                                  dxy=0.,dy=0.,rig_scale=0.)
        txt, _ = P.build_sagittal_xml(D.XML,xml_text=repaired)
        root=ET.fromstring(txt)
        # Match the paper's 3+3+3 independent sagittal coordinates.
        eq=root.find('equality')
        for side in 'rl':
            ET.SubElement(eq,'joint',name='lock_mtp_'+side,joint1='mtp_angle_'+side)
        for mesh in root.findall('./asset/mesh'):
            mesh.set('file',str(Path(D.XML).parent/mesh.get('file')))
        contact=root.find('contact')
        for pair in list(contact):
            if not any(x in pair.get('geom2','') for x in ('calcn_','toes_')):
                contact.remove(pair)
        MODEL=mujoco.MjModel.from_xml_string(ET.tostring(root,encoding='unicode'))
        MODEL.opt.timestep=.002
        MODEL.opt.integrator=mujoco.mjtIntegrator.mjINT_IMPLICITFAST
        REF=dict(np.load(D.BSOLVE,allow_pickle=True))
    return MODEL, REF


def project_followers(model,data):
    for i in range(model.neq):
        if model.eq_type[i]!=mujoco.mjtEq.mjEQ_JOINT:
            continue
        dep,ind=int(model.eq_obj1id[i]),int(model.eq_obj2id[i])
        if ind<0:
            continue
        da,ia=model.jnt_qposadr[dep],model.jnt_qposadr[ind]
        x=data.qpos[ia]-model.qpos0[ia]
        coefficients=model.eq_data[i,:5]
        data.qpos[da]=model.qpos0[da]+sum(coefficients[k]*x**k for k in range(5))


def imitation_weights():
    """Bounded equation-8 primitive fit, shared by the two sides."""
    global IMITATION
    if IMITATION is not None:
        return IMITATION
    from scipy.optimize import lsq_linear
    model,ref=model_and_reference()
    _,_,strikes=D.grf_onsets()
    t=ref['t']; ai={str(n):i for i,n in enumerate(ref['act_names'])}
    patterns=[]; targets=[]
    for side,s in (('L','l'),('R','r')):
        ons=strikes[side]
        period=np.median(np.diff(ons))
        phase=((t-ons[0])/period)%1.
        patterns.append(np.stack([C.primitives(p) for p in phase]))
        targets.append(np.stack([ref['acts'][:,[ai[n.format(s=s)] for n in C.MAP[m]]].mean(1)
                                 for m in C.MUSCLES],axis=1))
    x,y=np.concatenate(patterns),np.concatenate(targets)
    w=np.stack([lsq_linear(x,y[:,j],bounds=(-1.,1.)).x for j in range(9)])
    IMITATION=(np.concatenate([w,w]),dict(rmse=float(np.sqrt(np.mean((x@w.T-y)**2)))))
    return IMITATION


def simulate(params=None,duration=10.,support=0.,save=None):
    started=time.perf_counter()
    params=params or {}
    model,ref=model_and_reference()
    data=mujoco.MjData(model)
    ids,geoms,pelvis,tiltadr,txadr=P.index_maps(model,mujoco)
    group_ids=[ids[(m,s)] for s in ('L','R') for m in C.MUSCLES]
    assert all(min(group)>=0 for group in group_ids)
    gains={k:(v if np.isscalar(v) else v*.1) for k,v in C.default_gains().items()}
    for k in gains:
        if isinstance(gains[k],np.ndarray):
            gains[k]=np.array([params.get(f'{k}_{m}',gains[k][j]) for j,m in enumerate(C.MUSCLES)])
    W=np.array(params.get('W',imitation_weights()[0]))
    for j,m in enumerate(C.MUSCLES):
        W[j]*=params.get('ff_'+m,1.)
        W[j+9]*=params.get('ff_'+m,1.)
    index=int(params.get('initial_frame',0))%(len(ref['t'])-1)
    data.qpos[:]=ref['qpos'][index]
    for name in P.LOCK_JOINTS+['mtp_angle_r','mtp_angle_l']:
        jid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,name)
        data.qpos[model.jnt_qposadr[jid]]=model.qpos0[model.jnt_qposadr[jid]]
    project_followers(model,data)
    data.qpos[1]+=params.get('initial_height_offset',0.)
    data.qpos[tiltadr]+=params.get('initial_pitch_offset',0.)
    dtref=ref['t'][index+1]-ref['t'][index]
    mujoco.mj_differentiatePos(model,data.qvel,dtref,ref['qpos'][index],ref['qpos'][index+1])
    for name in P.LOCK_JOINTS+['mtp_angle_r','mtp_angle_l']:
        jid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,name)
        data.qvel[model.jnt_dofadr[jid]]=0.
    if 'initial_vx' in params:
        data.qvel[0]=params['initial_vx']
    data.qvel[1]+=params.get('initial_vz_offset',0.)
    for side in 'rl':
        for joint in ('hip_flexion','knee_angle','ankle_angle'):
            name=joint+'_'+side
            jid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,name)
            data.qpos[model.jnt_qposadr[jid]]+=params.get('initial_q_'+name,0.)
            data.qvel[model.jnt_dofadr[jid]]+=params.get('initial_v_'+name,0.)
    project_followers(model,data)
    data.act[:]=ref['acts'][index]
    mujoco.mj_forward(model,data)
    h0=float(data.subtree_com[0,2])
    weight=float(model.body_mass.sum()*-model.opt.gravity[2])
    ctrl=CompiledDiRusso(W,omega=params.get('omega',2*np.pi/1.23),
                        gamma=params.get('gamma',1.),gains=gains,
                        balance_gains={key:np.array([params.get(
                            f'balance_{key}_{muscle}',params.get('balance_'+key,default))
                            for muscle in C.BALANCE_MUSCLES])
                            for key,default in (('kp',8.),('kd',1.5),('target',0.))})
    _,_,onsets=D.grf_onsets()
    for si,side in enumerate(('L','R')):
        period=np.median(np.diff(onsets[side]))
        ctrl.phi[si]=((ref['t'][index]-onsets[side][0])/period)%1.
    n=int(duration/.002)
    qids=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,j+'_'+s)
          for s in 'rl' for j in ('hip_flexion','knee_angle','ankle_angle')]
    qadr=model.jnt_qposadr[qids]
    lengths=np.stack([model.actuator_lengthrange[g].mean(0) for g in group_ids])
    span=np.maximum(lengths[:,1]-lengths[:,0],1e-6)
    capacity=np.array([model.actuator_gainprm[g,2].sum() for g in group_ids])
    f6=np.zeros(6); sensors=np.zeros((2,4,9)); heel_prev=np.zeros(2)
    strikes=np.zeros(2,dtype=np.bool_)
    log=dict(t=np.zeros(n),q=np.zeros((n,6)),com=np.zeros((n,3)),
             heel=np.zeros((n,2)),ctrl=np.zeros((n,18)),support=np.zeros(n),
             phi=np.zeros((n,2)),tilt=np.zeros(n))
    fell=None
    numerical_failure=False
    warning_counts=np.zeros(8,dtype=int)
    for i in range(n):
        heel=np.zeros(2); toe=np.zeros(2)
        for k in range(data.ncon):
            c=data.contact[k]
            mujoco.mj_contactForce(model,data,k,f6)
            # MuJoCo contact-frame x is normal; transform to world vertical.
            fz=abs(c.frame[2]*f6[0]+c.frame[5]*f6[1]+c.frame[8]*f6[2])
            for si,s in enumerate(('L','R')):
                if geoms[('heel',s)] in (c.geom1,c.geom2):
                    heel[si]+=fz
                if geoms[('toe',s)] in (c.geom1,c.geom2):
                    toe[si]+=fz
        length=np.array([data.actuator_length[g].mean() for g in group_ids])
        velocity=np.array([data.actuator_velocity[g].mean() for g in group_ids])
        force=np.array([np.abs(data.actuator_force[g]).sum() for g in group_ids])
        ln=np.clip((length-lengths[:,0])/span,0.,2.).reshape(2,9)
        vn=np.clip(velocity/span,-5.,5.).reshape(2,9)
        fn=np.clip(force/np.maximum(capacity,1e-6),0.,2.).reshape(2,9)
        sensors[:,0]=(65./200.)*np.sqrt(np.maximum(vn,0.))
        sensors[:,1]=ln; sensors[:,2]=fn
        sensors[:,3]=np.clip((heel+toe)/weight,0.,1.)[:,None]
        strikes[:]=(heel_prev<=20.)&(heel>20.)
        heel_prev[:]=heel
        exc=ctrl.step_array(sensors,strikes,float(data.qpos[tiltadr]))
        data.ctrl[:]=0.
        for j,g in enumerate(group_ids):
            data.ctrl[g]=exc.reshape(-1)[j]
        data.xfrc_applied[:]=0.
        data.xfrc_applied[pelvis,2]=support*weight
        log['t'][i]=i*.002; log['q'][i]=data.qpos[qadr]
        log['com'][i]=data.subtree_com[0]; log['heel'][i]=heel
        log['ctrl'][i]=exc.reshape(-1); log['support'][i]=support*weight
        log['phi'][i]=ctrl.phi; log['tilt'][i]=data.qpos[tiltadr]
        if (not np.isfinite(data.qpos).all()) or data.subtree_com[0,2]<.8*h0:
            fell=i*.002
            break
        mujoco.mj_step(model,data)
        warning_counts=np.array(data.warning.number,dtype=int)
        if (np.any(warning_counts[[0,4,5,6]]>0) or
                not np.isfinite(data.qpos).all() or
                not np.isfinite(data.qvel).all()):
            numerical_failure=True
            fell=(i+1)*.002
            break
    log={k:v[:i+1] for k,v in log.items()}
    hits=[]
    for side in range(2):
        on=log['heel'][:,side]>20.
        raw=np.flatnonzero(on[1:]&~on[:-1])+1
        primary=[]
        for index in raw:
            if not primary or log['t'][index]-log['t'][primary[-1]]>.45:
                primary.append(int(index))
        hits.append(primary)
    phases=[]
    for a,b in zip(hits[1][:-1],hits[1][1:]):
        for left in hits[0]:
            if a<=left<b:
                phases.append((log['t'][left]-log['t'][a])/(log['t'][b]-log['t'][a]))
    elapsed=max(float(log['t'][-1]),.002)
    speed=float((log['com'][-1,0]-log['com'][0,0])/elapsed)
    phase=float(np.median(phases)) if phases else None
    # Continuous duration penalty gives a gradient among falling trials.
    cost=100.*(duration-elapsed)/duration+10.*max(1.10-speed,0.)+10.*max(speed-1.25,0.)
    cost+=.01*float(np.mean(log['ctrl']**2))
    if not np.isfinite(cost):
        cost=1e6
    if numerical_failure:
        cost=1e6
    result=dict(fell=fell is not None,fall_t=fell,elapsed=elapsed,
                cycL=len(hits[0]),cycR=len(hits[1]),phase=phase,speed=speed,
                cost=cost,wall_s=time.perf_counter()-started,
                support_fraction=support,initial_com_height=h0,
                numerical_failure=numerical_failure,warning_counts=warning_counts.tolist(),
                survive=(fell is None and min(map(len,hits))>=3 and
                         phase is not None and .35<=phase<=.65))
    if save:
        np.savez_compressed(HERE/(save+'.npz'),**log)
        (HERE/(save+'.json')).write_text(json.dumps(dict(params=params,metrics=result),indent=2))
    return result


if __name__=='__main__':
    for support in (0.,.7):
        result=simulate(support=support,save=f'dirusso_corrected_ws{round(100*support)}')
        print(json.dumps(result),flush=True)
