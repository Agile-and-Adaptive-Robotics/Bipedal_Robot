"""Isolated root-inclusive static equilibrium with unknown regional normal loads.

No model flags, caller state, force profiles, or existing solver are changed.
Linear contact-force feasibility is separately checked against actual forward
dynamics at the same pose; optimized reactions are never applied as fake forces.
"""
import numpy as np
import mujoco
from scipy.optimize import lsq_linear


def _scratch(model,data):
    d=mujoco.MjData(model);d.qpos[:]=data.qpos;d.qvel[:]=0.;d.act[:]=0.
    d.qfrc_applied[:]=data.qfrc_applied;d.xfrc_applied[:]=data.xfrc_applied
    return d


def _external(model,d):
    result=d.qfrc_applied.copy()
    for body in range(1,model.nbody):
        w=d.xfrc_applied[body]
        if np.any(w):mujoco.mj_applyFT(model,d,w[:3],w[3:],d.xipos[body],body,result)
    return result


def _support_regions(model,d):
    """Four regional lower-surface points, preserving actual body ownership."""
    pg=model.geom('ground-plane').id;origin=d.geom_xpos[pg];normal=d.geom_xmat[pg].reshape(3,3)[:,2].copy()
    if normal[2]<0:normal=-normal
    regions=[]
    for side in ('r','l'):
        heel=model.body('calcn_'+side).id;toe=model.body('toes_'+side).id
        start=d.xpos[heel];direction=d.xpos[toe]-start;den=max(float(direction@direction),1e-12)
        candidates=[[],[]]
        for g in range(model.ngeom):
            body=int(model.geom_bodyid[g])
            if body not in (heel,toe) or model.geom_type[g]!=mujoco.mjtGeom.mjGEOM_MESH:continue
            mesh=model.geom_dataid[g];first=model.mesh_vertadr[mesh];count=model.mesh_vertnum[mesh]
            vertices=model.mesh_vert[first:first+count]@d.geom_xmat[g].reshape(3,3).T+d.geom_xpos[g]
            fraction=(vertices-start)@direction/den
            for region in (0,1):
                keep=np.ones(count,dtype=bool) if body==toe and region==1 else np.zeros(count,dtype=bool) if body==toe else (fraction<.5 if region==0 else fraction>=.5)
                if not np.any(keep):continue
                points=vertices[keep];gaps=(points-origin)@normal;k=int(np.argmin(gaps))
                candidates[region].append((float(gaps[k]),points[k].copy(),body))
        for region in (0,1):
            if not candidates[region]:raise ValueError('No foot surface region found')
            gap,point,body=min(candidates[region],key=lambda x:x[0])
            regions.append(dict(label=('hindfoot_' if region==0 else 'forefoot_')+side,gap_m=gap,point=point,body=body))
    return regions,normal


def solve_contact_standing_activations(model,data,Fmax=None,preload=False,*,bilateral=True,support_fraction=None,
        contact_gap_tolerance_m=.003,prior_weight=.001,return_diagnostics=False):
    """Solve 23 independent DOF equilibrium, activations and four regional loads.

    Unknown normal loads are nonnegative and only allowed within the declared
    surface gap. They are optimization variables, NOT applied virtual forces.
    A supplied support_fraction replaces the scratch pelvis upward force; None
    means use existing data forces only. Caller model/data remain untouched.
    This static normal-only approximation omits unknown tangential friction.
    """
    if preload:raise ValueError('No artificial net-torque preload')
    if model.na!=model.nu:raise ValueError('Requires one activation state per muscle')
    from mechanics import CoupledArms,fill_projection
    from build_network import _group_weight
    from muscle_map import classify
    from params import POSTURE_OVERRIDE,W_POSTURE
    d=_scratch(model,data);ar=CoupledArms(model);ar.project(d)
    BW=float(model.body_mass.sum()*np.linalg.norm(model.opt.gravity))
    if BW<=0:raise ValueError('Gravity/body weight required')
    if support_fraction is not None:
        if not 0.<=support_fraction<=1.:raise ValueError('Support fraction must be in [0,1]')
        d.xfrc_applied[model.body('pelvis').id,2]=support_fraction*BW
    mujoco.mj_forward(model,d)
    dependent=set(ar.dependent.tolist());drivers=np.array([i for i in range(model.nv) if i not in dependent],dtype=np.int64)
    E=np.zeros((model.nv,len(drivers)))
    fill_projection(E,drivers,ar.dependent,ar.independent,ar.coefficients,ar.offsets,d.qpos,ar.qadr)
    moment=d.actuator_moment.copy();passive=d.actuator_force.copy()
    b=E.T@(d.qfrc_bias-d.qfrc_passive-_external(model,d)-moment.T@passive)
    regions,normal=_support_regions(model,d);C=np.zeros((model.nv,4))
    # C[:,k] is strided: use a contiguous MuJoCo output buffer.
    for k,r in enumerate(regions):
        column=np.zeros(model.nv);mujoco.mj_applyFT(model,d,normal,np.zeros(3),r['point'],r['body'],column);C[:,k]=column
    allowed=np.array([r['gap_m']<=contact_gap_tolerance_m for r in regions])
    d.act[:]=1.;mujoco.mj_forward(model,d);capacity=-(d.actuator_force-passive)
    if np.min(capacity)<-1e-7:raise ValueError('Non-pulling activation response')
    A=E.T@(-moment.T*np.maximum(capacity,0.)[None,:])
    names=[model.actuator(i).name for i in range(model.nu)];seen=set();groups=[]
    for i,name in enumerate(names):
        if name in seen:continue
        other=name[:-2]+('_l' if name.endswith('_r') else '_r') if name.endswith(('_r','_l')) else name
        ids=[i]+([names.index(other)] if other!=name and other in names else [])
        groups.append(ids);seen.update(names[j] for j in ids)
    S=np.zeros((model.nu,len(groups))) if bilateral else np.eye(model.nu)
    if bilateral:
        for k,ids in enumerate(groups):S[ids,k]=1.
    prior=np.array([min(.5,.6*_group_weight(classify(n),W_POSTURE,POSTURE_OVERRIDE)) for n in names])
    contact=E.T@C[:,allowed]*BW;design=np.column_stack([A@S,contact]);nact=S.shape[1];nc=int(allowed.sum())
    scales=[];driver_names=[]
    for i in drivers:
        j=int(model.dof_jntid[i]);kind=model.jnt_type[j];driver_names.append(model.joint(j).name)
        if kind not in (mujoco.mjtJoint.mjJNT_SLIDE,mujoco.mjtJoint.mjJNT_HINGE):raise ValueError('Scalar-joint model required')
        scales.append(1./(BW if kind==mujoco.mjtJoint.mjJNT_SLIDE else BW*.5))
    scales=np.asarray(scales);reg=np.zeros((model.nu+nc,nact+nc));reg[:model.nu,:nact]=prior_weight*S
    if nc:reg[model.nu:,nact:]=np.eye(nc)*prior_weight*.1
    remaining=max(0.,1.-float(support_fraction or 0.));rhs=np.r_[prior_weight*prior,np.full(nc,prior_weight*.1*remaining/max(nc,1))]
    fit=lsq_linear(np.vstack([scales[:,None]*design,reg]),np.r_[scales*b,rhs],bounds=(np.zeros(nact+nc),np.r_[np.ones(nact),np.full(nc,2.)]),method='bvls')
    if not fit.success:raise RuntimeError('Static contact solve failed: '+fit.message)
    x=S@fit.x[:nact];force=np.zeros(4);force[allowed]=fit.x[nact:]*BW
    residual=design@fit.x-b
    # Independently evaluate real contact solver. Never apply optimized reactions.
    d.act[:]=x;d.ctrl[:]=x;d.qvel[:]=0.;mujoco.mj_forward(model,d)
    actual=E.T@(d.qfrc_actuator+d.qfrc_passive+_external(model,d)+d.qfrc_constraint-d.qfrc_bias)
    roots=np.array([k for k,n in enumerate(driver_names) if n.startswith('pelvis_')]);internal=np.array([k for k,n in enumerate(driver_names) if not n.startswith('pelvis_')])
    diagnostics=dict(scope='Root-inclusive static normal-force optimization; real forward contact forces evaluated separately. No virtual reactions applied; no holding or walking claim.',
        independent_coordinates=driver_names,free_activation_parameters=nact,contact_force_N=force.tolist(),contact_force_bw=(force/BW).tolist(),
        regions=[dict(label=r['label'],gap_m=r['gap_m'],position_m=r['point'].tolist(),force_allowed=bool(a)) for r,a in zip(regions,allowed)],
        normalized_linear_residual=float(np.linalg.norm(scales*residual)),linear_residual_by_coordinate=dict(zip(driver_names,residual.tolist())),
        actual_normalized_forward_residual=float(np.linalg.norm(scales*actual)),actual_projected_net_force=dict(zip(driver_names,actual.tolist())),
        root_coordinate_accelerations={driver_names[k]:float(d.qacc[drivers[k]]) for k in roots},
        max_internal_coordinate_acceleration=float(np.max(abs(d.qacc[drivers[internal]]))),activation_max=float(x.max()),
        contacts_in_forward=int(d.ncon),contact_gap_tolerance_m=contact_gap_tolerance_m,prior_weight=prior_weight,
        body_weight_N=BW,support_fraction_argument=support_fraction,model_disableflags=int(model.opt.disableflags))
    return (x,diagnostics) if return_diagnostics else x
