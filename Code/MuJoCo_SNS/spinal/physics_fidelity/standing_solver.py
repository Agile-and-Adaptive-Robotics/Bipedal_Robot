"""Opt-in bounded equality-aware standing initializer; no imposed net preload.

Static fixed-pose approximation with frozen zero-activation contact forces.
Returns nu activations. Does not integrate or establish balance/walking.
"""
import numpy as np
import mujoco
from scipy.optimize import lsq_linear


def solve_standing_activations(model,data,Fmax=None,preload=False,*,bilateral=True,return_diagnostics=False):
    """Derive actual active capacity and passive tension from the current model.

    Fmax is accepted for runner-call compatibility but is not used. Preload
    must be False: the historical artificial torque preload is not reproduced.
    Current qfrc_applied/xfrc_applied are included; no support force is invented.
    Incoming activation/velocity are restored, then caches are recomputed.
    Optional diagnostics describe the frozen-contact solve, not a static proof.
    """
    if preload:raise ValueError('Corrected initializer has no artificial net-torque preload; pass preload=False')
    if model.na!=model.nu:raise ValueError('Prototype requires one activation state per muscle actuator')
    from mechanics import CoupledArms,JOINTS
    from build_network import _group_weight
    from muscle_map import classify
    from params import POSTURE_OVERRIDE,W_POSTURE
    saved_act=data.act.copy();saved_vel=data.qvel.copy()
    try:
        data.qvel[:]=0.;data.act[:]=0.;mujoco.mj_forward(model,data)
        ar=CoupledArms(model);ar.update(data);E=ar.E.copy();moment=data.actuator_moment.copy()
        passive=data.actuator_force.copy();external=data.qfrc_applied.copy()
        for body in range(1,model.nbody):
            wrench=data.xfrc_applied[body]
            if np.any(wrench):mujoco.mj_applyFT(model,data,wrench[:3],wrench[3:],data.xipos[body],body,external)
        demand=data.qfrc_bias-data.qfrc_constraint-data.qfrc_passive-external
        target=E.T@(demand-moment.T@passive)
        data.act[:]=1.;mujoco.mj_forward(model,data)
        capacity=-(data.actuator_force-passive)
        if np.min(capacity)<-1e-7:raise ValueError('Unexpected non-pulling muscle force response')
        A=E.T@(-moment.T*np.maximum(capacity,0.)[None,:])
        names=[model.actuator(i).name for i in range(model.nu)]
        prior=np.array([min(.5,.6*_group_weight(classify(a),W_POSTURE,POSTURE_OVERRIDE)) for a in names])
        groups=[];seen=set()
        for i,name in enumerate(names):
            if name in seen:continue
            other=name[:-2]+('_l' if name.endswith('_r') else '_r') if name.endswith(('_r','_l')) else name
            indices=[i]+([names.index(other)] if other!=name and other in names else [])
            for j in indices:seen.add(names[j])
            groups.append(indices)
        S=np.zeros((model.nu,len(groups))) if bilateral else np.eye(model.nu)
        if bilateral:
            for k,ids in enumerate(groups):S[ids,k]=1.
        fit=lsq_linear(np.vstack([.01*A@S,.25*S]),np.r_[.01*target,.25*prior],bounds=(0.,1.),method='bvls')
        if not fit.success:raise RuntimeError('Standing bounded solve did not converge: '+fit.message)
        activation=S@fit.x;residual=A@activation-target
        diagnostics=dict(scope='Fixed-pose, zero-velocity active-capacity/passive-force fit with frozen zero-activation contacts; not proof of balance.',
            bilateral=bool(bilateral),free_parameters=int(S.shape[1]),actuators=model.nu,
            frozen_constraint_residual_L2_Nm=float(np.linalg.norm(residual)),residual_by_joint=dict(zip(JOINTS,residual.tolist())),
            external_projected_force_Nm=(E.T@external).tolist(),incoming_external_wrench_present=bool(np.any(data.xfrc_applied)),
            incoming_generalized_external_force_present=bool(np.any(data.qfrc_applied)),solver_optimality=float(fit.optimality))
    finally:
        data.act[:]=saved_act;data.qvel[:]=saved_vel;mujoco.mj_forward(model,data)
    return (activation,diagnostics) if return_diagnostics else activation
