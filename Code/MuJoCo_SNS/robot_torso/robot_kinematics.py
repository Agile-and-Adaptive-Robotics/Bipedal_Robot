"""Constrained knee coordinates in the x-anterior, y-superior, z-right frame."""
import mujoco
import numpy as np

def knee_pose(model,data,side,angle,seed=None):
    joints=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,x+'_'+side+'_joint') for x in ('back','shank','front')]
    addresses=model.jnt_qposadr[joints];signs=model.jnt_axis[joints,2]
    assert np.allclose(np.abs(signs),1.,atol=1e-8)
    eq=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_EQUALITY,'knee_loop_'+side)
    first,second=model.eq_obj1id[eq],model.eq_obj2id[eq]
    def residual(q):
        data.qpos[addresses]=q;mujoco.mj_kinematics(model,data)
        a=data.xpos[first]+data.xmat[first].reshape(3,3)@model.eq_data[eq,:3]
        b=data.xpos[second]+data.xmat[second].reshape(3,3)@model.eq_data[eq,3:6]
        return np.r_[(a-b)[:2],q[:2]@signs[:2]-angle]
    q=np.zeros(3) if seed is None else np.array(seed).copy()
    for _ in range(40):
        error=residual(q)
        if np.linalg.norm(error)<1e-11:
            residual(q);mujoco.mj_forward(model,data);return q
        h=1e-7
        jacobian=np.column_stack([(residual(q+np.eye(3)[i]*h)-residual(q-np.eye(3)[i]*h))/(2*h) for i in range(3)])
        update=np.linalg.solve(jacobian,-error)
        for factor in (1.,.5,.25,.125,.0625,.03125):
            candidate=q+factor*update
            if np.linalg.norm(residual(candidate))<np.linalg.norm(error):q=candidate;break
        else:raise ValueError('Unreachable knee angle %.6f'%angle)
    raise ValueError('Knee solve did not converge')

def pose(model,dynamic,kind,side,angle):
    data=mujoco.MjData(model)
    if kind=='knee':knee_pose(model,data,side,angle)
    else:
        name='thigh_'+side+'_joint' if kind=='hip' else 'lumbar_ball'
        joint=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,name)
        tangent=np.zeros(model.nv);tangent[model.jnt_dofadr[joint]+2]=1.
        mujoco.mj_integratePos(model,data.qpos,tangent,angle)
    dynamic.update(data)
    return data

def camera_axes(position,target,up=np.array([0.,1.,0.])):
    forward=np.array(target)-position;forward/=np.linalg.norm(forward)
    right=np.cross(forward,up);right/=np.linalg.norm(right)
    actual_up=np.cross(right,forward)
    return np.column_stack((right,actual_up,-forward))
