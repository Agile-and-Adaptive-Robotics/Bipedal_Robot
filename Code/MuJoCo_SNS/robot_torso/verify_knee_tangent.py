"""Check that allocation retains all six independent hip rotations."""
import json
from run_neural_robot import *
model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_abduction_diagnostic_v2.xml'));data=mujoco.MjData(model)
for side in ('r','l'):knee_pose(model,data,side,-.35)
records=[]
for perturbation in (0.,.002):
    trial=mujoco.MjData(model);trial.qpos[:]=data.qpos
    for side in ('r','l'):
        jid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'back_'+side+'_joint')
        trial.qpos[model.jnt_qposadr[jid]]+=perturbation
    mujoco.mj_forward(model,trial)
    mask=trial.efc_type==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY)
    jac=trial.efc_J.reshape(trial.nefc,model.nv)[mask,6:]
    _,singular,vh=np.linalg.svd(jac,full_matrices=True);rank=int((singular>1e-8).sum());old=vh[rank:].T
    new=anatomical_knee_tangent(model,trial)
    assert new.shape==(15,11) and np.max(abs(new.T@new-np.eye(11)))<1e-12
    hips=[]
    for side in ('r','l'):
        jid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'thigh_'+side+'_joint');a=model.jnt_dofadr[jid]-6;hips.extend(range(a,a+3))
    assert np.allclose(np.diag(new@new.T)[hips],1.,atol=1e-12)
    if perturbation==0:assert np.max(abs(jac@new))<1e-10
    records.append(dict(hinge_perturbation_rad=perturbation,old_numerical_constraint_rank=rank,
        old_internal_tangent_dofs=old.shape[1],corrected_internal_tangent_dofs=new.shape[1],
        old_hip_projection_diagonal=np.diag(old@old.T)[hips].tolist(),
        corrected_hip_projection_diagonal=np.diag(new@new.T)[hips].tolist()))
(HERE/'knee_tangent_verification.json').write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
