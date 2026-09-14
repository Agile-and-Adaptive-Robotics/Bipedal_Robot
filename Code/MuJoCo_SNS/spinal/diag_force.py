"""One-frame diagnostic: how do the converted muscle actuators produce
force, and what does mj_inverse say at the leg dofs?"""
import numpy as np
import mujoco
import runner

model, data = runner.build_air_model() if hasattr(runner, "build_air_model") \
    else (None, None)
if model is None:
    m0 = mujoco.MjModel.from_xml_path(str(runner.MODEL))
    d0 = mujoco.MjData(m0)
    mujoco.mj_resetDataKeyframe(m0, d0, 0)
    model = runner.apply_harness(m0, d0, kxy=0, kz=0, ky=0, dxy=1.0, dy=1.0,
                                 no_ground=True, pin_rot=False)
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)

print(f"nu={model.nu} na={model.na} dyntype[0]={model.actuator_dyntype[0]} "
      f"gaintype[0]={model.actuator_gainprm[0]}")
print(f"gaintype id: {int(model.actuator_gaintype[0])}, "
      f"biastype: {int(model.actuator_biastype[0])}")
i_soleus = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, "soleus_r")
i_vas = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, "vas_lat_r")
for nm, i in (("soleus_r", i_soleus), ("vas_lat_r", i_vas)):
    print(f"{nm}: gainprm={model.actuator_gainprm[i]} "
          f"biasprm={model.actuator_biasprm[i]}")

for c in (0.0, 1.0):
    data.ctrl[:] = c
    mujoco.mj_forward(model, data)
    print(f"ctrl={c}: soleus_r force {data.actuator_force[i_soleus]:9.2f} "
          f"vas_lat_r {data.actuator_force[i_vas]:9.2f} "
          f"len soleus {data.actuator_length[i_soleus]:.4f}")
for a in (0.0, 1.0):
    data.act[:] = a
    mujoco.mj_forward(model, data)
    print(f"act={a}: soleus_r force {data.actuator_force[i_soleus]:9.2f} "
          f"vas_lat_r {data.actuator_force[i_vas]:9.2f}")

# inverse dynamics at the keyframe (standing), qvel=qacc=0, no GRF
data.ctrl[:] = 0.0
data.qvel[:] = 0.0
data.qacc[:] = 0.0
mujoco.mj_inverse(model, data)
jn = ("hip_flexion_r", "knee_angle_r", "ankle_angle_r", "lumbar_extension",
      "pelvis_tilt", "pelvis_ty")
for nm in jn:
    jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, nm)
    adr = model.jnt_dofadr[jid]
    print(f"qfrc_inverse[{nm:16s}] = {data.qfrc_inverse[adr]:9.3f}  "
          f"(constraint {data.qfrc_constraint[adr]:8.3f}, "
          f"passive {data.qfrc_passive[adr]:8.3f}, "
          f"bias {data.qfrc_bias[adr]:8.3f})")
