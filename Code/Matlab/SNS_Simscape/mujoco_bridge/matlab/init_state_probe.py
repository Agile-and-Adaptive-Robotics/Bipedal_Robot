"""Inspect initial-state choices for the simbridge model.

The blockset's initData() = mj_makeData only (no mj_resetData), so the plant
starts at qpos = zeros. Compare: qpos0, keyframe qpos, and what zero-qpos
means geometrically (pelvis world height, foot z).
"""
import mujoco
import numpy as np

m = mujoco.MjModel.from_xml_path(
    r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody"
    r"\mjc\gait2392_simbody\gait2392_simbody_cvt3_simbridge.xml")
d = mujoco.MjData(m)

jn = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_JOINT, i) for i in range(m.njnt)]
for i, n in enumerate(jn):
    if n and ("pelvis" in n or n.startswith("knee_angle")):
        adr = m.jnt_qposadr[i]
        print(f"{n:24s} qpos0={float(m.qpos0[adr]):+.3f} key={float(m.key_qpos[0, adr]):+.3f}")

print("nq", m.nq, "nkey", m.nkey)

# where does the torso body sit at qpos=0 vs qpos0?
bid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, "torso")
mujoco.mj_forward(m, d)
print(f"qpos=0    torso z = {d.xpos[bid, 2]:+.3f}  pelvis_z={d.xpos[1, 2]:+.3f}")
d.qpos[:] = m.qpos0
mujoco.mj_forward(m, d)
print(f"qpos=qpos0 torso z = {d.xpos[bid, 2]:+.3f}  pelvis_z={d.xpos[1, 2]:+.3f}")
d.qpos[:] = m.key_qpos[0]
mujoco.mj_forward(m, d)
print(f"qpos=key  torso z = {d.xpos[bid, 2]:+.3f}  pelvis_z={d.xpos[1, 2]:+.3f}")

# 1 s of zero-ctrl stepping from qpos=0: does the zero state survive?
d2 = mujoco.MjData(m)
t, nfe = 0, 0
for _ in range(200):
    mujoco.mj_step(m, d2)
    nfe += d2.ncon
    if not np.isfinite(d2.qpos).all():
        print("zero-state sim NaN at step", _)
        break
else:
    print(f"zero-state 1 s clean: ncon_mean={nfe/200:.0f} "
          f"pelvis_z_end={d2.xpos[1, 2]:+.3f} knee_r={d2.qpos[m.jnt_qposadr[jn.index('knee_angle_r')]]:+.3f}")
