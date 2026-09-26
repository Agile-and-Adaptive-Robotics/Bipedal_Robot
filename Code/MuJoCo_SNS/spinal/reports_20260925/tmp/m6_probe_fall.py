"""M6 probe: what happens in the first 0.6 s of the free-pelvis ground run?
Logs joints, pelvis z, plate forces, limit violations, qacc, contacts.
Run from w2l_mujoco with the myo env."""
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
HERE = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
sys.path.insert(0, HERE)

import numpy as np
import mujoco
import test_w2l_air as TW

m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "w2l_ground.xml"))
d = mujoco.MjData(m)
m.dof_damping[:] = 3.0
m.jnt_solimp[:] = np.array([0.9, 0.99, 0.001, 0.5, 2.0])
m.jnt_solref[:] = np.array([0.006, 1.0])
jadr = m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "root")]
d.qpos[jadr + 2] -= 0.026

names = ("hip_L", "knee_L", "ankle_L", "toe_L", "hip_R", "knee_R", "ankle_R", "toe_R")
qadr = {n: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, n)] for n in names}
root = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, "Root")
plates = {}
for g in ("foot_L_contact", "foot_R_contact", "toe_L_contact", "toe_R_contact"):
    plates[g] = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g)
c6 = np.zeros(6)

# ground clearance of each contact geom at t=0 (after drop)
mujoco.mj_forward(m, d)
print("== t=0 state after drop ==")
for g, gi in plates.items():
    print(f"  {g}: geom zpos {d.geom_xpos[gi][2]:+.4f}")
print(f"  pelvis z {d.xpos[root][2]:.4f}  COM z {d.subtree_com[0][2]:.4f}")
for n in names:
    jid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, n)
    lo, hi = m.jnt_range[jid]
    q = d.qpos[qadr[n]]
    flag = "  <-- OUTSIDE RANGE" if (m.jnt_limited[jid] and (q < lo - 1e-9 or q > hi + 1e-9)) else ""
    print(f"  {n}: q {np.degrees(q):+7.2f} deg  range [{np.degrees(lo):+.1f},{np.degrees(hi):+.1f}]{flag}")

print("\n== first 0.6 s, every 20 ms ==")
print("  t     pelZ   qaccN   hipL   kneeL  anklL  toeL   heelL  toeL   ncon")
for i in range(600):
    mujoco.mj_step(m, d)
    if i % 20 == 0:
        f = {"foot_L_contact": 0.0, "toe_L_contact": 0.0}
        for k in range(d.ncon):
            con = d.contact[k]
            for nm, gi in plates.items():
                if con.geom1 == gi or con.geom2 == gi:
                    mujoco.mj_contactForce(m, d, k, c6)
                    if nm in f:
                        f[nm] += abs(c6[0])
        print(f"  {i*0.001:5.2f} {d.xpos[root][2]:6.3f} {np.abs(d.qacc).max():7.1f} "
              f"{np.degrees(d.qpos[qadr['hip_L']]):+6.1f} {np.degrees(d.qpos[qadr['knee_L']]):+6.1f} "
              f"{np.degrees(d.qpos[qadr['ankle_L']]):+6.1f} {np.degrees(d.qpos[qadr['toe_L']]):+6.1f} "
              f"{f['foot_L_contact']:6.1f} {f['toe_L_contact']:6.1f} {d.ncon:4d}")
