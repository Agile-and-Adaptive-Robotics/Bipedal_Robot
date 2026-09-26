"""M3 probe v6: debug why knee/ankle tendon arms look wrong.
Perturb knee_L +10 deg, mj_forward, print tibia/foot body world pos + quat,
all tendon lengths, and site world positions of the knee chain.
"""
import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import numpy as np
import mujoco

AIR = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_air.xml"
m = mujoco.MjModel.from_xml_path(AIR)
kn = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "knee_L")
print("knee_L type:", m.jnt_type[kn], "qposadr:", m.jnt_qposadr[kn],
      "axis:", m.jnt_axis[kn], "body:", mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_BODY, m.jnt_bodyid[kn]))
for i in range(m.njnt):
    print(mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_JOINT, i),
          "type", m.jnt_type[i], "qadr", m.jnt_qposadr[i], "dofadr", m.jnt_dofadr[i])

for name in ("tibia_L", "foot_L", "toe_L"):
    b = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, name)
    print(name, "body_jntnum", m.body_jntnum[b], "body_jntadr", m.body_jntadr[b], "parent",
          mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_BODY, m.body_parentid[b]))

d = mujoco.MjData(m)
mujoco.mj_forward(m, d)
L0 = d.ten_length.copy()
P0 = {s: d.site_xpos[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_SITE, s)].copy()
      for s in ("site_tibia_L_b", "site_femur_L_b_low", "site_knee_L", "site_tibia_L", "site_foot_L_b")}
d.qpos[m.jnt_qposadr[kn]] += np.radians(10.0)
mujoco.mj_forward(m, d)
for s, p in P0.items():
    pn = d.site_xpos[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_SITE, s)]
    print(f"site {s}: moved {np.linalg.norm(pn - p)*1000:.2f} mm")
for i in range(m.ntendon):
    tn = mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_TENDON, i)
    print(f"tendon {tn:<16} dL = {(d.ten_length[i]-L0[i])*1000:+.3f} mm")
