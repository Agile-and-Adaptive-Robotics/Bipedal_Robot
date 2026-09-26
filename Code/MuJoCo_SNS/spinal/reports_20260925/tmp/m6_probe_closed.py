"""M6 probe: closed-loop latch stand, first 1.5 s — which joint folds first,
which muscles pull, are the contact encoders alive?"""
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
HERE = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
sys.path.insert(0, HERE)

import numpy as np
import mujoco
import build_w2l_aff_net as B

B.SPLIT.NAP["tau_max_h"] = 0.25
net = B.build(comm=1.0)

m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "w2l_ground.xml"))
d = mujoco.MjData(m)
m.dof_damping[:] = 3.0
m.jnt_solimp[:] = np.array([0.9, 0.99, 0.001, 0.5, 2.0])
m.jnt_solref[:] = np.array([0.006, 1.0])
jadr = m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "root")]
d.qpos[jadr + 2] -= 0.026
for anm in ("ankle_L", "ankle_R"):
    jid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, anm)
    m.jnt_range[jid] = np.array([np.radians(-20.0), np.radians(5.0)])

acts = ("hip_L_flx", "hip_L_ext", "knee_L_flx", "knee_L_ext",
        "ankle_L_flx", "ankle_L_ext", "hip_R_flx", "hip_R_ext",
        "knee_R_flx", "knee_R_ext", "ankle_R_flx", "ankle_R_ext")
aid = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a) for a in acts}
names = ("hip_L", "knee_L", "ankle_L", "hip_R", "knee_R", "ankle_R")
qadr = {n: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, n)]
        for n in names}
root = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, "Root")
plist = ("foot_L_contact", "toe_L_contact", "foot_R_contact", "toe_R_contact")
plates = {g: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g) for g in plist}
c6 = np.zeros(6)

u = net.make_inputs()
itE = [net.input_index(f"TONIC {s} RG ext") for s in ("L", "R")]
iheel = [net.input_index(f"PORT heel {s}") for s in ("L", "R")]
itoe = [net.input_index(f"PORT toe {s}") for s in ("L", "R")]
iib = [net.input_index(f"PORT Ib {s}") for s in ("L", "R")]
iE = {s: net.idx[f"{s} RG ext"] for s in ("L", "R")}

print("== closed-loop latch (te=3 both), every 100 ms ==")
print("  t   pelZ   hipL  kneeL anklL  hipR  kneeR anklR | "
      "ctrlL: hipE kneE ankE hipF kneF ankF | forceL[N]: same order | "
      "plates heelL toeL heelR toeR | RG-E L/R")
cf = {g: 0.0 for g in plates}
V = np.zeros(len(net.idx))
for i in range(1500):
    t = i * 0.001
    if i % 2 == 0:
        u[:] = 0.0
        u[itE[0]] = u[itE[1]] = 3.0
        for k in range(2):
            u[iheel[k]] = 4.0 * min(cf[plist[k]] / 50.0, 1.0)
            u[itoe[k]] = 4.0 * min(cf[plist[2 + k]] / 50.0, 1.0)
        u[iib[0]] = u[iib[1]] = 1.0
        V = net.step(u)
        for a, c in net.muscle_ctrl(V).items():
            d.ctrl[aid[a]] = min(c, 0.5)
    mujoco.mj_step(m, d)
    for g in cf:
        cf[g] = 0.0
    for k in range(d.ncon):
        con = d.contact[k]
        for nm, gi in plates.items():
            if con.geom1 == gi or con.geom2 == gi:
                mujoco.mj_contactForce(m, d, k, c6)
                cf[nm] += abs(float(c6[0]))
    if i % 100 == 0:
        qs = " ".join(f"{np.degrees(d.qpos[qadr[n]]):+5.1f}" for n in names)
        cs = " ".join(f"{d.ctrl[aid[a]]:.2f}" for a in
                      ("hip_L_ext", "knee_L_ext", "ankle_L_ext",
                       "hip_L_flx", "knee_L_flx", "ankle_L_flx"))
        fs = " ".join(f"{d.actuator_force[aid[a]]:5.0f}" for a in
                      ("hip_L_ext", "knee_L_ext", "ankle_L_ext",
                       "hip_L_flx", "knee_L_flx", "ankle_L_flx"))
        pf = " ".join(f"{cf[g]:5.0f}" for g in plist)
        print(f"  {t:4.1f} {d.xpos[root][2]:6.3f}  {qs} | {cs} | {fs} | "
              f"{pf} | {V[iE['L']]:.2f}/{V[iE['R']]:.2f}")
