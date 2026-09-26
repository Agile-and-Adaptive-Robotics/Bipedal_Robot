"""M3 probe v3: direction table in the AIR rig (root lifted +0.3 m so no
ground contact is possible), co-contraction 0.30, test actuator 1.0 vs 0.0.
Reads joint qpos + world foot/toe pose to settle anatomical signs.
Also writes the lifted XML once to w2l_mujoco\\w2l_air.xml (new file).
"""
import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import numpy as np
import mujoco

SRC = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_mjcf.xml"
AIR = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_air.xml"
LIFT = 0.30   # m; leg vertical reach ~0.95 m, pelvis at 0.993 -> feet can never reach z=0

txt = open(SRC, encoding="utf-8").read()
old = 'pos="-3.454 0 0.99298"'
assert old in txt, "root pos pattern not found"
new = ('pos="-3.454 0 1.29298"')
hdr = ("<!-- AIR VARIANT for milestone 3 (w2l_mujoco/test_w2l_air.py): identical to\n"
       "     w2l_mjcf.xml except the welded Root is lifted +0.30 m so that NO leg\n"
       "     pose can reach the ground plane (leg vertical reach ~0.95 m). Keeps the\n"
       "     M1 artifact untouched and makes 'no ground contact' provable (ncon==0). -->\n")
open(AIR, "w", encoding="utf-8").write(hdr + txt.replace(old, new))
print("wrote", AIR)

JNTS = ["hip_L", "knee_L", "ankle_L", "toe_L", "hip_R", "knee_R", "ankle_R", "toe_R"]
m = mujoco.MjModel.from_xml_path(AIR)
qadr = {j: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)] for j in JNTS}
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(m.nu)]
g_foot_L = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "foot_L")
g_toe_L = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "toe_L")
g_foot_R = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "foot_R")
g_toe_R = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "toe_R")
g_fc_L = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "foot_L_contact")
g_tc_L = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "toe_L_contact")
g_fc_R = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "foot_R_contact")
g_tc_R = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "toe_R_contact")


def run(ctrl_test, a_i, steps=600):
    d = mujoco.MjData(m)
    d.ctrl[:] = 0.30
    d.ctrl[a_i] = ctrl_test
    maxncon = 0
    for i in range(steps):
        mujoco.mj_step(m, d)
        maxncon = max(maxncon, d.ncon)
    return (np.degrees([d.qpos[qadr[j]] for j in JNTS]),
            d.geom_xpos[g_foot_L][0], d.geom_xpos[g_toe_L][2],
            d.geom_xpos[g_foot_R][0], d.geom_xpos[g_toe_R][2],
            d.geom_xpos[g_fc_L][2], d.geom_xpos[g_tc_L][2],
            d.geom_xpos[g_fc_R][2], d.geom_xpos[g_tc_R][2], maxncon)


q0, fx0, tz0, fxR0, tzR0, z1, z2, z3, z4, nc0 = run(0.0, 0)
print(f"baseline (all 0.30): qpos(deg) " +
      " ".join(f"{j}={v:+.1f}" for j, v in zip(JNTS, q0)) +
      f" | footL x={fx0:+.3f} toeL z={tz0:+.3f} | footR x={fxR0:+.3f} toeR z={tzR0:+.3f}"
      f" | contact plate z L=({z1:+.3f},{z2:+.3f}) R=({z3:+.3f},{z4:+.3f}) ncon={nc0}")

print("\n== test 1.0 minus 0.0 (all others 0.30), air rig, t=0.6 s ==")
print("actuator".ljust(14) + "".join(j.rjust(9) for j in JNTS) +
      "   footLx   toeLz   footRx   toeRz   ncon")
for a_i, a in enumerate(acts):
    q1, fx1, tz1, fxR1, tzR1, *_z, nc1 = run(1.0, a_i)
    qd = q1 - q0
    print(a.ljust(14) + "".join(f"{v:+9.2f}" for v in qd) +
          f"  {fx1-fx0:+8.3f} {tz1-tz0:+8.3f} {fxR1-fxR0:+8.3f} {tzR1-tzR0:+8.3f}   {nc1}")
