"""M3 probe v4: actuator direction table with GRAVITY OFF (air rig).
No collapse, no contact - the pure kinematic direction each actuator drives
its joint, plus world foot/toe motion for anatomical labeling.
"""
import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import numpy as np
import mujoco

AIR = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_air.xml"
JNTS = ["hip_L", "knee_L", "ankle_L", "toe_L", "hip_R", "knee_R", "ankle_R", "toe_R"]
m = mujoco.MjModel.from_xml_path(AIR)
m.opt.gravity[:] = 0.0
qadr = {j: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)] for j in JNTS}
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(m.nu)]
g = {n: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, n) for n in
     ("foot_L", "toe_L", "foot_R", "toe_R")}


def run(ctrl_test, a_i, steps=400):
    d = mujoco.MjData(m)
    d.ctrl[:] = 0.30
    d.ctrl[a_i] = ctrl_test
    nc = 0
    for i in range(steps):
        mujoco.mj_step(m, d)
        nc = max(nc, d.ncon)
    return (np.degrees([d.qpos[qadr[j]] for j in JNTS]),
            d.geom_xpos[g["foot_L"]][0], d.geom_xpos[g["toe_L"]][2],
            d.geom_xpos[g["foot_R"]][0], d.geom_xpos[g["toe_R"]][2], nc)


q0, fx0, tz0, fxR0, tzR0, nc0 = run(0.0, 0)
print(f"baseline: qpos(deg) " + " ".join(f"{j}={v:+.2f}" for j, v in zip(JNTS, q0)) +
      f" | footL x={fx0:+.4f} toeL z={tz0:+.4f} | footR x={fxR0:+.4f} toeR z={tzR0:+.4f} ncon={nc0}")
print("\n== test 1.0 minus 0.0, gravity OFF, others 0.30, t=0.4 s ==")
print("actuator".ljust(14) + "".join(j.rjust(9) for j in JNTS) +
      "   footLx   toeLz   footRx   toeRz   ncon")
for a_i, a in enumerate(acts):
    q1, fx1, tz1, fxR1, tzR1, nc1 = run(1.0, a_i)
    print(a.ljust(14) + "".join(f"{v:+9.2f}" for v in (q1 - q0)) +
          f"  {fx1-fx0:+8.3f} {tz1-tz0:+8.3f} {fxR1-fxR0:+8.3f} {tzR1-tzR0:+8.3f}   {nc1}")
