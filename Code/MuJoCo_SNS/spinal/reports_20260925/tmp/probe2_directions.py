"""M3 probe v2: per-actuator direction DIFFS under co-contraction posture.
All actuators held at 0.30; test actuator stepped to 1.0 vs 0.0 (two runs);
report the trajectory difference at t=0.4 s + foot/toe world x (facing).
"""
import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import numpy as np
import mujoco

MJCF = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_mjcf.xml"
JNTS = ["hip_L", "knee_L", "ankle_L", "toe_L", "hip_R", "knee_R", "ankle_R", "toe_R"]

m = mujoco.MjModel.from_xml_path(MJCF)
qadr = {j: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)] for j in JNTS}
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(m.nu)]

d = mujoco.MjData(m)
mujoco.mj_forward(m, d)
gname = lambda i: mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_GEOM, i)
for i in range(m.ngeom):
    n = gname(i)
    if n in ("foot_L_contact", "toe_L_contact", "foot_L", "toe_L"):
        print(f"geom {n}: x={d.geom_xpos[i][0]:+.4f} z={d.geom_xpos[i][2]:+.4f}")
for a in ("hip_L_flx", "hip_L_ext"):
    i = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_TENDON, "t_" + a)
    print(f"tendon t_{a} rest length {d.ten_length[i]:.4f} (range midpoint is FLV peak)")


def run(ctrl_test, a_i, steps=400):
    d = mujoco.MjData(m)
    d.ctrl[:] = 0.30
    d.ctrl[a_i] = ctrl_test
    q = np.zeros((steps, len(JNTS)))
    for i in range(steps):
        mujoco.mj_step(m, d)
        for k, j in enumerate(JNTS):
            q[i, k] = d.qpos[qadr[j]]
    return np.degrees(q)


print("\n== direction diffs (deg, test1.0 minus test0.0, all others 0.30), t=0.4s ==")
hdr = "actuator".ljust(14) + "".join(j.rjust(9) for j in JNTS)
print(hdr)
for a_i, a in enumerate(acts):
    q1 = run(1.0, a_i)
    q0 = run(0.0, a_i)
    dd = q1[-1] - q0[-1]
    print(a.ljust(14) + "".join(f"{v:+9.2f}" for v in dd))
