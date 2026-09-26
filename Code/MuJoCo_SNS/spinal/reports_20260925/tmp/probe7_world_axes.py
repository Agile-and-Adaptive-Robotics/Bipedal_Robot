"""M3 probe v7: WORLD axes of all 8 hinges (R_body @ axis_local), and where
the foot actually moves under a knee/hip/ankle rotation. Checks the M1
transport: M1 report section 4b claims all four joints per leg are sagittal
(lateral AL-z) hinges.
"""
import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import numpy as np
import mujoco

AIR = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_air.xml"
m = mujoco.MjModel.from_xml_path(AIR)
d = mujoco.MjData(m)
mujoco.mj_forward(m, d)

print("joint   world axis (MJ)          world axis (AL remap back: AL=(x,-y... ) n/a)   body")
for i in range(m.njnt):
    n = mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_JOINT, i)
    b = m.jnt_bodyid[i]
    R = d.xmat[b].reshape(3, 3)
    ax_w = R @ m.jnt_axis[i]
    print(f"{n:<8} ({ax_w[0]:+.4f}, {ax_w[1]:+.4f}, {ax_w[2]:+.4f})   body={mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_BODY, b)}")

# displacement direction of toe geom per +10 deg of each L joint
g_toe = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "toe_L")
g_foot = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, "foot_L")
for j in ("hip_L", "knee_L", "ankle_L", "toe_L"):
    d2 = mujoco.MjData(m)
    mujoco.mj_forward(m, d2)
    p0t = d2.geom_xpos[g_toe].copy(); p0f = d2.geom_xpos[g_foot].copy()
    d2.qpos[m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)]] += np.radians(10.0)
    mujoco.mj_forward(m, d2)
    dt = d2.geom_xpos[g_toe] - p0t
    dfl = d2.geom_xpos[g_foot] - p0f
    print(f"+10deg {j}: toe dxyz=({dt[0]:+.3f},{dt[1]:+.3f},{dt[2]:+.3f}) m  foot dxyz=({dfl[0]:+.3f},{dfl[1]:+.3f},{dfl[2]:+.3f})")
