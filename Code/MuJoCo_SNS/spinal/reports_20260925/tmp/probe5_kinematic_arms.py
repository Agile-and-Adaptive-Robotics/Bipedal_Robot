"""M3 probe v5: kinematic moment-arm table. qpos perturbed +1 deg per joint,
mj_forward only (no dynamics), record every tendon's length change (mm).
dl < 0 = muscle SHORTENS when joint moves +  => that muscle's contraction
drives the joint toward -qpos. Pure kinematics; settles flexion signs.
"""
import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import numpy as np
import mujoco

AIR = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_air.xml"
JNTS = ["hip_L", "knee_L", "ankle_L", "toe_L", "hip_R", "knee_R", "ankle_R", "toe_R"]
m = mujoco.MjModel.from_xml_path(AIR)
qadr = {j: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)] for j in JNTS}
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(m.nu)]
ntn = m.ntendon
tname = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_TENDON, i) for i in range(ntn)]

d = mujoco.MjData(m)
mujoco.mj_forward(m, d)
L0 = d.ten_length.copy()

print("rest tendon lengths (m):")
for i in range(ntn):
    if i < 12:
        rng = (m.actuator_lengthrange[i][0], m.actuator_lengthrange[i][1])
        mid = 0.5 * (rng[0] + rng[1])
        frac = (L0[i] - rng[0]) / (rng[1] - rng[0])
        print(f"  {tname[i]:<16} L={L0[i]:.4f}  FLV-normalized={frac:.3f}")

DD = 1.0  # deg
print(f"\n dl [mm] per +{DD:.0f} deg of joint angle (mj_forward only)")
print("joint".ljust(9) + "".join(a.replace('_ext', '/E').replace('_flx', '/F').rjust(13) for a in acts))
for j in JNTS:
    d = mujoco.MjData(m)
    d.qpos[qadr[j]] += np.radians(DD)
    mujoco.mj_forward(m, d)
    dl = (d.ten_length - L0) * 1000.0
    print(j.ljust(9) + "".join(f"{v:+13.2f}" for v in dl[:12]))
