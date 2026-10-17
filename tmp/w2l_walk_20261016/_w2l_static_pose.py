r"""Static geometry probe: where are the heel/toe contact plates at the
spawn pose (qpos0) and at the shipped keyframe (rest pose)? Explains the
L-heel-never-loads finding from the 2026-10-16 drive sweep."""
import sys

sys.path.insert(0, r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco")
import mujoco                # noqa: E402
import numpy as np           # noqa: E402

M = mujoco.MjModel.from_xml_path(
    r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_mjcf_fixed.xml")
d = mujoco.MjData(M)

GEOMS = ("foot_L_contact", "foot_R_contact", "toe_L_contact", "toe_R_contact")


def report(tag):
    mujoco.mj_forward(M, d)
    print(f"== {tag} ==")
    root = M.jnt_qposadr[mujoco.mj_name2id(M, mujoco.mjtObj.mjOBJ_JOINT, "root")]
    print(f"  root freejoint qpos: {np.round(d.qpos[root:root+7], 4)}")
    for g in GEOMS:
        gid = mujoco.mj_name2id(M, mujoco.mjtObj.mjOBJ_GEOM, g)
        print(f"  {g:16s} world z = {d.geom_xpos[gid][2]:+.4f} m")
    for j in ("ankle_L", "ankle_R", "knee_L", "knee_R", "hip_L", "hip_R"):
        jid = mujoco.mj_name2id(M, mujoco.mjtObj.mjOBJ_JOINT, j)
        qadr = M.jnt_qposadr[jid]
        print(f"  {j:8s} q = {np.degrees(d.qpos[qadr]):+7.2f} deg "
              f"range [{np.degrees(M.jnt_range[jid][0]):+.1f}, "
              f"{np.degrees(M.jnt_range[jid][1]):+.1f}]")
    print(f"  ncon = {d.ncon}")
    for k in range(d.ncon):
        c = d.contact[k]
        n1 = mujoco.mj_id2name(M, mujoco.mjtObj.mjOBJ_GEOM, c.geom1)
        n2 = mujoco.mj_id2name(M, mujoco.mjtObj.mjOBJ_GEOM, c.geom2)
        f = np.zeros(6)
        mujoco.mj_contactForce(M, d, k, f)
        print(f"    contact {n1} / {n2}: normal force {abs(f[0]):.1f} N")


d.qpos[:] = M.qpos0
report("qpos0 (spawn pose used by test_w2l_ground.py)")
print("nkey =", M.nkey)
if M.nkey:
    mujoco.mj_resetDataKeyframe(M, d, 0)
    report("keyframe 0 (shipped rest pose)")
