"""M3 probe: per-actuator joint-direction table + ground-contact check on the
shipped M1 body (w2l_mjcf.xml). Read-only w.r.t. the artifact.

For each of the 12 actuators: hold ctrl=1.0 for 1.5 s with the pelvis as
shipped (no joint -> welded to world), all other ctrl 0; record dqpos of the
6 leg hinges. Also: 2 s all-zero passive run -> which joints move, min geom
height, max ncon (can anything touch the ground plane?).
"""
import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import numpy as np
import mujoco

MJCF = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_mjcf.xml"
JNTS = ["hip_L", "knee_L", "ankle_L", "toe_L", "hip_R", "knee_R", "ankle_R", "toe_R"]

m = mujoco.MjModel.from_xml_path(MJCF)
print(f"nq={m.nq} nu={m.nu} njnt={m.njnt} nbody={m.nbody}")
print("joints:", [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_JOINT, i) for i in range(m.njnt)])
root_jnt = m.body_jntadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, "Root")]
print("Root body jntadr:", root_jnt, "(0/-1 -> welded)")
qadr = {j: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)] for j in JNTS}
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(m.nu)]

# rest pose geom heights
d = mujoco.MjData(m)
mujoco.mj_forward(m, d)
gl = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_GEOM, i) for i in range(m.ngeom)]
zs = {gl[i]: float(d.geom_xpos[i][2]) for i in range(m.ngeom)}
print("rest geom z (m):", {k: round(v, 4) for k, v in sorted(zs.items(), key=lambda kv: kv[1])})
print("rest qpos (deg):", {j: round(float(np.degrees(d.qpos[qadr[j]])), 2) for j in JNTS})
print("rest actuator lengths:", [round(float(d.actuator_length[i]), 4) for i in range(m.nu)])

print("\n== passive 2 s (all ctrl 0) ==")
d = mujoco.MjData(m)
maxncon = 0
for i in range(2000):
    mujoco.mj_step(m, d)
    maxncon = max(maxncon, d.ncon)
    if not np.isfinite(d.qpos).all():
        print("NON-FINITE at", i); break
print("after 2 s qpos (deg):", {j: round(float(np.degrees(d.qpos[qadr[j]])), 2) for j in JNTS})
print("max ncon over passive run:", maxncon)

print("\n== per-actuator probe: ctrl=1.0 for 1.5 s, others 0, fresh model each ==")
hdr = "actuator".ljust(14) + "".join(j.rjust(9) for j in JNTS)
print(hdr)
for a_i, a in enumerate(acts):
    d = mujoco.MjData(m)
    mujoco.mj_forward(m, d)
    q0 = np.array([d.qpos[qadr[j]] for j in JNTS])
    d.ctrl[:] = 0.0
    d.ctrl[a_i] = 1.0
    ncon = 0
    for i in range(1500):
        mujoco.mj_step(m, d)
        ncon = max(ncon, d.ncon)
    q1 = np.array([d.qpos[qadr[j]] for j in JNTS])
    dq = np.degrees(q1 - q0)
    print(a.ljust(14) + "".join(f"{v:+9.2f}" for v in dq) + f"   ncon={ncon}")
