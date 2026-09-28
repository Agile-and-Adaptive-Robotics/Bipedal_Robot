"""Focused diagnosis of the Li closed-loop retry on w2l_mjcf_fixed.xml
(2026-09-26; run1 = reports_20260925/logs/test_li_stepping_fixed_run1_20s.log:
body stays up 0.63-0.99 m but 0 stance episodes for 20 s, ankles tens of
degrees outside their [-20,-5] deg limit).

Same loop as test_li_stepping.py (defaults, no retune), plus:
  - what actually contacts what (all geom pairs with normal force > 1e-6)
  - which contacts carry the weight (support audit)
  - joint-limit violation census (which joint, when it first blows through)
  - is the neural side alive (MV voltages, ctrl, CPG/SN voltages)
  - where the feet are (world z of the 4 monitored contact bodies)
Run:  D:/Anaconda/envs/myo/python.exe diag_li_fixed.py
"""
from __future__ import annotations

import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"D:\Anaconda\envs\myo")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
W2L = os.path.join(HERE, "..", "..", "w2l_mujoco")
sys.path.insert(0, W2L)

import mujoco  # noqa: E402
from test_li_stepping import DT_NET, N_SUB  # noqa: E402
from build_li_net import build  # noqa: E402

MJCF = os.path.join(W2L, "w2l_mjcf_fixed.xml")
DUR = 20.0
# verbatim copies of test_li_stepping.main()'s runtime constants (they are
# function locals there, not module constants)
KICK, KICK_T = 3.0, 0.2
SETTLE_T = 0.5
JOINT_DAMP = 1.5
SETTLE = {"knee_L_ext": 0.35, "knee_R_ext": 0.35,
          "hip_L_ext": 0.25, "hip_R_ext": 0.25,
          "ankle_L_ext": 0.15, "ankle_R_ext": 0.15}
JNTS = ("hip_L", "knee_L", "ankle_L", "hip_R", "knee_R", "ankle_R")
MON = ("foot_L_contact", "toe_L_contact", "foot_R_contact", "toe_R_contact")


def main() -> int:
    net = build(dt=DT_NET)
    m = mujoco.MjModel.from_xml_path(MJCF)
    d = mujoco.MjData(m)
    m.dof_damping[:] = JOINT_DAMP

    act_ids = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
               for a in net.muscle_outputs}
    jnt_ids = {j: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j) for j in JNTS}
    qadr = {j: m.jnt_qposadr[i] for j, i in jnt_ids.items()}
    g_ids = [mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g) for g in MON]
    body_ids = {g: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, g) for g in MON}
    gname = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_GEOM, i) for i in range(m.ngeom)]

    nsteps = int(DUR / 0.001)
    cols = {}
    def L(key, i, val):
        cols.setdefault(key, np.zeros(nsteps))[i] = val

    u = net.make_inputs()
    pair_force = {}          # (geom1, geom2) -> sum normal force over run
    pair_n = {}              # (geom1, geom2) -> steps in contact
    first_viol = {}          # joint -> (t, deg beyond limit) first breach > 1 deg
    viol_max = {}            # joint -> max deg beyond limit
    mv_rows, ctrl_rows = [], []

    for i in range(nsteps):
        t = i * 0.001
        c = np.zeros(4)
        # all-pairs contact audit this step
        for k in range(d.ncon):
            con = d.contact[k]
            f6 = np.zeros(6)
            mujoco.mj_contactForce(m, d, k, f6)
            if f6[0] > 1e-6:
                key = (gname[con.geom1], gname[con.geom2])
                pair_force[key] = pair_force.get(key, 0.0) + f6[0]
                pair_n[key] = pair_n.get(key, 0) + 1
            for j, gid in enumerate(g_ids):
                if con.geom1 == gid or con.geom2 == gid:
                    f6b = np.zeros(6)
                    mujoco.mj_contactForce(m, d, k, f6b)
                    if f6b[0] > 1e-6:
                        c[j] = 1.0
        hipL = np.degrees(d.qpos[qadr["hip_L"]])
        hipR = np.degrees(d.qpos[qadr["hip_R"]])
        u[:] = 0.0
        net.set_all_tonics(u, extra={"R tonic drive": KICK} if t < KICK_T else None)
        net.set_port(u, "L_foot ground contact", net.contact_current(None, c[0]))
        net.set_port(u, "L_toe ground contact", net.contact_current(None, c[1]))
        net.set_port(u, "R_foot ground contact", net.contact_current(None, c[2]))
        net.set_port(u, "R_toe ground contact", net.contact_current(None, c[3]))
        net.set_port(u, "L_hip middle", net.hip_current(None, hipL))
        net.set_port(u, "R_hip middle", net.hip_current(None, hipR))
        for _ in range(N_SUB):
            V = net.step(u)
        ctrl = net.muscle_ctrl(V)
        if t < SETTLE_T:
            for a, v in SETTLE.items():
                ctrl[a] = max(ctrl[a], v)
        for a, val in ctrl.items():
            d.ctrl[act_ids[a]] = min(val, 0.9)
        mujoco.mj_step(m, d)

        L("pz", i, d.qpos[2]); L("px", i, d.qpos[0]); L("py", i, d.qpos[1])
        for j, jn in enumerate(JNTS):
            deg = np.degrees(d.qpos[qadr[jn]])
            L(jn, i, deg)
            lo, hi = np.degrees(m.jnt_range[jnt_ids[jn]])
            over = max(lo - deg, deg - hi, 0.0)
            if over > 1.0:
                viol_max[jn] = max(viol_max.get(jn, 0.0), over)
                if jn not in first_viol:
                    first_viol[jn] = (t, deg, lo, hi)
        for j, gn in enumerate(MON):
            L(gn + "_z", i, d.xpos[body_ids[gn]][2])
        for a, aid in act_ids.items():
            L("F_" + a, i, d.actuator_force[aid])
            L("c_" + a, i, d.ctrl[aid])
        L("ncon", i, d.ncon)
        L("cL", i, c[0]); L("ctL", i, c[1]); L("cR", i, c[2]); L("ctR", i, c[3])
        if not np.all(np.isfinite(d.qpos)):
            print(f"[FAIL] qpos non-finite at t={t:.3f}")
            return 1
        if i % int(nsteps / 8) == 0:
            mv_rows.append((t, {a: V[net.idx[net.muscle_outputs[a]]] for a in net.muscle_outputs}))
            ctrl_rows.append((t, dict(ctrl)))

    t = np.arange(nsteps) * 0.001
    print(f"== diag_li_fixed: {DUR:.0f} s on {os.path.basename(MJCF)} ==")
    print(f"pelvis z: start {cols['pz'][0]:.3f}  min {cols['pz'].min():.3f} (t={t[np.argmin(cols['pz'])]:.2f}s)"
          f"  end {cols['pz'][-1]:.3f}   |  x end {cols['px'][-1]:.3f}")

    print("\n-- joint limits (fixed xml) vs realized --")
    for jn in JNTS:
        lo, hi = np.degrees(m.jnt_range[jnt_ids[jn]])
        q = cols[jn]
        fv = first_viol.get(jn)
        fvs = (f"first >1deg breach t={fv[0]:.2f}s (q={fv[1]:+.1f}, range [{fv[2]:+.1f},{fv[3]:+.1f}])"
               if fv else "never breached")
        print(f"  {jn:<7} range [{lo:+.1f},{hi:+.1f}]  realized {q.min():+.1f}..{q.max():+.1f}"
              f"  max violation {viol_max.get(jn, 0.0):.1f} deg  {fvs}")

    print("\n-- contacts (normal force > 1e-6), whole run --")
    if not pair_force:
        print("  NONE. d.ncon max =", int(cols['ncon'].max()))
    else:
        for key, f in sorted(pair_force.items(), key=lambda kv: -kv[1]):
            print(f"  {key[0]:<16} x {key[1]:<16}  sumF {f:12.1f} N  steps {pair_n[key]} ({1000*pair_n[key]/nsteps:.1f} s)")
    mon_s = sum(pair_force.get((g, g2), 0.0) + pair_force.get((g2, g), 0.0)
                for g in MON for g2 in [x for pair in pair_force for x in pair] if g2 != g)
    print(f"  monitored-geom total normal force: {mon_s:.1f} N·s·1e-3 units (per-step sums)")
    print(f"  monitored contact flags ON fraction: L_foot {cols['cL'].mean():.3f}"
          f"  L_toe {cols['ctL'].mean():.3f}  R_foot {cols['cR'].mean():.3f}  R_toe {cols['ctR'].mean():.3f}")

    print("\n-- foot contact-body world z (m), sampled --")
    for k in range(0, nsteps, nsteps // 10):
        zs = [f"{cols[g + '_z'][k]:+.3f}" for g in MON]
        print(f"  t={t[k]:5.1f}s  " + "  ".join(f"{g.split('_')[0][:2]}{g.split('_')[1][0]}:{z}" for g, z in zip(MON, zs)))

    print("\n-- neural side (is it alive?) --")
    for tt, mvs in mv_rows[:3] + mv_rows[-3:]:
        top = sorted(mvs.items(), key=lambda kv: -kv[1])[:4]
        print(f"  t={tt:5.1f}s  top MV: " + ", ".join(f"{a.replace('_ext','E').replace('_flx','F')}={v:.1f}mV" for a, v in top))
    for lab in ("L_CPG stance", "R_CPG stance", "L_foot ground contact", "R_foot ground contact"):
        v = None
    # full-run voltage probe: rerun-free approximation via stored ctrl rows is
    # not enough; report ctrl instead (V rows sampled above are representative)
    for tt, cc in ctrl_rows[:2] + ctrl_rows[-2:]:
        nz = {a: round(v, 2) for a, v in cc.items() if v > 0.02}
        print(f"  t={tt:5.1f}s  ctrl>0.02: {nz if nz else 'NONE (open loop - no neural drive reaches muscles)'}")
    for a in ("knee_L_ext", "knee_R_ext", "ankle_L_ext", "ankle_R_ext",
              "hip_L_ext", "hip_R_ext"):
        F = cols["F_" + a]
        C = cols["c_" + a]
        print(f"  {a:<13} ctrl {C.min():.2f}..{C.max():.2f}   actuator force {F.min():8.1f}..{F.max():8.1f} N"
              f"  (last 5 s mean {F[-5000:].mean():7.1f} N)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
