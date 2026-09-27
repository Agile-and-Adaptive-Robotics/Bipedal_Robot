"""SINGLE-CAUSE A/B for the 2026-09-26 Li closed-loop retry diagnosis.

Suspected cause (from diag_li_fixed.py): w2l_mjcf_fixed.xml ships the Root
WITHOUT a freejoint (make_w2l_mjcf.py:49 says "pelvis must be FREE in MuJoCo"
but the emitted Root body has no joint; validate_body.py:52 enshrines
"8 hinge joints" as the gate). The pelvis is therefore bolted at z=0.993 m,
the feet hover 2.6-3.1 cm above the floor, and Li's CONTACT-driven CPG can
never see heel/toe contact -> structurally open loop.

A/B: identical closed loop (defaults, no gain touched) on a TEMP COPY of the
fixed xml with ONLY <freejoint name="root"/> added to the Root body.
Model of record is NOT modified. Q: do feet reach ground, do contacts
register, does the body stand or collapse?
Run:  D:/Anaconda/envs/myo/python.exe diag_li_freejoint_ab.py
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

SRC = os.path.join(W2L, "w2l_mjcf_fixed.xml")
PROBE = os.path.join(HERE, "li_free_probe.xml")
DUR = 20.0
KICK, KICK_T = 3.0, 0.2
SETTLE_T = 0.5
JOINT_DAMP = 1.5
SETTLE = {"knee_L_ext": 0.35, "knee_R_ext": 0.35,
          "hip_L_ext": 0.25, "hip_R_ext": 0.25,
          "ankle_L_ext": 0.15, "ankle_R_ext": 0.15}
JNTS = ("hip_L", "knee_L", "ankle_L", "hip_R", "knee_R", "ankle_R")
MON = ("foot_L_contact", "toe_L_contact", "foot_R_contact", "toe_R_contact")


def run(mjcf: str, label: str) -> None:
    net = build(dt=DT_NET)
    m = mujoco.MjModel.from_xml_path(mjcf)
    d = mujoco.MjData(m)
    m.dof_damping[:] = JOINT_DAMP
    act_ids = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
               for a in net.muscle_outputs}
    qadr = {j: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)]
            for j in JNTS}
    root_z = 2 if m.njnt == 9 else None
    if root_z is None:
        print(f"[{label}] model has {m.njnt} joints (welded) - qpos[2] is ankle_L")
    g_ids = [mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g) for g in MON]
    gname = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_GEOM, i) for i in range(m.ngeom)]

    nsteps = int(DUR / 0.001)
    pz = np.zeros(nsteps); px = np.zeros(nsteps)
    cflag = np.zeros((nsteps, 4))
    viol = {j: 0.0 for j in JNTS}
    qend = {}
    pair_force = {}
    footz = {g: np.zeros(nsteps) for g in MON}
    bid = {g: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, g) for g in MON}

    for i in range(nsteps):
        t = i * 0.001
        c = np.zeros(4)
        for k in range(d.ncon):
            con = d.contact[k]
            f6 = np.zeros(6)
            mujoco.mj_contactForce(m, d, k, f6)
            if f6[0] > 1e-6:
                key = (gname[con.geom1], gname[con.geom2])
                pair_force[key] = pair_force.get(key, 0.0) + f6[0]
                for j, gid in enumerate(g_ids):
                    if con.geom1 == gid or con.geom2 == gid:
                        c[j] = 1.0
        hipL = np.degrees(d.qpos[qadr["hip_L"]])
        hipR = np.degrees(d.qpos[qadr["hip_R"]])
        u = net.make_inputs()
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
        pz[i] = d.qpos[root_z] if root_z is not None else 0.99298
        px[i] = d.qpos[0] if root_z is not None else -3.454
        cflag[i] = c
        for g in MON:
            footz[g][i] = d.xpos[bid[g]][2]
        for jn in JNTS:
            deg = np.degrees(d.qpos[qadr[jn]])
            lo, hi = np.degrees(m.jnt_range[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, jn)])
            viol[jn] = max(viol[jn], max(lo - deg, deg - hi, 0.0))
        if not np.all(np.isfinite(d.qpos)):
            print(f"[{label}] qpos NON-FINITE at t={t:.3f}")
            return
        if i in (0, 500, 1000, 2000, 5000, 10000, nsteps - 1):
            qend[t] = {jn: round(float(np.degrees(d.qpos[qadr[jn]])), 1) for jn in JNTS}

    print(f"\n== [{label}] {DUR:.0f}s ==")
    print(f"  pelvis z: min {pz.min():+.3f} (t={np.argmin(pz)*0.001:.2f}s)  end {pz[-1]:+.3f}   "
          f"pelvis x: start {px[0]:+.3f} end {px[-1]:+.3f}")
    nf = int(cflag.sum())
    print(f"  monitored-geom contact steps: {nf} ({1000*nf/nsteps:.1f} s); "
          f"flags ON frac L_foot {cflag[:,0].mean():.3f} L_toe {cflag[:,1].mean():.3f} "
          f"R_foot {cflag[:,2].mean():.3f} R_toe {cflag[:,3].mean():.3f}")
    print("  all-pairs contacts (top 6 by sumF):")
    for key, f in sorted(pair_force.items(), key=lambda kv: -kv[1])[:6]:
        print(f"    {key[0]:<16} x {key[1]:<16} sumF {f:11.1f}")
    print("  max joint-limit violation (deg):",
          {k: round(v, 1) for k, v in viol.items()})
    print("  pose census (deg):")
    for tt, qq in qend.items():
        print(f"    t={tt:5.1f}s {qq}")
    print("  foot-body z at end (m):", {g: round(float(footz[g][-1]), 3) for g in MON})


def main() -> int:
    txt = open(SRC, encoding="utf-8").read()
    tag = '<body name="Root"'
    i0 = txt.index(tag)
    i1 = txt.index(">", i0) + 1
    assert "<freejoint" not in txt[:i1][-400:], "Root already has a joint?"
    open(PROBE, "w", encoding="utf-8").write(
        txt[:i1] + '\n      <freejoint name="root"/>  <!-- A/B probe only -->' + txt[i1:])
    print("wrote probe model:", PROBE)
    run(SRC, "A: shipped (welded Root)")
    run(PROBE, "B: +freejoint (pelvis free)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
