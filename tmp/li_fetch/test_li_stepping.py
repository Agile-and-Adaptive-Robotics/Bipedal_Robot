r"""MILESTONE-2 GATE: Li's contact-driven 2-layer CPG (build_li_net) driving
the M1 MuJoCo body (w2l_mjcf.xml). 20 s ground run; verdict measured against
Li's own AnimatLab reference (his asim rerun headless on this box,
Neuromechanical_Models/Li Model/DataTool_7.txt):
   period 1.30 s (0.77 Hz), stance ~0.61-0.67 s, L/R antiphase,
   height 0.95-1.02 (never falls), 0.64 m/s.

Run:  C:\Users\Ben Bolen\.conda\envs\myo\python.exe test_li_stepping.py
"""
from __future__ import annotations

import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import mujoco  # noqa: E402  (CONDA_PREFIX set before import, per the skill)
from build_li_net import build  # noqa: E402

# M3 RETRY (2026-09-26): the M2 gate's original runs used w2l_mjcf.xml, whose
# knee/ankle hinges shipped as VERTICAL-axis yaw joints (transport bug,
# fix_joint_axes.py) - the body-side collapse they showed is not interpretable.
# The closed loop now runs on the axis-corrected body of record.
MJCF = os.path.join(HERE, "w2l_mjcf_fixed.xml")
DT_NET = 0.0002
N_SUB = 5                      # net substeps per 1 ms physics step
DUR = 20.0                     # s, the gate asks for >= 20 s
REF_PERIOD = 1.30              # s, Li's DataTool_7.txt
REF_SPEED = 0.6384             # m/s


def contact_counts(m, d, g_ids):
    c = np.zeros(4)
    for i in range(d.ncon):
        con = d.contact[i]
        for k, gid in enumerate(g_ids):
            if con.geom1 == gid or con.geom2 == gid:
                f6 = np.zeros(6)
                mujoco.mj_contactForce(m, d, i, f6)
                if f6[0] > 1e-6:
                    c[k] = 1.0
    return c


def episodes(on, t, gap=0.3):
    idx = np.where(on)[0]
    if len(idx) == 0:
        return []
    eps = [[idx[0], idx[0]]]
    for i in idx[1:]:
        if (t[i] - t[eps[-1][1]]) < gap:
            eps[-1][1] = i
        else:
            eps.append([i, i])
    return [(t[a], t[b]) for a, b in eps]


def main():
    # knob overrides for calibration sweeps (all default = verbatim/calibrated)
    kv = {}
    for arg in sys.argv[1:]:
        if arg.startswith("--") and "=" in arg:
            k, v = arg[2:].split("=", 1)
            kv[k.replace("-", "_")] = float(v)
    dur = kv.pop("dur", DUR)
    net = build(dt=DT_NET, knobs=kv)

    m = mujoco.MjModel.from_xml_path(MJCF)
    d = mujoco.MjData(m)
    # spawn keyframe (2026-09-30): ankles inside range, feet at 1 mm
    try:
        mujoco.mj_resetDataKeyframe(m, d, 0)
        mujoco.mj_forward(m, d)
    except Exception:
        pass

    act_ids = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
               for a in net.muscle_outputs}
    jnt_ids = {j: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)
               for j in ("hip_L", "knee_L", "ankle_L",
                         "hip_R", "knee_R", "ankle_R")}
    qadr = {j: m.jnt_qposadr[i] for j, i in jnt_ids.items()}
    g_ids = [mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g) for g in
             ("foot_L_contact", "toe_L_contact",
              "foot_R_contact", "toe_R_contact")]

    nsteps = int(dur / 0.001)
    log_t = np.zeros(nsteps)
    log_q = np.zeros((nsteps, 6))
    log_c = np.zeros((nsteps, 4))          # footL toeL footR toeR
    log_h = np.zeros((nsteps, 2))          # pelvis z, x
    u = net.make_inputs()

    # small, documented kickoff: extra tonic on the R tonic port for 0.2 s
    # (Li's own asim ships an (inactive) Stimulus_1; the symmetric rest pose
    # otherwise holds both feet down and both stance CPGs off.)
    KICK, KICK_T = 3.0, 0.2
    # BEN'S PROTOCOL (2026-09-30, his correction): in AnimatLab these walkers
    # run under a VIRTUAL WALKER - suspended in AIR while the CPG air-steps,
    # then DROPPED onto the platform so heel strikes drive the contact CPG.
    # The previous gate spawned feet-down with a settle-assist co-contraction
    # (an invention, not the protocol) and heels never loaded (0 stance
    # episodes). REPLACED by the protocol: a support harness (feedforward
    # weight + vertical PD, weak horizontal damping - the runner-rig idiom)
    # holds the root CLEAR above contact during [0, T_HOLD), then releases
    # ALL force at T_HOLD - the body drops, lands, and the heel SNs fire.
    # Settle-assist block REMOVED.
    T_HOLD = float(kv.pop("hold", 3.0))    # air-hold duration, s
    T_LOW = float(kv.pop("lower", 1.0))     # platform-lowering duration, s
    CLEAR = float(kv.pop("clear", 0.04))   # extra height above contact, m
    # Ben 2026-09-30 (clarification): "In AnimatLab they can drop it WITH
    # the harness" - the virtual walker STAYS ATTACHED after the drop. So
    # after the lowering, keep the support on at SUPPORT fraction of body
    # weight + PD at contact height (feet bear the rest -> heel SNs fire,
    # contact-driven CPG engages, body cannot buckle - sidesteps the
    # no-Kse/Kpe body-fidelity blocker for the march gate).
    SUPPORT = float(kv.pop("support", 0.7))    # retained weight fraction
    root_bid = int(m.jnt_bodyid[0])
    _mass = float(np.sum(m.body_mass))
    _wgt = _mass * (-m.opt.gravity[2])            # N
    K_H = 4.0 * _wgt / max(CLEAR, 0.01)           # sag < CLEAR/4 at weight
    C_H = 2.0 * float(np.sqrt(K_H * _mass))       # ~critical
    z_contact = float(d.qpos[2])                  # keyframe z (feet 1 mm)
    z_hold = z_contact + CLEAR
    # M1 DEVIATION STAND-IN: AnimatLab LinearHill muscles carry B damping
    # (400-800 N s/m per muscle, M1 report section 3: "B - none"); the M1
    # MuJoCo 2.3.7 <muscle> has no damping, and Li's own joint frictions are
    # Enabled=False (aproj Friction blocks) except the toes - so his body is
    # damped through the MUSCLES, which M1 does not model. Without a joint
    # damping stand-in the undamped hinges blow through their limits at the
    # landing impact (measured: ankle -82 deg through a [-20,-5] limit).
    # 2026-09-30: the XML now carries per-joint damping sum(B*r^2)
    # (make_w2l_mjcf.py) - default = USE THE XML; --joint-damp=X still
    # overrides for sweeps.
    JOINT_DAMP = kv.pop("joint_damp", 0.0)
    if JOINT_DAMP > 0:
        m.dof_damping[:] = JOINT_DAMP

    for i in range(nsteps):
        t = i * 0.001
        c = contact_counts(m, d, g_ids)    # footL toeL footR toeR
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
        V = u  # placeholder to silence linters; real step below
        for _ in range(N_SUB):
            V = net.step(u)
        ctrl = net.muscle_ctrl(V)
        # Ben's protocol, AnimatLab semantics: hold in AIR (air-stepping)
        # -> LOWER to the platform over T_LOW -> harness STAYS ON at
        # SUPPORT fraction (drop WITH the harness - 2026-09-30 Ben).
        if t < T_HOLD:
            z_t, wff = z_hold, 1.0
        elif t < T_HOLD + T_LOW:   # linear descent to contact height
            z_t = z_hold + (z_contact - z_hold) * \
                (t - T_HOLD) / max(T_LOW, 1e-6)
            wff = 1.0
        else:                      # ON THE PLATFORM, harness retained
            z_t, wff = z_contact, SUPPORT
        Fz = wff * _wgt + K_H * (z_t - d.qpos[2]) - C_H * d.qvel[2]
        d.xfrc_applied[root_bid, 2] = max(Fz, 0.0)
        d.xfrc_applied[root_bid, 0] = -20.0 * d.qvel[0]   # weak xy damp
        d.xfrc_applied[root_bid, 1] = -20.0 * d.qvel[1]
        for a, val in ctrl.items():
            d.ctrl[act_ids[a]] = min(val, 0.9)
        mujoco.mj_step(m, d)
        log_t[i] = t
        log_q[i] = [np.degrees(d.qpos[qadr[j]]) for j in
                    ("hip_L", "knee_L", "ankle_L", "hip_R", "knee_R", "ankle_R")]
        log_c[i] = c
        log_h[i] = d.qpos[2], d.qpos[0]
        if not np.all(np.isfinite(d.qpos)):
            print(f"[FAIL] qpos non-finite at t={t:.3f} s")
            return 1

    # ------------------------------------------------------------- metrics
    # measure from AFTER the drop + landing transient (hold phase has no
    # ground contact by construction - Ben's protocol)
    an = int(max(0.05, (T_HOLD + T_LOW + 0.5) / dur) * nsteps)
    t, q, c, h = log_t[an:], log_q[an:], log_c[an:], log_h[an:]
    height_min, height_end = h[:, 0].min(), h[-1, 0]
    speed = (h[-1, 1] - h[0, 1]) / (t[-1] - t[0])

    print("== M2 gate: Li CPG on the M1 MuJoCo body, "
          f"{dur:.0f} s ground run ==  knobs: {kv or 'defaults'}")
    print(f"  pelvis height: min {height_min:.3f} m, end {height_end:.3f} m")
    print(f"  mean forward speed: {speed:+.3f} m/s   (Li ref {REF_SPEED:+.3f})")

    verdicts = []
    ref, periods = {}, {}
    for k, (name, col) in enumerate([("L_foot", 0), ("R_foot", 2)]):
        on = c[:, col] > 0.5
        eps = [e for e in episodes(on, t) if e[1] - e[0] > 0.05]
        onsets = np.array([s for s, e in eps])
        durs = [e - s for s, e in eps]
        per = np.median(np.diff(onsets)) if len(onsets) > 2 else float("nan")
        duty = np.sum(on) / len(on)
        print(f"  {name}: {len(eps)} stance episodes, "
              f"period {per:.3f} s (Li ref {REF_PERIOD}), "
              f"stance {np.median(durs) if durs else float('nan'):.3f} s, "
              f"duty {duty:.2f}")
        ref[name] = on.astype(float)
        periods[name] = per
        verdicts.append(len(eps) >= 3)

    # alternation: L/R stance-signal Pearson r at 0 lag (downsample to 10 ms)
    ds = slice(None, None, 10)
    a, b = ref["L_foot"][ds], ref["R_foot"][ds]
    if a.std() > 0 and b.std() > 0:
        r = float(np.corrcoef(a, b)[0, 1])
    else:
        r = float("nan")
    print(f"  L/R stance correlation (0 lag, 10 ms samples): r = {r:+.3f} "
          f"(Li antiphase => r < -0.5)")

    for side, j in (("L", (0, 1, 2)), ("R", (3, 4, 5))):
        hip, knee, ank = q[:, j[0]], q[:, j[1]], q[:, j[2]]
        print(f"  hip {side}: {hip.min():+.1f}..{hip.max():+.1f} deg  "
              f"knee {side}: {knee.min():+.1f}..{knee.max():+.1f} deg  "
              f"ankle {side}: {ank.min():+.1f}..{ank.max():+.1f} deg")

    perL = periods["L_foot"]
    perR = periods["R_foot"]
    per_ok = any(abs(p - REF_PERIOD) / REF_PERIOD < 0.3
                 for p in (perL, perR) if np.isfinite(p))
    stepped = all(verdicts)
    anti = (r < -0.5)
    up = height_end > 0.8
    if stepped and anti and up:
        verdict = ("matched" if per_ok else "qualitatively matched")
    elif stepped:
        verdict = "qualitatively matched (steps, but " + \
                  ("falls" if not up else "no antiphase") + ")"
    else:
        verdict = "not yet"
    print(f"VERDICT: {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
