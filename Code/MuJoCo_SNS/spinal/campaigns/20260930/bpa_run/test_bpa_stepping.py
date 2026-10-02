r"""BPA-ACTUATION GATE (campaign 2026-09-30/10-01): the split-RG W2L net
(build_w2l_split_net, M4 corrected-V3 commissurals) driving the w2l body
whose 12 Hill muscles are REPLACED by Ben's BPA model (bpa_muscle.py) on
the same tendon routes (make_bpa_walker.py -> w2l_mjcf_bpa.xml).

Protocol = Ben's virtual-walker drop (verbatim from test_li_stepping.py):
  [0, T_HOLD)      harness holds the root CLEAR of contact while the CPG
                   air-steps (feedforward weight + vertical PD, weak xy damp)
  [T_HOLD, +T_LOW) linear descent to contact height
  after            harness RETAINED at SUPPORT fraction of body weight
                   ("drop it WITH the harness", Ben 2026-09-30)

GATE (the campaign ask):
  air   : >= 3 air-stepping bursts WITH knee swing on BOTH legs
          (knee-flexion excursions, test_w2l_air swing_peaks convention)
  stand : retained-support stand where feet bear load
          (height held > 0.8 m + per-foot mean normal force > 20 N over
           the retained phase)

Pressure ceiling 620 kPa: activation a in [0,1] -> a*620 kPa; --cap
(default 0.5, the air-gate ctrl cap) keeps peaks <= cap*620 kPa.

Run (easteregg2):  D:\Anaconda\envs\myo\python.exe test_bpa_stepping.py [knobs]
  knobs: --dur=20 --hold=6 --te=3 --tf=4 --tau=0.25 --cap=0.5
         --stand-te=3 --stand-tf=4  (post-drop drive; set tf=0 for a
         latched-extensor quiet stand instead of marching)
         --dia=20 --count=1 --kmax-frac=0.25 --tendon=0.04 --bias=0.02
         --damp-bpa=0 --support=0.7 --lower=1.0 --clear=0.04
"""
from __future__ import annotations

import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
# easteregg2 env (remote); harmless no-op elsewhere
os.environ.setdefault("CONDA_PREFIX", r"D:\Anaconda\envs\myo")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SPINAL = os.path.dirname(HERE)
ROOT = os.path.dirname(SPINAL)              # Code\MuJoCo_SNS
for p in (HERE, SPINAL, ROOT):
    if p not in sys.path:
        sys.path.insert(0, p)

import mujoco  # noqa: E402

import build_w2l_split_net as B  # noqa: E402
from bpa_muscle import BPAMuscle, maxBPAforce, tendon_spring_rate  # noqa: E402
from bpa_mujoco import BPAMuscleSystem  # noqa: E402

MJCF = os.path.join(HERE, "w2l_mjcf_bpa.xml")
DT_PHY = 0.001
NSUB = 2                       # physics steps per 2 ms net step (net dt=DT)

KV = dict(dur=20.0, hold=6.0, lower=1.0, clear=0.04, support=0.7,
          te=3.0, tf=4.0, tau=0.25, cap=0.5, comm=1.0,
          stand_te=None, stand_tf=None,
          dia=20.0, count=1.0, kmax_frac=0.25, tendon=0.04, bias=0.02,
          damp_bpa=0.0, stiff=1.0, servo_rate=0.01)
for arg in sys.argv[1:]:
    if arg.startswith("--") and "=" in arg:
        k, v = arg[2:].split("=", 1)
        KV[k.replace("-", "_")] = float(v)
if KV["stand_te"] is None:
    KV["stand_te"] = KV["te"]
if KV["stand_tf"] is None:
    KV["stand_tf"] = KV["tf"]

JNTS = ("hip_L", "knee_L", "ankle_L", "hip_R", "knee_R", "ankle_R")
PRESSURE_MAX = 620.0


def swing_peaks(sig: np.ndarray, min_gap_ms: int = 400):
    """Onsets of sustained flexion excursions (test_w2l_air convention)."""
    thr = sig.mean() + 0.3 * sig.std()
    on = sig > thr
    starts = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    if len(starts) == 0:
        return starts
    keep = [starts[0]]
    for s in starts[1:]:
        if s - keep[-1] > min_gap_ms:
            keep.append(s)
    return np.array(keep)


def foot_forces(m, d, g_ids):
    """Per-geom vertical contact normal force (N) + stance flags."""
    f = np.zeros(4)
    for i in range(d.ncon):
        con = d.contact[i]
        for k, gid in enumerate(g_ids):
            if con.geom1 == gid or con.geom2 == gid:
                f6 = np.zeros(6)
                mujoco.mj_contactForce(m, d, i, f6)
                f[k] += f6[0]
    return f


def main() -> int:
    # ------------------------------------------------------------------ net
    B.NAP["tau_max_h"] = KV["tau"]
    tpl = os.path.join(HERE, "connectome_templates.json")
    net = B.build(comm=KV["comm"],
                  template_path=tpl if os.path.exists(tpl) else None)

    # ------------------------------------------------------------------ body
    m = mujoco.MjModel.from_xml_path(MJCF)
    d = mujoco.MjData(m)
    mujoco.mj_resetDataKeyframe(m, d, 0)
    mujoco.mj_forward(m, d)

    act_ids = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
               for a in net.muscle_outputs}
    assert len(act_ids) == 12, f"actuator count {len(act_ids)} != 12"
    # stiff joint limits (test_w2l_air.py stand-in for missing muscle
    # damping; the BPA force is bounded but the BPA has NO length stop of
    # its own, so the LIMIT must not yield - measured: without this the
    # ankles blow through their range to +115 deg and the routes stretch
    # past the sampled max)
    if KV["stiff"] > 0.0:
        m.jnt_solimp[:] = np.array([0.9, 0.99, 0.001, 0.5, 2.0])
        m.jnt_solref[:] = np.array([0.006, 1.0])
    jnt_ids = {j: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)
               for j in JNTS}
    qadr = {j: m.jnt_qposadr[i] for j, i in jnt_ids.items()}
    g_ids = [mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g) for g in
             ("foot_L_contact", "toe_L_contact",
              "foot_R_contact", "toe_R_contact")]

    # ------------------------------------------------------------ BPA sizing
    # One BPA muscle-tendon unit per existing route. SIZING LAW: rest =
    # MAX route length over the side's joint ranges - tendon, so the BPA
    # strain stays in [0,1] for every reachable pose (zero strain at the
    # joint extreme where the route is longest). festo4() grows
    # exponentially for NEGATIVE strain (stretch) - sizing to the max
    # route keeps forces bounded by p*Fmax by construction (measured: a
    # spawn-sized rest put knee forces at 13 kN and blew the sim up).
    ranges = {}
    for side in ("L", "R"):
        for j in ("hip", "knee", "ankle"):
            jid = jnt_ids[f"{j}_{side}"]
            adr = m.jnt_qposadr[jid]
            ranges[(side, j)] = (m.jnt_range[jid][0], m.jnt_range[jid][1],
                                 adr)
    n_grid = 7
    max_len = {a: 0.0 for a in act_ids}
    min_len = {a: 1e9 for a in act_ids}
    q0 = d.qpos.copy()
    for side in ("L", "R"):
        axes = ("hip", "knee", "ankle")
        for iH in range(n_grid):
            for iK in range(n_grid):
                for iA in range(n_grid):
                    for j, ix in zip(axes, (iH, iK, iA)):
                        lo, hi, adr = ranges[(side, j)]
                        d.qpos[adr] = lo + (hi - lo) * ix / (n_grid - 1)
                    mujoco.mj_forward(m, d)
                    for a, aid in act_ids.items():
                        if f"_{side}_" in a:
                            L = float(d.actuator_length[aid])
                            max_len[a] = max(max_len[a], L)
                            min_len[a] = min(min_len[a], L)
    d.qpos[:] = q0
    mujoco.mj_resetDataKeyframe(m, d, 0)
    mujoco.mj_forward(m, d)

    muscles = {}
    print("== BPA units (rest = max sampled route - tendon; 7^3 pose grid "
          "per side over joint ranges) ==")
    print(f"  {'name':<14}{'route min..max':>20}{'spawn':>8}{'rest':>8}"
          f"{'Fmax':>8}{'kmax_len':>9}{'rel_spawn':>10}")
    for a, aid in act_ids.items():
        tln = KV["tendon"]
        rest = max_len[a] - tln + KV["bias"]
        kmf = KV["kmax_frac"]
        mus = BPAMuscle(
            name=a, diameter=KV["dia"], resting_length=rest,
            kmax_length=rest * (1.0 - kmf), tendon_length=tln,
            bpa_count=int(KV["count"]), pressure_max_kpa=PRESSURE_MAX,
            extra_damping=KV["damp_bpa"])
        muscles[a] = mus
        rel0 = (rest - (float(d.actuator_length[aid]) - tln)) / rest / kmf
        print(f"  {a:<14}{min_len[a]:10.4f}..{max_len[a]:.4f}"
              f"{float(d.actuator_length[aid]):8.4f}{rest:8.4f}"
              f"{mus.bpa_count * mus.fmax:8.0f}{mus.kmax_length:9.4f}"
              f"{rel0:10.3f}")
    system = BPAMuscleSystem(muscles)
    system.attach(m)
    mujoco.set_mjcb_control(system.cb)

    # ------------------------------------------------------------- harness
    T_HOLD, T_LOW, CLEAR, SUPPORT = (KV["hold"], KV["lower"],
                                     KV["clear"], KV["support"])
    root_bid = int(m.jnt_bodyid[0])
    _mass = float(np.sum(m.body_mass))
    _wgt = _mass * (-m.opt.gravity[2])
    K_H = 4.0 * _wgt / max(CLEAR, 0.01)
    C_H = 2.0 * float(np.sqrt(K_H * _mass))
    z_contact = float(d.qpos[2])
    z_hold = z_contact + CLEAR
    # RETAINED-PHASE FORCE SERVO (documented deviation from the open-loop
    # "lower to z_contact"): the spawn keyframe has feet 1 mm clear at the
    # SPAWN pose, but the march/stand pose differs (hips swing to
    # extension, knees cycle), so holding z_contact leaves the feet in the
    # air (measured: 0 N over 11.5 s). The AnimatLab virtual walker
    # descends until the feet BEAR on the platform - replicated here by
    # servoing the PD target down/up at SERVO_RATE until the feet carry
    # (1-SUPPORT)*weight, then holding that height.
    LOAD_TGT = (1.0 - SUPPORT) * _wgt
    z_stand = z_contact
    load_prev = 0.0
    SERVO_RATE = KV["servo_rate"]
    print(f"== harness: mass {_mass:.1f} kg, weight {_wgt:.0f} N, "
          f"z_contact {z_contact:.4f}, z_hold {z_hold:.4f}, "
          f"load target {LOAD_TGT:.0f} N ==")

    # ------------------------------------------------------------------ run
    dur = KV["dur"]
    nsteps = int(dur / DT_PHY)
    log_t = np.zeros(nsteps)
    log_q = np.zeros((nsteps, 6))
    log_f = np.zeros((nsteps, 4))          # vertical normal force per geom
    log_h = np.zeros((nsteps, 2))          # pelvis z, x
    log_p = np.zeros((nsteps, 12))         # pressure per muscle (kPa)
    log_rgE = np.zeros(nsteps)
    u = net.make_inputs()
    iE = net.input_index("TONIC L RG ext")
    iF = net.input_index("TONIC L RG flx")
    iEr = net.input_index("TONIC R RG ext")
    iFr = net.input_index("TONIC R RG flx")
    iS = net.input_index("Stimulus_1")
    iS2 = net.input_index("Stimulus_2")

    for i in range(nsteps):
        t = i * DT_PHY
        if i % NSUB == 0:
            te_now = KV["te"] if t < T_HOLD + T_LOW else KV["stand_te"]
            tf_now = KV["tf"] if t < T_HOLD + T_LOW else KV["stand_tf"]
            u[iE], u[iEr] = te_now, te_now
            u[iF], u[iFr] = tf_now, tf_now
            u[iS] = 10.0 if t < 0.01 else 0.0     # verbatim kickoff
            u[iS2] = 10.0 if t < 0.01 else 0.0    # antiphase kickoff
            V = net.step(u)
            log_rgE[i] = V[net.idx["L RG ext"]]
            for a, c in net.muscle_ctrl(V).items():
                system.set_activation(a, min(c, KV["cap"]))
        else:
            log_rgE[i] = log_rgE[i - 1]

        if t < T_HOLD:
            z_t, wff = z_hold, 1.0
        elif t < T_HOLD + T_LOW:
            z_t = z_hold + (z_contact - z_hold) * (t - T_HOLD) / max(T_LOW, 1e-6)
            wff = 1.0
        else:
            z_t, wff = z_stand, SUPPORT
            if load_prev < LOAD_TGT * 0.9:
                z_stand = max(z_stand - SERVO_RATE * DT_PHY, 0.55)
            elif load_prev > LOAD_TGT * 1.15:
                z_stand = min(z_stand + SERVO_RATE * DT_PHY, z_contact)
        Fz = wff * _wgt + K_H * (z_t - d.qpos[2]) - C_H * d.qvel[2]
        d.xfrc_applied[root_bid, 2] = max(Fz, 0.0)
        d.xfrc_applied[root_bid, 0] = -20.0 * d.qvel[0]
        d.xfrc_applied[root_bid, 1] = -20.0 * d.qvel[1]

        mujoco.mj_step(m, d)
        log_t[i] = t
        log_q[i] = [np.degrees(d.qpos[qadr[j]]) for j in JNTS]
        log_f[i] = foot_forces(m, d, g_ids)
        load_prev = float(log_f[i, 0] + log_f[i, 1] + log_f[i, 2] + log_f[i, 3])
        log_h[i] = d.qpos[2], d.qpos[0]
        log_p[i] = [muscles[a].pressure_kpa for a in sorted(muscles)]
        if not np.all(np.isfinite(d.qpos)):
            print(f"[FAIL] qpos non-finite at t={t:.3f} s")
            return 1

    mujoco.set_mjcb_control(None)

    # ============================================================== metrics
    t_air = (log_t >= 1.0) & (log_t < T_HOLD)      # drop kickoff transient
    t_ret = log_t >= (T_HOLD + T_LOW + 0.5)        # after landing transient

    print(f"== BPA actuation gate: split-RG W2L + BPA muscles, {dur:.0f} s ==")
    print(f"   knobs: te={KV['te']} tf={KV['tf']} tau={KV['tau']} "
          f"cap={KV['cap']} dia={KV['dia']:.0f} count={int(KV['count'])} "
          f"kmax_frac={KV['kmax_frac']} tendon={KV['tendon']} "
          f"bias={KV['bias']} damp_bpa={KV['damp_bpa']} "
          f"support={SUPPORT} stand_te={KV['stand_te']} "
          f"stand_tf={KV['stand_tf']}")

    # --- air phase: knee swing bursts
    flex_air = -log_q[t_air][:, (1, 4)]            # knee flexion-positive
    bursts = {}
    for side, col in (("L", 0), ("R", 1)):
        pk = swing_peaks(flex_air[:, col])
        bursts[side] = pk
        rng = flex_air[:, col].max() - flex_air[:, col].min()
        per = float(np.diff(pk).mean() * DT_PHY) if len(pk) >= 3 else float("nan")
        print(f"   AIR knee {side}: {len(pk)} flexion excursions, "
              f"period {per:.3f} s, flexion range {rng:.1f} deg "
              f"[{flex_air[:, col].min():+.1f},{flex_air[:, col].max():+.1f}]")
    hip_air = -log_q[t_air][:, (0, 3)]
    for side, col in (("L", 0), ("R", 1)):
        rng = hip_air[:, col].max() - hip_air[:, col].min()
        print(f"   AIR hip {side}: flexion range {rng:.1f} deg")
    rg = log_rgE[t_air]
    thr = 0.5 * rg.max()
    on = rg > thr
    st = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    keep = [s for i2, s in enumerate(st)
            if i2 == 0 or (s - st[i2 - 1]) * DT_PHY > 0.4]
    rg_per = (float(np.diff(keep).mean() * DT_PHY)
              if len(keep) >= 3 else float("nan"))
    print(f"   AIR neural: L RG-E bursts {len(keep)}, period {rg_per:.3f} s, "
          f"max {rg.max():.2f} mV")

    # --- retained phase: stand + load
    def episodes(on, tt, gap=0.3):
        idx = np.where(on)[0]
        if len(idx) == 0:
            return []
        eps = [[idx[0], idx[0]]]
        for i2 in idx[1:]:
            if (tt[i2] - tt[eps[-1][1]]) < gap:
                eps[-1][1] = i2
            else:
                eps.append([i2, i2])
        return [(tt[a], tt[b]) for a, b in eps if tt[b] - tt[a] > 0.05]

    h = log_h[t_ret]
    q = log_q[t_ret]
    f = log_f[t_ret]
    p = log_p[t_ret]
    tt_ret = log_t[t_ret]
    footL = f[:, 0] + f[:, 1]
    footR = f[:, 2] + f[:, 3]
    print(f"   RETAINED phase ({tt_ret[-1] - tt_ret[0]:.1f} s of it): "
          f"pelvis z min {h[:, 0].min():.3f} m, end {h[-1, 0]:.3f} m")
    print(f"   RETAINED vertical load: mean L {footL.mean():.1f} N / "
          f"R {footR.mean():.1f} N; PEAK L {footL.max():.0f} N / "
          f"R {footR.max():.0f} N; "
          f"stance episodes L {len(episodes(footL > 5.0, tt_ret))} / "
          f"R {len(episodes(footR > 5.0, tt_ret))} "
          f"(harness carries {SUPPORT * 100:.0f}% of {_wgt:.0f} N)")
    for side, j in (("L", (0, 1, 2)), ("R", (3, 4, 5))):
        hip, knee, ank = q[:, j[0]], q[:, j[1]], q[:, j[2]]
        print(f"   RETAINED joints {side}: hip {hip.min():+.1f}..{hip.max():+.1f} "
              f"knee {knee.min():+.1f}..{knee.max():+.1f} "
              f"ankle {ank.min():+.1f}..{ank.max():+.1f} deg")
    print(f"   BPA pressure: peak {p.max():.0f} kPa (ceiling "
          f"{PRESSURE_MAX:.0f}), mean-of-active {p[p > 1].mean():.0f} kPa")
    fmax_now = {a: muscles[a].force(float(d.actuator_length[act_ids[a]]))
                for a in sorted(muscles)}
    print("   BPA force at final pose (N): "
          + " ".join(f"{a}={fmax_now[a]:.0f}" for a in sorted(fmax_now)))

    # --- verdicts
    air_ok = all(len(bursts[s]) >= 3 for s in ("L", "R")) and all(
        flex_air[:, c].max() - flex_air[:, c].min() >= 10.0 for c in (0, 1))
    up = h[-1, 0] > 0.8 and h[:, 0].min() > 0.8
    loaded = (footL.mean() > 20.0 and footR.mean() > 20.0) or \
        (len(episodes(footL > 5.0, tt_ret)) >= 3 and
         len(episodes(footR > 5.0, tt_ret)) >= 3 and
         max(footL.max(), footR.max()) > 50.0)
    finite = np.all(np.isfinite(log_q))
    for k, v in (("finite", finite),
                 ("air_bursts_ge3_both_knees", air_ok),
                 ("stand_height_held", up),
                 ("feet_bear_load", loaded),
                 ("pressure_le_620kPa", p.max() <= PRESSURE_MAX + 1e-9)):
        print(f"   [{'PASS' if v else 'FAIL'}] {k}")
    if finite and air_ok and up and loaded:
        print("VERDICT: BPA ACTUATION PASS "
              "(air-stepping + retained-support stand on BPA muscles)")
        return 0
    print("VERDICT: BPA ACTUATION not yet")
    return 1


if __name__ == "__main__":
    sys.exit(main())
