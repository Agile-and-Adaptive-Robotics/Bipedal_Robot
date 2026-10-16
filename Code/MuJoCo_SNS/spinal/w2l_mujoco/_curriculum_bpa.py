"""BPA-WALKER curriculum, non-spiking (2026-10-03 overnight campaign).

Ben's ask: "the Robot with BPA actuated by spiking and non-spiking
models. All of it tuned." This module tunes the NON-SPIKING split-RG W2L
controller (build_w2l_split_net) driving the BPA-actuated w2l body
(w2l_mjcf_bpa.xml, gate-passing 2026-09-30: campaigns/20260930/
bpa_actuation.md) under the virtual-walker protocol (air-step held
clear -> descent -> retained-support march).

Protocol + BPA sizing law adapted VERBATIM from test_bpa_stepping.py
(2026-09-30, gate PASS) into a parameterized run_bpa(net, p) so the
SPIKING twin (_curriculum_bpa_spiking.py) can reuse the exact same
plant loop with a different net. Gains in the template are Ben's
literature-derived starting points, NOT fixed - they are searched here
at the controller level (te/tf/tau/comm) and the actuator-interface
level (cap/kmax_frac/tendon/count), per Ben's 2026-10-03 ruling.

Score (maximize; graded gate not pass/fail):
  air  : 3*(bursts_L+bursts_R) + 0.4*(knee_range_L+R) + period
         regularity bonus; RHYTHM GATE (<3 bursts either knee -> ~0)
  stand: height held + 4*(stance episodes L+R) + load evenness +
         pressure headroom; collapse / non-finite -> sentinels below.

Usage: python _curriculum_bpa.py [n_trials]
  (run from w2l_mujoco/ with the myo python)
"""
from __future__ import annotations

import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
SPINAL = os.path.dirname(HERE)
ROOT = os.path.dirname(SPINAL)
for p in (HERE, SPINAL, ROOT):
    if p not in sys.path:
        sys.path.insert(0, p)
os.chdir(HERE)                     # relative MJCF/template paths
DB = f"sqlite:///{SPINAL}/optuna_bpa.db"
STUDY = "curr_bpa_ns_20261003"

import numpy as np                 # noqa: E402
import optuna                      # noqa: E402
import mujoco                      # noqa: E402

import build_w2l_split_net as B    # noqa: E402
from bpa_muscle import BPAMuscle   # noqa: E402
from bpa_mujoco import BPAMuscleSystem  # noqa: E402

MJCF = os.path.join(HERE, "w2l_mjcf_bpa.xml")
TPL = os.path.join(SPINAL, "connectome_templates.json")
DT_PHY = 0.001
NSUB = 2
JNTS = ("hip_L", "knee_L", "ankle_L", "hip_R", "knee_R", "ankle_R")
PRESSURE_MAX = 620.0

# gate-passing config (campaigns/20260930/bpa_actuation.md) = the seed
SEED = dict(te=3.0, tf=4.0, tau=0.25, cap=0.6, comm=1.0,
            stand_te=3.0, stand_tf=4.0, kmax_frac=0.6, tendon=0.04,
            bias=0.02, count=2, support=0.7, dia=20.0, damp_bpa=0.0)


def swing_peaks(sig: np.ndarray, min_gap_ms: int = 400) -> np.ndarray:
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
    f = np.zeros(4)
    for i in range(d.ncon):
        con = d.contact[i]
        for k, gid in enumerate(g_ids):
            if con.geom1 == gid or con.geom2 == gid:
                f6 = np.zeros(6)
                mujoco.mj_contactForce(m, d, i, f6)
                f[k] += f6[0]
    return f


def run_bpa(net, p, save_npz=None, verbose=False):
    """One protocol run. p keys: te tf tau cap stand_te stand_tf
    kmax_frac tendon bias count support dia damp_bpa dur hold lower.
    Returns a metrics dict (empty dict on non-finite blowup)."""
    dur = float(p.get("dur", 20.0))
    T_HOLD = float(p.get("hold", 6.0))
    T_LOW = float(p.get("lower", 2.0))
    SUPPORT = float(p["support"])
    m = mujoco.MjModel.from_xml_path(MJCF)
    d = mujoco.MjData(m)
    mujoco.mj_resetDataKeyframe(m, d, 0)
    mujoco.mj_forward(m, d)
    act_ids = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
               for a in net.muscle_outputs}
    m.jnt_solimp[:] = np.array([0.9, 0.99, 0.001, 0.5, 2.0])
    m.jnt_solref[:] = np.array([0.006, 1.0])
    jnt_ids = {j: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)
               for j in JNTS}
    qadr = {j: m.jnt_qposadr[i] for j, i in jnt_ids.items()}
    g_ids = [mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g) for g in
             ("foot_L_contact", "toe_L_contact",
              "foot_R_contact", "toe_R_contact")]
    ranges = {}
    for side in ("L", "R"):
        for j in ("hip", "knee", "ankle"):
            jid = jnt_ids[f"{j}_{side}"]
            ranges[(side, j)] = (m.jnt_range[jid][0], m.jnt_range[jid][1],
                                 m.jnt_qposadr[jid])
    n_grid = 7
    max_len = {a: 0.0 for a in act_ids}
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
                            max_len[a] = max(max_len[a],
                                             float(d.actuator_length[aid]))
    d.qpos[:] = q0
    mujoco.mj_resetDataKeyframe(m, d, 0)
    mujoco.mj_forward(m, d)

    muscles = {}
    for a, aid in act_ids.items():
        tln = float(p["tendon"])
        rest = max_len[a] - tln + float(p["bias"])
        mus = BPAMuscle(
            name=a, diameter=float(p["dia"]),
            resting_length=rest,
            kmax_length=rest * (1.0 - float(p["kmax_frac"])),
            tendon_length=tln, bpa_count=int(p["count"]),
            pressure_max_kpa=PRESSURE_MAX,
            extra_damping=float(p.get("damp_bpa", 0.0)))
        muscles[a] = mus
    system = BPAMuscleSystem(muscles)
    system.attach(m)
    mujoco.set_mjcb_control(system.cb)

    root_bid = int(m.jnt_bodyid[0])
    _mass = float(np.sum(m.body_mass))
    _wgt = _mass * (-m.opt.gravity[2])
    K_H = 4.0 * _wgt / max(0.04, 0.01)
    C_H = 2.0 * float(np.sqrt(K_H * _mass))
    z_contact = float(d.qpos[2])
    z_hold = z_contact + 0.04
    LOAD_TGT = (1.0 - SUPPORT) * _wgt
    z_stand = z_contact
    load_prev = 0.0
    SERVO_RATE = 0.01

    nsteps = int(dur / DT_PHY)
    log_t = np.zeros(nsteps)
    log_q = np.zeros((nsteps, 6))
    log_f = np.zeros((nsteps, 4))
    log_h = np.zeros((nsteps, 2))
    log_p = np.zeros((nsteps, 12))
    log_rgE = np.zeros(nsteps)
    log_V = []
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
            te_now = float(p["te"]) if t < T_HOLD + T_LOW \
                else float(p["stand_te"])
            tf_now = float(p["tf"]) if t < T_HOLD + T_LOW \
                else float(p["stand_tf"])
            u[iE], u[iEr] = te_now, te_now
            u[iF], u[iFr] = tf_now, tf_now
            u[iS] = 10.0 if t < 0.01 else 0.0
            u[iS2] = 10.0 if t < 0.01 else 0.0
            V = net.step(u)
            log_rgE[i] = V[net.idx["L RG ext"]]
            for a, c in net.muscle_ctrl(V).items():
                system.set_activation(a, min(c, float(p["cap"])))
        else:
            log_rgE[i] = log_rgE[i - 1]
        if t < T_HOLD:
            z_t, wff = z_hold, 1.0
        elif t < T_HOLD + T_LOW:
            z_t = z_hold + (z_contact - z_hold) * (t - T_HOLD) / T_LOW
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
        load_prev = float(log_f[i].sum())
        log_h[i] = d.qpos[2], d.qpos[0]
        log_p[i] = [muscles[a].pressure_kpa for a in sorted(muscles)]
        if not np.all(np.isfinite(d.qpos)):
            mujoco.set_mjcb_control(None)
            return {"finite": False}
    mujoco.set_mjcb_control(None)

    t_air = (log_t >= 1.0) & (log_t < T_HOLD)
    t_ret = log_t >= (T_HOLD + T_LOW + 0.5)
    flex_air = -log_q[t_air][:, (1, 4)]
    bursts, brng, pers = {}, {}, {}
    for side, col in (("L", 0), ("R", 1)):
        pk = swing_peaks(flex_air[:, col])
        bursts[side] = len(pk)
        brng[side] = float(flex_air[:, col].max()
                           - flex_air[:, col].min())
        pers[side] = (float(np.diff(pk).mean() * DT_PHY)
                      if len(pk) >= 3 else float("nan"))
    rg = log_rgE[t_air]
    on = rg > 0.5 * max(rg.max(), 1e-9)
    st = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    keep = [s for i2, s in enumerate(st)
            if i2 == 0 or (s - st[i2 - 1]) * DT_PHY > 0.4]
    rg_per = (float(np.diff(keep).mean() * DT_PHY)
              if len(keep) >= 3 else float("nan"))
    rg_per_std = (float(np.std(np.diff(keep)) * DT_PHY)
                  if len(keep) >= 4 else float("nan"))

    def episodes(onv, tt, gap=0.3):
        idx = np.where(onv)[0]
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
    f = log_f[t_ret]
    tt_ret = log_t[t_ret]
    footL = f[:, 0] + f[:, 1]
    footR = f[:, 2] + f[:, 3]
    met = {
        "finite": True,
        "bursts_L": bursts["L"], "bursts_R": bursts["R"],
        "knee_rng_L": brng["L"], "knee_rng_R": brng["R"],
        "hip_rng_L": float(-log_q[t_air][:, 0].max()
                           - -log_q[t_air][:, 0].min()),
        "rg_bursts": len(keep), "rg_period": rg_per,
        "rg_period_std": rg_per_std,
        "h_min": float(h[:, 0].min()), "h_end": float(h[-1, 0]),
        "epi_L": len(episodes(footL > 5.0, tt_ret)),
        "epi_R": len(episodes(footR > 5.0, tt_ret)),
        "load_L": float(footL.mean()), "load_R": float(footR.mean()),
        "p_peak": float(log_p.max()),
        "up": bool(h[-1, 0] > 0.8 and h[:, 0].min() > 0.8),
    }
    if save_npz:
        np.savez_compressed(save_npz, t=log_t, q=log_q, f=log_f,
                            h=log_h, p=log_p, rgE=log_rgE,
                            metrics=json.dumps(met))
    if verbose:
        print(json.dumps(met, indent=1))
    return met


def score_bpa(met):
    if not met.get("finite", False):
        return -200.0
    # rhythm gate: a non-stepping config must not win (2026-09-20 law)
    if met["bursts_L"] < 3 or met["bursts_R"] < 3:
        return -10.0 + 0.5 * (met["knee_rng_L"] + met["knee_rng_R"])
    air = (3.0 * (met["bursts_L"] + met["bursts_R"])
           + 0.4 * (met["knee_rng_L"] + met["knee_rng_R"]))
    if np.isfinite(met["rg_period_std"]) and met["rg_period"] > 0:
        air += 10.0 * max(0.0, 1.0 - met["rg_period_std"]
                          / max(met["rg_period"], 1e-6))
    if not met["up"]:
        ret = -50.0                       # collapsed during retained
    else:
        ret = 30.0
        ret += 4.0 * (met["epi_L"] + met["epi_R"])
        sym = abs(met["load_L"] - met["load_R"]) / \
            max(met["load_L"] + met["load_R"], 1e-6)
        ret -= 20.0 * sym
        if met["load_L"] < 5.0 or met["load_R"] < 5.0:
            ret -= 30.0                   # one foot never bears
    if met["p_peak"] > PRESSURE_MAX:
        ret -= 100.0                      # pressure ceiling breach
    return air + ret


def _obj(trial):
    p = dict(
        te=trial.suggest_float("te", 2.0, 5.5),
        tf=trial.suggest_float("tf", 2.0, 5.5),
        tau=trial.suggest_float("tau", 0.15, 0.45),
        cap=trial.suggest_float("cap", 0.4, 0.8),
        comm=trial.suggest_float("comm", 0.5, 1.5),
        stand_te=trial.suggest_float("stand_te", 1.0, 6.0),
        stand_tf=trial.suggest_float("stand_tf", 0.0, 4.0),
        support=trial.suggest_float("support", 0.5, 0.85),
        kmax_frac=trial.suggest_float("kmax_frac", 0.35, 0.7),
        tendon=trial.suggest_float("tendon", 0.02, 0.06),
        count=trial.suggest_int("count", 1, 3),
        dia=20.0, bias=0.02, damp_bpa=0.0,
    )
    B.NAP["tau_max_h"] = p["tau"]
    net = B.build(comm=p["comm"], template_path=TPL)
    met = run_bpa(net, p)
    return score_bpa(met)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    st = optuna.create_study(direction="maximize", storage=DB,
                             study_name=STUDY, load_if_exists=True,
                             sampler=optuna.samplers.TPESampler(
                                 seed=31, n_startup_trials=8))
    if len(st.trials) == 0:
        st.enqueue_trial({k: v for k, v in SEED.items()
                          if k not in ("dia", "bias", "damp_bpa")})
        print("seeded", json.dumps(SEED), flush=True)
    st.optimize(_obj, n_trials=n, gc_after_trial=True)
    b = st.best_trial
    print(f"== BPA ns best {b.value:.3f} (trial {b.number})")
    print(json.dumps(b.params, indent=1))
    with open(os.path.join(SPINAL, "curriculum_bpa_ns.json"), "w",
              encoding="utf-8") as f:
        json.dump({"score": b.value, "params": b.params,
                   "trial": b.number, "study": STUDY, "db": DB,
                   "fixed": {k: SEED[k] for k in
                             ("dia", "bias", "damp_bpa")}}, f, indent=2)


if __name__ == "__main__":
    main()
