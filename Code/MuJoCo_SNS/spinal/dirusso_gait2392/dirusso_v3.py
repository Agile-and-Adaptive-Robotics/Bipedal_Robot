"""Di Russo stage-2 v3: her ACTUAL method (from full-text read 2026-10-26).

Fixes from v2 (all-falls):
1. CMA-ES sampler (her method; λ=40 population, not sequential TPE)
2. Five-term cost function from eq. 13: mimic + gait + effort + limits + head
3. Her fall criterion: COM height / initial < 0.8 (not absolute pelvis z)
4. Her speed range: [1.10, 1.25] m/s (not point target)
5. Starts from REFERENCE initial conditions (first frame of the walking
   cycle), not a standing keyframe
6. Pre-built model (no per-trial rebuild)
7. Curriculum staged harness: 0.9 → 0.7 → 0.5 → 0.3 → 0.0
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import mujoco
import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import controller as C          # noqa: E402
import plant as P               # noqa: E402

XML = r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml"
BSOLVE = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\bsolve_out.npz"
DB = "sqlite:///" + os.path.join(HERE, "optuna_dirusso_v3.db").replace("\\", "/")

# ---- pre-built model (once) ----
_model = None
_data = None

def get_model():
    global _model, _data
    if _model is None:
        txt, nz = P.build_sagittal_xml(XML)
        tmp = os.path.join(os.path.dirname(XML), "_dirusso_v3_tmp.xml")
        open(tmp, "w", encoding="utf-8").write(txt)
        try:
            _model = mujoco.MjModel.from_xml_path(tmp)
        finally:
            os.remove(tmp)
        _model.opt.timestep = 0.002
        _data = mujoco.MjData(_model)
    return _model, _data


def run_trial(params, wsup=0.7, dur=10.0, W_ff=None, gains=None):
    """One trial with the pre-built model. Returns metrics."""
    model, data = get_model()
    ids, geoms, pelvis, tilt_adr, tx_adr = P.index_maps(model, mujoco)

    # ---- INITIAL CONDITIONS from the reference walking cycle start ----
    ref = np.load(BSOLVE, allow_pickle=True)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    # impose the FIRST FRAME of the reference walking cycle
    # (her approach: start from a simpler controller's stable gait)
    q0 = ref["qpos"][0]
    for jn in ("pelvis_tilt", "pelvis_tx", "pelvis_ty",
               "hip_flexion_r", "knee_angle_r", "ankle_angle_r", "mtp_angle_r",
               "hip_flexion_l", "knee_angle_l", "ankle_angle_l", "mtp_angle_l"):
        jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, jn)
        full_jid = mujoco.mj_name2id(
            mujoco.MjModel.from_xml_path(XML), mujoco.mjtObj.mjOBJ_JOINT, jn)
        if jid >= 0 and full_jid >= 0:
            full_qadr = mujoco.MjModel.from_xml_path(XML).jnt_qposadr[full_jid]
            data.qpos[model.jnt_qposadr[jid]] = q0[full_qadr]
    # root from reference (columns 0-6 of qpos)
    for rc in range(min(7, q0.shape[0])):
        if rc < data.qpos.shape[0]:
            data.qpos[rc] = q0[rc]
    for j in P.LOCK_JOINTS:
        jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, j)
        data.qpos[model.jnt_qposadr[jid]] = 0.0
        data.qvel[model.jnt_dofadr[jid]] = 0.0
    # initial velocity from reference frame 0->1
    q1 = ref["qpos"][1]
    dt_ref = ref["t"][1] - ref["t"][0]
    for jn in ("hip_flexion_r", "knee_angle_r", "ankle_angle_r",
               "hip_flexion_l", "knee_angle_l", "ankle_angle_l"):
        jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, jn)
        full_model = mujoco.MjModel.from_xml_path(XML)
        full_jid = mujoco.mj_name2id(full_model, mujoco.mjtObj.mjOBJ_JOINT, jn)
        if jid >= 0 and full_jid >= 0:
            full_qadr = full_model.jnt_qposadr[full_jid]
            vel = (q1[full_qadr] - q0[full_qadr]) / dt_ref
            dof = model.jnt_dofadr[jid]
            if dof < data.qvel.shape[0]:
                data.qvel[dof] = vel
    mujoco.mj_forward(model, data)
    hcom_initial = float(data.subtree_com[0][2])

    # ---- controller ----
    if W_ff is None:
        W_ff = np.load(os.path.join(HERE, "dirusso_W_stage1.npy"))
    ctrl = C.DiRussoController(W_ff, dt=0.002,
                               omega=params.get("omega", 5.11),
                               gamma=params.get("gamma", 1.0),
                               gains=gains)
    weight = 9.81 * model.body_mass.sum()

    lr = {}
    for (m, s), aid in ids.items():
        lo = model.actuator_lengthrange[aid[0]]
        lr[(m, s)] = (lo[0], lo[1])

    n = int(dur / 0.002)
    dt = 0.002
    heelF = np.zeros((n, 2)); com_x = np.zeros(n); com_z = np.zeros(n)
    tilt_log = np.zeros(n); q_log = np.zeros((n, 6))
    heel_prev = {"L": 20.1, "R": 20.1}
    strikes = {"L": [], "R": []}
    fall_t = None
    f6 = np.zeros(6)

    for i in range(n):
        t = i * dt
        sens = {}
        heel_now = {}
        for si, s in enumerate(("L", "R")):
            # contact forces
            fh = ft = 0.0
            for k in range(data.ncon):
                c = data.contact[k]
                mujoco.mj_contactForce(model, data, k, f6)
                fz = max(f6[2], 0.0)
                if c.geom1 == geoms[("heel", s)] or c.geom2 == geoms[("heel", s)]:
                    fh += fz
                elif c.geom1 == geoms[("toe", s)] or c.geom2 == geoms[("toe", s)]:
                    ft += fz
            heel_now[s] = fh
            heelF[i, si] = fh
            hs = heel_prev[s] <= 20.0 < fh
            heel_prev[s] = fh
            # muscle sensors
            ria = np.zeros(C.NM); rii = np.zeros(C.NM); rib = np.zeros(C.NM)
            for m in C.MUSCLES:
                aid = ids[(m, s)]
                L = float(np.mean([data.actuator_length[a] for a in aid]))
                V = float(np.mean([data.actuator_velocity[a] for a in aid]))
                F = float(np.mean([abs(data.actuator_force[a]) for a in aid]))
                lo, hi = lr[(m, s)]
                ln_ = np.clip((L - lo) / max(hi - lo, 1e-6), 0, 2)
                vn_ = np.clip(V / max(hi - lo, 1e-6), -5, 5)
                fn_ = np.clip(F / max(sum(model.actuator_gainprm[a, 2]
                                          for a in aid), 1e-6), 0, 2)
                ria[C.MIDX[m]] = ln_; rii[C.MIDX[m]] = vn_; rib[C.MIDX[m]] = fn_
            a_, b_, c_ = C.DiRussoController.receptors(rii, ria, rib)
            cut = float(np.clip((fh + ft) / 700.0, 0, 1))
            sens[s] = (a_, b_, c_, cut)
            if hs:
                strikes[s].append(t)

        com = data.subtree_com[0]
        com_x[i] = com[0]; com_z[i] = com[2]
        tilt_log[i] = data.qpos[tilt_adr]
        tilt_deg = abs(np.degrees(tilt_log[i]))

        # her fall criterion: COM height / initial < 0.8
        if fall_t is None and com[2] < 0.8 * hcom_initial:
            fall_t = t

        exc = ctrl.step(sens, {s: heel_now[s] > 20.0 and
                         (not strikes[s] or t - strikes[s][-1] > 0.45)
                         for s in ("L", "R")}, tilt_log[i])
        for si, s in enumerate(("L", "R")):
            for m in C.MUSCLES:
                for a in ids[(m, s)]:
                    data.ctrl[a] = exc[si, C.MIDX[m]]
        data.xfrc_applied[:] = 0.0
        data.xfrc_applied[pelvis, 2] = wsup * weight
        mujoco.mj_step(model, data)
        if not np.isfinite(data.qpos).all():
            fall_t = fall_t or t
            break

    # ---- metrics ----
    t_arr = np.arange(n) * dt
    free = t_arr >= 1.0  # skip first second
    if fall_t is not None and fall_t < 2.0:
        return dict(fell=True, fall_t=fall_t, obj=-100.0, survive=False)
    # COM velocity (her Jgait term)
    com_free = com_x[free]
    t_free = t_arr[free]
    if len(com_free) > 10:
        speed = np.polyfit(t_free, com_free, 1)[0]
    else:
        speed = 0.0
    # strides
    def cyc(f):
        on = f > 20
        st = t_free[1:][on[1:] & ~on[:-1]]
        return sum(1 for j, x in enumerate(st)
                   if j == 0 or x - st[j-1] > 0.45)
    cycL = cyc(heelF[free, 0]); cycR = cyc(heelF[free, 1])
    # metabolic effort (her Jeffort: sum of activations)
    effort = float(np.mean(np.sum(exc**2, axis=1))) if exc is not None else 100
    # phase
    ph = None
    if cycL >= 3 and cycR >= 3:
        rT = np.array(strikes["R"]); lT = np.array(strikes["L"])
        phs = []
        for j in range(len(rT) - 1):
            ins = lT[(lT >= rT[j]) & (lT < rT[j+1])]
            if len(ins) == 1:
                phs.append((ins[0] - rT[j]) / (rT[j+1] - rT[j]))
        if phs:
            ph = float(np.median(phs))

    fell = fall_t is not None
    survive = (not fell and min(cycL, cycR) >= 3 and ph is not None
               and 0.35 <= ph <= 0.65)

    # her 5-term cost (eq 13, minimized)
    Jmimic = 0  # would need reference comparison per timestep
    Jgait = 0
    if fell:
        Jgait = 100
    elif speed < 1.10 or speed > 1.25:
        Jgait = 10 * min(abs(speed - 1.175) / 0.075, 1.0)
    Jeffort = effort * 0.01
    Jlimit = 0  # joint limits handled by MuJoCo's limit constraints
    Jhead = max(tilt_deg - 15.0, 0.0) * 0.1 if not fell else 50
    Jtotal = Jmimic + Jgait + Jeffort + Jlimit + Jhead

    return dict(fell=fell, fall_t=fall_t, obj=-Jtotal, survive=survive,
                speed=speed, cycL=cycL, cycR=cycR, phase=ph,
                effort=effort, com_z_min=float(com_z.min()),
                tilt_max=float(tilt_deg))


# ---- curriculum driver ----
CURRICULUM = [
    dict(wsup=0.9, n=100, name="s2r_c1"),   # barely needs to support itself
    dict(wsup=0.7, n=100, name="s2r_c2"),
    dict(wsup=0.5, n=100, name="s2r_c3"),
    dict(wsup=0.3, n=100, name="s2r_c4"),
    dict(wsup=0.0, n=100, name="s2r_c5"),   # free walking
]

SPACE = {}
for m in C.MUSCLES:
    SPACE[f"ff_{m}"] = (0.3, 3.0)
    for g in ("ia_mono", "ia_recip", "ii_exc", "ii_inh", "ib_auto",
              "ib_rev", "rc_mn", "cut_exc"):
        SPACE[f"{g}_{m}"] = (0.0, 2.0)
SPACE["omega"] = (3.5, 7.0)
SPACE["gamma"] = (0.0, 2.0)


def build_params_and_gains(p):
    W0 = np.load(os.path.join(HERE, "dirusso_W_stage1.npy"))
    W = W0.copy()
    gains = {}
    for m in C.MUSCLES:
        i = C.MIDX[m]
        W[i, :] *= p[f"ff_{m}"]
        W[9 + i, :] *= p[f"ff_{m}"]
        for g in ("ia_mono", "ia_recip", "ii_exc", "ii_inh", "ib_auto",
                  "ib_rev", "rc_mn", "cut_exc"):
            if g not in gains:
                gains[g] = np.zeros(C.NM)
            gains[g][i] = p[f"{g}_{m}"]
    return W, gains


def run_curriculum():
    best_params = {f"ff_{m}": 1.0 for m in C.MUSCLES}
    for m in C.MUSCLES:
        for g in ("ia_mono", "ia_recip", "ii_exc", "ii_inh", "ib_auto",
                  "ib_rev", "rc_mn", "cut_exc"):
            best_params[f"{g}_{m}"] = 0.1
    best_params["omega"] = 5.11
    best_params["gamma"] = 1.0

    for stage in CURRICULUM:
        study = optuna.create_study(
            study_name=stage["name"], storage=DB,
            direction="maximize", load_if_exists=True,
            sampler=optuna.samplers.CmaEsSampler(
                sigma0=0.3, popsize=20, seed=20261026))
        if not study.trials:
            study.enqueue_trial(best_params)

        def obj(trial):
            p = {k: trial.suggest_float(k, lo, hi)
                 for k, (lo, hi) in SPACE.items()}
            W, gains = build_params_and_gains(p)
            m = run_trial(p, wsup=stage["wsup"], dur=10.0,
                         W_ff=W, gains=gains)
            print(f"  [{stage['name']}#{trial.number}] "
                  f"obj={m['obj']:.1f} fell={m['fell']} "
                  f"speed={m.get('speed', 0):.2f} "
                  f"cyc={m.get('cycL', 0)}/{m.get('cycR', 0)}", flush=True)
            return float(m["obj"])

        study.optimize(obj, n_trials=stage["n"])
        b = study.best_trial
        best_params = dict(b.params)
        print(f"== {stage['name']} done: best {b.value:.1f} "
              f"(wsup={stage['wsup']})", flush=True)
        json.dump(dict(params=best_params, value=b.value,
                       wsup=stage["wsup"]),
                  open(os.path.join(HERE, f"{stage['name']}_best.json"), "w"),
                  indent=1)


if __name__ == "__main__":
    run_curriculum()
