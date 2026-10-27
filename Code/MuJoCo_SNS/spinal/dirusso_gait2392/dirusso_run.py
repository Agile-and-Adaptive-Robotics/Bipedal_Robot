r"""Di Russo on gait2392 — driver (D2, 2026-10-25).

Stage 1 (her eq. 12 imitation): least-squares fit of the signed per-muscle
W against the reference walk's back-solved activations (bsolve_out.npz,
121 frames along subject01 IK), phased per side by measured GRF onsets —
exactly her "imitate a simulated behavior" step, our reference being the
SO/NNLS solution along the human walk.

Closed loop: sagittal-locked plant (plant.py), controller.py at the 2 ms
physics step, heel-strike phase reset from calcaneus normal force (>20 N
rising edge), vertical weight-support harness (wsup; her sim is unassisted
2D — wsup is the bring-up/wean knob, labeled adaptation), MuJoCo muscle
activation dynamics driven by MN excitation via ctrl.

Usage:
  python dirusso_run.py [--dur=30] [--wsup=0.7] [--omega=5.7] [--out=...]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import mujoco

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import controller as C          # noqa: E402
import plant as P               # noqa: E402

XML = r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml"
BSOLVE = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\bsolve_out.npz"
GRF_MOT = r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\subject01_walk1_grf.mot"


def grf_onsets():
    """Per-side heel-strike times from the measured GRF (rising >20 N)."""
    rows, names, header_done = [], [], False
    for ln in open(GRF_MOT, encoding="utf-8", errors="replace"):
        if ln.startswith("endheader"):
            header_done = True
            continue
        if not header_done or ln.startswith("time") or not ln.strip():
            if ln.startswith("time"):
                names = ln.split()
            continue
        rows.append([float(x) for x in ln.split()])
    a = np.array(rows)
    col = {n: i for i, n in enumerate(names)}
    t = a[:, 0]
    out = {}
    for side, key in (("R", "ground_force_vy"), ("L", "1_ground_force_vy")):
        f = a[:, col[key]]
        on = f > 20.0
        strikes = t[1:][on[1:] & ~on[:-1]]
        # merge double-taps closer than 0.45 s (kine_ref convention)
        merged = []
        for x in strikes:
            if not merged or (x - merged[-1]) > 0.45:
                merged.append(x)
        out[side] = np.array(merged)
    return t[0], t[-1], out


def stage1_fit():
    """Signed W [2*NM,5] by least squares vs reference group activations."""
    d = np.load(BSOLVE, allow_pickle=True)
    t, acts, anames = d["t"], d["acts"], list(d["act_names"])
    ai = {n: i for i, n in enumerate(anames)}
    t0, t1, onsets = grf_onsets()
    ref = np.zeros((len(t), 2, C.NM))
    for si, s in enumerate(("r", "l")):
        for m in C.MUSCLES:
            cols = [ai[n.format(s=s)] for n in C.MAP[m]]
            ref[:, si, C.MIDX[m]] = acts[:, cols].mean(axis=1)
    W = np.zeros((2 * C.NM, 5))
    fit_r2 = np.zeros((2, C.NM))
    for si, side in enumerate(("L", "R")):
        ons = onsets[side]
        for k in range(len(ons) - 1):
            seg = (t >= ons[k]) & (t < ons[k + 1])
            T = ons[k + 1] - ons[k]
            phi = ((t[seg] - ons[k]) / T) % 1.0
            Pm = np.stack([C.primitives(p) for p in phi])
            Wside, *_ = np.linalg.lstsq(Pm, ref[seg, si, :], rcond=None)
            W[si * C.NM:(si + 1) * C.NM] = Wside.T
            pred = Pm @ Wside
            for m in range(C.NM):
                ss = np.sum((ref[seg, si, m] - ref[seg, si, m].mean()) ** 2)
                fit_r2[si, m] = 1 - np.sum(
                    (ref[seg, si, m] - pred[:, m]) ** 2) / max(ss, 1e-9)
    return W, fit_r2


def load_model():
    """Patch + compile with the original's relative Geometry/ assets."""
    txt, n_zero = P.build_sagittal_xml(XML)
    tmp = os.path.join(os.path.dirname(XML), "_dirusso_sagittal_tmp.xml")
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(txt)
    try:
        model = mujoco.MjModel.from_xml_path(tmp)
    finally:
        os.remove(tmp)
    model.opt.timestep = 0.002   # frozen stack contract (cvt3 ships 0.005)
    print(f"[plant] sagittal locks: {len(P.LOCK_JOINTS)} joints; "
          f"non-her actuators zeroed: {n_zero}; "
          f"her actuators kept: {len(P.HER_ACTUATORS)}")
    return model


def run(W, dur=30.0, wsup=0.7, omega=2 * np.pi / 1.1, gamma=1.0, gscale=0.1, gains_in=None,
        out="dirusso_run.npz", replay_s=2.0):
    """One trial. If replay_s > 0, the first replay_s seconds impose the
    reference walking kinematics (bsolve qpos) to build up neural state
    before releasing into free simulation (Di Russo's bootstrap)."""
    model = load_model()
    data = mujoco.MjData(model)
    ids, geoms, pelvis, tilt_adr, tx_adr = P.index_maps(model, mujoco)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    for j in P.LOCK_JOINTS:
        jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, j)
        data.qpos[model.jnt_qposadr[jid]] = 0.0
        data.qvel[model.jnt_dofadr[jid]] = 0.0
    mujoco.mj_forward(model, data)

    # ---- REPLAY REFERENCE (bsolve kinematics for the sagittal joints) ----
    ref = np.load(BSOLVE, allow_pickle=True)
    ref_t, ref_q = ref["t"], ref["qpos"]
    ref_acts, ref_names = ref["acts"], list(ref["act_names"])
    # find the qpos column indices for the sagittal joints in the FULL model
    SAGITTAL = ("pelvis_tilt", "pelvis_tx", "pelvis_ty",
                "hip_flexion_r", "knee_angle_r", "ankle_angle_r", "mtp_angle_r",
                "hip_flexion_l", "knee_angle_l", "ankle_angle_l", "mtp_angle_l")
    # the bsolve qpos is from the full (unlocked) model; find cols by joint
    # order: root(7) then joints in model order — the first 7 are the root
    # freejoint (3 trans + 4 quat). Columns 7+ map to model joints.
    # From the earlier check: col 8=hip_flexion_r, 11=knee_angle_r, etc.
    # We use a simpler approach: look up by joint name.
    full_model = mujoco.MjModel.from_xml_path(XML)
    ref_cols = {}
    for jn in SAGITTAL:
        jid = mujoco.mj_name2id(full_model, mujoco.mjtObj.mjOBJ_JOINT, jn)
        if jid >= 0:
            ref_cols[jn] = full_model.jnt_qposadr[jid]
    # also the root free-joint qpos (columns 0-6)
    root_jadr = full_model.jnt_qposadr[mujoco.mj_name2id(
        full_model, mujoco.mjtObj.mjOBJ_JOINT, "root")] \
        if mujoco.mj_name2id(full_model, mujoco.mjtObj.mjOBJ_JOINT, "root") >= 0 \
        else 0

    gains = {k: (v if np.isscalar(v) else v * gscale)
             for k, v in C.default_gains().items()}
    gains["ff"] = 1.0
    if gains_in is not None:
        gains.update(gains_in)
    ctrl = C.DiRussoController(W, dt=model.opt.timestep,
                               omega=omega, gamma=gamma, gains=gains)
    g = model.opt.gravity[2] * -1.0
    weight = -g * model.body_mass.sum()
    dt = model.opt.timestep
    n = int(dur / dt)
    lr = {}
    for (m, s), aid in ids.items():
        lo = model.actuator_lengthrange[aid[0]]
        lr[(m, s)] = (lo[0], lo[1])
    log = {k: np.zeros(n) for k in
           ("t", "hL", "hR", "tL", "tR", "z", "x", "tilt", "cycleL", "cycleR")}
    heel_prev = {"L": 20.1, "R": 20.1}
    strikes = {"L": [], "R": []}

    n_replay = int(replay_s / dt) if replay_s > 0 else 0
    for i in range(n):
        t = i * dt

        # ---- REPLAY PHASE: impose reference kinematics ----
        if i < n_replay:
            # find the reference frame (loop the 2-second reference)
            ref_time = (t % (ref_t[-1] - ref_t[0] + ref_t[1] - ref_t[0]))
            idx = np.searchsorted(ref_t, ref_time + ref_t[0]) - 1
            idx = np.clip(idx, 0, len(ref_t) - 2)
            frac = (ref_time + ref_t[0] - ref_t[idx]) / \
                (ref_t[idx + 1] - ref_t[idx])
            # interpolate reference qpos
            for jn, col in ref_cols.items():
                # find the corresponding joint in the sagittal model
                sjid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, jn)
                if sjid >= 0 and jn not in P.LOCK_JOINTS:
                    target = ref_q[idx, col] + frac * \
                        (ref_q[idx + 1, col] - ref_q[idx, col])
                    data.qpos[model.jnt_qposadr[sjid]] = target
            # root position from reference
            for rc in range(min(7, ref_q.shape[1])):
                if rc < data.qpos.shape[0]:
                    target = ref_q[idx, rc] + frac * \
                        (ref_q[idx + 1, rc] - ref_q[idx, rc])
                    data.qpos[rc] = target
            # compute velocity from the reference (finite difference)
            # CRITICAL FIX 2026-10-26: releasing with zero velocity makes
            # the walker fall from a dead stop; it needs walking momentum
            ref_dt = ref_t[idx + 1] - ref_t[idx]
            for jn, col in ref_cols.items():
                sjid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, jn)
                if sjid >= 0 and jn not in P.LOCK_JOINTS:
                    qadr = model.jnt_qposadr[sjid]
                    dofadr = model.jnt_dofadr[sjid]
                    vel = (ref_q[idx + 1, col] - ref_q[idx, col]) / ref_dt
                    if dofadr < data.qvel.shape[0]:
                        data.qvel[dofadr] = vel
            # root velocity (columns 0-6: 3 translational + angular via
            # quaternion; approximate with translational finite diff)
            root_dof = model.jnt_dofadr[mujoco.mj_name2id(
                model, mujoco.mjtObj.mjOBJ_JOINT, "root")] \
                if mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT,
                                     "root") >= 0 else 0
            for rc in range(min(3, ref_q.shape[1])):
                dof = root_dof + rc
                if dof < data.qvel.shape[0]:
                    vel = (ref_q[idx + 1, rc] - ref_q[idx, rc]) / ref_dt
                    data.qvel[dof] = vel
            mujoco.mj_forward(model, data)

        # ---- SENSORS (same for replay and free phases) ----
        sens = {}
        heel_now = {}
        for s in ("L", "R"):
            fh = mujoco_contacts(data, model, geoms[("heel", s)])
            ft = mujoco_contacts(data, model, geoms[("toe", s)])
            heel_now[s] = fh
            log[f"h{s}"][i], log[f"t{s}"][i] = fh, ft
            hs = heel_prev[s] <= 20.0 < fh
            heel_prev[s] = fh
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
                ria[C.MIDX[m]], rii[C.MIDX[m]], rib[C.MIDX[m]] = ln_, vn_, fn_
            a_, b_, c_ = C.DiRussoController.receptors(rii, ria, rib)
            cut = float(np.clip((fh + ft) / 700.0, 0, 1))
            sens[s] = (a_, b_, c_, cut)
            if hs:
                strikes[s].append(t)
            log[f"cycle{s}"][i] = len(strikes[s])

        tilt = data.qpos[tilt_adr]
        exc = ctrl.step(sens, {s: heel_now[s] > 20.0 and
                         (not strikes[s] or t - strikes[s][-1] > 0.45)
                         for s in ("L", "R")}, tilt)
        for si, s in enumerate(("L", "R")):
            for m in C.MUSCLES:
                for a in ids[(m, s)]:
                    data.ctrl[a] = exc[si, C.MIDX[m]]

        # ---- HARNESS ----
        data.xfrc_applied[:] = 0.0
        if i < n_replay:
            data.xfrc_applied[pelvis, 2] = 1.0 * weight  # full support during replay
        else:
            t_free = t - replay_s
            data.xfrc_applied[pelvis, 2] = (wsup if t_free > 0.6 else 1.0) * weight

        log["t"][i], log["z"][i], log["x"][i], log["tilt"][i] = \
            t, data.xpos[pelvis, 2], data.qpos[tx_adr], tilt

        # ---- STEP ----
        if i < n_replay:
            # during replay: step physics (lets actuator activation build
            # up from ctrl), THEN snap qpos back to the reference (so the
            # kinematics stay on-track while the internal actuator state
            # evolves). Without mj_step during replay, data.act stays at 0
            # and the muscles have no force at the release transition.
            qpos_ref = data.qpos.copy()
            mujoco.mj_step(model, data)
            data.qpos[:] = qpos_ref
            data.qvel[:] = 0.0
            mujoco.mj_forward(model, data)
        else:
            mujoco.mj_step(model, data)
            if not np.isfinite(data.qpos).all():
                print(f"NON-FINITE at t={t:.2f}")
                break

    np.savez(os.path.join(HERE, out), W=W, **log)
    return log, strikes


def mujoco_contacts(data, model, geom_id):
    tot = 0.0
    f = np.zeros(6)
    for k in range(data.ncon):
        c = data.contact[k]
        if c.geom1 == geom_id or c.geom2 == geom_id:
            mujoco.mj_contactForce(model, data, k, f)
            tot += max(f[2], 0.0)   # contact-frame z = normal
    return tot


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dur", type=float, default=30.0)
    ap.add_argument("--wsup", type=float, default=0.7)
    ap.add_argument("--omega", type=float, default=2 * np.pi / 1.1)
    ap.add_argument("--fit-only", action="store_true")
    a = ap.parse_args()
    W, r2 = stage1_fit()
    print("stage-1 imitation fit R^2 (mean per muscle):")
    for si, s in enumerate(("L", "R")):
        print("  " + s + ": " + " ".join(
            f"{m}={r2[si, C.MIDX[m]]:.2f}" for m in C.MUSCLES))
    np.save(os.path.join(HERE, "dirusso_W_stage1.npy"), W)
    if a.fit_only:
        sys.exit(0)
    log, strikes = run(W, dur=a.dur, wsup=a.wsup, omega=a.omega)
    z, x = log["z"], log["x"]
    seg = slice(len(z) // 3, None)
    slope = np.polyfit(log["t"][seg], x[seg], 1)[0]
    print(f"\nclosed loop: dur {a.dur}s wsup {a.wsup}")
    print(f"  pelvis z min/mean: {z.min():.3f}/{z.mean():.3f} m")
    print(f"  forward speed (last 2/3): {slope:+.3f} m/s")
    print(f"  heel strikes L/R: {len(strikes['L'])}/{len(strikes['R'])}")
    print(f"  tilt range: {np.degrees(log['tilt']).min():.1f}.."
          f"{np.degrees(log['tilt']).max():.1f} deg")
