"""Fast trial runner for the Di Russo gait2392 MuJoCo adaptation.

Speed layers (the 2026-10-26/27 set Ben asked for):
  1. NUMBA neural step: controller_numba.CompiledDiRusso, parity-verified
     against controller.py by ../reports_20261027_codex/test_dirusso_equations.py
     (max output difference 3.34e-16 over 2000 steps, ~450x isolated).
  2. Vectorized sensor extraction: one actuator_length/velocity/force read
     per step and precomputed per-group index arrays, replacing the
     per-muscle Python loops in dirusso_v3.run_trial.
  3. One-time model loads: dirusso_v3.run_trial reloaded the FULL XML via
     MjModel.from_xml_path once per joint name inside the initial-condition
     block (17+ parses per trial); loaded once here.

NOT included on purpose: the sparse neural/physics backend
(w2l_cpg/fast_eval.SparseForward) -- the 2026-10-27 recovery report shows it
fails full-trajectory parity on the W2L stack (45.56 deg max joint
deviation). Only parity-safe layers are used here.

Parity contract: verify_fast_parity.py runs _reference_trial (a line-faithful
port of dirusso_v3.run_trial) and run_trial_fast with identical inputs and
asserts trajectory agreement; it also cross-checks the port against
dirusso_v3.run_trial metrics.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import controller as C          # noqa: E402
import plant as P               # noqa: E402
from controller_numba import CompiledDiRusso  # noqa: E402

XML = r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml"
BSOLVE = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\bsolve_out.npz"

IC_JOINTS = ("pelvis_tilt", "pelvis_tx", "pelvis_ty",
             "hip_flexion_r", "knee_angle_r", "ankle_angle_r", "mtp_angle_r",
             "hip_flexion_l", "knee_angle_l", "ankle_angle_l", "mtp_angle_l")
VEL_JOINTS = ("hip_flexion_r", "knee_angle_r", "ankle_angle_r",
              "hip_flexion_l", "knee_angle_l", "ankle_angle_l")

# ---------------------------------------------------------------------------
# Environment (built once per process; workers each build their own)
# ---------------------------------------------------------------------------

class FastEnv:
    """Prebuilt model + every index/constant array the trial loop needs."""

    def __init__(self):
        import mujoco
        self.mujoco = mujoco
        txt, _ = P.build_sagittal_xml(XML)
        tmp = os.path.join(os.path.dirname(XML), "_dirusso_fast_tmp.xml")
        with open(tmp, "w", encoding="utf-8") as fh:
            fh.write(txt)
        try:
            self.model = mujoco.MjModel.from_xml_path(tmp)
        finally:
            os.remove(tmp)
        self.model.opt.timestep = 0.002
        self.data = mujoco.MjData(self.model)
        self.full_model = mujoco.MjModel.from_xml_path(XML)  # for qpos addresses
        self.ref = np.load(BSOLVE, allow_pickle=True)

        ids, geoms, pelvis, tilt_adr, tx_adr = P.index_maps(self.model, mujoco)
        self.ids, self.geoms = ids, geoms
        self.pelvis, self.tilt_adr, self.tx_adr = pelvis, tilt_adr, tx_adr

        # per (side, muscle) index arrays in the SAME order as ids[(m, side)]
        self.grp_idx = np.full((2, C.NM), -1, dtype=np.int64)  # 1-elem groups only
        self.grp_slices = []        # (si, mi, idx_array)
        for si, s in enumerate(("L", "R")):
            for m in C.MUSCLES:
                idx = np.asarray(ids[(m, s)], dtype=np.int64)
                self.grp_slices.append((si, C.MIDX[m], idx))
        # ctrl write order: (si, m) then actuator ids -> flat index list
        flat = []
        flat_mn = []   # mn.ravel() index for each ctrl slot (group-shared value)
        for si, s in enumerate(("L", "R")):
            for m in C.MUSCLES:
                for a in ids[(m, s)]:
                    flat.append(a)
                    flat_mn.append(si * C.NM + C.MIDX[m])
        self.ctrl_flat = np.asarray(flat, dtype=np.int64)
        self.ctrl_mn_idx = np.asarray(flat_mn, dtype=np.int64)

        # length ranges + force normalization per group (constant)
        self.lo = np.zeros((2, C.NM))
        self.span = np.zeros((2, C.NM))
        self.fmax = np.zeros((2, C.NM))
        for si, s in enumerate(("L", "R")):
            for m in C.MUSCLES:
                aid = ids[(m, s)]
                lr = self.model.actuator_lengthrange[aid[0]]
                self.lo[si, C.MIDX[m]] = lr[0]
                self.span[si, C.MIDX[m]] = max(lr[1] - lr[0], 1e-6)
                self.fmax[si, C.MIDX[m]] = max(
                    sum(self.model.actuator_gainprm[a, 2] for a in aid), 1e-6)

        self.heel_gid = np.array([geoms[("heel", "L")], geoms[("heel", "R")]])
        self.toe_gid = np.array([geoms[("toe", "L")], geoms[("toe", "R")]])
        self.weight = 9.81 * self.model.body_mass.sum()

        # IC address maps (computed once; v3 recomputed these per call with a
        # fresh from_xml_path parse per joint name)
        self.ic_qadr = []
        for jn in IC_JOINTS:
            jid = self.mujoco.mj_name2id(self.model, self.mujoco.mjtObj.mjOBJ_JOINT, jn)
            fjid = self.mujoco.mj_name2id(self.full_model, self.mujoco.mjtObj.mjOBJ_JOINT, jn)
            self.ic_qadr.append((jid, self.full_model.jnt_qposadr[fjid]
                                 if fjid >= 0 else -1))
        self.ic_dof = []
        for jn in VEL_JOINTS:
            jid = self.mujoco.mj_name2id(self.model, self.mujoco.mjtObj.mjOBJ_JOINT, jn)
            fjid = self.mujoco.mj_name2id(self.full_model, self.mujoco.mjtObj.mjOBJ_JOINT, jn)
            self.ic_dof.append((jid, self.full_model.jnt_qposadr[fjid]
                                if fjid >= 0 else -1))
        self.lock_jids = [self.mujoco.mj_name2id(
            self.model, self.mujoco.mjtObj.mjOBJ_JOINT, j) for j in P.LOCK_JOINTS]

    def set_reference_ic(self):
        """Same values as dirusso_v3.run_trial's IC block (first ref frame)."""
        mjc, model, data = self.mujoco, self.model, self.data
        mjc.mj_resetDataKeyframe(model, data, 0)
        q0 = self.ref["qpos"][0]
        for jid, full_qadr in self.ic_qadr:
            if jid >= 0 and full_qadr >= 0:
                data.qpos[model.jnt_qposadr[jid]] = q0[full_qadr]
        for rc in range(min(7, q0.shape[0])):
            if rc < data.qpos.shape[0]:
                data.qpos[rc] = q0[rc]
        for jid in self.lock_jids:
            data.qpos[model.jnt_qposadr[jid]] = 0.0
            data.qvel[model.jnt_dofadr[jid]] = 0.0
        q1 = self.ref["qpos"][1]
        dt_ref = self.ref["t"][1] - self.ref["t"][0]
        for jid, full_qadr in self.ic_dof:
            if jid >= 0 and full_qadr >= 0:
                vel = (q1[full_qadr] - q0[full_qadr]) / dt_ref
                dof = model.jnt_dofadr[jid]
                if dof < data.qvel.shape[0]:
                    data.qvel[dof] = vel
        mjc.mj_forward(model, data)
        return float(data.subtree_com[0][2])


_ENV = None


def get_env() -> FastEnv:
    global _ENV
    if _ENV is None:
        _ENV = FastEnv()
    return _ENV


# ---------------------------------------------------------------------------
# Fast trial
# ---------------------------------------------------------------------------

def run_trial_fast(params, wsup=0.7, dur=10.0, W_ff=None, gains=None,
                    balance_gains=None, record=False, disable_strikes=False):
    """Same semantics as dirusso_v3.run_trial; NUMBA + vectorized sensors."""
    env = get_env()
    mjc, model, data = env.mujoco, env.model, env.data
    hcom_initial = env.set_reference_ic()

    if W_ff is None:
        W_ff = np.load(os.path.join(HERE, "dirusso_W_stage1.npy"))
    ctrl = CompiledDiRusso(W_ff,
                           omega=params.get("omega", 5.11),
                           gamma=params.get("gamma", 1.0),
                           gains=gains, balance_gains=balance_gains)

    n = int(dur / 0.002)
    dt = 0.002
    heelF = np.zeros((n, 2)); com_x = np.zeros(n); com_z = np.zeros(n)
    tilt_log = np.zeros(n); q_log = np.zeros((n, 6))
    heel_prev = np.array([20.1, 20.1])
    strikes_T = [[], []]
    fall_t = None
    f6 = np.zeros(6)
    sens4 = np.zeros((2, 4, 9))
    strike_arr = np.zeros(2, dtype=np.bool_)
    last_mn = np.zeros((2, C.NM))

    grp = env.grp_slices
    for i in range(n):
        t = i * dt
        al = data.actuator_length
        av = data.actuator_velocity
        af = data.actuator_force
        for si in range(2):
            heel = toe = 0.0
            for k in range(data.ncon):
                c = data.contact[k]
                mjc.mj_contactForce(model, data, k, f6)
                fz = max(f6[2], 0.0)
                if c.geom1 == env.heel_gid[si] or c.geom2 == env.heel_gid[si]:
                    heel += fz
                elif c.geom1 == env.toe_gid[si] or c.geom2 == env.toe_gid[si]:
                    toe += fz
            heelF[i, si] = heel
            hs = heel_prev[si] <= 20.0 < heel
            heel_prev[si] = heel
            # muscle sensors (same math and order as the reference loop)
            ln_ = np.empty(C.NM); vn_ = np.empty(C.NM); fn_ = np.empty(C.NM)
            for gsi, mi, idx in grp:
                if gsi != si:
                    continue
                L = float(np.mean(al[idx]))
                V = float(np.mean(av[idx]))
                F = float(np.mean(np.abs(af[idx])))
                ln_[mi] = min(max((L - env.lo[si, mi]) / env.span[si, mi], 0.0), 2.0)
                vn_[mi] = min(max(V / env.span[si, mi], -5.0), 5.0)
                fn_[mi] = min(max(F / env.fmax[si, mi], 0.0), 2.0)
            # receptors(velocity, length, force) -- the same argument order
            # as dirusso_v3's receptors(rii, ria, rib)
            a_, b_, c_ = C.DiRussoController.receptors(vn_, ln_, fn_)
            cut = float(min(max((heel + toe) / 700.0, 0.0), 1.0))
            sens4[si, 0] = a_; sens4[si, 1] = b_; sens4[si, 2] = c_
            sens4[si, 3, 0] = cut
            if hs:
                strikes_T[si].append(t)

        com = data.subtree_com[0]
        com_x[i] = com[0]; com_z[i] = com[2]
        tilt_log[i] = data.qpos[env.tilt_adr]
        tilt_deg = abs(np.degrees(tilt_log[i]))
        q_log[i] = (data.qpos[env.tilt_adr], data.qpos[env.tx_adr], com[2],
                    data.qpos[env.tilt_adr] * 0 + 0, tilt_log[i], com[0])
        if fall_t is None and com[2] < 0.8 * hcom_initial:
            fall_t = t

        # heel-strike reset with the same 0.45 s refractory
        for si in range(2):
            strike_arr[si] = (heelF[i, si] > 20.0 and
                              (not strikes_T[si] or t - strikes_T[si][-1] > 0.45))
        if disable_strikes:
            strike_arr[:] = False
        last_mn = ctrl.step_array(sens4, strike_arr, tilt_log[i]).copy()
        data.ctrl[env.ctrl_flat] = last_mn.reshape(-1)[env.ctrl_mn_idx]
        data.xfrc_applied[:] = 0.0
        data.xfrc_applied[env.pelvis, 2] = wsup * env.weight
        mjc.mj_step(model, data)
        if not np.isfinite(data.qpos).all():
            fall_t = fall_t if fall_t is not None else t
            break

    metrics = _metrics(n, dt, heelF, com_x, com_z, tilt_deg, fall_t,
                       strikes_T, last_mn)
    if record:
        metrics["trace"] = dict(heelF=heelF[:i + 1], com_x=com_x[:i + 1],
                                com_z=com_z[:i + 1], tilt=tilt_log[:i + 1],
                                steps=i + 1)
    return metrics


def _metrics(n, dt, heelF, com_x, com_z, tilt_deg, fall_t, strikes_T, exc):
    """Identical to the dirusso_v3 metrics block."""
    t_arr = np.arange(n) * dt
    free = t_arr >= 1.0
    if fall_t is not None and fall_t < 2.0:
        return dict(fell=True, fall_t=fall_t, obj=-100.0, survive=False)
    com_free = com_x[free]
    t_free = t_arr[free]
    speed = (np.polyfit(t_free, com_free, 1)[0]
             if len(com_free) > 10 else 0.0)

    def cyc(f):
        on = f > 20
        st = t_free[1:][on[1:] & ~on[:-1]]
        return sum(1 for j, x in enumerate(st)
                   if j == 0 or x - st[j - 1] > 0.45)
    cycL = cyc(heelF[free, 0]); cycR = cyc(heelF[free, 1])
    effort = float(np.mean(np.sum(exc ** 2, axis=1))) if exc is not None else 100
    ph = None
    if cycL >= 3 and cycR >= 3:
        rT = np.array(strikes_T[1]); lT = np.array(strikes_T[0])
        phs = []
        for j in range(len(rT) - 1):
            ins = lT[(lT >= rT[j]) & (lT < rT[j + 1])]
            if len(ins) == 1:
                phs.append((ins[0] - rT[j]) / (rT[j + 1] - rT[j]))
        if phs:
            ph = float(np.median(phs))
    fell = fall_t is not None
    survive = (not fell and min(cycL, cycR) >= 3 and ph is not None
               and 0.35 <= ph <= 0.65)
    Jgait = 0
    if fell:
        Jgait = 100
    elif speed < 1.10 or speed > 1.25:
        Jgait = 10 * min(abs(speed - 1.175) / 0.075, 1.0)
    Jeffort = effort * 0.01
    Jhead = max(tilt_deg - 15.0, 0.0) * 0.1 if not fell else 50
    Jtotal = Jeffort + Jgait + Jhead
    return dict(fell=fell, fall_t=fall_t, obj=-Jtotal, survive=survive,
                speed=speed, cycL=cycL, cycR=cycR, phase=ph,
                effort=effort, com_z_min=float(com_z.min()),
                tilt_max=float(tilt_deg))


# ---------------------------------------------------------------------------
# Reference port (line-faithful copy of dirusso_v3.run_trial, for parity)
# ---------------------------------------------------------------------------

def _reference_trial(params, wsup=0.7, dur=10.0, W_ff=None, gains=None,
                     record=False, disable_strikes=False):
    import dirusso_v3 as V3
    mujoco = env = None
    import mujoco as _mj
    mujoco = _mj
    model, data = V3.get_model()
    ids, geoms, pelvis, tilt_adr, tx_adr = P.index_maps(model, mujoco)
    ref = np.load(BSOLVE, allow_pickle=True)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    q0 = ref["qpos"][0]
    for jn in IC_JOINTS:
        jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, jn)
        full_jid = mujoco.mj_name2id(mujoco.MjModel.from_xml_path(XML),
                                     mujoco.mjtObj.mjOBJ_JOINT, jn)
        if jid >= 0 and full_jid >= 0:
            full_qadr = mujoco.MjModel.from_xml_path(XML).jnt_qposadr[full_jid]
            data.qpos[model.jnt_qposadr[jid]] = q0[full_qadr]
    for rc in range(min(7, q0.shape[0])):
        if rc < data.qpos.shape[0]:
            data.qpos[rc] = q0[rc]
    for j in P.LOCK_JOINTS:
        jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, j)
        data.qpos[model.jnt_qposadr[jid]] = 0.0
        data.qvel[model.jnt_dofadr[jid]] = 0.0
    q1 = ref["qpos"][1]
    dt_ref = ref["t"][1] - ref["t"][0]
    for jn in VEL_JOINTS:
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
    tilt_log = np.zeros(n)
    heel_prev = {"L": 20.1, "R": 20.1}
    strikes = {"L": [], "R": []}
    fall_t = None
    f6 = np.zeros(6)
    exc = None

    for i in range(n):
        t = i * dt
        sens = {}
        heel_now = {}
        for si, s in enumerate(("L", "R")):
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
        if fall_t is None and com[2] < 0.8 * hcom_initial:
            fall_t = t

        strike_sig = {s: heel_now[s] > 20.0 and
                      (not strikes[s] or t - strikes[s][-1] > 0.45)
                      for s in ("L", "R")}
        if disable_strikes:
            strike_sig = {"L": False, "R": False}
        exc = ctrl.step(sens, strike_sig, tilt_log[i])
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

    metrics = _metrics(n, dt, heelF, com_x, com_z, tilt_deg, fall_t,
                       (strikes["L"], strikes["R"]), exc)
    if record:
        metrics["trace"] = dict(heelF=heelF[:i + 1], com_x=com_x[:i + 1],
                                com_z=com_z[:i + 1], tilt=tilt_log[:i + 1],
                                steps=i + 1)
    return metrics


# ---------------------------------------------------------------------------
# Parallel campaign (multiprocessing; no optuna search is launched here)
# ---------------------------------------------------------------------------

def _worker_init():
    global _ENV
    _ENV = FastEnv()


def _worker_eval(job):
    params, wsup, dur, W_ff, gains = job
    return run_trial_fast(params, wsup=wsup, dur=dur, W_ff=W_ff, gains=gains)


def evaluate_batch(jobs, n_workers=None):
    """Evaluate trial jobs in parallel fresh processes.

    Each worker builds its own model/env once (module-level FastEnv), so
    there is no cross-process MuJoCo sharing.
    """
    import multiprocessing as mp
    n_workers = n_workers or min(10, os.cpu_count() or 4)
    ctx = mp.get_context("spawn")
    with ctx.Pool(n_workers, initializer=_worker_init) as pool:
        return pool.map(_worker_eval, jobs)


if __name__ == "__main__":
    print("import-only module; run verify_fast_parity.py")
