"""Probe: do frontal sign flips actually change the ID torques?"""
import numpy as np
import mujoco
from scipy.signal import savgol_filter

import bsolve_ik as B
import runner

t_ik, ik_names, ik_vals = B.read_mot(B.IK_MOT)
ik = {n: ik_vals[:, i] for i, n in enumerate(ik_names)}
t_g, g_names, g_vals = B.read_mot(B.GRF_MOT)

model, data = B.build_air_model()
T = len(t_ik)
Fr, Pr, Tr = B.grf_at(t_ik, t_g, g_names, g_vals, pref="")
Fl, Pl, Tl = B.grf_at(t_ik, t_g, g_names, g_vals, pref="1_")
rows = np.array([model.jnt_dofadr[
    mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, j)]
    for j in B.FIT_JOINTS])
cal_r = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "calcn_r")
cal_l = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "calcn_l")

dt = float(np.mean(np.diff(t_ik)))
win = min(11, T - (1 - T % 2))
jacp = np.zeros((3, model.nv))
jacr = np.zeros((3, model.nv))

def run(sgn, label):
    qd = B.frames_qpos(model, data, sgn, ik, ik_names, t_ik)
    v = savgol_filter(qd, win, 3, deriv=1, delta=dt, axis=0)
    a = savgol_filter(qd, win, 3, deriv=2, delta=dt, axis=0)
    med = {j: [] for j in ("hip_adduction_r", "subtalar_angle_r",
                           "knee_angle_r", "hip_flexion_r")}
    n_med = 0
    for k in range(0, T, 4):
        B.set_frame(model, data, sgn, {n: ik[n][k] for n in ik_names})
        data.qpos[:] = qd[k]
        data.qvel[:] = v[k]
        data.qacc[:] = a[k]
        data.act[:] = 0.0
        mujoco.mj_inverse(model, data)
        tau = data.qfrc_inverse.copy()
        for (F, P, Tq, bid) in ((Fr, Pr, Tr, cal_r), (Fl, Pl, Tl, cal_l)):
            mujoco.mj_jac(model, data, jacp, jacr, P[k], bid)
            tau -= jacp.T @ F[k] + jacr.T @ Tq[k]
        for j in med:
            adr = model.jnt_dofadr[
                mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, j)]
            med[j].append(abs(tau[adr]))
        n_med += 1
    print(f"{label}: " + " ".join(
        f"{j}={np.median(m):7.1f}" for j, m in med.items()))

    # geometry probe: stance-foot CoP distance to its ankle axis, frame 60
    k = 60
    B.set_frame(model, data, sgn, {n: ik[n][k] for n in ik_names})
    mujoco.mj_forward(model, data)
    tal = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "talus_r")
    cal = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "calcn_r")
    print(f"   t={t_ik[k]:.2f} talus_r pos {np.round(data.xpos[tal], 3)} "
          f"calcn_r {np.round(data.xpos[cal], 3)} cop_r {np.round(Pr[k], 3)} "
          f"F_r {np.round(Fr[k], 0)}")

base = {n: 1.0 for n in B.ROT_DRIVERS}
run(base, "all+1        ")
flip = dict(base)
flip.update(hip_adduction_r=-1.0, hip_adduction_l=-1.0)
run(flip, "hipadd -1/-1 ")
flip2 = dict(base)
flip2.update(subtalar_angle_r=-1.0, subtalar_angle_l=-1.0)
run(flip2, "subtal -1/-1 ")
