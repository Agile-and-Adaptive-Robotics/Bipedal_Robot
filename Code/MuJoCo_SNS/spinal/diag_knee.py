"""Probe: where does the knee-row ID torque come from?
Decompose qfrc_inverse[knee] into: static gravity, +GRF, +velocity,
+acceleration, and the constraint-force share, per gait phase."""
import numpy as np
import mujoco

import bsolve_ik as B

t_ik, ik_names, ik_vals = B.read_mot(B.IK_MOT)
ik = {n: ik_vals[:, i] for i, n in enumerate(ik_names)}
t_g, g_names, g_vals = B.read_mot(B.GRF_MOT)
model, data = B.build_air_model(id_mode=True)
T = len(t_ik)
Fr, Pr, Tr = B.grf_at(t_ik, t_g, g_names, g_vals, pref="")
Fl, Pl, Tl = B.grf_at(t_ik, t_g, g_names, g_vals, pref="1_")
signs = {n: 1.0 for n in B.ROT_DRIVERS}
signs.update(hip_adduction_l=-1.0, subtalar_angle_l=-1.0,
             hip_rotation_r=-1.0, hip_rotation_l=1.0)

qd = B.frames_qpos(model, data, signs, ik, ik_names, t_ik)
dt = float(np.mean(np.diff(t_ik)))
qvel, qacc = B.smooth_derivs(qd, dt)

knee = model.jnt_dofadr[
    mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, "knee_angle_r")]
cal_r = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "calcn_r")
cal_l = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "calcn_l")
jacp = np.zeros((3, model.nv))
jacr = np.zeros((3, model.nv))

def id_at(k, v, a, with_grf):
    B.set_frame(model, data, signs, {n: ik[n][k] for n in ik_names})
    data.qpos[:] = qd[k]
    data.qvel[:] = v
    data.qacc[:] = a
    data.act[:] = 0.0
    mujoco.mj_inverse(model, data)
    tau = data.qfrc_inverse.copy()
    g_tau = data.qfrc_constraint.copy() if False else np.zeros(model.nv)
    if with_grf:
        for (F, P, Tq, bid) in ((Fr, Pr, Tr, cal_r), (Fl, Pl, Tl, cal_l)):
            mujoco.mj_jac(model, data, jacp, jacr, P[k], bid)
            tau -= jacp.T @ F[k] + jacr.T @ Tq[k]
    return tau

print(f"{'t':>5s} {'Fz_r':>6s} {'static':>8s} {'+dyn':>8s} {'constraint':>10s} "
      f"{'knee qacc':>9s}")
for k in range(0, T, 10):
    tau_static = id_at(k, 0.0, 0.0, True)             # gravity+GRF only
    tau_full = id_at(k, qvel[k], qacc[k], True)       # everything
    dyn = tau_full[knee] - tau_static[knee]
    print(f"{t_ik[k]:5.2f} {Fr[k, 2]:6.0f} {tau_static[knee]:8.1f} "
          f"{tau_full[knee]:8.1f} {dyn:10.1f} {qacc[k, knee]:9.1f}")
