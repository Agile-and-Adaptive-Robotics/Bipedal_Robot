"""Locate the first NaN/instability: which DoF diverges, under what ctrl."""
import sys
import numpy as np
import mujoco

import build_network as bn
import runner as R
from params import DT, E_HI


def main():
    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    model = R.apply_harness(model, data)
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    mujoco.mj_forward(model, data)
    acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(model.nu)]
    net = bn.build(acts, dt=DT)
    aid = {n: i for i, n in enumerate(acts)}
    u = net.make_inputs()
    u[net.input_index("POSTURE")] = 1.0

    for k in range(500):
        t = k * DT
        drive, posture = R.drive_posture(t, 4.0)
        for i, a in enumerate(acts):
            u[net.input_index("POST_" + a)] = 0.0
        u[net.input_index("DRIVE")] = drive
        u[net.input_index("POSTURE")] = posture
        v = net.step(u)
        for i, a in enumerate(acts):
            data.ctrl[aid[a]] = np.clip(v[net.idx[net.mn_names[a]]] / E_HI, 0.0, 1.0)
        if t < R.WARMUP:
            data.qvel[:] = 0.0
            mujoco.mj_forward(model, data)
        else:
            mujoco.mj_step(model, data)
        if not np.all(np.isfinite(data.qvel)):
            print(f"NaN at t={t:.3f}")
            bad = np.flatnonzero(~np.isfinite(data.qvel))
            # map dofs back to joints
            names = []
            for b in bad:
                for j in range(model.njnt):
                    da, nd = model.jnt_dofadr[j], 1
                    if da <= b < da + nd:
                        names.append(mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, j))
            print("bad dofs -> joints:", names[:8])
            top = np.argsort(data.qvel[np.isfinite(data.qvel)])[::-1][:5]
            print("largest qvel before:", top)
            break
        # report any dof velocity > 50 rad/s as precursor
        if k % 20 == 0 and t > R.WARMUP:
            w = np.abs(data.qvel)
            j = int(np.argmax(w))
            for jj in range(model.njnt):
                if model.jnt_dofadr[jj] <= j < model.jnt_dofadr[jj] + 1:
                    nm = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, jj)
                    if w[j] > 20:
                        print(f"t={t:.2f} fast: {nm} |qvel|={w[j]:.1f}")
                    break


if __name__ == "__main__":
    main()
