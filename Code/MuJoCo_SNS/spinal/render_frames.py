"""Render offscreen snapshots of the deafferented air-stepping run at
chosen gait-cycle phases (Ben's dissertation figure).

Replays the exact runner air configuration headless, records qpos at every
step, finds RG-E rise-to-rise gait cycles from the saved npz, and renders
one tiled PNG of the model at the requested cycle phases.

Usage: python render_frames.py [out.png] [--phases 0.08,0.33,0.58,0.83]
"""
from __future__ import annotations

import sys

import numpy as np

import runner as R
import mujoco
from params import AFF, DT, E_HI, MOD, SCHEDULE

OUT = "render_walk_phases.png"


def main(argv):
    phases = [0.08, 0.33, 0.58, 0.83]
    out = OUT
    args = list(argv)
    while args:
        a = args.pop(0)
        if a == "--phases":
            phases = [float(x) for x in args.pop(0).split(",")]
        elif not a.startswith("--"):
            out = a

    # ---- build the air-stepping configuration (same as runner's
    # deafferented air path) and replay, recording qpos every step
    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    key_pose = R.capture_pose(model)

    model_solve = R.apply_harness(model, data, kxy=1500.0)
    data_solve = mujoco.MjData(model_solve)
    R.seed_pose(model_solve, data_solve, key_pose)
    mujoco.mj_forward(model_solve, data_solve)

    model = R.apply_harness(model, data, kxy=1500.0, no_ground=True,
                            pin_rot=True)
    data = mujoco.MjData(model)
    R.seed_pose(model, data, key_pose)
    mujoco.mj_forward(model, data)

    acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
            for i in range(model.nu)]
    import build_network as bn
    net = bn.build(acts, dt=DT, interleg=True)
    aid = {n: i for i, n in enumerate(acts)}

    Lrange = model.actuator_lengthrange
    Lmid = Lrange.mean(axis=1)
    Lhalf = np.maximum((Lrange[:, 1] - Lrange[:, 0]) / 2, 1e-3)
    Fmax = np.maximum(model.actuator_gainprm[:, 2], 5.0)
    act_stand = R.solve_standing_activations(model_solve, data_solve, Fmax)

    iport = {p: net.input_index(p) for p in net.inputs}
    walk_drive = 1.5
    WARMUP = R.WARMUP
    dur = SCHEDULE["stand2"][1]
    nsteps = int(dur / DT)
    qposes = np.zeros((nsteps, model.nq))
    u = net.make_inputs()
    for k in range(nsteps):
        t = k * DT
        drive, posture = R.drive_posture(t, walk_drive)
        for i, a in enumerate(acts):
            u[iport[f"Ia_{a}"]] = 0.0        # deafferented
            u[iport[f"II_{a}"]] = 0.0
            u[iport[f"Ib_{a}"]] = 0.0
        u[iport["BAL_PF"]] = 0.0
        u[iport["BAL_DF"]] = 0.0
        u[iport["DRIVE"]] = drive
        u[iport["POSTURE"]] = posture
        stand_frac = 0.05 + 0.95 * (1.0 - min(drive / walk_drive, 1.0))
        for i, a in enumerate(acts):
            u[iport[f"POST_{a}"]] = 6.0 * act_stand[i] * stand_frac
        v = net.step(u)
        for i, a in enumerate(acts):
            data.ctrl[aid[a]] = np.clip(v[net.idx[net.mn_names[a]]] / E_HI,
                                        0.0, 1.0)
        if t < WARMUP:
            data.qvel[:] = 0.0
            mujoco.mj_forward(model, data)
        else:
            mujoco.mj_step(model, data)
        qposes[k] = data.qpos

    # ---- gait cycles from the runner's own npz (RG_E_r rises)
    dnp = np.load("spinal_run.npz", allow_pickle=True)
    tn, neuro = dnp["t"], dnp["neuro"]
    rge = neuro[:, 2]                     # RG_E_r column
    walk = tn > WARMUP + 0.5
    thr = 0.5 * np.max(rge[walk])
    on = (rge > thr) & walk
    rises = np.flatnonzero(np.diff(on.astype(int)) == 1) + 1
    rises = [r for r in rises if tn[r] > WARMUP + 0.5]
    if len(rises) < 3:
        print(f"only {len(rises)} cycles detected - not rendering")
        return
    c0, c1 = int(rises[1]), int(rises[2])   # one clean mid-run cycle
    cyc_t = tn[c0:c1]
    print(f"gait cycle: t={cyc_t[0]:.2f}..{cyc_t[-1]:.2f} s "
          f"({cyc_t[-1] - cyc_t[0]:.2f} s)")

    # ---- render the requested phases
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, len(phases), figsize=(2.9 * len(phases), 4.4))
    try:
        # enlarge the offscreen framebuffer if the model allows it
        try:
            model.vis.global_.offwidth = 700
            model.vis.global_.offheight = 900
        except Exception:
            pass
        renderer = mujoco.Renderer(model, height=880, width=660)
    except Exception as e:
        print(f"Renderer unavailable: {e}")
        return
    cam = mujoco.MjvCamera()
    cam.lookat[:] = [0.0, 0.0, 0.72]
    cam.distance = 2.8
    cam.azimuth = 90
    cam.elevation = -5
    # hide the ground plane (air-stepping) - move it out of the visible
    # geom groups so the print figure is model-only
    gid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM, "ground-plane")
    if gid >= 0:
        model.geom_group[gid] = 5

    def to_white(pix):
        """Black MuJoCo backdrop -> white for print."""
        dark = pix.sum(axis=2) < 40
        pix = pix.copy()
        pix[dark] = (255, 255, 255)
        return pix

    for ax, ph in zip(np.atleast_1d(axs), phases):
        t_target = cyc_t[0] + ph * (cyc_t[-1] - cyc_t[0])
        i = int(np.argmin(np.abs(tn - t_target)))
        data.qpos[:] = qposes[i]
        # recompute POSITION-DEPENDENT fields - mj_kinematics alone moves
        # the bones but NOT the tendon paths the muscle cylinders are
        # drawn from (muscles appeared detached from the skeleton)
        mujoco.mj_kinematics(model, data)
        mujoco.mj_comPos(model, data)
        try:
            mujoco.mj_fwdPosition(model, data)   # includes mj_tendon
        except Exception:
            mujoco.mj_forward(model, data)
        renderer.update_scene(data, camera=cam)
        renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SKYBOX] = 0
        ax.imshow(to_white(renderer.render()))
        ax.set_title(f"gait cycle {100 * ph:.0f}%", fontsize=10)
        ax.axis("off")
    fig.suptitle("Spinal-CPG air-stepping (deafferented, trunk clamped)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    print(f"saved {out}")


if __name__ == "__main__":
    main(sys.argv[1:])
