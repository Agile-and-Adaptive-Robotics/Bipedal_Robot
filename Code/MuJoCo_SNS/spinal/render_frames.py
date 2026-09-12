"""Render offscreen snapshots of the LAST spinal_run.npz at chosen
gait-cycle phases. Reads the recorded full joint state (qfull), so the
rendered poses are EXACTLY the run's - muscles follow bones, config
(air vs ground) detected from the npz.

Usage: python render_frames.py [out.png] [--phases 0.08,0.33,0.58,0.83]
"""
from __future__ import annotations

import sys

import numpy as np

import runner as R
import mujoco
from params import DT

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

    d = np.load("spinal_run.npz", allow_pickle=True)
    t, qfull, neuro = d["t"], d["qfull"], d["neuro"]
    cfg = str(d["cfg"]) if "cfg" in d else "air"

    # ---- model with the same repairs the run used
    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    key_pose = R.capture_pose(model)
    data0 = mujoco.MjData(model)
    model = R.apply_harness(model, data0, kxy=1500.0,
                            no_ground=(cfg == "air"), pin_rot=(cfg == "air"))
    data = mujoco.MjData(model)

    # ---- gait cycles from RG_E_r rises during the walk window
    drive = neuro[:, 0]
    walk = drive > 0.8 * max(np.max(drive), 1e-6)
    rge = neuro[:, 2]
    thr = 0.5 * np.max(rge[walk])
    on = (rge > thr) & walk
    rises = np.flatnonzero(np.diff(on.astype(int)) == 1) + 1
    if len(rises) < 3:
        print(f"only {len(rises)} cycles detected - not rendering")
        return
    c0, c1 = int(rises[1]), int(rises[2])
    cyc_t = t[c0:c1]
    print(f"gait cycle: t={cyc_t[0]:.2f}..{cyc_t[-1]:.2f} s "
          f"({cyc_t[-1] - cyc_t[0]:.2f} s) [{cfg}]")

    # ---- render
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, len(phases), figsize=(2.9 * len(phases), 4.4))
    try:
        model.vis.global_.offwidth = 700
        model.vis.global_.offheight = 900
    except Exception:
        pass
    renderer = mujoco.Renderer(model, height=880, width=660)
    cam = mujoco.MjvCamera()
    cam.lookat[:] = [0.0, 0.0, 0.72]
    cam.distance = 2.8
    cam.azimuth = 90
    cam.elevation = -5
    if cfg == "air":
        gid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM,
                                "ground-plane")
        if gid >= 0:
            model.geom_group[gid] = 5          # hide ground for air steps

    def to_white(pix):
        dark = pix.sum(axis=2) < 40
        pix = pix.copy()
        pix[dark] = (255, 255, 255)
        return pix

    for ax, ph in zip(np.atleast_1d(axs), phases):
        t_target = cyc_t[0] + ph * (cyc_t[-1] - cyc_t[0])
        i = int(np.argmin(np.abs(t - t_target)))
        R.seed_pose(model, data, key_pose)
        data.qpos[:] = qfull[i]
        mujoco.mj_kinematics(model, data)
        mujoco.mj_comPos(model, data)
        try:
            mujoco.mj_fwdPosition(model, data)  # tendon paths for muscles
        except Exception:
            mujoco.mj_forward(model, data)
        renderer.update_scene(data, camera=cam)
        renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SKYBOX] = 0
        ax.imshow(to_white(renderer.render()))
        ax.set_title(f"gait cycle {100 * ph:.0f}%", fontsize=10)
        ax.axis("off")
    title = ("Spinal-CPG ground walking (afferent feedback, semi-supported)"
             if cfg == "ground" else
             "Spinal-CPG air-stepping (deafferented, trunk clamped)")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    print(f"saved {out}")


if __name__ == "__main__":
    main(sys.argv[1:])
