"""Render the deafferented air-stepping run to an animated GIF.

Replays the exact runner air configuration headless, renders frames with
mujoco.Renderer during the walk window, and encodes with Pillow (no
ffmpeg needed). Output: walk_air.gif (30 fps, half-size frames).

Usage: python render_video.py [out.gif] [--fps 30] [--t0 4.6] [--t1 15.5]
"""
from __future__ import annotations

import sys

import numpy as np

import runner as R
import mujoco
from params import AFF, DT, E_HI, SCHEDULE


def main(argv):
    out = "walk_air.gif"
    fps = 30.0
    t0, t1 = 4.6, 15.6          # walk window (drive on -> ramp down)
    args = list(argv)
    while args:
        a = args.pop(0)
        if a == "--fps":
            fps = float(args.pop(0))
        elif a == "--t0":
            t0 = float(args.pop(0))
        elif a == "--t1":
            t1 = float(args.pop(0))
        elif not a.startswith("--"):
            out = a

    # ---- build + replay (identical to render_frames.py)
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
    nsteps = int(SCHEDULE["stand2"][1] / DT)

    gid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM, "ground-plane")
    if gid >= 0:
        model.geom_group[gid] = 5          # hide ground (air-stepping)
    try:
        model.vis.global_.offwidth = 700
        model.vis.global_.offheight = 900
    except Exception:
        pass
    renderer = mujoco.Renderer(model, height=640, width=480)
    cam = mujoco.MjvCamera()
    cam.lookat[:] = [0.0, 0.0, 0.72]
    cam.distance = 2.8
    cam.azimuth = 90
    cam.elevation = -5
    renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SKYBOX] = 0

    from PIL import Image
    frame_every = max(1, int(round((1.0 / fps) / DT)))
    frames = []
    u = net.make_inputs()
    print(f"replaying {nsteps} steps, capturing every {frame_every} "
          f"(~{fps:.0f} fps)...", flush=True)
    for k in range(nsteps):
        t = k * DT
        drive, posture = R.drive_posture(t, walk_drive)
        for i, a in enumerate(acts):
            u[iport[f"Ia_{a}"]] = 0.0
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
        if t0 <= t <= t1 and k % frame_every == 0:
            mujoco.mj_kinematics(model, data)
            mujoco.mj_comPos(model, data)
            renderer.update_scene(data, camera=cam)
            pix = renderer.render()
            dark = pix.sum(axis=2) < 40
            pix = pix.copy()
            pix[dark] = (255, 255, 255)    # print-white backdrop
            frames.append(Image.fromarray(pix).resize((240, 320)))
        if k % 1000 == 0:
            print(f"  t={t:5.1f}s frames={len(frames)}", flush=True)

    try:
        renderer.close()   # not present in mujoco 2.3.x; GC handles it
    except Exception:
        pass
    if not frames:
        print("no frames captured - check t0/t1")
        return
    frames[0].save(out, save_all=True, append_images=frames[1:],
                   duration=int(1000 / fps), loop=0)
    print(f"saved {out} ({len(frames)} frames, {fps:.0f} fps, "
          f"{len(frames) / fps:.1f} s)")


if __name__ == "__main__":
    main(sys.argv[1:])
