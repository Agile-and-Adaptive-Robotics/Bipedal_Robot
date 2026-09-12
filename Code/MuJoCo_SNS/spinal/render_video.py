"""Render the LAST spinal_run.npz to an animated GIF (30 fps default).

Reads the recorded full joint state (qfull) - the rendered motion is
exactly the simulated run, config (air vs ground) auto-detected.
Encodes with Pillow (no ffmpeg needed).

Usage: python render_video.py [out.gif] [--fps 30]
"""
from __future__ import annotations

import sys

import numpy as np

import runner as R
import mujoco
from params import DT


def main(argv):
    out = "walk_air.gif"
    fps = 30.0
    args = list(argv)
    while args:
        a = args.pop(0)
        if a == "--fps":
            fps = float(args.pop(0))
        elif not a.startswith("--"):
            out = a

    d = np.load("spinal_run.npz", allow_pickle=True)
    t, qfull, neuro = d["t"], d["qfull"], d["neuro"]
    cfg = str(d["cfg"]) if "cfg" in d else "air"
    drive = neuro[:, 0]
    walk = drive > 0.8 * max(np.max(drive), 1e-6)
    idx = np.flatnonzero(walk)
    if len(idx) < 100:
        print("walk window too short")
        return
    i0, i1 = int(idx[0]), int(idx[-1])
    print(f"walk window: t={t[i0]:.2f}..{t[i1]:.2f} s [{cfg}]")

    model = mujoco.MjModel.from_xml_path(str(R.MODEL))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)
    key_pose = R.capture_pose(model)
    data0 = mujoco.MjData(model)
    model = R.apply_harness(model, data0, kxy=1500.0,
                            no_ground=(cfg == "air"), pin_rot=(cfg == "air"))
    data = mujoco.MjData(model)

    if cfg == "air":
        gid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM,
                                "ground-plane")
        if gid >= 0:
            model.geom_group[gid] = 5
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
    for k in range(i0, i1 + 1):
        if k % frame_every:
            continue
        R.seed_pose(model, data, key_pose)
        data.qpos[:] = qfull[k]
        mujoco.mj_kinematics(model, data)
        mujoco.mj_comPos(model, data)
        try:
            mujoco.mj_fwdPosition(model, data)
        except Exception:
            mujoco.mj_forward(model, data)
        renderer.update_scene(data, camera=cam)
        pix = renderer.render()
        dark = pix.sum(axis=2) < 40
        pix = pix.copy()
        pix[dark] = (255, 255, 255)
        frames.append(Image.fromarray(pix).resize((240, 320)))
    try:
        renderer.close()
    except Exception:
        pass
    if not frames:
        print("no frames captured")
        return
    frames[0].save(out, save_all=True, append_images=frames[1:],
                   duration=int(1000 / fps), loop=0)
    print(f"saved {out} ({len(frames)} frames, {fps:.0f} fps, "
          f"{len(frames) / fps:.1f} s)")


if __name__ == "__main__":
    main(sys.argv[1:])
