"""FIGURES task (2026-09-25): cheap offscreen GIF of the w2lvar stage-5
winner's ground walk.

Replays the saved qfull (qpos) from figs_npz_w2lvar.npz into the raw cvt3
model (visual geoms only - harness/patch elements are irrelevant to the
bones) and renders offscreen with mujoco.Renderer. Walk window 2-13 s
(runner --eval schedule), 10 fps.

Usage: python make_gif.py [probe]
  probe  render ONE frame, report timing, no gif write (capability check)
"""
import os
import sys
import time

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, SPINAL)
os.chdir(SPINAL)

import mujoco  # noqa: E402

MODEL_XML = (r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim"
             r"\Gait2392_Robotbody\mjc\gait2392_simbody"
             r"\gait2392_simbody_cvt3.xml")
NPZ = os.path.join("reports_20260925", "figs", "tmp",
                   "figs_npz_w2lvar.npz")
OUT = os.path.join("reports_20260925", "figs", "w2lvar_s5_walk.gif")
WALK_T0, WALK_T1 = 2.0, 13.0
FPS = 10
H, W = 360, 480

probe = "probe" in sys.argv
t0 = time.time()
model = mujoco.MjModel.from_xml_path(MODEL_XML)
data = mujoco.MjData(model)
print(f"[gif] model loaded in {time.time()-t0:.1f} s "
      f"(nq={model.nq})", flush=True)

z = np.load(NPZ, allow_pickle=True)
t, qfull = z["t"], z["qfull"]
print(f"[gif] npz t {t[0]:.2f}..{t[-1]:.2f} qfull {qfull.shape}")

renderer = mujoco.Renderer(model, height=H, width=W)
cam = mujoco.MjvCamera()
cam.azimuth = 90.0      # side view (world x = travel, z = up)
cam.elevation = -8.0
cam.distance = 2.4
com = z["com"]
k = int(np.searchsorted(t, WALK_T0 + 0.5))
data.qpos[:] = qfull[k]
mujoco.mj_forward(model, data)
cam.lookat[:] = [com[k, 0], 0.0, 0.85]
renderer.update_scene(data, camera=cam)
t0 = time.time()
img = renderer.render()
dt1 = time.time() - t0
print(f"[gif] first frame render {dt1:.2f} s, shape {img.shape}")
if probe:
    sys.exit(0)

from PIL import Image  # noqa: E402

step = 1.0 / FPS
ts = np.arange(WALK_T0, WALK_T1, step)
frames = []
t0 = time.time()
t_start_render = time.time()
for i, tt in enumerate(ts):
    k = int(np.searchsorted(t, tt))
    data.qpos[:] = qfull[k]
    mujoco.mj_forward(model, data)
    cam.lookat[:] = [com[k, 0], 0.0, 0.85]
    renderer.update_scene(data, camera=cam)
    frames.append(Image.fromarray(renderer.render()).resize(
        (W, H), Image.BILINEAR))
    if i == 9:
        rate = (time.time() - t_start_render) / 10.0
        est = rate * len(ts)
        print(f"[gif] ~{rate:.2f} s/frame, est total {est:.0f} s",
              flush=True)
        if est > 300:
            print("[gif] ABORT: too slow, skipping gif")
            sys.exit(2)
render_s = time.time() - t_start_render
frames[0].save(OUT, save_all=True, append_images=frames[1:],
               duration=int(1000 / FPS), loop=0)
print(f"[gif] wrote {OUT} ({len(frames)} frames, render {render_s:.0f} s, "
      f"{os.path.getsize(OUT)/1e6:.1f} MB)")
