"""Render a walking animation (GIF) by replaying a runner npz through
the MuJoCo offscreen renderer (the intended use of log_qfull - "full
state for renderers").

Usage: python _render_gif.py [npz] [out.gif] [t0] [t1]
Default: spinal_run.npz -> curr_ground_walk.gif, walk window 2-14 s,
side view, 0.5x speed, 480x360.
"""
import io
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import mujoco
from PIL import Image

import runner as R

HERE = Path(__file__).parent
DISS = (HERE.parents[2] / "Documentation" / "Reports and Papers"
        / "Dissertation" / "CPG_airstepping_figs")

npz = sys.argv[1] if len(sys.argv) > 1 else "spinal_run.npz"
out = Path(sys.argv[2] if len(sys.argv) > 2 else "curr_ground_walk.gif")
t0 = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
t1 = float(sys.argv[4]) if len(sys.argv) > 4 else 14.0

z = np.load(npz, allow_pickle=True)
t, qfull = z["t"], z["qfull"]
z.close()
model = mujoco.MjModel.from_xml_path(str(R.MODEL))
data = mujoco.MjData(model)

STEP = 20          # sim steps per output frame (0.04 s at DT 0.002)
FPS = 25           # playback fps -> 1.0x realtime at STEP/FPS = 0.04 s
idx = np.flatnonzero((t >= t0) & (t <= t1))[::STEP]
print(f"frames: {len(idx)} ({t[idx[0]]:.1f}-{t[idx[-1]]:.1f} s, "
      f"0.04 s/frame @ {FPS} fps = 1.0x)")

renderer = mujoco.Renderer(model, height=360, width=640)
cam = mujoco.MjvCamera()
cam.lookat[:] = (0.0, 0.0, 0.9)
cam.distance = 3.0
cam.azimuth = 90.0     # side view
cam.elevation = -8.0

frames = []
for i in idx:
    data.qpos[:] = qfull[i, :model.nq]
    data.qvel[:] = 0.0
    mujoco.mj_forward(model, data)
    renderer.update_scene(data, camera=cam)
    frames.append(Image.fromarray(renderer.render()))
try:
    renderer.close()   # mujoco >= 3
except AttributeError:
    pass               # mujoco 2.3.7: no close(), GC handles it

out_local = HERE / out.name
frames[0].save(out_local, save_all=True, append_images=frames[1:],
               duration=int(1000 / FPS), loop=0)
frames[0].save(DISS / out.name, save_all=True, append_images=frames[1:],
               duration=int(1000 / FPS), loop=0)
alt = (f"Animation - {out.name}: side view (X anterior to the left, "
       f"Y up) of the ground-walk run replayed from simulation state, "
       f"{t0:.0f}-{t1:.0f} s window at 1x speed. The right leg cycles "
       f"through stance and swing; the left leg remains planted "
       f"(the current limitation).")
out_local.with_suffix(".alt.txt").write_text(alt, encoding="utf-8")
(DISS / out.name).with_suffix(".alt.txt").write_text(
    alt, encoding="utf-8")
print(f"saved {out_local} + copy in Dissertation CPG_airstepping_figs "
      f"({len(frames)} frames)")
