"""Smoke test: can mujoco 2.3.7 do offscreen rendering on this box?
Tries mujoco.Renderer; on failure prints the exception and returns 1.
Also dumps body names so the matplotlib fallback skeleton can use them.
"""
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")

HERE = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
FIGS = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\figs"
sys.path.insert(0, HERE)

import mujoco  # noqa: E402
import numpy as np  # noqa: E402

import test_w2l_air as TW  # noqa: E402  (parses OUR argv; defaults are fine)
TW.KV["lift"] = 0.30
TW.write_air_xml()

m = mujoco.MjModel.from_xml_path(TW.AIR)
d = mujoco.MjData(m)
print("nbody", m.nbody, "nq", m.nq, "nu", m.nu)
for i in range(m.nbody):
    print("body", i, mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_BODY, i))
mujoco.mj_forward(m, d)
com = d.subtree_com[0]
print("subtree COM:", com)

ok = False
try:
    r = mujoco.Renderer(m, height=420, width=360)
    cam = mujoco.MjvCamera()
    cam.lookat[:] = [com[0], com[1], 0.95]
    cam.azimuth = 90.0
    cam.elevation = -5.0
    cam.distance = 2.6
    r.update_scene(d, camera=cam)
    rgb = r.render()
    from PIL import Image
    Image.fromarray(rgb).save(os.path.join(FIGS, "tools", "_smoke_render.png"))
    print("RENDERER OK  rgb", rgb.shape, "nonzero px", int((rgb > 0).any(axis=2).sum()))
    r.close()
    ok = True
except Exception as e:  # noqa: BLE001
    print("RENDERER FAIL:", type(e).__name__, e)

sys.exit(0 if ok else 1)
