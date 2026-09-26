# -*- coding: utf-8 -*-
"""End-to-end pose check: mj body world pos (t=0) vs aproj FK (remapped); t=0 tendon lengths."""
import os
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import sys, json
import mujoco, numpy as np

HERE = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
sys.path.insert(0, HERE)
import make_w2l_mjcf as M

m = M.parse_aproj()
model = mujoco.MjModel.from_xml_path(os.path.join(HERE, "w2l_mjcf.xml"))
data = mujoco.MjData(model)
mujoco.mj_forward(model, data)

print("== world pose readback: aproj FK -> remap vs MuJoCo xpos (t=0) ==")
worst = 0.0
for b in m["bodies"]:
    if b.type != "Box":
        continue
    mj_name = b.name
    bid = model.body(mj_name).id
    xpos = data.xpos[bid]
    expect = M.al2mj_vec(b.p_world)
    dev = float(np.abs(np.array(xpos) - np.array(expect)).max())
    worst = max(worst, dev)
    print("  %-16s aproj->mj=(%.4f, %.4f, %.4f)  mj=(%.4f, %.4f, %.4f)  dev=%.2e"
          % (b.name, *expect, *xpos, dev))
print("  worst pose dev = %.3e m" % worst)

print("\n== t=0 tendon lengths (m) vs muscle rest expectations ==")
mus = {x["name"]: x for x in m["muscles"]}
for t in range(model.ntendon):
    tname = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_TENDON, t)
    L = float(data.ten_length[t])
    extra = ""
    key = tname[2:] if tname.startswith("t_") else tname
    if key in mus:
        extra = " L0=%.3f (L0-Lw=%.3f, L0+Lw=%.3f)" % (
            mus[key]["resting_length"],
            mus[key]["resting_length"] - mus[key]["lwidth"],
            mus[key]["resting_length"] + mus[key]["lwidth"])
    print("  %-16s len=%.4f%s" % (tname, L, extra))

print("\n== muscle path lengths from FK vs RestingLength (from source) ==")
for mus_i in m["muscles"]:
    pts = [m["by_id"][aid].site_pos_world_al for aid in mus_i["attach_ids"]]
    L = sum(M.vnorm(M.vsub(pts[i + 1], pts[i])) for i in range(len(pts) - 1))
    print("  %-12s path=%.3f  L0=%.3f  ratio=%.2f  maxT=%.0f N  Kse=%.0f  Kpe=%.0f  B=%.0f"
          % (mus_i["name"], L, mus_i["resting_length"], L / mus_i["resting_length"],
             mus_i["max_tension"], mus_i["kse"], mus_i["kpe"], mus_i["B"]))
