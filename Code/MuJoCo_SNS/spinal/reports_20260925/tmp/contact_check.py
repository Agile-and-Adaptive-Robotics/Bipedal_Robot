# -*- coding: utf-8 -*-
"""Which geoms touch at settled pose?"""
import os
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import mujoco, numpy as np
HERE = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
model = mujoco.MjModel.from_xml_path(os.path.join(HERE, "w2l_mjcf.xml"))
data = mujoco.MjData(model)
for i in range(2000):
    mujoco.mj_step(model, data)
print("t=%.2f ncon=%d" % (data.time, data.ncon))
for c in range(data.ncon):
    con = data.contact[c]
    g1 = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_GEOM, con.geom1)
    g2 = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_GEOM, con.geom2)
    print("  contact: %-16s <-> %-16s dist=%.4f" % (g1, g2, con.dist))
