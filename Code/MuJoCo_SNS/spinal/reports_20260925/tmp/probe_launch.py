# -*- coding: utf-8 -*-
"""Diagnose the launch: tendon forces + foot height in first 100 ms; spring damping A/B."""
import os
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import mujoco, numpy as np

XML = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_mjcf.xml"
model = mujoco.MjModel.from_xml_path(XML)
data = mujoco.MjData(model)
mujoco.mj_forward(model, data)
print("t=0: tendon len:", np.round(data.ten_length, 4))
print("t=0: spring tendon forces:", np.round(data.ten_length[12:], 4))
fz0 = data.xpos[model.body("foot_L").id][2]
print("t=0: foot_L z =", fz0)
for i in range(100):
    mujoco.mj_step(model, data)
    if i % 10 == 9:
        fl = data.xpos[model.body("foot_L").id][2]
        print("  t=%.3f foot z=%.4f  springF=%s  max|qvel|=%.3g" %
              (data.time, fl, np.round(data.ten_length[12:] * 16000, 2), np.abs(data.qvel).max()))
# estimate spring-driven energy check done; report max qacc too
print("final max|qacc|:", np.abs(data.qacc).max())
