# -*- coding: utf-8 -*-
"""A/B the toe-spring damping: does c=20000 cause the QACC spike? (2 s passive drop)"""
import os
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import mujoco, numpy as np, io, contextlib

XML = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco\w2l_mjcf.xml"

def run(damp_val, tag):
    tree = open(XML).read()
    if damp_val is not None:
        import re
        tree = re.sub(r'(<spatial name="t_Body_11"[^/]*?)damping="[0-9.]+"',
                      r'\1damping="%g"' % damp_val, tree)
        tree = re.sub(r'(<spatial name="t_toe_R_spring"[^/]*?)damping="[0-9.]+"',
                      r'\1damping="%g"' % damp_val, tree)
    model = mujoco.MjModel.from_xml_string(tree)
    data = mujoco.MjData(model)
    max_qacc = 0.0
    for i in range(2000):
        mujoco.mj_step(model, data)
        max_qacc = max(max_qacc, float(np.abs(data.qacc).max()))
        assert np.isfinite(data.qpos).all()
    fz = min(data.xpos[model.body("foot_L").id][2], data.xpos[model.body("foot_R").id][2])
    print("%-14s max|qacc|=%.4g  final foot z=%.4f  toe angle(deg)=%.3f" %
          (tag, max_qacc, fz, np.degrees(data.qpos[-1]) if model.nq else 0))

run(None, "source c=20000")
run(800.0, "c=800")
run(400.0, "c=400")
run(172.0, "c=172 (crit)")
run(0.0, "c=0")
