# Pure-Python MuJoCo 2.3.7 stepping time, same machine, small models.
import glob
import time
import mujoco

BASE = mujoco.__file__.rsplit("\\", 1)[0].rsplit("/", 1)[0]
cands = (glob.glob(BASE + "/models/*.xml")
         + glob.glob("D:/Github/Bipedal_Robot/Code/Matlab/SNS_Simscape/"
                     "mujoco_bridge/blockset/mujoco-simulink-blockset-main/"
                     "lib/win64/mujoco/model/**/*.xml", recursive=True))
picks = [c for c in cands if "cartpole" in c] or cands[:1]
picks += [c for c in cands if "humanoid" in c][:1]
for p in picks:
    m = mujoco.MjModel.from_xml_path(p)
    d = mujoco.MjData(m)
    mujoco.mj_resetData(m, d)
    t0 = time.perf_counter()
    for i in range(2000):
        mujoco.mj_step(m, d)
    el = time.perf_counter() - t0
    print("PYTHON-MUJOCO %s: nq=%d 2000 steps in %.3f s (%.0f steps/s)"
          % (p.rsplit("/", 1)[-1].rsplit("\\", 1)[-1], m.nq, el, 2000 / el))
