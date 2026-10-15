import io
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import mujoco
import runner

z = np.load("spinal_run.npz", allow_pickle=True)
q = np.asarray(z["q"], float)
model, data = runner.build_air_model() if hasattr(runner, "build_air_model") \
    else (None, None)
m = runner.MODEL and __import__("mujoco", fromlist=["MjModel"])
import mujoco as mj
mm = mj.MjModel.from_xml_path(str(runner.MODEL))
dof_names = []
for i in range(mm.njnt):
    nm = mj.mj_id2name(mm, mj.mjtObj.mjOBJ_JOINT, i)
    adr = mm.jnt_dofadr[i]
    typ = mm.jnt_type[i]
    ndof = {0: 0, 1: 1, 2: 1, 3: 3}[int(typ)]
    for k in range(ndof):
        dof_names.append(f"{nm}.{k}")
for col in range(min(12, q.shape[1])):
    v = np.degrees(q[:, col])
    print(f"qpos[{col:2d}] {dof_names[col] if col < len(dof_names) else '?':<26s} "
          f"{v.min():9.1f}..{v.max():9.1f}")
