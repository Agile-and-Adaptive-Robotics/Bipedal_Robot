import os
os.environ.pop("AARL_NET", None)
import sys
sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
import numpy as np
import mujoco
import build_network as bn
from pathlib import Path
HERE = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
xml = (HERE.parents[2] / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
       / "mjc" / "gait2392_simbody" / "gait2392_simbody_cvt3.xml")
m = mujoco.MjModel.from_xml_path(str(xml))
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(m.nu)]
net = bn.build(acts, dt=0.002, interleg=True)
c = net.compiled
g = c.g_max_non
print("g_max_non shape:", g.shape, "nonzero:", int(np.count_nonzero(g)))
print("nonzero rows:", int(np.count_nonzero((g != 0).any(axis=1))))
print("net.net.connections type:", type(net.net.connections),
      "len:", len(net.net.connections) if hasattr(net.net.connections, "__len__") else "n/a")
