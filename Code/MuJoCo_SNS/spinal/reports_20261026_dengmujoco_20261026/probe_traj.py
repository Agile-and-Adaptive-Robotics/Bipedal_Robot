# Probe: V[0] trajectory from step_neural vs simulate, first 0.3 s.
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent))
import deng_cpg_ode as D                     # noqa: E402
from deng_mujoco_runner import step_neural, NSUB_NEU, DT_NET  # noqa: E402

print("NSUB_NEU =", NSUB_NEU, " DT_NET =", DT_NET)
ref, _ = D.simulate("fixed", t_end=300.0, dt=0.1)
print("simulate   V_ext @ 0/50/100/200/300 ms:",
      np.round(ref[[0, 500, 1000, 2000, 3000], 0], 2))

V = D.VREST.copy()
h = np.ones(14) * 2 / 3
h[[12, 13]] = 0.0
traj = [V[0]]
for k in range(150):
    V, h = step_neural(V, h, np.zeros(14))
    traj.append(V[0])
print("step_neural V_ext @ 0/50/100/200/300 ms:",
      np.round(np.array(traj)[[0, 25, 50, 100, 150]], 2))
