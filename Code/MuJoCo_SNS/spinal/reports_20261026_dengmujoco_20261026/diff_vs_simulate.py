# Differ: D.simulate vs the runner's step_neural, identical kick, 1 s.
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent))
import deng_cpg_ode as D                # noqa: E402
from deng_mujoco_runner import step_neural  # noqa: E402

ref, _ = D.simulate("fixed", t_end=1000.0, dt=0.1)     # kick at 100 ms
print("simulate   : V_ext max %.2f mV (min %.2f)" %
      (ref[:, 0].max(), ref[:, 0].min()))

V = D.VREST.copy()
h = np.ones(len(D.NAMES)) * (1 / (1 + D.KH * np.exp(D.SH * (D.EH - D.VREST_HC))))
h[[12, 13]] = 0.0
mx, mn = -1e9, 1e9
for k in range(int(1.0 / 0.002)):
    t = k * 0.002
    iapp = np.zeros(len(D.NAMES))
    if 0.1 <= t < 0.101:
        iapp[0] = 10.0
    V, h = step_neural(V, h, iapp)
    mx = max(mx, V[0]); mn = min(mn, V[0])
print("step_neural: V_ext max %.2f mV (min %.2f)" % (mx, mn))
# one-substep micro diff: same state, one 0.1 ms step, compare dV
V1 = D.VREST.copy()
h1 = np.ones(len(D.NAMES)) * 2 / 3
iapp = np.zeros(len(D.NAMES)); iapp[0] = 10.0
V2 = V1.copy(); h2 = h1.copy()
# module formula inline (single 0.1 ms)
sat = np.clip((V1[D.src] - D.elo) / (D.ehi - D.elo), 0, 1)
isyn = np.zeros(14); np.add.at(isyn, D.dst, D.g * sat * (D.esyn - V1[D.dst]))
m_inf = 1 / (1 + D.KM * np.exp(D.SM * (D.EM - V1)))
i_na = np.zeros(14); na = list(D.NA)
i_na[na] = D.G_NA * m_inf[na] * h1[na] ** D.KH * (D.E_NA - V1[na])
dV_ref = (D.GM * (D.VREST - V1) + iapp + isyn + i_na) * 0.1 / D.CM
print("inline dV[0:2]:", dV_ref[0:2], " | i_na[0]:", i_na[0], "m_inf[0]:",
      m_inf[0], "h[0]:", h1[0])
