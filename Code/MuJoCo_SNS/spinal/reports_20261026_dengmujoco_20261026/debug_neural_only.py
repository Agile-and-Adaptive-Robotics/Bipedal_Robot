# Debug: neural-only run of the runner's step_neural (no MuJoCo).
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from deng_mujoco_runner import step_neural  # noqa: E402
import deng_cpg_ode as D                    # noqa: E402

V = D.VREST.copy()
h = np.ones(len(D.NAMES)) * (1 / (1 + D.KH * np.exp(D.SH * (D.EH - D.VREST_HC))))
h[[12, 13]] = 0.0
for k in range(int(12.0 / 0.002)):
    t = k * 0.002
    iapp = np.zeros(len(D.NAMES))
    if 0.1 <= t < 0.101:
        iapp[0] = 10.0
    V, h = step_neural(V, h, iapp)
print("V[0] (RG ext) range: %.2f .. %.2f mV" % (V.min() * 0 + min(V[0], 0), V[0]))
# track extrema over run instead: redo with logging
V = D.VREST.copy()
h = np.ones(len(D.NAMES)) * (1 / (1 + D.KH * np.exp(D.SH * (D.EH - D.VREST_HC))))
h[[12, 13]] = 0.0
lo, hi = 1e9, -1e9
for k in range(int(12.0 / 0.002)):
    t = k * 0.002
    iapp = np.zeros(len(D.NAMES))
    if 0.1 <= t < 0.101:
        iapp[0] = 10.0
    V, h = step_neural(V, h, iapp)
    lo = min(lo, V[0]); hi = max(hi, V[0])
print("RG ext over 12 s: %.2f .. %.2f mV (bursting needs swing > ~30 mV)"
      % (lo, hi))
