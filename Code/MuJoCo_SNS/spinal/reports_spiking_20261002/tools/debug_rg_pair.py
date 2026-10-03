import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_spiking_20261002\tools")
import numpy as np
from tune_rg_pair import build_rig, DT
m = build_rig(0.55, 5.32, 1.23, 0.005, 0.020, 0.35, 0.8)
n = int(6.0 / DT)
Vt = np.zeros((n+1, 5)); Th = np.zeros((n+1, 5)); SP = np.zeros((n+1, 5))
for k in range(n):
    m([2.5])
    Vt[k+1] = m.V; Th[k+1] = m.theta; SP[k+1] = m.spikes
names = ["DRIVE","RG_E","RG_F","InE","InF"]
for i in range(5):
    sp = np.flatnonzero(SP[:,i] == -1)
    print(f"{names[i]:6s} V[{Vt[:,i].min():7.2f},{Vt[:,i].max():7.2f}] "
          f"thr[{Th[:,i].min():7.2f},{Th[:,i].max():7.2f}] spikes={len(sp)} "
          f"first_at={sp[0]*DT if len(sp) else '-':>6}")
# window around t=2s
i0 = int(1.8/DT); i1 = int(2.4/DT)
print("t=1.8..2.4 s every 10 ms: V / theta")
for k in range(i0, i1, int(0.01/DT)):
    print(f"{k*DT:5.2f} " + " ".join(f"{Vt[k,i]:7.2f}/{Th[k,i]:6.2f}" for i in (1,2,3,4)))
