import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_spiking_20261002\tools")
import numpy as np
import tune_rg_pair as T
m = T.build_rig(0.70, T._CAL["k_s2s_exc"], T._CAL["k_s2s_inh"], 0.005, 0.020, 0.35, 0.8)
n = int(6.0 / T.DT)
Vt = np.zeros((n+1, 5)); Th = np.zeros((n+1, 5)); SP = np.zeros((n+1, 5)); Gf = []
for k in range(n):
    m([2.5])
    Vt[k+1] = m.V; Th[k+1] = m.theta; SP[k+1] = m.spikes
names = ["DRIVE","RG_E","RG_F","InE","InF"]
for i in range(5):
    sp = np.flatnonzero(SP[:,i] == -1)
    print(f"{names[i]:6s} V[{Vt[:,i].min():7.2f},{Vt[:,i].max():7.2f}] spikes={len(sp)}")
i0 = int(2.0/T.DT); i1 = int(2.3/T.DT)
for k in range(i0, i1, int(0.01/T.DT)):
    print(f"{k*T.DT:5.2f} E={Vt[k,1]:7.2f}/{Th[k,1]:6.2f} F={Vt[k,2]:7.2f}/{Th[k,2]:6.2f} InE={Vt[k,3]:7.2f}")
