import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_spiking_20261002\tools")
import numpy as np
import tune_rg_pair as T
m = T.build_rig(0.25, T._CAL["k_s2s_exc"], T._CAL["k_s2s_inh"], 0.005, 0.020, 0.35, 1.5)
n = int(6.0 / T.DT)
Vt = np.zeros((n+1, 5)); SPt = np.zeros((n+1, 5))
gF = np.zeros(n+1)   # total inhibitory g on F (InE->F)
for k in range(n):
    m([2.5])
    Vt[k+1] = m.V; SPt[k+1] = m.spikes
    gF[k+1] = m.g_spike[2, 3]  # g_spike[dest, src]: F<-InE
names = ["DRIVE","RG_E","RG_F","InE","InF"]
for i in range(5):
    sp = np.flatnonzero(SPt[:,i] == -1)
    print(f"{names[i]:6s} V[{Vt[:,i].min():7.2f},{Vt[:,i].max():7.2f}] spikes={len(sp)}")
i0 = int(2.0/T.DT); i1 = int(2.6/T.DT)
spE = np.flatnonzero(SPt[:,1] == -1)
spInE = np.flatnonzero(SPt[:,3] == -1)
print("E spike times 2.0-2.6:", np.round(spE[(spE> i0)&(spE<i1)]*T.DT, 3))
print("InE spike times:", np.round(spInE[(spInE>i0)&(spInE<i1)]*T.DT, 3))
for k in range(i0, i1, int(0.02/T.DT)):
    print(f"t={k*T.DT:5.2f} E={Vt[k,1]:7.2f} th={m.theta[1]:6.2f} F={Vt[k,2]:7.2f} gF={gF[k]:7.3f}")
