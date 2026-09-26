import numpy as np
from pathlib import Path
from scipy.signal import find_peaks
HERE = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
log = np.load(HERE / "reports_20260925" / "logs" / "syn6_net_smoke.npy")
watch = (["RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l"]
         + [f"PF_S{k}_r" for k in range(1, 7)]
         + [f"PF_S{k}_l" for k in range(1, 7)])
T = 12.0
dt = 0.002
tt = np.arange(log.shape[0]) * dt
msk = tt >= 2.0
for sd in ("r", "l"):
    P = np.array([log[msk, watch.index(f"PF_S{k}_{sd}")] for k in range(1, 7)])
    c = np.corrcoef(P)
    print(f"side {sd} pairwise corr (S1..S6):")
    print(np.round(c, 3))
    # peak phase within the RG cycle (period 0.904 s)
    rg_e = log[msk, watch.index(f"RG_E_{sd}")]
    pk, _ = find_peaks(rg_e, prominence=0.3)
    period = np.mean(np.diff(tt[pk]))
    phases = []
    for k in range(6):
        p, _ = find_peaks(P[k], prominence=0.2)
        if len(p) == 0:
            phases.append(float("nan"))
            continue
        ph = (tt[p] % period) / period * 100
        # circular mean
        ang = np.deg2rad(ph * 3.6)
        phases.append(float(np.rad2deg(np.angle(np.mean(np.exp(1j * ang)))) / 3.6 % 100))
    print(f"  mean peak phase %cycle: {np.round(phases, 1).tolist()}")
