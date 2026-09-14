"""Deng 2019 two-layer CPG - direct ODE reference (controls tau_h exactly).

Same circuit/parameters as deng_cpg_numpy.py (Nourse 2023 Fig 6A, Tables
A4-A7) but integrated directly so the Na h-gate time constant can be FIXED
at tau_h.max = 350 ms (Animatlab LinearHill interpretation) vs the
sns_toolbox voltage-dependent tau (collapses at depolarized V - suspected
reason the toolbox realization quenches and the Animatlab port latches).

Synapses (all): Isyn = g*clip((Vpre-Elo)/(Ehi-Elo),0,1)*(Esyn-Vpost) [nA]
Neurons: Cm dV/dt [ms] = Gm*(Vr-V) + Iapp + Isyn (+ I_Na on HCs)
I_Na = GNa*m*h^Kh*(ENa-V);  m = 1/(1+Km*exp(Sm*(Em-V)))  (instantaneous)
h: dh/dt = (h_inf - h)/tau_h,  h_inf = 1/(1+Kh*exp(Sh*(Eh-V)))
"""
from pathlib import Path

import numpy as np
from scipy.io import savemat
from scipy.signal import find_peaks

HERE = Path(__file__).parent

CM, GM = 5.0, 1.0
VREST_HC, VREST_MN = -60.0, -100.0
G_NA, E_NA = 1.5, 50.0
SM, EM, KM = 0.2, -40.0, 1.0
SH, EH, KH = -0.6, -60.0, 0.5
TAU_H = 350.0                       # ms, fixed (Animatlab interpretation)

# neuron order and types
NAMES = ["RG_HC_ext", "RG_HC_flx", "RG_IN_ext", "RG_IN_flx",
         "PF_hip_HC_ext", "PF_hip_HC_flx", "PF_ka_HC_ext", "PF_ka_HC_flx",
         "PF_hip_IN_ext", "PF_hip_IN_flx", "PF_ka_IN_ext", "PF_ka_IN_flx",
         "MN_hip_ext", "MN_hip_flx"]
NA = {0, 1, 4, 5, 6, 7}                       # persistent-Na HCs
VREST = np.full(len(NAMES), VREST_HC)
VREST[12] = VREST[13] = VREST_MN

# synapses: (src, dst, g, Esyn, Elo, Ehi)
W = (-60.0, -25.0)
S = []
def add(a, b, g, e, w):
    S.append((NAMES.index(a), NAMES.index(b), g, e, w[0], w[1]))
# RG mutual inhibition (IN-laminated)
add("RG_HC_ext", "RG_IN_ext", 2.749, -40.0, W)
add("RG_IN_ext", "RG_HC_flx", 2.749, -70.0, W)
add("RG_HC_flx", "RG_IN_flx", 2.749, -40.0, W)
add("RG_IN_flx", "RG_HC_ext", 2.749, -70.0, W)
# RG -> PF same-half weak exc
for tag in ("hip", "ka"):
    add("RG_HC_ext", f"PF_{tag}_HC_ext", 0.1, -40.0, (-60.0, -40.0))
    add("RG_HC_flx", f"PF_{tag}_HC_flx", 0.1, -40.0, (-60.0, -40.0))
# PF pairs
for tag in ("hip", "ka"):
    add(f"PF_{tag}_HC_ext", f"PF_{tag}_IN_ext", 2.749, -40.0, W)
    add(f"PF_{tag}_IN_ext", f"PF_{tag}_HC_flx", 2.749, -70.0, W)
    add(f"PF_{tag}_HC_flx", f"PF_{tag}_IN_flx", 2.749, -40.0, W)
    add(f"PF_{tag}_IN_flx", f"PF_{tag}_HC_ext", 2.749, -70.0, W)
# hip MNs
add("PF_hip_HC_ext", "MN_hip_ext", 2.565, -10.0, (-60.0, -50.0))
add("PF_hip_HC_flx", "MN_hip_flx", 3.632, -10.0, (-60.0, -50.0))
S = np.array(S)
src, dst, g, esyn, elo, ehi = (S[:, 0].astype(int), S[:, 1].astype(int),
                               S[:, 2], S[:, 3], S[:, 4], S[:, 5])

def simulate(tau_h_mode="fixed", t_end=20.0, dt=0.1, i_pulse=10.0,
             t_pulse=100.0, pulse_ms=1.0):
    n = int(round(t_end / dt))
    V = VREST.copy()
    h = np.ones(len(NAMES)) * (1 / (1 + KH * np.exp(SH * (EH - VREST_HC))))
    h[[12, 13]] = 0.0
    sat = np.zeros(len(S))
    out = np.zeros((n + 1, len(NAMES)))
    hh = np.zeros(n + 1)
    k0, k1 = int(t_pulse / dt), int((t_pulse + pulse_ms) / dt)
    for k in range(n):
        sat = np.clip((V[src] - elo) / (ehi - elo), 0.0, 1.0)
        isyn = np.zeros(len(NAMES))
        np.add.at(isyn, dst, g * sat * (esyn - V[dst]))
        m_inf = 1 / (1 + KM * np.exp(SM * (EM - V)))
        h_inf = 1 / (1 + KH * np.exp(SH * (EH - V)))
        if tau_h_mode == "fixed":
            tau = np.full(len(NAMES), TAU_H)
        else:   # sns_toolbox form
            x = KH * np.exp(SH * (EH - V))
            tau = TAU_H * h_inf * np.sqrt(np.maximum(x, 1e-300))
        i_na = np.zeros(len(NAMES))
        i_na[list(NA)] = (G_NA * m_inf[list(NA)] * h[list(NA)] ** KH
                          * (E_NA - V[list(NA)]))
        iapp = np.zeros(len(NAMES))
        if k0 <= k < k1:
            iapp[0] = i_pulse
        dV = (GM * (VREST - V) + iapp + isyn + i_na) * dt / CM
        V = V + dV
        h = h + dt * (h_inf - h) / tau
        out[k + 1] = V
        hh[k + 1] = h[0]
    return out, hh

for mode in ("fixed", "toolbox"):
    out, hh = simulate(mode, t_end=20000.0)
    tt = np.arange(out.shape[0]) * 0.1
    v = out[:, 0]
    pk, _ = find_peaks(v, prominence=1.0)
    per = np.diff(pk)
    late = out[tt > 15000]
    print(f"tau_h {mode:8s}: HC_ext [{v.min():.1f},{v.max():.1f}] mV, peaks "
          f"{len(pk)}, period {per.mean() if len(per) else np.nan:.0f} ms, "
          f"late swing {late[:, 0].max() - late[:, 0].min():.2f} mV")

# fixed-tau run is the reference; save it
out, hh = simulate("fixed", t_end=20000.0)
tt = np.arange(out.shape[0]) * 0.1
savemat(Path(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\results")
        / "deng_cpg_ref.mat",
        {"t": tt, "traces": out, "names": NAMES, "h_RG_HC_ext": hh,
         "pulse": [0.1, 10.0, 1.0], "tau_h_mode": "fixed"})
print("saved deng_cpg_ref.mat (fixed-tau run)")
