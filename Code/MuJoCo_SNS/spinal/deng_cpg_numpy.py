"""Deng 2019 / Nourse 2023 two-layer CPG in SNS-Toolbox (numpy reference).

Exact circuit per Nourse 2023 Fig 6A + Tables A4-A7 (text: _nourse2023.txt):

  RG layer:  HC_ext, HC_flx  = persistent-Na half-centers (Table A5)
             Ext_IN, Flx_IN  = plain non-spiking INs (Table A4)
    HC_ext -> Ext_IN  (2.749 uS, Esyn -40 mV, window -60..-25 mV)
    Ext_IN -> HC_flx  (2.749 uS, Esyn -70 mV, window -60..-25 mV)
    (mirror: HC_flx -> Flx_IN -> HC_ext)          <- IN-laminated mutual inhibition
  PF layer: hip pair + knee/ankle pair, same construction (8 neurons,
    16 synapses); RG HC -> same-half PF HC, weak exc 0.1 uS (window -60..-40).
  MN layer (hip only, demo): PF -> MNhip ext 2.565 / flx 3.632 uS
    (Esyn -10 mV, window -60..-50 mV); MN Vrest -100 mV.

Neurons: Cm 5 nF, Gm 1 uS; HC Vrest -60 mV, MN Vrest -100 mV (Table A4).
Na channel (Table A5): GNa 1.5 uS, ENa 50 mV, Sm 0.2 Sh -0.6, Km 1 Kh 0.5,
Em -40 mV, Eh -60 mV, tau_h.max 350 ms.
Kick: single 10 nA 1 ms pulse into HC_ext (Bolen port) at t = 0.1 s.
dt = 0.1 ms (Table A7). Run 20 s; report peak structure per window.
Saves deng_cpg_ref.mat for the Simulink side.
"""
from pathlib import Path

import numpy as np
from scipy.io import savemat
from scipy.signal import find_peaks

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import (NonSpikingNeuron,
                                 NonSpikingNeuronWithPersistentSodiumChannel)
from sns_toolbox.networks import Network

HERE = Path(__file__).parent

# ---------------- parameters (Nourse 2023 Tables A4-A6) --------------------
CM, GM = 5.0, 1.0                     # nF, uS
VREST_HC, VREST_MN = -60.0, -100.0    # mV
G_NA, E_NA = 1.5, 50.0
SM, EM, KM = 0.2, -40.0, 1.0
SH, EH, KH, TAU_H = -0.6, -60.0, 0.5, 350.0

def hc_neuron():
    return NonSpikingNeuronWithPersistentSodiumChannel(
        g_ion=np.array([G_NA]), e_ion=np.array([E_NA]),
        k_m=np.array([KM]), slope_m=np.array([SM]), e_m=np.array([EM]),
        k_h=np.array([KH]), slope_h=np.array([SH]), e_h=np.array([EH]),
        tau_max_h=np.array([TAU_H]),
        membrane_capacitance=CM, membrane_conductance=GM,
        resting_potential=VREST_HC, bias=0.0)

def in_neuron():
    return NonSpikingNeuron(membrane_capacitance=CM, membrane_conductance=GM,
                            resting_potential=VREST_HC, bias=0.0)

def mn_neuron():
    return NonSpikingNeuron(membrane_capacitance=CM, membrane_conductance=GM,
                            resting_potential=VREST_MN, bias=0.0)

def syn(g, esyn, elo, ehi):
    return NonSpikingSynapse(max_conductance=g, reversal_potential=esyn,
                             e_lo=elo, e_hi=ehi)

W = (-60.0, -25.0)    # HC<->IN synapse window
net = Network(name="deng cpg")
for nm in ("RG_HC_ext", "RG_HC_flx"):
    net.add_neuron(hc_neuron(), name=nm)
for nm in ("RG_IN_ext", "RG_IN_flx"):
    net.add_neuron(in_neuron(), name=nm)
for tag in ("hip", "ka"):
    for nm in (f"PF_{tag}_HC_ext", f"PF_{tag}_HC_flx"):
        net.add_neuron(hc_neuron(), name=nm)
    for nm in (f"PF_{tag}_IN_ext", f"PF_{tag}_IN_flx"):
        net.add_neuron(in_neuron(), name=nm)
for nm in ("MN_hip_ext", "MN_hip_flx"):
    net.add_neuron(mn_neuron(), name=nm)
net.add_input("RG_HC_ext")            # the single stimulus port

# RG mutual inhibition, IN-laminated
net.add_connection(syn(2.749, -40.0, *W), "RG_HC_ext", "RG_IN_ext")
net.add_connection(syn(2.749, -70.0, *W), "RG_IN_ext", "RG_HC_flx")
net.add_connection(syn(2.749, -40.0, *W), "RG_HC_flx", "RG_IN_flx")
net.add_connection(syn(2.749, -70.0, *W), "RG_IN_flx", "RG_HC_ext")
# RG -> PF (same half), weak exc
for tag in ("hip", "ka"):
    net.add_connection(syn(0.1, -40.0, -60.0, -40.0),
                       "RG_HC_ext", f"PF_{tag}_HC_ext")
    net.add_connection(syn(0.1, -40.0, -60.0, -40.0),
                       "RG_HC_flx", f"PF_{tag}_HC_flx")
# PF pairs, same IN-laminated construction
for tag in ("hip", "ka"):
    net.add_connection(syn(2.749, -40.0, *W), f"PF_{tag}_HC_ext", f"PF_{tag}_IN_ext")
    net.add_connection(syn(2.749, -70.0, *W), f"PF_{tag}_IN_ext", f"PF_{tag}_HC_flx")
    net.add_connection(syn(2.749, -40.0, *W), f"PF_{tag}_HC_flx", f"PF_{tag}_IN_flx")
    net.add_connection(syn(2.749, -70.0, *W), f"PF_{tag}_IN_flx", f"PF_{tag}_HC_ext")
# hip MNs
net.add_connection(syn(2.565, -10.0, -60.0, -50.0), "PF_hip_HC_ext", "MN_hip_ext")
net.add_connection(syn(3.632, -10.0, -60.0, -50.0), "PF_hip_HC_flx", "MN_hip_flx")

DT = 1e-4                              # 0.1 ms (Table A7)
c = net.compile(dt=DT, backend="numpy")
names = ["RG_HC_ext", "RG_HC_flx", "RG_IN_ext", "RG_IN_flx",
         "PF_hip_HC_ext", "PF_hip_HC_flx", "PF_ka_HC_ext", "PF_ka_HC_flx",
         "PF_hip_IN_ext", "PF_hip_IN_flx", "PF_ka_IN_ext", "PF_ka_IN_flx",
         "MN_hip_ext", "MN_hip_flx"]
idx = [c.net_params and None] if False else None  # (names index directly below)

# neurons were added in the order constructed above
order = ["RG_HC_ext", "RG_HC_flx", "RG_IN_ext", "RG_IN_flx",
         "PF_hip_HC_ext", "PF_hip_HC_flx", "PF_ka_HC_ext", "PF_ka_HC_flx",
         "PF_hip_IN_ext", "PF_hip_IN_flx", "PF_ka_IN_ext", "PF_ka_IN_flx",
         "MN_hip_ext", "MN_hip_flx"]

T_END, T_PULSE, I_PULSE = 20.0, 0.1, 10.0     # s, s, nA
n_steps = int(round(T_END / DT))
p0, p1 = int(T_PULSE / DT), int((T_PULSE + 1e-3) / DT)   # 10 nA, 1 ms
tr = np.zeros((n_steps + 1, len(order)))
for k in range(n_steps):
    u = [0.0] * 1                      # single external input
    iapp = I_PULSE if p0 <= k < p1 else 0.0
    c.forward([iapp])
    tr[k + 1] = c.V

tt = np.arange(n_steps + 1) * DT
print(f"ran {T_END} s; HC_ext range [{tr[:, 0].min():.1f}, {tr[:, 0].max():.1f}] mV")
for j, nm in enumerate(order[:8]):
    v = tr[:, j]
    pk, _ = find_peaks(v, prominence=0.3)
    per = np.diff(pk * DT)
    print(f"{nm:16s} peaks n={len(pk):3d}  period mean {per.mean() if len(per) else np.nan:.3f} s "
          f"(late swing {v[tt > 15].max() - v[tt > 15].min():.2f} mV)")

savemat(HERE / "SNS_Simscape_link" / "deng_cpg_ref.mat" if False
        else Path(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\results") / "deng_cpg_ref.mat",
        {"t": tt, "traces": tr, "names": order,
         "pulse": [T_PULSE, I_PULSE, 1e-3]})
print("saved deng_cpg_ref.mat")
