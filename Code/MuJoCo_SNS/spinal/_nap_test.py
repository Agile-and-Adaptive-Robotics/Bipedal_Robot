"""Deng-faithful RG with the REAL toolbox persistent-Na neurons — v2.
Adds: asymmetric 200 ms kick (symmetry break, like the runner's standing
phase), final-voltage diagnostics, e_m in the sweep. Targets: plateau
V ~ 2-4 mV, self-sustained alternation 30 s, no ADAP."""
import io
import sys

import numpy as np
from scipy.signal import find_peaks

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import (NonSpikingNeuron,
                                 NonSpikingNeuronWithPersistentSodiumChannel)
from sns_toolbox.networks import Network

E_HI = 5.0
DT = 0.002
T_END = 30.0


def _neu(tau):
    return NonSpikingNeuron(membrane_capacitance=float(tau),
                            membrane_conductance=1.0,
                            resting_potential=0.0, bias=0.0)


def _syn(g, exc):
    return NonSpikingSynapse(max_conductance=float(g),
                             reversal_potential=E_HI if exc else -E_HI,
                             e_lo=0.0, e_hi=E_HI)


def run(g_na, s_m, e_m, s_h, e_h, tau_h, g_inh=4.0, g_w=0.4, drive=2.5,
        cm=0.05):
    n = Network()
    nap = lambda: NonSpikingNeuronWithPersistentSodiumChannel(
        membrane_capacitance=cm, membrane_conductance=1.0,
        resting_potential=0.0, bias=0.0,
        g_ion=np.array([g_na]), e_ion=np.array([50.0]),
        k_m=np.array([1.0]), slope_m=np.array([s_m]),
        e_m=np.array([e_m]),
        k_h=np.array([1.0]), slope_h=np.array([s_h]),
        e_h=np.array([e_h]), tau_max_h=np.array([tau_h]))
    n.add_neuron(nap(), name="RG_E")
    n.add_neuron(nap(), name="RG_F")
    n.add_neuron(_neu(0.05), name="InE")
    n.add_neuron(_neu(0.05), name="InF")
    n.add_input("RG_E", name="DRIVE_E")
    n.add_input("RG_F", name="DRIVE_F")
    n.add_connection(_syn(g_inh, True), "RG_E", "InE")
    n.add_connection(_syn(g_inh, False), "InE", "RG_F")
    n.add_connection(_syn(g_inh, True), "RG_F", "InF")
    n.add_connection(_syn(g_inh, False), "InF", "RG_E")
    n.add_connection(_syn(g_w, True), "RG_E", "RG_F")
    n.add_connection(_syn(g_w, True), "RG_F", "RG_E")
    net = n.compile(backend="numpy", dt=DT)
    n_steps = int(round(T_END / DT))
    kick = int(round(0.2 / DT))
    v = np.zeros(n_steps + 1)
    for k in range(n_steps):
        u = [drive + (2.0 if k < kick else 0.0), drive]
        net.forward(u)
        V = net.V
        v[k + 1] = V[0] - V[1]
    m = v[int(5.0 / DT):]
    pk, _ = find_peaks(m, prominence=0.3)
    cyc = len(pk)
    per = float(np.median(np.diff(pk))) * DT if cyc >= 2 else float("nan")
    amp = float(m.max() - m.min())
    return cyc, per, amp, float(V[0]), float(V[1])


cases = [
    # g_na s_m  e_m  s_h   e_h  tau_h g_w drive  (NaP-sustained plateau:
    # drive SMALL, g_na large enough that h-inactivation ends the burst)
    (2.5, 0.8, 1.5, -2.0, 2.5, 0.50, 0.4, 0.5),
    (2.5, 0.8, 1.5, -2.0, 2.5, 1.00, 0.4, 0.5),
    (4.0, 0.8, 1.5, -2.0, 2.5, 1.00, 0.4, 0.5),
    (4.0, 0.8, 2.0, -2.0, 3.0, 1.00, 0.4, 0.5),
    (4.0, 0.8, 2.0, -2.0, 3.0, 1.00, 0.4, 1.0),
    (4.0, 0.8, 2.0, -2.0, 3.0, 0.50, 0.4, 1.0),
    (6.0, 0.8, 2.0, -2.0, 3.0, 1.00, 0.4, 1.0),
    (6.0, 0.8, 2.0, -2.0, 3.5, 1.00, 0.4, 1.0),
    (6.0, 0.8, 2.5, -2.0, 3.5, 1.00, 0.4, 1.0),
    (6.0, 0.8, 2.5, -2.0, 3.5, 1.50, 0.4, 1.0),
    (6.0, 0.8, 2.5, -2.0, 3.5, 1.50, 0.8, 1.0),
    (6.0, 0.8, 2.5, -2.0, 3.5, 1.50, 0.4, 1.5),
]
print(f"{'g_na':>5s} {'s_m':>4s} {'e_m':>4s} {'s_h':>5s} {'e_h':>4s} "
      f"{'tau_h':>5s} {'g_w':>4s} {'drv':>4s} {'cyc':>4s} {'period':>7s} "
      f"{'amp':>6s} {'VE':>5s} {'VF':>5s}")
for g_na, s_m, e_m, s_h, e_h, tau_h, g_w, drive in cases:
    cyc, per, amp, ve, vf = run(g_na, s_m, e_m, s_h, e_h, tau_h,
                                g_w=g_w, drive=drive)
    print(f"{g_na:5.2f} {s_m:4.1f} {e_m:4.1f} {s_h:5.1f} {e_h:4.1f} "
          f"{tau_h:5.2f} {g_w:4.1f} {drive:4.1f} {cyc:4d} {per:7.2f} "
          f"{amp:6.2f} {ve:5.2f} {vf:5.2f}")
