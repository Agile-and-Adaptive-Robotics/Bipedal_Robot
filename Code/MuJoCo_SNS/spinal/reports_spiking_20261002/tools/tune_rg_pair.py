"""Sweep the spiking RG half-center pair (RG-E/RG-F + InE/InF + graded
DRIVE) over the mirror's own cell parameters to find the configuration
that ALTERNATES at the analog winner's period (0.40 s at the 'best'
tables, DRIVE 2.5 - measured by check_selfsustain_spiking.py).

This is the v1 calibration loop SPIKING_MIRROR_PLAN.md anticipates
("Parameter-mapping rules (v1 starting points, to be calibrated)"): the
synapse scales come from calibrate_spiking.py; this probe calibrates the
RG cell's spike-frequency-adaptation constants (thr_inc) and the DRIVE
rate-map gain (k_ns) against the analog rhythm, on a 5-neuron rig that
steps in milliseconds.

Usage: D:\\Anaconda\\envs\\myo\\python.exe tune_rg_pair.py
"""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np
from scipy.signal import find_peaks

from sns_toolbox.connections import NonSpikingSynapse, SpikingSynapse
from sns_toolbox.neurons import NonSpikingNeuron, SpikingNeuron
from sns_toolbox.networks import Network

DT = 0.0005
E_HI = 5.0
T_END = 20.0
DRIVE_I = 2.5          # nA into the analog DRIVE cell (winner operating pt)
G = dict(rg_mutual_inh=4.0, descend_to_rg_e=1.517816924128661,
         descend_to_rg_f=1.4)   # 'best' tables


def rg_cell(tau_theta, thr_inc, thr=-50.0):
    return SpikingNeuron(
        threshold_time_constant=tau_theta,
        threshold_initial_value=thr, threshold_proportionality_constant=0.0,
        threshold_leak_rate=1.0, threshold_increment=thr_inc,
        threshold_floor=-48.0, reset_potential=-60.0,
        membrane_capacitance=0.05, membrane_conductance=1.0,
        resting_potential=-70.0, bias=16.0)


def relay():
    return SpikingNeuron(
        threshold_time_constant=0.10, threshold_initial_value=-50.0,
        threshold_proportionality_constant=0.0, threshold_leak_rate=1.0,
        threshold_increment=0.0, threshold_floor=-50.0,
        reset_potential=-60.0, membrane_capacitance=0.05,
        membrane_conductance=1.0, resting_potential=-70.0, bias=16.0)


def build_rig(k_ns, k_s2s_exc, k_s2s_inh, tau_exc, tau_inh, tau_theta,
              thr_inc):
    net = Network(name="rg pair")
    net.add_neuron(NonSpikingNeuron(
        membrane_capacitance=0.10, membrane_conductance=1.0,
        resting_potential=0.0, bias=0.0), name="DRIVE")
    net.add_input("DRIVE")
    net.add_neuron(rg_cell(tau_theta, thr_inc), name="RG_E")
    net.add_neuron(rg_cell(tau_theta, thr_inc), name="RG_F")
    net.add_neuron(relay(), name="InE")
    net.add_neuron(relay(), name="InF")

    def ns2s(g, exc):
        return NonSpikingSynapse(max_conductance=g * k_ns,
                                 reversal_potential=0.0 if exc else -70.0,
                                 e_lo=0.0, e_hi=E_HI)

    def ss(g, exc):
        k = k_s2s_exc if exc else k_s2s_inh
        g_inc = g * k
        return SpikingSynapse(max_conductance=8 * g_inc,
                              reversal_potential=0.0 if exc else -70.0,
                              time_constant=tau_exc if exc else tau_inh,
                              transmission_delay=0,
                              conductance_increment=g_inc)

    net.add_connection(ns2s(G["descend_to_rg_e"], True), "DRIVE", "RG_E")
    net.add_connection(ns2s(G["descend_to_rg_f"], True), "DRIVE", "RG_F")
    net.add_connection(ss(G["rg_mutual_inh"], True), "RG_E", "InE")
    net.add_connection(ss(G["rg_mutual_inh"], False), "InE", "RG_F")
    net.add_connection(ss(G["rg_mutual_inh"], True), "RG_F", "InF")
    net.add_connection(ss(G["rg_mutual_inh"], False), "InF", "RG_E")
    return net.compile(backend="numpy", dt=DT)


import json as _json
import pathlib as _pl
_CAL = _json.loads((_pl.Path(r'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal') / 'spiking_calibration.json').read_text(encoding='utf-8'))


def run2(k_ns, tau_theta, thr_inc, k_exc=None, k_inh=None):
    k_exc = _CAL['k_s2s_exc'] if k_exc is None else k_exc
    k_inh = _CAL['k_s2s_inh'] if k_inh is None else k_inh
    m = build_rig(k_ns, k_exc, k_inh, 0.005, 0.020, tau_theta, thr_inc)
    n = int(T_END / DT)
    st_e = np.zeros(n + 1); st_f = np.zeros(n + 1)
    for k in range(n):
        m([DRIVE_I])
        st_e[k + 1] = st_e[k] + (m.spikes[1] == -1)
        st_f[k + 1] = st_f[k] + (m.spikes[2] == -1)
    tt = np.arange(n + 1) * DT
    # burst onset = rises in the smoothed spike-count derivative
    de = np.diff(st_e, prepend=0.0)
    win = int(0.05 / DT)
    env = np.convolve(de, np.ones(win) / win, mode="same")
    pk, _ = find_peaks(env, height=0.25 * max(env.max(), 1e-9), distance=int(0.1 / DT))
    per = float(np.mean(np.diff(tt[pk]))) if len(pk) >= 3 else float("nan")
    both = st_e[-1] > 20 and st_f[-1] > 20
    return per, st_e[-1] / T_END, st_f[-1] / T_END, len(pk)


def main():
    print(f"target: analog period 0.399 s at DRIVE {DRIVE_I}")
    print(f"{'k_ns':>6s} {'tau_th':>7s} {'inc':>5s} {'period':>8s} "
          f"{'fE':>6s} {'fF':>6s} {'bursts':>6s} {'alt':>4s}")
    best = None
    for k_ns in (0.25, 0.35, 0.45):
        for tau_theta in (0.15, 0.25, 0.35):
            for thr_inc in (1.5, 2.0, 2.5, 3.0, 4.0):
                per, fe, ff, nb = run2(k_ns, tau_theta, thr_inc)
                alt = fe > 3 and ff > 3
                print(f"{k_ns:6.2f} {tau_theta:7.2f} {thr_inc:5.2f} "
                      f"{per:8.3f} {fe:6.1f} {ff:6.1f} {nb:6d} "
                      f"{str(alt):>4s}")
                if alt and np.isfinite(per) and \
                        abs(per - 0.399) <= 0.08:
                    if best is None or abs(per - 0.399) < abs(best[0] - 0.399):
                        best = (per, k_ns, tau_theta, thr_inc, fe, ff)
    print("BEST:", best)


if __name__ == "__main__":
    main()
