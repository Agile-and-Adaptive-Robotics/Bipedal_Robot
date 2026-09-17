"""Optimize the toolbox NaP half-center parameters to hit the human gait
period (~1.2 s) with self-sustained alternation, NO ADAP.

Search: g_na, e_m, s_h, e_h, tau_max_h  (s_m=0.8, g_inh=4, g_w=0.4,
drive=1.0, cm=0.05 fixed). Objective: |period-1.2| with heavy penalties
for no-oscillation / blowup / tiny amplitude. Prints every improvement;
saves the winner to nap_rg_params.json.
"""
import io
import json
import sys

import numpy as np
from scipy.optimize import differential_evolution

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import (NonSpikingNeuron,
                                 NonSpikingNeuronWithPersistentSodiumChannel)
from sns_toolbox.networks import Network

E_HI, DT, T_END = 5.0, 0.002, 24.0
G_INH, G_W, DRIVE, CM, S_M = 4.0, 0.4, 1.0, 0.05, 0.8
TARGET_PER = 1.2


def _neu(tau):
    return NonSpikingNeuron(membrane_capacitance=float(tau),
                            membrane_conductance=1.0,
                            resting_potential=0.0, bias=0.0)


def _syn(g, exc):
    return NonSpikingSynapse(max_conductance=float(g),
                             reversal_potential=E_HI if exc else -E_HI,
                             e_lo=0.0, e_hi=E_HI)


def simulate(x):
    g_na, e_m, s_h, e_h, tau_h = x
    n = Network()
    nap = lambda: NonSpikingNeuronWithPersistentSodiumChannel(
        membrane_capacitance=CM, membrane_conductance=1.0,
        resting_potential=0.0, bias=0.0,
        g_ion=np.array([g_na]), e_ion=np.array([50.0]),
        k_m=np.array([1.0]), slope_m=np.array([S_M]),
        e_m=np.array([e_m]),
        k_h=np.array([1.0]), slope_h=np.array([s_h]),
        e_h=np.array([e_h]), tau_max_h=np.array([tau_h]))
    n.add_neuron(nap(), name="RG_E")
    n.add_neuron(nap(), name="RG_F")
    n.add_neuron(_neu(0.05), name="InE")
    n.add_neuron(_neu(0.05), name="InF")
    n.add_input("RG_E", name="DRIVE_E")
    n.add_input("RG_F", name="DRIVE_F")
    n.add_connection(_syn(G_INH, True), "RG_E", "InE")
    n.add_connection(_syn(G_INH, False), "InE", "RG_F")
    n.add_connection(_syn(G_INH, True), "RG_F", "InF")
    n.add_connection(_syn(G_INH, False), "InF", "RG_E")
    n.add_connection(_syn(G_W, True), "RG_E", "RG_F")
    n.add_connection(_syn(G_W, True), "RG_F", "RG_E")
    net = n.compile(backend="numpy", dt=DT)
    n_steps = int(round(T_END / DT))
    kick = int(round(0.2 / DT))
    v = np.zeros(n_steps + 1)
    vfull = np.zeros((n_steps + 1, 2))
    for k in range(n_steps):
        net.forward([DRIVE + (2.0 if k < kick else 0.0), DRIVE])
        V = net.V
        v[k + 1] = V[0] - V[1]
        vfull[k + 1] = (V[0], V[1])
    return v, vfull


def objective(x):
    try:
        v, vfull = simulate(x)
    except Exception:
        return 1e6
    if not np.all(np.isfinite(vfull)):
        return 1e6
    m = v[int(4.0 / DT):]
    if not np.all(np.isfinite(m)):
        return 1e6
    lo, hi = float(vfull.min()), float(vfull.max())
    if lo < -5.5 or hi > 7.0:            # left the operating range
        return 1e5 + abs(min(lo + 5.5, 0.0)) + abs(max(hi - 7.0, 0.0))
    from scipy.signal import find_peaks
    pk, _ = find_peaks(m, prominence=0.3)
    if len(pk) < 4:
        return 1e4                        # no sustained rhythm
    per = float(np.median(np.diff(pk))) * DT
    amp = float(m.max() - m.min())
    score = abs(per - TARGET_PER)
    if amp < 2.0:
        score += (2.0 - amp)              # want usable swing
    return score


bounds = [(1.0, 10.0),      # g_na
          (0.5, 3.5),       # e_m
          (-4.0, -0.5),     # s_h
          (1.0, 4.5),       # e_h
          (0.3, 4.0)]       # tau_max_h

best = {"score": float("inf")}
nit = [0]


def cb(intermediate_result=None, **kw):
    x, f = kw.get("x"), kw.get("fun")
    if f is not None and f < best["score"]:
        best.update(score=float(f), x=list(map(float, x)))
        print(f"[{nit[0]:4d}] score {f:.4f} x="
              f"{np.round(x, 3).tolist()}", flush=True)


res = differential_evolution(objective, bounds, maxiter=60, popsize=10,
                             tol=1e-3, seed=7, polish=True,
                             updating="deferred", workers=-1, callback=cb)
print("\n== winner ==")
print("score", res.fun, "x", res.x.tolist())
v, vfull = simulate(res.x)
from scipy.signal import find_peaks
m = v[int(4.0 / DT):]
pk, _ = find_peaks(m, prominence=0.3)
per = float(np.median(np.diff(pk))) * DT if len(pk) >= 2 else float("nan")
print(f"cycles {len(pk)} period {per:.3f} s amp {m.max()-m.min():.2f} mV "
      f"V range {vfull.min():.2f}..{vfull.max():.2f}")
with open("nap_rg_params.json", "w", encoding="utf-8") as f:
    json.dump({"g_na": res.x[0], "e_m": res.x[1], "s_h": res.x[2],
               "e_h": res.x[3], "tau_max_h": res.x[4], "s_m": S_M,
               "g_inh": G_INH, "g_w": G_W, "drive": DRIVE, "cm": CM,
               "period_s": per, "score": float(res.fun)}, f, indent=2)
print("saved nap_rg_params.json")
