"""Verify the Simulink session's tau_h(V) collapse claim for the
persistent-Na h-gate as parameterized in Nourse 2023 Table A5:
  z_inf(V) = 1 / (1 + K*exp(S*(E - V)))
  tau_z(V) = tau_max * z_inf(V) * sqrt(K*exp(S*(E - V)))
h-gate: S = -0.6, E = -60 mV, tau_max = 350 ms. K is the one free
parameter (SNS-Toolbox's default K for tau functions) - sweep it.
Deng's HCs operate around Vrest -60 mV bursting toward ~ -40 mV.
"""
import io
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

S, E, TAU_MAX = -0.6, -60.0, 350.0


def tau_v(v, k):
    z_inf = 1.0 / (1.0 + k * np.exp(S * (E - v)))
    return TAU_MAX * z_inf * np.sqrt(k * np.exp(S * (E - v)))


print("tau_h(V) [ms] for K candidates; V from -70 (hyperpol) to 0 mV")
vs = np.arange(-70.0, 1.0, 10.0)
print("V[mV]   " + "  ".join(f"K={k:<6g}" for k in (0.01, 0.1, 1.0)))
for v in vs:
    row = "  ".join(f"{tau_v(v, k):8.2f}" for k in (0.01, 0.1, 1.0))
    print(f"{v:5.0f}   {row}")
# claim: collapses to ~0.1 ms at depolarized voltages
print("\nminimum over V in [-40, 0] for each K:")
for k in (0.01, 0.1, 1.0):
    vv = np.linspace(-40.0, 0.0, 200)
    t = tau_v(vv, k)
    print(f"  K={k:<6g}: min tau_h = {t.min():.3f} ms at V={vv[t.argmin()]:.0f}")
print("\nreference (Simulink/Deng/Animatlab): FIXED tau_h = 350 ms")
