"""2-neuron units-mapping reference (SNS-Toolbox numpy side).

Circuit: constant 5 nA into neuron A (tau 0.2 s); A -E-> B (tau 0.5 s),
synapse g = 0.3 uS, Esyn = +8 mV, e_lo/e_hi = 0/5 mV — the exact parameter
family used by the spinal network. Steady states (analytic):
  V_A(inf) = I/G = 5 mV          (Sat_A = 1)
  V_B(inf) = g*(E - V_B)/G  ->  V_B = 0.3*8/1.3 = 1.8462 mV
Saves t, V_A, V_B to SNS_Simscape/results/units_ref_2n.mat for the Simulink
side of the units test.
"""
from pathlib import Path

import numpy as np
from scipy.io import savemat
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.networks import Network

import build_network as bn          # reuse _neu/_syn: identical conventions

here = Path(__file__).resolve().parent
out_dir = here.parents[1] / "Matlab" / "SNS_Simscape" / "results"
out_dir.mkdir(exist_ok=True)

TAU_A, TAU_B, G_SYN, I_EXT = 0.2, 0.5, 0.3, 5.0

net = Network(name="units2n")
net.add_neuron(bn._neu(TAU_A), name="A")
net.add_neuron(bn._neu(TAU_B), name="B")
net.add_input("A")
net.add_connection(bn._syn(G_SYN, exc=True), "A", "B")
c = net.compile(dt=1e-4, backend="numpy")

t_end, dt = 3.0, 1e-4
n_steps = int(round(t_end / dt))
va = np.zeros(n_steps + 1)
vb = np.zeros(n_steps + 1)
for k in range(n_steps):
    c.forward([I_EXT])
    va[k + 1] = c.V[0]
    vb[k + 1] = c.V[1]
t = np.arange(n_steps + 1) * dt

va_ss, vb_ss = I_EXT / 1.0, G_SYN * 8.0 / (1.0 + G_SYN)
ok = abs(va[-1] - va_ss) < 1e-3 and abs(vb[-1] - vb_ss) < 5e-3  # V_B still settling
print(f"numpy: V_A end {va[-1]:.5f} (analytic {va_ss:.5f}), "
      f"V_B end {vb[-1]:.5f} (analytic {vb_ss:.5f})  "
      f"{'OK' if ok else 'MISMATCH'}")
if not ok:
    raise SystemExit(1)

savemat(out_dir / "units_ref_2n.mat",
        {"t": t, "va": va, "vb": vb,
         "params": {"TAU_A": TAU_A, "TAU_B": TAU_B, "G_SYN": G_SYN,
                    "I_EXT": I_EXT, "ESYN": bn.E_REV_EXC,
                    "E_LO": bn.SYN_E_LO, "E_HI": bn.SYN_E_HI}})
print(f"wrote {out_dir / 'units_ref_2n.mat'}")
