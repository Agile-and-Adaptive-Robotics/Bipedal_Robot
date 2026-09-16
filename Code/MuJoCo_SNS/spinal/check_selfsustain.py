"""Self-sustainment check: does the tuned network oscillate indefinitely at
constant DRIVE = 2.5 nA, or does it settle to a fixed point (as the Simulink
E2 run did)? Runs 20 s, reports RG_E-RG_F antiphase peak times in windows.
"""
from pathlib import Path

import mujoco
import numpy as np
from scipy.signal import find_peaks

import build_network as bn
import params
from draw_circuit import effective_tables

HERE = Path(__file__).parent
t = effective_tables("best")
params.G.clear(); params.G.update(t["G"])
params.TAU.clear(); params.TAU.update(t["TAU"])
for ph in params.W_PF_MN:
    params.W_PF_MN[ph].clear(); params.W_PF_MN[ph].update(t["W"][ph])
params.W_POSTURE.clear(); params.W_POSTURE.update(t["WPOST"])

m = mujoco.MjModel.from_xml_path(
    str(HERE.parents[2] / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
        / "mjc" / "gait2392_simbody" / "gait2392_simbody_cvt3.xml"))
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(m.nu)]
net = bn.build(acts, dt=params.DT, interleg=True)
u = net.make_inputs()
u[net.input_index("DRIVE")] = 2.5

T_END = 20.0
n = int(round(T_END / params.DT))
v = np.zeros(n + 1)
for k in range(n):
    V = net.step(u)
    v[k + 1] = V[net.idx["RG_E_r"]] - V[net.idx["RG_F_r"]]
tt = np.arange(n + 1) * params.DT

pk, _ = find_peaks(v, prominence=0.3)
print("antiphase peaks t(s):", np.round(pk * params.DT, 2))
for w0 in (0.0, 5.0, 10.0, 15.0):
    msk = tt >= w0
    print(f"window {w0:.0f}-{T_END:.0f} s: swing {v[msk].max() - v[msk].min():.3f} mV")
