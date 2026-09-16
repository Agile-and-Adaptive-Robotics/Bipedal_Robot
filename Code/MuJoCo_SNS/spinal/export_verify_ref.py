"""Full-network verification reference for the generated Simulink model.

Rebuilds the tuned network (same composite as the export), drives it with
u = 0 except DRIVE = 2.5 nA for 2 s (1000 steps at DT), and saves the final
membrane-potential vector (410) + traces of a few key cells to
SNS_Simscape/results/verify_ref.mat. The Simulink model must land on the
same state (same wiring, same units) within integrator tolerance.
"""
from pathlib import Path

import mujoco
import numpy as np
from scipy.io import savemat

import build_network as bn
import params
from draw_circuit import effective_tables

HERE = Path(__file__).resolve().parent
MODEL = (HERE.parents[2] / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
         / "mjc" / "gait2392_simbody" / "gait2392_simbody_cvt3.xml")

t = effective_tables("best")
params.G.clear(); params.G.update(t["G"])
params.TAU.clear(); params.TAU.update(t["TAU"])
for ph in params.W_PF_MN:
    params.W_PF_MN[ph].clear(); params.W_PF_MN[ph].update(t["W"][ph])
params.W_POSTURE.clear(); params.W_POSTURE.update(t["WPOST"])

m = mujoco.MjModel.from_xml_path(str(MODEL))
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(m.nu)]
net = bn.build(acts, dt=params.DT, interleg=True)

u = net.make_inputs()
u[net.input_index("DRIVE")] = 2.5

# TWO references: the production dt = 2 ms Euler step the runner uses, and a
# fine dt = 0.1 ms integration (integrator-truth). The Simulink model (ode45)
# should match the FINE one to solver tolerance; the 2 ms one carries Euler
# phase error of its own (the relaxation oscillator's switching instants).
nets = {}
for tag, dt in (("coarse", params.DT), ("fine", 1e-4)):
    m2 = mujoco.MjModel.from_xml_path(str(MODEL))
    acts2 = [mujoco.mj_id2name(m2, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
             for i in range(m2.nu)]
    nets[tag] = bn.build(acts2, dt=dt, interleg=True)

names = list(net.idx)
trace_names = ["DRIVE", "RG_E_r", "RG_F_r", "ADAP_E_r", "ADAP_F_r",
               "PF_E1_r", "PF_F1_r", "MN_vas_med_r", "MN_glut_max1_r"]
tidx = [names.index(n) for n in trace_names]

finals, traces, mids = {}, {}, {}
for tag, nw in nets.items():
    n_steps = int(round(2.0 / (params.DT if tag == "coarse" else 1e-4)))
    mid_step = int(round(0.3 / (params.DT if tag == "coarse" else 1e-4)))
    tr = np.zeros((n_steps + 1, len(trace_names)))
    for k in range(n_steps):
        V = nw.step(u)
        tr[k + 1] = V[tidx]
        if k + 1 == mid_step:
            mids[tag] = V.copy()
    finals[tag] = nw.compiled.V.copy()
    traces[tag] = tr
    rg = names.index("RG_E_r")
    print(f"{tag}: V range final [{finals[tag].min():.3f}, "
          f"{finals[tag].max():.3f}] mV; RG_E_r end {finals[tag][rg]:.4f} mV")

dv = np.abs(finals["fine"] - finals["coarse"])
rg = names.index("RG_E_r")
print(f"fine-vs-coarse (integrator drift): max {dv.max():.3f} mV at "
      f"{names[int(np.argmax(dv))]}, RG_E_r {abs(finals['fine'][rg] - finals['coarse'][rg]):.3f} mV")

out_dir = HERE.parents[1] / "Matlab" / "SNS_Simscape" / "results"
savemat(out_dir / "verify_ref.mat",
        {"names": names, "V_final": finals["fine"],
         "V_final_coarse": finals["coarse"],
         "V_mid_coarse": mids["coarse"],
         "traces": traces["fine"], "traces_coarse": traces["coarse"],
         "trace_names": trace_names,
         "t": np.arange(int(2.0 / 1e-4) + 1) * 1e-4,
         "t_coarse": np.arange(int(2.0 / params.DT) + 1) * params.DT,
         "drive_nA": 2.5})
print(f"wrote {out_dir / 'verify_ref.mat'}")
