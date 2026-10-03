import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
os.chdir(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, os.getcwd())
assert "AARL_NET" not in os.environ
import numpy as np, mujoco, time
import build_network_spiking as bs
m = mujoco.MjModel.from_xml_path(str(
    __import__("pathlib").Path(os.getcwd()).parents[2] / "Solid_Models" /
    "OpenSim" / "Gait2392_Robotbody" / "mjc" / "gait2392_simbody" /
    "gait2392_simbody_cvt3.xml"))
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(m.nu)]
t0=time.time()
net = bs.build(acts, interleg=True)
print(f"built: {len(net.idx)} neurons, {len(net.inputs)} inputs, "
      f"{len(net.net.connections)} synapses in {time.time()-t0:.1f}s; "
      f"n_sub={net.n_sub}")
u = net.make_inputs()
u[net.input_index("DRIVE")] = 2.5
v0 = None
t0=time.time()
for k in range(1000):
    V = net.step(u)
print(f"2 s sim in {time.time()-t0:.2f}s; finite={np.all(np.isfinite(V))}")
print(f"RG_E_r readout level {V[net.idx['RG_E_r']]:.3f} mV; raw spikes "
      f"{int(net.spike_counts[net.ridx['RG_E_r']])} in 2 s; "
      f"MN max {float(np.max([V[net.idx[n]] for n in net.mn_names.values()])):.3f} mV")
print(f"DRIVE V {V[net.idx['DRIVE']]:.3f}")
