"""GOAL 4 figure data — replicate the GATE 2 constant-DRIVE recipe of
tools/check_selfsustain_spiking.py VERBATIM (same effective_tables("best"),
same DRIVE 2.5 nA, same 20 s, same networks) but additionally DUMP the
traces needed for the dissertation rhythm figure:

  goal4_gate2_traces.npz:
    t            time (s)
    v_ns         non-spiking V[RG_E_r] - V[RG_F_r] (mV, readout)
    cnt_e/cnt_f  cumulative spike counts of the RAW spiking RG_E_r / RG_F_r

The script prints the same summary lines as the gate script; the figure is
only valid if they match logs/gate2_rhythm2.log
(non-spiking 0.399 s, spiking 0.621 s).

Usage: D:\\Anaconda\\envs\\myo\\python.exe tools\\goal4_dump_gate2_traces.py
"""
import io
import os
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

HERE = Path(__file__).resolve().parent
SPINAL = HERE.parents[1]
os.chdir(SPINAL)
sys.path.insert(0, str(SPINAL))
assert "AARL_NET" not in os.environ

import numpy as np  # noqa: E402
from scipy.signal import find_peaks  # noqa: E402

import mujoco  # noqa: E402

import params  # noqa: E402
import build_network as bn  # noqa: E402
import build_network_spiking as bs  # noqa: E402
from draw_circuit import effective_tables  # noqa: E402

DRIVE = 2.5
T_END = 20.0

t = effective_tables("best")
params.G.clear(); params.G.update(t["G"])
params.TAU.clear(); params.TAU.update(t["TAU"])
for ph in params.W_PF_MN:
    params.W_PF_MN[ph].clear(); params.W_PF_MN[ph].update(t["W"][ph])
params.W_POSTURE.clear(); params.W_POSTURE.update(t["WPOST"])

m = mujoco.MjModel.from_xml_path(
    str(SPINAL.parents[2] / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
       / "mjc" / "gait2392_simbody" / "gait2392_simbody_cvt3.xml"))
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(m.nu)]

# ---- non-spiking (identical to the gate script) -------------------------
net = bn.build(acts, dt=params.DT, interleg=True)
u = net.make_inputs()
u[net.input_index("DRIVE")] = DRIVE
n = int(round(T_END / params.DT))
v = np.zeros(n + 1)
for k in range(n):
    V = net.step(u)
    v[k + 1] = V[net.idx["RG_E_r"]] - V[net.idx["RG_F_r"]]
tt = np.arange(n + 1) * params.DT
pk, _ = find_peaks(v, prominence=0.3)
tail = v[-int(round(10.0 / params.DT)):]
pk10, _ = find_peaks(tail, prominence=0.3)
periods = np.diff(tt[pk]) if len(pk) >= 2 else np.array([])
print(f"[non-spiking] peaks {len(pk)}, last-10s {len(pk10)}; "
      f"period {periods.mean():.3f} +/- {periods.std():.3f} s "
      f"({1.0/max(periods.mean(),1e-9):.2f} Hz)")
v_ns, T_ns = v, float(periods.mean())

# ---- spiking mirror (identical to the gate script) ----------------------
net = bs.build(acts, dt=params.DT, interleg=True)
u = net.make_inputs()
u[net.input_index("DRIVE")] = DRIVE
n = int(round(T_END / params.DT))
v = np.zeros(n + 1)
cnt = np.zeros((n + 1, 2), dtype=float)
c_e, c_f = net.ridx["RG_E_r"], net.ridx["RG_F_r"]
base = net.spike_counts.copy()
for k in range(n):
    V = net.step(u)
    v[k + 1] = V[net.idx["RG_E_r"]] - V[net.idx["RG_F_r"]]
    cnt[k + 1, 0] = net.spike_counts[c_e] - base[c_e]
    cnt[k + 1, 1] = net.spike_counts[c_f] - base[c_f]
tt = np.arange(n + 1) * params.DT
win = int(round(0.100 / params.DT))
env = np.convolve(np.diff(cnt[:, 0], prepend=0.0)
                  - np.diff(cnt[:, 1], prepend=0.0),
                  np.ones(win) / win, mode="same")
env *= 1000.0 * params.DT / 0.100
pk, _ = find_peaks(env, prominence=1.0, distance=int(0.15 / params.DT))
tail = env[-int(round(10.0 / params.DT)):]
pk10, _ = find_peaks(tail, prominence=1.0, distance=int(0.15 / params.DT))
periods = np.diff(tt[pk]) if len(pk) >= 2 else np.array([])
print(f"[spiking] E-F spike-envelope peaks {len(pk)}, last-10s {len(pk10)}; "
      f"period {periods.mean():.3f} +/- {periods.std():.3f} s "
      f"({1.0/max(periods.mean(),1e-9):.2f} Hz)")
sc = net.spike_counts
print("  raw RG spike rates [Hz]: " + " ".join(
    f"{c}={sc[net.ridx[c]] / T_END:.1f}"
    for c in ("RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l")))
T_sp = float(periods.mean())
print(f"period: non-spiking {T_ns:.6f} s, spiking {T_sp:.6f} s "
      f"({100 * (T_sp - T_ns) / T_ns:+.1f}%)")

np.savez_compressed(HERE / "goal4_gate2_traces.npz",
                    t=tt, v_ns=v_ns, cnt_e=cnt[:, 0], cnt_f=cnt[:, 1])
print(f"saved {HERE / 'goal4_gate2_traces.npz'}")
