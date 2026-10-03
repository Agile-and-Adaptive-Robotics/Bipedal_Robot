"""Full-net sweep of the RG SFA increment (thr_inc) for the rhythm gate:
the isolated RG pair (tune_rg_pair.py) runs ~0.34 s at thr_inc 4.0, but
the FULL network (commissurals, both legs) lengthens the cycle.  This
sweep measures the true E-F spike-envelope period of the full spiking
network at several thr_inc values and writes the winner into
spiking_calibration.json (key rg_thr_inc).

Usage: D:\\Anaconda\\envs\\myo\\python.exe sweep_thrinc_fullnet.py
"""
import io
import json
import os
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8",
                              errors="replace")

SPINAL = Path(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
os.chdir(SPINAL)
sys.path.insert(0, str(SPINAL))
assert "AARL_NET" not in os.environ

import numpy as np
from scipy.signal import find_peaks

import mujoco

import params
import build_network_spiking as bs
from draw_circuit import effective_tables

T = effective_tables("best")
params.G.clear(); params.G.update(T["G"])
params.TAU.clear(); params.TAU.update(T["TAU"])
for ph in params.W_PF_MN:
    params.W_PF_MN[ph].clear(); params.W_PF_MN[ph].update(T["W"][ph])
params.W_POSTURE.clear(); params.W_POSTURE.update(T["WPOST"])

m = mujoco.MjModel.from_xml_path(
    str(SPINAL.parents[2] / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
        / "mjc" / "gait2392_simbody" / "gait2392_simbody_cvt3.xml"))
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(m.nu)]

DRIVE, T_END = 2.5, 12.0
TARGET = 0.399

best = None
import itertools
for k_ns, inc in ((0.45, 4.0), (0.55, 4.0), (0.55, 3.0)):
    bs.CAL["rg_thr_inc"] = float(inc)
    bs.CAL["k_ns"] = float(k_ns)
    net = bs.build(acts, dt=params.DT, interleg=True)
    u = net.make_inputs()
    u[net.input_index("DRIVE")] = DRIVE
    n = int(round(T_END / params.DT))
    cnt = np.zeros((n + 1, 2))
    c_e, c_f = net.ridx["RG_E_r"], net.ridx["RG_F_r"]
    base = net.spike_counts.copy()
    for k in range(n):
        net.step(u)
        cnt[k + 1, 0] = net.spike_counts[c_e] - base[c_e]
        cnt[k + 1, 1] = net.spike_counts[c_f] - base[c_f]
    tt = np.arange(n + 1) * params.DT
    win = int(round(0.100 / params.DT))
    env = np.convolve(np.diff(cnt[:, 0], prepend=0.0)
                      - np.diff(cnt[:, 1], prepend=0.0),
                      np.ones(win) / win, mode="same")
    pk, _ = find_peaks(env, prominence=1.0, distance=int(0.15 / params.DT))
    per = float(np.mean(np.diff(tt[pk]))) if len(pk) >= 3 else float("nan")
    fe = float(cnt[-1, 0] / T_END); ff = float(cnt[-1, 1] / T_END)
    alt = fe > 1.5 and ff > 1.5 and np.isfinite(per)
    print(f"k_ns {k_ns} thr_inc {inc:4.1f}: period {per:6.3f} s (target {TARGET}), "
          f"E {fe:4.1f} Hz F {ff:4.1f} Hz peaks {len(pk)} alt={alt}",
          flush=True)
    if alt and abs(per - TARGET) < 0.20 * TARGET:
        if best is None or abs(per - TARGET) < abs(best[1] - TARGET):
            best = (inc, per, k_ns)

print("WINNER:", best)
if best is not None:
    p = SPINAL / "spiking_calibration.json"
    d = json.loads(p.read_text(encoding="utf-8"))
    d["rg_thr_inc"] = float(best[0])
    d["k_ns"] = float(best[2])
    d["rg_sweep_note"] = (
        f"k_ns 0.35 + rg_thr_inc {best[0]} set by the FULL-NET sweep "
        "(tools/sweep_thrinc_fullnet.py; period target 0.399 s = the "
        "analog 'best' build at DRIVE 2.5; measured "
        f"{best[1]:.3f} s)")
    p.write_text(json.dumps(d, indent=1), "utf-8")
    print("json updated")
