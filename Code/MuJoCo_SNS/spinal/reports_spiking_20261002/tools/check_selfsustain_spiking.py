"""GATE 2 (SPIKING_MIRROR_PLAN.md): constant-DRIVE self-sustained rhythm.

Mirrors check_selfsustain.py: same winner tables (draw_circuit
.effective_tables('best')), same model, constant DRIVE = 2.5 nA, 20 s,
interleg on.  Builds BOTH the non-spiking network and the spiking mirror
and compares:
  - RG alternation period (readout-tap envelope for the spiking build;
    NaP membrane for the analog build)
  - E-duty of the RG_E_r envelope
  - raw within-burst firing rate of the spiking RG cells
PASS = both builds sustain the rhythm (>= 3 antiphase peaks in the last
10 s) AND the spiking period is within 20% of the non-spiking period
(the plan's tolerance).

Usage: D:\\Anaconda\\envs\\myo\\python.exe check_selfsustain_spiking.py
"""
import io
import os
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

SPINAL = Path(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
os.chdir(SPINAL)
sys.path.insert(0, str(SPINAL))
assert "AARL_NET" not in os.environ

import numpy as np
from scipy.signal import find_peaks

import mujoco

import params
import build_network as bn
import build_network_spiking as bs
from draw_circuit import effective_tables

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


def run_net(net, label, count_spikes=False):
    u = net.make_inputs()
    u[net.input_index("DRIVE")] = DRIVE
    n = int(round(T_END / params.DT))
    v = np.zeros(n + 1)
    cnt = None
    if count_spikes:
        # cumulative spike-count series of the RAW RG cells -> the burst
        # envelope (the readout taps ripple at within-phase spike rate
        # and contaminate peak finding; the 100 ms spike envelope is the
        # honest phase signal for a spiking half-center)
        cnt = np.zeros((n + 1, 2), dtype=float)
        c_e, c_f = net.ridx["RG_E_r"], net.ridx["RG_F_r"]
        base = net.spike_counts.copy()
    for k in range(n):
        V = net.step(u)
        v[k + 1] = V[net.idx["RG_E_r"]] - V[net.idx["RG_F_r"]]
        if count_spikes:
            cnt[k + 1, 0] = net.spike_counts[c_e] - base[c_e]
            cnt[k + 1, 1] = net.spike_counts[c_f] - base[c_f]
    tt = np.arange(n + 1) * params.DT
    if count_spikes:
        win = int(round(0.100 / params.DT))
        env = np.convolve(np.diff(cnt[:, 0], prepend=0.0)
                          - np.diff(cnt[:, 1], prepend=0.0),
                          np.ones(win) / win, mode="same")
        env *= 1000.0 * params.DT / 0.100   # -> spikes/s scale
        pk, _ = find_peaks(env, prominence=1.0, distance=int(0.15 / params.DT))
        tail = env[-int(round(10.0 / params.DT)):]
        pk10, _ = find_peaks(tail, prominence=1.0,
                             distance=int(0.15 / params.DT))
        periods = np.diff(tt[pk]) if len(pk) >= 2 else np.array([])
        print(f"[{label}] E-F spike-envelope peaks {len(pk)}, last-10s "
              f"{len(pk10)}; period "
              f"{periods.mean():.3f} +/- {periods.std():.3f} s "
              f"({1.0 / max(periods.mean(), 1e-9):.2f} Hz)")
    else:
        env = v
        pk, _ = find_peaks(env, prominence=0.3)
        tail = env[-int(round(10.0 / params.DT)):]
        pk10, _ = find_peaks(tail, prominence=0.3)
        periods = np.diff(tt[pk]) if len(pk) >= 2 else np.array([])
        print(f"[{label}] peaks total {len(pk)}, last-10s {len(pk10)}; "
              f"period {periods.mean():.3f} +/- {periods.std():.3f} s "
              f"({1.0 / max(periods.mean(), 1e-9):.2f} Hz) "
              f"swing {v.max() - v.min():.2f}")
    for w0 in (0.0, 5.0, 10.0, 15.0):
        msk = tt >= w0
        print(f"  window {w0:.0f}-{T_END:.0f}s: swing "
              f"{v[msk].max() - v[msk].min():.3f} mV")
    if count_spikes:
        sc = net.spike_counts
        print("  raw RG spike rates [Hz]: " + " ".join(
            f"{c}={sc[net.ridx[c]] / T_END:.1f}"
            for c in ("RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l")))
    return (len(pk10), float(periods.mean()) if len(pk) >= 2 else None,
            v)


net_ns = bn.build(acts, dt=params.DT, interleg=True)
n_ns, T_ns, v_ns = run_net(net_ns, "non-spiking")

net_sp = bs.build(acts, dt=params.DT, interleg=True)
n_sp, T_sp, v_sp = run_net(net_sp, "spiking", count_spikes=True)

ok_rhythm = n_ns >= 3 and n_sp >= 3
ok_period = (T_ns is not None and T_sp is not None and
             abs(T_sp - T_ns) / T_ns <= 0.20)
print(f"\nperiod: non-spiking {T_ns} s, spiking {T_sp} s "
      f"({None if not (T_ns and T_sp) else 100 * (T_sp - T_ns) / T_ns:+.1f}%)")
print(f"rhythm-sustained={ok_rhythm} period-within-20%={ok_period}")
print("GATE2 RHYTHM MIRROR:", "PASS" if (ok_rhythm and ok_period)
      else "FAIL")
sys.exit(0 if (ok_rhythm and ok_period) else 1)
