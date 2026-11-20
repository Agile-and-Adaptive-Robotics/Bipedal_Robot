"""Standalone calibration of the spinal network's rhythm layer (no MuJoCo).

Builds the full network against the synthetic 92-muscle list, drives it with
constant DRIVE (walking) + POSTURE (standing), and reports:
  - RG-E/RG-F oscillation period and duty per leg,
  - left-right antiphase quality,
  - PF burst window onsets within the half-cycle (phase staggering).

Usage:
    python check_rhythm.py [drive_nA] [duration_s]
Writes check_rhythm.png next to this file.
"""
import sys
from pathlib import Path

import numpy as np

import build_network as bn
from muscle_map import _GROUPS_BY_NAME

HERE = Path(__file__).parent


def actuator_names():
    names = []
    for side in ("r", "l"):
        for base in _GROUPS_BY_NAME:
            names.append(f"{base}_{side}")
    return names


def main():
    drive_nA = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
    dur = float(sys.argv[2]) if len(sys.argv) > 2 else 12.0
    dt = bn.DT

    net = bn.build(actuator_names(), dt=dt)
    print(f"neurons: {len(net.idx)}  inputs: {len(net.inputs)}")

    u = net.make_inputs()
    u[net.input_index("DRIVE")] = drive_nA
    u[net.input_index("POSTURE")] = 1.0

    nsteps = int(dur / dt)
    watch = ["RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l",
             "PF_E1_r", "PF_E2_r", "PF_F1_r", "PF_F2_r"]
    cols = {w: np.zeros(nsteps) for w in watch}

    for k in range(nsteps):
        v = net.step(u)
        for w in watch:
            cols[w][k] = v[net.idx[w]]

    t = np.arange(nsteps) * dt
    thr = 0.5 * max(cols["RG_E_r"].max(), 1e-6)

    def bursts(sig, level):
        on = sig > level
        starts = np.flatnonzero(on[1:] & ~on[:-1]) + 1
        ends = np.flatnonzero(~on[1:] & on[:-1]) + 1
        return starts, ends

    sE, eE = bursts(cols["RG_E_r"], thr)
    sF, eF = bursts(cols["RG_F_r"], thr)
    print(f"RG threshold: {thr:.2f} mV | RG_E_r max {cols['RG_E_r'].max():.2f} "
          f"RG_F_r max {cols['RG_F_r'].max():.2f}")
    if len(sE) > 1:
        per = np.diff(t[sE])
        duty = (t[eE[:len(per)]] - t[sE[:len(per)]]) / per if len(eE) else np.nan
        print(f"RG-E periods: mean {per.mean():.3f} s  {np.round(per, 3)}")
        print(f"duty (E fraction): {np.round(duty, 2)}")
    else:
        print("RG-E: no repeated bursts detected -> NOT oscillating")
    if len(sF) > 0 and len(sE) > 0:
        lag = t[sF[:3]] - t[sE[0]]
        print(f"RG-F burst starts relative to first RG-E start: {np.round(lag, 3)}")
        alt = cols["RG_E_l"][sE[0]:] , cols["RG_E_r"][sE[0]:]
        same = np.corrcoef(alt[0][:2000], alt[1][:2000])[0, 1] if nsteps > 2000 else np.nan
        print(f"left-right RG-E correlation (want < 0 for antiphase): {same:.2f}")

    # PF staggering: onset lag of each PF burst within the RG-E cycle
    for phase in ("E1", "E2", "F1", "F2"):
        s, _ = bursts(cols[f"PF_{phase}_r"], 0.5 * max(cols[f"PF_{phase}_r"].max(), 1e-6))
        if len(s):
            print(f"PF_{phase}_r first onsets: {np.round(t[s[:4]], 2)}")

    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(len(watch), 1, figsize=(9, 1.4 * len(watch)),
                               sharex=True)
        for a, w in zip(ax, watch):
            a.plot(t, cols[w], lw=0.8)
            a.set_ylabel(w, fontsize=7)
            a.grid(True, alpha=0.3)
        ax[-1].set_xlabel("t [s]")
        fig.suptitle(f"Rhythm check, DRIVE={drive_nA} nA")
        fig.tight_layout()
        out = HERE / "check_rhythm.png"
        fig.savefig(out, dpi=130)
        print(f"plot: {out}")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
