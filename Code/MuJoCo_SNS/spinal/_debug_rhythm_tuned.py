"""Debug gate (2026-09-20): does the TUNED network still self-sustain at
constant DRIVE on current code?  Loads the production tables exactly the
way the curriculum does (fitted baseline + v10 multipliers + the stage-2
air winner = the current best rhythm config), then runs network-only,
constant DRIVE, 20 s, and reports RG-E/RG-F antiphase swing per window.

PASS = last-5 s swing >= max(0.1 mV, 1/2 first-window swing).
Usage: python _debug_rhythm_tuned.py [drive] [duration_s]
"""
import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np
from scipy.signal import find_peaks

import _curriculum as C
import build_network as bn
import optuna_walk_v10 as OW
import params

HERE = Path(__file__).parent

drive = float(sys.argv[1]) if len(sys.argv) > 1 else None
T_END = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0

# --- production tables, exactly as _curriculum.set_stage does it ----------
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
stage2 = json.loads((HERE / "curriculum_stage2.json")
                    .read_text(encoding="utf-8"))["params"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
C.set_stage(2, {**mul, **stage2})
if drive is None:
    drive = float(stage2["drive"])
print(f"config: fitted baseline + v10 multipliers + stage2 winner "
      f"(drive {drive:.3f}, renshaw {params.G['renshaw']})")

acts = [f"{b}_{s}" for s in ("r", "l") for b in
        __import__("muscle_map")._GROUPS_BY_NAME]
net = bn.build(acts, dt=params.DT, interleg=True)
print(f"neurons: {len(net.idx)}  inputs: {len(net.inputs)}")

u = net.make_inputs()
u[net.input_index("DRIVE")] = drive
u[net.input_index("POSTURE")] = 1.0

n = int(round(T_END / params.DT))
v = np.zeros(n + 1)
ve = np.zeros(n + 1)
vf = np.zeros(n + 1)
for k in range(n):
    V = net.step(u)
    v[k + 1] = V[net.idx["RG_E_r"]] - V[net.idx["RG_F_r"]]
    ve[k + 1] = V[net.idx["RG_E_r"]]
    vf[k + 1] = V[net.idx["RG_F_r"]]
tt = np.arange(n + 1) * params.DT

pk, _ = find_peaks(v, prominence=0.3)
print("antiphase peak times (s):", np.round(pk * params.DT, 2)[:12])
if len(pk) > 2:
    per = np.diff(pk * params.DT)
    print(f"period: mean {per.mean():.3f} s  {np.round(per, 3)}")
sw0 = None
for w0 in (0.0, 5.0, 10.0, 15.0):
    if w0 >= T_END:
        break
    msk = tt >= w0
    sw = float(v[msk].max() - v[msk].min())
    if sw0 is None:
        sw0 = sw
    print(f"window {w0:.0f}-{T_END:.0f} s: swing {sw:.3f} mV "
          f"(E {ve[msk].max():.2f} F {vf[msk].max():.2f})")
last = tt >= max(0.0, T_END - 5.0)
sw_last = float(v[last].max() - v[last].min())
ok = sw_last >= max(0.1, 0.5 * (sw0 or 0.0))
print(f"SELF-SUSTAIN: {'PASS' if ok else 'FAIL'} "
      f"(last-5s swing {sw_last:.3f} mV vs first-window {sw0:.3f} mV)")
