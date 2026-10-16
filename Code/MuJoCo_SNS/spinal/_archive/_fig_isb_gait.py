"""ISB/OpenSim-convention gait figures, BOTH legs (Ben 2026-09-20:
"I want axes in the ISB standard used by OpenSim... gait cycles for
both legs better matched the OpenSim (mean, amplitude, phase,
periodicity)").

Fig 1  curr_gait_isb.png — mean gait cycles, 3 joints x 2 legs:
       OpenSim reference (black) vs sim (right leg red, left leg blue).
       OpenSim joint conventions: hip_flexion + = flexion;
       knee_angle NEGATIVE = flexion; ankle_angle + = dorsiflexion.
       Global frame of the underlying data is the OpenSim/ISB frame
       (X anterior, Y up, Z lateral-right) preserved by the converter.
Fig 2  curr_gait_contact.png — walk-window time series: hip/knee/ankle
       both legs + per-foot vertical contact force (N), stance shading.

Reads spinal_run.npz (run any winner first, e.g. _diag_stage3.py).
Cycles/duty use the LOGGED heel+toe contact forces (kine_ref v2).
Usage: python _fig_isb_gait.py [out_prefix]
"""
import io
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import kine_ref as KR

HERE = Path(__file__).parent
FIGS = (HERE.parents[2] / "Documentation" / "Reports and Papers"
        / "Dissertation" / "CPG_airstepping_figs")

z = np.load("spinal_run.npz", allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
contact = z["contact"] if "contact" in z.files else None
kj = [str(s) for s in z["key_joints"]]
# EVAL runs use the runner's 16-s schedule: walk=(2.0, 13.0) — runner.py
# rewrites SCHEDULE for --eval (its walk[0] is 2.0, not params' 5.0)
WALK = 2.0
m = t >= WALK

ref = KR.ref_cached()
legs = {s: KR.sim_side(t, q, neuro, WALK, s, contact) for s in ("r", "l")}

prefix = sys.argv[1] if len(sys.argv) > 1 else "curr_gait"
plt.rcParams.update({"font.size": 10, "axes.grid": True,
                     "grid.alpha": 0.3})
JOINTS = ("hip", "knee", "ankle")
LABEL = {"hip": "Hip flexion [$\\deg$]  (+ flexion)",
         "knee": "Knee angle [$\\deg$]  ($-$ flexion)",
         "ankle": "Ankle angle [$\\deg$]  (+ dorsiflexion)"}
SIM_C = {"r": "#d62728", "l": "#1f77b4"}

fig, ax = plt.subplots(3, 2, figsize=(9.5, 8.5), sharex=True)
for col, side in enumerate(("r", "l")):
    mean = legs[side][0]
    for row, j in enumerate(JOINTS):
        a = ax[row, col]
        a.plot(KR.GRID, ref[side][j], "k-", lw=2.0, label="OpenSim ref")
        if mean is not None:
            a.plot(KR.GRID, mean[j], color=SIM_C[side], lw=1.8,
                   label="SNS sim")
        a.set_ylabel(LABEL[j], fontsize=9)
        if row == 0:
            a.set_title(f"{'Right' if side == 'r' else 'Left'} leg",
                        fontsize=11, fontweight="bold")
            a.legend(fontsize=8, loc="upper right")
        if row == 2:
            a.set_xlabel("gait cycle [%]  (0 = loading response)")
d = {s: legs[s] for s in ("r", "l")}
txt = (f"OpenSim/ISB conventions: hip +flexion, knee $-$flexion, "
       f"ankle +dorsi; global X anterior, Y up, Z right.  "
       f"duty r/l {d['r'][2]:.2f}/{d['l'][2]:.2f} "
       f"(ref {ref['duty_r']:.2f}/{ref['duty_l']:.2f})")
fig.text(0.5, 0.005, txt, ha="center", fontsize=8)
fig.suptitle("Mean gait cycles vs OpenSim reference — both legs", y=0.995)
fig.tight_layout(rect=(0, 0.015, 1, 0.98))
fig.savefig(HERE / f"{prefix}_gait_isb.png", dpi=200)
fig.savefig(FIGS / f"{prefix}_gait_isb.png", dpi=200)
plt.close(fig)

# ---- Fig 2: time series + contact ----
fig2, ax2 = plt.subplots(4, 1, figsize=(10, 8.5), sharex=True)
for row, j in enumerate(JOINTS):
    a = ax2[row]
    for side in ("r", "l"):
        i = kj.index(f"{j if j != 'hip' else 'hip'}_flexion_{side}"
                     if j == "hip" else f"{j}_angle_{side}")
        a.plot(t[m], q[m, i], color=SIM_C[side], lw=1.2,
               label=f"{side.upper()}")
    a.set_ylabel(LABEL[j], fontsize=9)
    a.legend(fontsize=8, ncol=2, loc="upper right")
a = ax2[3]
if contact is not None:
    a.plot(t[m], contact[m, 0], color=SIM_C["r"], lw=1.0, label="R")
    a.plot(t[m], contact[m, 1], color=SIM_C["l"], lw=1.0, label="L")
    a.axhline(KR.LOAD_N, color="k", ls=":", lw=0.8)
    a.set_ylabel("foot contact force [N]")
    a.legend(fontsize=8, ncol=2, loc="upper right")
a.set_xlabel("t [s]  (walk window)")
fig2.suptitle("Walk-window joint angles and per-foot contact — "
              "both legs", y=0.995)
fig2.tight_layout(rect=(0, 0, 1, 0.98))
fig2.savefig(HERE / f"{prefix}_gait_contact.png", dpi=200)
fig2.savefig(FIGS / f"{prefix}_gait_contact.png", dpi=200)
plt.close(fig2)

for s in ("r", "l"):
    mean, per, duty_c, cfrac = legs[s][0], legs[s][1], legs[s][2], legs[s][3]
    T = np.mean(per) if per else float("nan")
    print(f"leg {s}: cycles {len(per)} duty {duty_c:.2f} T {T:.2f}s "
          f"contact frac {cfrac:.2f}")
print(f"saved {prefix}_gait_isb.png + {prefix}_gait_contact.png "
      f"(here + Dissertation figs)")
